"""
Step 5b — Ipopt feasibility solve.

`build_model!`/`build_model_full!` builds a square system of levels-form
equality constraints (see PLAN.md for the %-change-to-levels conversion) with
no meaningful objective (a CGE benchmark/shock solve is a feasibility problem,
not an optimization), so `solve_model!` sets a constant zero objective and
tunes Ipopt's NLP solver settings for a large, sparse, square system before
calling `JuMP.optimize!`.

Diagnostics default to `print_level=5` (per-iteration convergence table only).
Higher levels are dangerous at this model's scale: empirically `print_level=8`
already dumps full per-iteration vectors (`new vars[i]`, `curr_c[i]`,
`curr_x[1][i]`, `final y_c/y_d/z_L/z_U` for all ~2.4M variables) — not just
`>=11` as Ipopt's docs imply — ballooning the log to gigabytes in minutes.
Keep the default at 5 unless a specific solve genuinely needs the extra detail.
Ipopt's own iteration log is mirrored to `output_file` on disk;
Ipopt writes that file directly (flushed per iteration) instead of going
through Julia's stdout, which is often block-buffered when redirected to a
file — so `output_file` is the one that actually shows live progress while a
long solve is still running.
"""
function solve_model!(m::JuMP.Model;
                       max_iter::Int=3000,
                       tol::Float64=1.0e-6,
                       constr_viol_tol::Float64=1.0e-6,
                       mu_strategy::String="adaptive",
                       nlp_scaling_method::String="gradient-based",
                       bound_push::Float64=1.0e-2,
                       bound_frac::Float64=1.0e-2,
                       print_level::Int=5,
                       output_file::String="ipopt_log.txt")
    JuMP.set_attribute(m, "max_iter", max_iter)
    JuMP.set_attribute(m, "tol", tol)
    JuMP.set_attribute(m, "constr_viol_tol", constr_viol_tol)
    JuMP.set_attribute(m, "mu_strategy", mu_strategy)
    # ── NLP scaling. Ipopt's default "gradient-based" rescales each constraint
    #    so its max gradient element ≈ nlp_scaling_max_gradient (100). For this
    #    levels model a handful of constraints have tiny gradients (near-index
    #    terms), so gradient-based scaling AMPLIFIES their residual by up to
    #    ~1e8×: a benchmark start whose true worst residual is 1.8e-4 is
    #    reported by Ipopt as inf_pr≈2.3e4, and the badly-scaled rows then stall
    #    the solve (inf_du stuck ~1e9, oscillating inf_pr). Passing "none" lets
    #    Ipopt see the true (feasible) infeasibility so it converges directly.
    JuMP.set_attribute(m, "nlp_scaling_method", nlp_scaling_method)
    # Keep the benchmark start point undisturbed: default bound_push/bound_frac
    # (1e-2) shove every variable 1% off its bounds before iter 0, needlessly
    # perturbing an already-feasible start; callers may pass a smaller value.
    JuMP.set_attribute(m, "bound_push", bound_push)
    JuMP.set_attribute(m, "bound_frac", bound_frac)
    JuMP.set_attribute(m, "print_level", print_level)
    JuMP.set_attribute(m, "print_user_options", "yes")
    JuMP.set_attribute(m, "print_timing_statistics", "yes")
    JuMP.set_attribute(m, "output_file", output_file)

    if objective_sense(m) == MOI.FEASIBILITY_SENSE
        @objective(m, Min, 0.0)
    end

    JuMP.optimize!(m)

    status = JuMP.termination_status(m)
    ok = status in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED)
    return (status=status, solved=ok, primal_status=JuMP.primal_status(m))
end
