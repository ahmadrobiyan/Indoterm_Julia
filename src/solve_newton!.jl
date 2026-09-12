"""
Step 5c — square-system Newton solve (the model's real solver).

**Why not Ipopt.** A GEMPACK-style CGE closure produces a *square* system of
nonlinear equalities: once the closure fixes the exogenous variables, the number
of free variables equals the number of equations and there is nothing left to
optimize. Ipopt is an interior-point **optimizer**; handed a square system it has
zero degrees of freedom, which is degenerate for it — empirically it hangs in the
very first KKT factorization (56 minutes at `iter 0`, 5.6 GB, no iteration output,
after printing `Too few degrees of freedom (n_x, n_c)`). A square `F(x) = 0` wants
a Newton solver, which is what this file provides. `solve_model!.jl` (Ipopt) is
kept for reference and for any future genuinely-under-determined variant.

**Scaling is mandatory, not a refinement.** This model's Jacobian rows span ~22
orders of magnitude (≈2.7e-10 to 1.7e12): value identities carry ~1e6
coefficients while price equations are O(1). Without equilibration an absolute
residual of 1e-4 is meaningless float noise on one row and a real error on
another, `lu` reports `SingularException` on a perfectly good matrix, and
SuiteSparseQR's rank estimate is garbage (its tolerance keys off the largest
column norm). Every row is therefore divided by its largest Jacobian entry at the
start point, and residual/step/convergence all use those units — so `tol` is
effectively a RELATIVE tolerance.

**Squaring.** After the equation-level fixes in `build_equations.jl` the model is
naturally square, but two degenerate cases are still handled defensively because
they are silent and fatal:
  * *orphan variables* — free variables appearing in no constraint at all are
    pinned at their benchmark. Today there are none. Until 2026-07-31 there was
    exactly one, `natfhou`, and this defensive step is precisely what hid the
    fact that its equation (TERM.TAB:2053) had never been written: pinning an
    orphan and implementing its equation give the same answer for every closure
    that leaves the variable endogenous. `coalprice.CMF` swaps it exogenous, the
    orphan vanished, and the model came out one variable short of square. **An
    orphan reported here is a missing equation until proven otherwise** — do not
    treat the pin as a fix;
  * *dead rows* — constraints all of whose variables are fixed, i.e. `0 == 0`
    (today none). They are satisfied but contribute zero Jacobian rows, so they
    are deleted. Note Ipopt's `dependency_detector` does **not** remove these.
If the system is still not square afterwards, we raise with a diagnostic instead
of silently solving the wrong problem.
"""

using SparseArrays
using LinearAlgebra

# ── CGNR: Conjugate Gradient on Normal Residuals ──────────────────────────
# Solves  min_d ||J*d - b||² + μ||d||²  using CGNR applied to the augmented
# system [J; sqrt(μ)I] * d ≈ [b; 0].  μ = 0 gives the standard minimum-norm
# least-squares solution.  μ > 0 is Levenberg-Marquardt / Tikhonov
# regularization: it bounds the singular values of the normal-equation matrix
# away from zero, making CGNR converge fast on ill-conditioned Jacobians.
# Uses pre-allocated work vectors to avoid GC overhead.
function cgnr!(d::Vector{Float64}, J::SparseMatrixCSC{Float64,Int},
               b::AbstractVector{Float64};
               mu::Float64=0.0, maxit::Int=200, tol::Float64=1e-10,
               work_r::Vector{Float64}=zeros(size(J,1)),
               work_z::Vector{Float64}=zeros(size(J,2)),
               work_p::Vector{Float64}=zeros(size(J,2)),
               work_rd::Vector{Float64}=zeros(size(J,2)),
               work_qb::Vector{Float64}=zeros(size(J,1)))
    n = size(J, 2); m = size(J, 1)
    Jt = J'
    fill!(d, 0.0)
    fill!(work_rd, 0.0)
    copyto!(work_r, b)          # r_b = b (d=0)
    mul!(work_z, Jt, work_r)    # z = J'*r_b (gradient of ||J*d - b||²)
    copyto!(work_p, work_z)     # p = z
    γ = dot(work_z, work_z)
    γ0 = γ
    r_d = work_rd
    q_b = work_qb
    for _ in 1:maxit
        mul!(q_b, J, work_p)    # q_b = J*p
        qn = dot(q_b, q_b)
        if mu > 0
            @inbounds for i in 1:n
                qn += mu * work_p[i] * work_p[i]
            end
        end
        α = γ / qn
        @inbounds for i in 1:n
            d[i] += α * work_p[i]
            r_d[i] -= α * sqrt(mu) * work_p[i]
        end
        @inbounds for i in 1:m
            work_r[i] -= α * q_b[i]
        end
        # z = J'*r_b + sqrt(μ)*r_d
        mul!(work_z, Jt, work_r)
        if mu > 0
            s = sqrt(mu)
            @inbounds for i in 1:n
                work_z[i] += s * r_d[i]
            end
        end
        γnext = dot(work_z, work_z)
        if sqrt(γnext) < tol * sqrt(γ0); break; end
        β = γnext / γ
        @inbounds for i in 1:n
            work_p[i] = work_z[i] + β * work_p[i]
        end
        γ = γnext
    end
    return d, sqrt(γ) / sqrt(γ0)
end

"""
Result of [`solve_newton!`](@ref).

`status` is `:converged`, `:maxit`, `:no_progress`, or `:linear_solve_failed`.
`residual` is the **scaled** (relative) `||F||_inf`; `raw_residual` is in the
model's own units and is expected to look large (~1e-4) purely because of the
1e6-magnitude value identities — judge convergence on `residual`.
`values` maps each variable name to its solved level, shaped like `vars`.
"""
struct NewtonResult
    status::Symbol
    solved::Bool
    iters::Int
    residual::Float64
    raw_residual::Float64
    values::Dict{String,Any}
    n_free::Int
    n_eqs::Int
end

# Collect the free-variable ids referenced by a constraint function.
function _collect_free!(buf::Vector{Int}, f, freeid::Dict{Int64,Int})
    if f isa VariableRef
        id = get(freeid, f.index.value, 0); id != 0 && push!(buf, id)
    elseif f isa GenericAffExpr
        for (v, _) in f.terms
            id = get(freeid, v.index.value, 0); id != 0 && push!(buf, id)
        end
    elseif f isa GenericQuadExpr
        for (v, _) in f.aff.terms
            id = get(freeid, v.index.value, 0); id != 0 && push!(buf, id)
        end
        for (p, _) in f.terms
            id = get(freeid, p.a.index.value, 0); id != 0 && push!(buf, id)
            id = get(freeid, p.b.index.value, 0); id != 0 && push!(buf, id)
        end
    elseif f isa GenericNonlinearExpr
        for a in f.args; _collect_free!(buf, a, freeid); end
    elseif f isa AbstractArray
        for a in f; _collect_free!(buf, a, freeid); end
    end
    return buf
end

_rhs_of(s) = s isa MOI.EqualTo ? s.value :
             s isa MOI.LessThan ? s.upper :
             s isa MOI.GreaterThan ? s.lower : 0.0

"""
    solve_newton!(m, vars; maxit=25, tol=1e-8, verbose=true)

Solve the square levels system `F(x) = 0` defined by `m` (built by
`build_model_full!` and closed by `initialize_model!`) with a damped sparse
Newton method, and return a [`NewtonResult`](@ref).

`tol` is a **relative** tolerance on the equilibrated residual. ~1e-8 is the
realistic double-precision floor for ~90k equations built from ~1e6-magnitude
data: below it the residual is dominated by floating-point evaluation error, the
line search can find no improving step, and a stricter `tol` mislabels an
already-converged solve as `:no_progress`.

Mutates `m`: pins orphan variables and deletes dead rows (see the module
docstring), and writes the solution back as each variable's start value so a
subsequent solve warm-starts from it.
"""


function solve_newton!(m::JuMP.Model, vars::Dict{String,Any};
                       maxit::Int=25, tol::Float64=1e-8, verbose::Bool=true,
                       diagnose_columns::Bool=false, dry_run::Bool=false,
                       linsolve::Symbol=:lu)
    log(msg) = verbose && println(msg)

    # ── index the free variables ────────────────────────────────────────────
    freeid = Dict{Int64,Int}(); freeref = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value
        haskey(freeid, iv) && continue
        push!(freeref, vr); freeid[iv] = length(freeref)
    end

    # ── scan constraints: find dead rows, and which free vars are referenced ─
    seen = falses(length(freeref))
    deadcons = Any[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            buf = _collect_free!(Int[], JuMP.constraint_object(con).func, freeid)
            if isempty(buf)
                push!(deadcons, con)              # all-fixed `0 == 0` row
            else
                for id in buf; seen[id] = true; end
            end
        end
    end

    # ── square the system defensively ───────────────────────────────────────
    orphans = [freeref[i] for i in eachindex(freeref) if !seen[i]]
    for vr in orphans
        sv = JuMP.start_value(vr)
        JuMP.fix(vr, sv === nothing ? 0.0 : sv; force=true)
    end
    for con in deadcons; JuMP.delete(m, con); end
    if !isempty(orphans) || !isempty(deadcons)
        log("  squaring: pinned $(length(orphans)) orphan var(s) " *
            "$(isempty(orphans) ? "" : "[" * join(name.(orphans[1:min(end,5)]), ", ") * "]") " *
            ", deleted $(length(deadcons)) dead row(s)")
    end

    # ── the Newton unknowns are whatever is still free ──────────────────────
    # fam_nowfree parallels nowfree: the vars-dict family key per free column.
    # (JuMP.name(vr) would give "xtradmar[1,...]" — the indexed name, which the
    # Schur plan cannot match against family keys. A wrong namespace here fails
    # SILENTLY: empty S, plan nothing, every iteration degrades to CGNR.)
    nowfree = VariableRef[]
    fam_nowfree = String[]
    for (nm, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        if JuMP.is_fixed(vr); continue; end
        push!(nowfree, vr); push!(fam_nowfree, String(nm))
    end
    N = length(nowfree)

    # ── nonlinear evaluator over ALL variables (fixed ones enter as constants;
    #    passing only free vars errors, since fixed vars appear in expressions) ─
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

    nlmodel = MOI.Nonlinear.Model()
    rhs = Float64[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            r = _rhs_of(co.set)
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
            push!(rhs, r)
        end
    end
    NC = length(rhs)
    log("  system: $NC equations, $N unknowns" * (NC == N ? " (square)" : " ⚠ NOT SQUARE"))
    if NC != N
        error("""
              Newton needs a square system but got $NC equations vs $N free variables
              (difference $(N - NC)). More free variables than equations means the
              closure is under-determined — localize with a bipartite matching over
              (constraints × free vars); unmatched variable families name the culprit.
              Fewer means over-determination — usually redundant/dead rows.
              See PLAN.md "Gate #5c RESOLVED" for the diagnosis procedure.
              """)
    end

    evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(evaluator, [:Grad, :Jac])

    x = Vector{Float64}(undef, length(allv))
    for (k, v) in enumerate(allv)
        x[k] = JuMP.is_fixed(v) ? JuMP.fix_value(v) :
               (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))
    end
    freecols = [pos[JuMP.index(v).value] for v in nowfree]
    lb = [JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf for v in nowfree]

    st = MOI.jacobian_structure(evaluator)
    jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
    Jval = zeros(length(jrows_all))
    col2free = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols); col2free[c] = i; end
    keep = [col2free[c] != 0 for c in jcols_all]
    jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

    gbuf = zeros(NC)
    function residual!(g, xv)
        MOI.eval_constraint(evaluator, g, xv)
        @inbounds for i in eachindex(g); g[i] -= rhs[i]; end
        return g
    end

    # ── row equilibration (computed ONCE, then frozen) ──────────────────────
    # Jacobian row scales span many orders of magnitude, so the residual has to
    # be measured in scaled units. But the scaling must be FROZEN: it is the
    # yardstick that ‖F‖ is measured against, and the line search, the
    # convergence test and the iteration-to-iteration comparison all compare
    # readings taken at different points. Recomputing it each iteration changed
    # the length of the ruler between readings, so a step could look like
    # progress purely because its rows had been rescaled. Conditioning of the
    # linear solve does not depend on this choice — Ruiz re-equilibrates the
    # scaled Jacobian from scratch every iteration anyway.
    function compute_rowscale!(rs, xv)
        MOI.eval_constraint_jacobian(evaluator, Jval, xv)
        v0 = Jval[keep]
        fill!(rs, 1.0)
        mx = zeros(Float64, NC)
        @inbounds for k in eachindex(jr)
            a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
        end
        @inbounds for i in 1:NC; rs[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
        return rs
    end
    rowscale = compute_rowscale!(ones(NC), x)
    sres(xv, buf, rs) = (residual!(buf, xv); buf ./ rs)

    # ── Ruiz two-sided equilibration ────────────────────────────────────────
    # Row-only scaling can make mixed-derivative rows (e.g. value identities that
    # contain both a 1e6 flow coefficient and a 1.0 shifter coefficient) hide
    # weak columns after scaling. Ruiz scales rows and columns iteratively so
    # every row and column has infinity-norm 1. This removes the artificial
    # rank deficiency seen in the row-only scaled Jacobian.
    function ruiz_scales(J0::SparseMatrixCSC{Float64,Int};
                         maxiter::Int=20, tol::Float64=1e-6)
        m, n = size(J0)
        Dr = ones(Float64, m); Dc = ones(Float64, n)
        Jwork = copy(J0)
        for it in 1:maxiter
            # row scaling
            rowmax = zeros(Float64, m)
            @inbounds for k in 1:nnz(Jwork)
                i = Jwork.rowval[k]
                a = abs(Jwork.nzval[k])
                a > rowmax[i] && (rowmax[i] = a)
            end
            @inbounds for i in 1:m
                s = rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0
                Dr[i] *= s
            end
            @inbounds for k in 1:nnz(Jwork)
                i = Jwork.rowval[k]
                Jwork.nzval[k] /= (rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0)
            end
            # column scaling
            colmax = zeros(Float64, n)
            @inbounds for j in 1:n
                for k in Jwork.colptr[j]:(Jwork.colptr[j+1]-1)
                    a = abs(Jwork.nzval[k])
                    a > colmax[j] && (colmax[j] = a)
                end
            end
            @inbounds for j in 1:n
                s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0
                Dc[j] *= s
            end
            @inbounds for j in 1:n
                s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0
                for k in Jwork.colptr[j]:(Jwork.colptr[j+1]-1)
                    Jwork.nzval[k] /= s
                end
            end
            # convergence check
            dev = 0.0
            @inbounds for i in 1:m
                s = rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0
                dev = max(dev, abs(s - 1.0))
            end
            @inbounds for j in 1:n
                s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0
                dev = max(dev, abs(s - 1.0))
            end
            dev < tol && break
        end
        return Dr, Dc
    end

    gs = sres(x, gbuf, rowscale)
    log("  initial scaled ||F||_inf = $(norm(gs, Inf))")

    if dry_run || diagnose_columns
        # Evaluate Jacobian at current point and compute column norms
        MOI.eval_constraint_jacobian(evaluator, Jval, x)
        vals = Jval[keep]
        J = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
        colnrm = [norm(J[:, j]) for j in 1:N]
        nfree = count(colnrm .< 1e-10)
        if diagnose_columns
            weak = findall(colnrm .< 1e-10)
            pc = Dict{String,Int}()
            for wi in weak
                pre = split(JuMP.name(nowfree[wi]), "[")[1]
                pc[pre] = get(pc, pre, 0) + 1
            end
            println("  [diagnose] weak cols (<1e-10): $(nfree) / $N")
            for (pre, cnt) in sort(collect(pc), by=x->-x[2])
                println("    $(rpad(pre,18)) $cnt")
            end
            for wi in weak[1:min(5, length(weak))]
                println("    $(JuMP.name(nowfree[wi]))  colnrm=$(round(colnrm[wi], sigdigits=2))")
            end
        end
        if dry_run
            return NewtonResult(:dry_run, false, 0, norm(gs, Inf), norm(gbuf, Inf),
                               Dict{String,Any}(), N, NC)
        end
    end

    status = :maxit; iters = 0
    # A benchmark solve starts AT the solution by construction, so check before
    # stepping. Taking a step anyway amplifies the ~1e-9 float-noise residual
    # through the Jacobian into a ~1e-3 drift off the benchmark: the line search
    # accepts on the 2-norm while the inf-norm — the actual convergence
    # criterion — gets worse.
    if norm(gs, Inf) < tol
        log("  start point already satisfies F(x)=0 to tol — no step taken")
        status = :converged
    end
    # Pre-allocate CGNR work vectors (reused across iterations)
    cgnr_r = zeros(NC)
    cgnr_z = zeros(N)
    cgnr_p = zeros(N)
    cgnr_d = zeros(N)
    cgnr_rd = zeros(N)
    cgnr_qb = zeros(NC)

    mu = 0.1  # Levenberg-Marquardt regularization (adapted per iteration)

    # ── merit function and non-monotone history ─────────────────────────────
    # Steps are ACCEPTED on ½‖F‖₂², not on ‖F‖∞. With 83k rows a single
    # high-curvature row can veto a step that improves all the others: measured
    # at the 25x6 shock, `gret[i,d] = log(pcap+τ) − log(pinvitot+τ)` for
    # Livestock/BaliNusa picks up ½(0.046)² ≈ 1.1e-3 of pure log-curvature at
    # α=1, more than the entire starting residual of 3.0e-4, because `xcap` and
    # `xlnd` are exogenous in the base closure and the fixed factors absorb the
    # whole adjustment in their price. A monotone ‖F‖∞ test therefore rejected
    # α=1 and settled on α=0.25, which then shrank every iteration — linear
    # convergence with a decaying rate, i.e. the observed stall. The 2-norm
    # averages over rows, so one such row cannot veto. ‖F‖∞ remains the
    # CONVERGENCE criterion; only the ACCEPTANCE criterion changes.
    merit(g) = 0.5 * sum(abs2, g)
    NM_WINDOW = 5                # Grippo non-monotone window
    mem = Float64[merit(gs)]     # merits of accepted iterates
    Delta = Inf                  # trust radius on ‖d‖∞ (first LU step is trusted)
    # A non-monotone search is allowed to wander uphill, so the best point seen
    # must be remembered explicitly or it can be walked away from and lost.
    x_best = copy(x); gs_best = copy(gs); inf_best = norm(gs, Inf)
    for it in 1:(status == :converged ? 0 : maxit)
        iters = it
        MOI.eval_constraint_jacobian(evaluator, Jval, x)
        vals = Jval[keep]     # rowscale is deliberately NOT recomputed — see above
        Jbase = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)

        # Two-sided Ruiz equilibration on top of row scaling.  Dr/Dc are the
        # additional row/col multipliers that make every row and column of
        # J_s = Dr * Jbase * Dc have infinity-norm 1.
        Dr_ruiz, Dc_ruiz = ruiz_scales(Jbase)
        J_s = sparse(jr, jc,
                     [vals[k] / (rowscale[jr[k]] * Dr_ruiz[jr[k]] * Dc_ruiz[jc[k]])
                      for k in eachindex(jr)], NC, N)
        colnrm = [norm(J_s[:, j]) for j in 1:N]
        # right-hand side in Ruiz-scaled units: gs is already divided by rowscale
        rhs_s = gs ./ Dr_ruiz

        f0_inf = norm(gs, Inf)
        f0_m   = merit(gs)
        rel_cap = 0.10
        local d = zeros(N)
        accepted = false; alpha = 0.0; cgnr_rel = Inf

        # ── Direct sparse LU, with CGNR/Levenberg-Marquardt as fallback ─
        # `linsolve=:lu` (default): exact historical behaviour. `linsolve=:schur`
        # routes the first-try direction through `schur_linsolve.jl` (exact
        # Jacobian-level condensation + refinement); the `exact` lin_rel check
        # below validates either path identically, and any Schur failure falls
        # into the same CGNR fallback.
        # Measured at 25x6 (83,541 square, κ(J_s) ≈ 2.6e10): lu(J_s) costs 1.8s
        # and returns ‖J_s·dy + rhs_s‖/‖rhs_s‖ ≈ 4e-10, whereas CGNR at the 200
        # iterations budgeted below returns max|dy| = 0.025 against a true 1221 —
        # CGNR's convergence rate goes as κ², so it was never computing the Newton
        # direction the line search was then evaluating. LU is tried first each
        # iteration; CGNR with escalating mu remains the fallback for a failed
        # factorization or for scales where the fill-in is prohibitive (see the
        # region-scaling measurements: LU fill-in is the wall beyond ~20 regions,
        # which is what the Excerpt 49 condensation is for).
        F_lu = nothing
        if linsolve == :lu
            try
                F_lu = lu(J_s)
            catch e
                log("  it $it: lu(J_s) failed — $(sprint(showerror, e)); using CGNR")
            end
        end

        # Start with current mu. If the step does not improve the residual,
        # increase mu (more gradient-like, shorter step) and retry.
        for lmtry in 1:5
            t_solve = time()
            local dy_try
            # `exact` marks a verified Newton direction: it earns the full step
            # (no magnitude cap), since capping an exact direction destroys the
            # very property the line search relies on.
            exact = false
            if (F_lu !== nothing || linsolve == :schur) && lmtry == 1
                try
                    if linsolve == :schur
                        dy_schur = schur_linsolve(J_s, -rhs_s, fam_nowfree; verbose=false)
                        if dy_schur === nothing
                            # LOUD: a requested Schur path that never engages is a
                            # silent CGNR downgrade (bit-identical stalls). Never log().
                            println("  it $it: Schur path unavailable — falling back to CGNR " *
                                    "(S match failed; check family namespace)")
                            flush(stdout)
                            error("Schur path unavailable for this iteration")
                        end
                        dy_try = dy_schur
                    else
                        dy_try = F_lu \ (-rhs_s)
                    end
                    lin_rel = norm(J_s * dy_try .+ rhs_s) / max(norm(rhs_s), eps())
                    cgnr_rel = lin_rel
                    exact = isfinite(lin_rel) && lin_rel < 1e-6 && all(isfinite, dy_try)
                    exact || log("  it $it: LU step rejected (lin_rel=$(round(lin_rel;sigdigits=3))), using CGNR")
                    # Trust region. An accurately-solved direction is still not a
                    # usable one once J goes near-singular: the step is then
                    # dominated by the near-null mode and the line search can only
                    # scale it, never redirect it. Demote such a step to the
                    # damped path, which changes the direction.
                    if exact
                        dmag = maximum(abs, dy_try ./ Dc_ruiz)
                        if dmag > Delta
                            exact = false
                            log("  it $it: LU step |d|=$(round(dmag;sigdigits=4)) exceeds trust radius " *
                                "$(round(Delta;sigdigits=4)) — using damped CGNR")
                        end
                    end
                catch e
                    dy_try = zeros(N); cgnr_rel = Inf
                    log("  it $it: LU solve failed — $(sprint(showerror, e))")
                end
            end
            if !exact
                try
                    cgnr_rel = cgnr!(cgnr_d, J_s, -rhs_s; mu=mu, maxit=200, tol=1e-10,
                                     work_r=cgnr_r, work_z=cgnr_z, work_p=cgnr_p,
                                     work_rd=cgnr_rd, work_qb=cgnr_qb)[2]
                    dy_try = copy(cgnr_d)
                    @inbounds for i in 1:N; colnrm[i] < 1e-12 && (dy_try[i] = 0.0); end
                catch e
                    dy_try = zeros(N); cgnr_rel = Inf
                    log("  it $it lmtry=$lmtry: CGNR failed — $(sprint(showerror, e))")
                end
            end
            if any(!isfinite, dy_try); dy_try = zeros(N); exact = false; end

            # Convert direction back to original-variable units.
            # ruiz_scales returns Dr/Dc as DIVISORS: J_s = Dr⁻¹·Jbase·Dc⁻¹ (see the
            # construction at line ~412 and rhs_s = gs./Dr_ruiz above). CGNR solves
            #     J_s·dy = -gs./Dr  ⟺  Dr⁻¹·Jbase·Dc⁻¹·dy = -Dr⁻¹·gs
            #                       ⟺  Jbase·(dy./Dc) = -gs,
            # so the step in original units is dy ./ Dc, NOT dy .* Dc. Multiplying
            # here corrupted every step by a factor of Dc² per variable, producing a
            # direction that was not a descent direction at all: the line search then
            # saw ‖F(x+αd)‖ ≈ α‖J·d‖ grow linearly in α and rejected every α. It was
            # invisible at the benchmark only because a start point that already
            # satisfies F(x)=0 returns above without ever taking a step.
            d_try = dy_try ./ Dc_ruiz

            # ── Cap and freeze ─────────────────────────────────────────
            # Only for an INEXACT direction. The caps below are heuristics whose
            # job is to keep a poorly-solved direction from doing damage; applied
            # to a verified Newton step they turn it into a non-Newton one and the
            # line search then sees ‖F(x+αd)‖ grow linearly in α. Measured: the
            # additive floor of 10.0 was clipping delVGDPEXP/delXGDPEXP/delGDPINC
            # steps of 300–1221 (these are value-denominated deltas in Rp bn whose
            # benchmark is 0, so `max(0.10|x|, 10.0)` bears no relation to their
            # scale) — 157 variables were being clipped, the worst by 122×.
            # Domain safety is handled by fraction-to-boundary below, which is
            # retained unconditionally.
            if !exact
                @inbounds for i in 1:N
                    absx = abs(x[freecols[i]])
                    if isfinite(lb[i]) && lb[i] >= 1e-6 - 1e-9
                        # ratio/index variable: 10% relative cap
                        maxstep = rel_cap * max(absx, 1e-3)
                    else
                        # additive flow/shifter: allow O(1) absolute moves (e.g. log-shifters
                        # that start at 0 and may need to move several units); still relative
                        # for large genuine flows.
                        maxstep = max(rel_cap * absx, 10.0)
                    end
                    d_try[i] > maxstep && (d_try[i] = maxstep)
                    d_try[i] < -maxstep && (d_try[i] = -maxstep)
                end
            end
            # Same reasoning as the magnitude cap: zeroing components of a
            # VERIFIED Newton direction destroys the very property the line
            # search relies on. (Measured inert at 25x6 — Ruiz leaves every
            # column with infinity-norm 1 — but it is a silent hazard.)
            if !exact
                @inbounds for i in 1:N
                    if colnrm[i] < 1e-3; d_try[i] = 0.0; end
                end
            end
            a = 1.0
            @inbounds for i in 1:N
                if isfinite(lb[i]) && d_try[i] < -1e-10
                    gap = x[freecols[i]] - lb[i]
                    if gap > 1e-8
                        cap = 0.99 * gap / (-d_try[i])
                        cap > 0 && cap < a && (a = cap)
                    else; d_try[i] = 0.0; end
                end
            end
            n_frozen = count(colnrm .< 1e-3)
            if lmtry == 1 || verbose
                log("  it $it lmtry=$lmtry mu=$(round(mu;sigdigits=3)): $(exact ? "LU" : "CGNR") $(round(time()-t_solve;digits=1))s  max|d|=$(round(maximum(abs, d_try);sigdigits=4))  a0=$(round(a;sigdigits=3))  frozen=$n_frozen  lin_rel=$(round(cgnr_rel;sigdigits=3))")
            end

            # ── Line search (on the merit, ½‖F‖₂²) ─────────────────────
            xt = copy(x); tbuf = zeros(NC)
            best_a = 0.0; best_m = Inf; best_g = f0_inf; best_gt = gs
            for ls in 1:20
                t_eval = time()
                @inbounds for i in 1:N; xt[freecols[i]] = x[freecols[i]] + a * d_try[i]; end
                gt = sres(xt, tbuf, rowscale)
                mt = merit(gt)
                if isfinite(mt) && mt < best_m
                    best_m = mt; best_a = a; best_g = norm(gt, Inf)
                    best_gt = copy(gt); copyto!(d, d_try)
                end
                if verbose && it <= 2
                    println("      ls $ls a=$(round(a;sigdigits=3)) merit=$(round(mt;sigdigits=4)) ginf=$(round(norm(gt,Inf);sigdigits=4)) teval=$(round(time()-t_eval;digits=3))s")
                    flush(stdout)
                end
                a *= 0.5
                a < 1e-6 && break
            end

            # ── Accept / adapt mu ──────────────────────────────────────
            # Non-monotone (Grippo) reference: the WORST merit in the recent
            # window, not the immediately preceding one. Newton's answer to a
            # curvature overshoot is to correct it at the NEXT iteration, which
            # a strictly monotone test never permits — it just shrinks α.
            nm_ref = maximum(mem)
            if best_a > 0 && isfinite(best_m) && best_m < nm_ref
                @inbounds for i in 1:N
                    x[freecols[i]] = x[freecols[i]] + best_a * d[i]
                end
                gs = best_gt
                accepted = true; alpha = best_a
                push!(mem, best_m)
                length(mem) > NM_WINDOW && popfirst!(mem)
                # Regularization is adapted on the QUALITY of the accepted step,
                # not merely on acceptance. Accepting α≈1e-6 for a 1e-4 relative
                # merit change is not evidence the direction is trustworthy —
                # treating it as such drove mu down to 1e-8 (i.e. no damping at
                # all) exactly as the Jacobian was becoming near-singular.
                if best_a > 0.5 && best_m < 0.9 * f0_m
                    mu = max(1e-8, mu * 0.1)
                elseif best_a < 0.1
                    mu = min(1e6, mu * 10.0)
                end
                # Trust radius keyed on the merit reduction achieved, NOT on alpha.
                # An earlier version grew Delta only when alpha >= 0.99; since the
                # damped path essentially never returns alpha = 1, Delta ratcheted
                # monotonically down to ~0.01 against an LU step of ~1650 and the
                # exact path was locked out for good (gate2: 28 consecutive
                # demotions, merit flat at 2.5e-9). A trust region must be able to
                # EXPAND after a good step or it is just a decaying step cap.
                moved = best_a * maximum(abs, d)
                red = (f0_m - best_m) / max(f0_m, eps())
                Delta = if red > 0.1
                    max(2 * moved, 2 * Delta)      # good step: allow the region to grow back
                elseif red > 0.01
                    max(moved, Delta)              # acceptable: hold
                else
                    max(moved, 1e-3)               # poor: contract to what we actually used
                end
                if best_g < inf_best
                    inf_best = best_g; copyto!(x_best, x); gs_best = copy(gs)
                end
                break
            else
                # No α on this ray decreases the merit even against the
                # non-monotone reference. Make the next solve more gradient-like
                # and retry (lmtry ≥ 2 always takes the damped CGNR path).
                #
                # This applies to a VERIFIED exact direction too. An earlier
                # version stopped here instead, on the argument that CGNR cannot
                # beat an exact solve. That argument was measured at the
                # BENCHMARK start point, where J is well conditioned and the LU
                # step is right; it does not transfer to a late iterate. The
                # 0.9997 gate log shows the exact step degenerating —
                # max|d| 1221 → 916800 with lin_rel 3.2e-10 → 1.3e-6, the
                # signature of a near-singular J — while iteration 7, the one
                # iteration that fell through to damped CGNR (mu = 1e-7),
                # produced max|d| = 1.171, α = 0.55 and the largest merit drop
                # of the whole run (4.99e-7 → 1.04e-7). Near a singular point a
                # damped step beats the exact one precisely because it is damped.
                if exact
                    log("  it $it: exact Newton ray gives no merit decrease " *
                        "(best=$(round(best_m;sigdigits=4)) vs ref=$(round(nm_ref;sigdigits=4)))" *
                        " — falling back to damped CGNR")
                end
                mu = min(1e6, mu * 10.0)
                if lmtry < 5
                    log("  it $it lmtry=$lmtry: no descent, increasing mu → $(round(mu;sigdigits=3))")
                end
            end
        end

        log("  it $it: scaled ||F||_inf = $(norm(gs, Inf))  merit = $(round(merit(gs);sigdigits=4))  alpha = $alpha")

        if norm(gs, Inf) < tol
            status = :converged; break
        elseif !accepted
            status = norm(gs, Inf) < 10tol ? :converged : :no_progress
            break
        end
    end

    # A non-monotone search may end on an uphill excursion. Return the best
    # iterate actually seen, never the last one visited.
    if status != :converged && inf_best < norm(gs, Inf)
        log("  restoring best iterate: ||F||_inf $(norm(gs, Inf)) → $inf_best")
        copyto!(x, x_best); gs = gs_best
        norm(gs, Inf) < tol && (status = :converged)
    end

    residual!(gbuf, x)
    values = Dict{String,Any}()
    for (nm, v) in vars
        if v isa AbstractArray
            out = Array{Float64}(undef, size(v))
            for idx in eachindex(v); out[idx] = x[pos[JuMP.index(v[idx]).value]]; end
            values[nm] = out
        else
            values[nm] = x[pos[JuMP.index(v).value]]
        end
    end
    # Warm-start any subsequent solve from this solution.
    for (nm, v) in vars, (k, vr) in enumerate(v isa AbstractArray ? vec(v) : [v])
        JuMP.is_fixed(vr) || JuMP.set_start_value(vr, (values[nm] isa AbstractArray ?
                                                       vec(values[nm])[k] : values[nm]))
    end

    return NewtonResult(status, status == :converged, iters,
                        norm(gs, Inf), norm(gbuf, Inf), values, N, NC)
end

"""
    diagnose_jacobian(m, vars; threshold=1e-10)

Evaluate the Jacobian at the current point and print statistics about weak
columns (columns whose 2-norm < `threshold`). Groups them by variable-family
prefix. Returns `(total_free, n_weak, prefix_counts)`.
"""
function diagnose_jacobian(m::JuMP.Model, vars::Dict{String,Any};
                           threshold::Float64=1e-10)

    # Identify free variables
    freeid = Dict{Int64,Int}(); freeref = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value; haskey(freeid, iv) && continue
        push!(freeref, vr); freeid[iv] = length(freeref)
    end
    N = length(freeref)

    # Find dead rows and orphans for a square system
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
    freecols = [pos[JuMP.index(v).value] for v in freeref]
    seen = falses(N)
    deadcons = Any[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            buf = _collect_free!(Int[], JuMP.constraint_object(con).func, freeid)
            if isempty(buf)
                push!(deadcons, con)
            else
                for id in buf; seen[id] = true; end
            end
        end
    end
    orphans = [freeref[i] for i in eachindex(freeref) if !seen[i]]
    for vr in orphans
        sv = JuMP.start_value(vr)
        JuMP.fix(vr, sv === nothing ? 0.0 : sv; force=true)
    end
    for con in deadcons; JuMP.delete(m, con); end

    # Rebuild free-variable list after squaring
    nowfree = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) || push!(nowfree, vr)
    end
    N = length(nowfree)

    # Build evaluator
    nlmodel = MOI.Nonlinear.Model()
    rhs = Float64[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            r = _rhs_of(co.set)
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
            push!(rhs, r)
        end
    end
    NC = length(rhs)
    evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(evaluator, [:Grad, :Jac])

    # Evaluate Jacobian
    x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
    st = MOI.jacobian_structure(evaluator)
    jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
    Jval_raw = zeros(length(jrows_all))
    col2free = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols); col2free[c] = i; end
    keep = [col2free[c] != 0 for c in jcols_all]
    jc = [col2free[c] for c in jcols_all[keep]]
    MOI.eval_constraint_jacobian(evaluator, Jval_raw, x)
    vals = Jval_raw[keep]

    # Column norms in O(nnz)
    colnrm = zeros(N)
    @inbounds for k in eachindex(vals); colnrm[jc[k]] += vals[k]^2; end
    colnrm .= sqrt.(colnrm)
    weak = findall(colnrm .< threshold)

    println("  [diagnose] free vars: $N   weak cols (<$threshold): $(length(weak))")
    pc = Dict{String,Int}()
    for wi in weak
        pre = split(JuMP.name(nowfree[wi]), "[")[1]
        pc[pre] = get(pc, pre, 0) + 1
    end
    for (pre, cnt) in sort(collect(pc), by=x->-x[2])
        println("    $(rpad(pre,18)) $cnt")
    end
    if !isempty(weak)
        println("  Sample:")
        for wi in weak[1:min(5, length(weak))]
            println("    $(JuMP.name(nowfree[wi]))  colnrm=$(round(colnrm[wi], sigdigits=2))")
        end
    end
    return N, length(weak), pc
end

"""
    euler_continuation!(m, vars, shock_vr, target_final, nsteps; tol=1e-6)

Solve a single shock using the Euler continuation method:
1. `nsteps` Euler steps (1 Newton iteration each), reusing the evaluator.
2. Final Newton cleanup (up to 25 iterations).

The evaluator is built ONCE, avoiding the ~8s rebuild overhead per step.
Returns a [`NewtonResult`](@ref) with the final solution.
"""
function euler_continuation!(m::JuMP.Model, vars::Dict{String,Any},
                             shock_vr::JuMP.VariableRef, target_final::Float64,
                             nsteps::Int; tol::Float64=1e-6, verbose::Bool=true)
    log(msg) = verbose && println(msg)

    freeid = Dict{Int64,Int}(); freeref = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value; haskey(freeid, iv) && continue
        push!(freeref, vr); freeid[iv] = length(freeref)
    end
    seen = falses(length(freeref))
    deadcons = Any[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            buf = _collect_free!(Int[], JuMP.constraint_object(con).func, freeid)
            if isempty(buf); push!(deadcons, con)
            else; for id in buf; seen[id] = true; end; end
        end
    end
    orphans = [freeref[i] for i in eachindex(freeref) if !seen[i]]
    for vr in orphans
        sv = JuMP.start_value(vr)
        JuMP.fix(vr, sv === nothing ? 0.0 : sv; force=true)
    end
    for con in deadcons; JuMP.delete(m, con); end
    if !isempty(orphans) || !isempty(deadcons)
        log("  squaring: pinned $(length(orphans)) orphan(s), deleted $(length(deadcons)) dead row(s)")
    end

    nowfree = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) || push!(nowfree, vr)
    end
    N = length(nowfree)
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
    freecols = [pos[JuMP.index(v).value] for v in nowfree]
    lb = [JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf for v in nowfree]

    nlmodel = MOI.Nonlinear.Model()
    rhs = Float64[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            r = _rhs_of(co.set)
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
            push!(rhs, r)
        end
    end
    NC = length(rhs)
    evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(evaluator, [:Grad, :Jac])
    log("  system: $NC equations, $N unknowns")

    st = MOI.jacobian_structure(evaluator)
    jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
    Jval = zeros(length(jrows_all))
    col2free = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols); col2free[c] = i; end
    keep = [col2free[c] != 0 for c in jcols_all]
    jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

    # x vector (includes fixed vars — updated each iteration)
    x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
    shock_idx = pos[JuMP.index(shock_vr).value]

    # Row scaling (computed at the initial point)
    MOI.eval_constraint_jacobian(evaluator, Jval, x)
    rowscale = ones(Float64, NC)
    let v0 = Jval[keep], mx = zeros(NC)
        @inbounds for k in eachindex(jr)
            a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
        end
        @inbounds for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
    end
    gbuf = zeros(NC)
    sres(xv, buf) = (residual!(buf, xv, evaluator, rhs); buf ./ rowscale)

    gs = sres(x, gbuf)
    log("  initial ||F|| = $(norm(gs, Inf))")

    # ═══ Euler steps (Jacobian re-evaluated each step, with inner Newton corrections) ═══
    sh0 = JuMP.fix_value(shock_vr)
    inner_tol = 1e-7
    iters_total = nsteps
    for s in 1:nsteps
        cumul = sh0 + (target_final - sh0) * s / nsteps
        JuMP.fix(shock_vr, cumul; force=true)
        x[shock_idx] = cumul
        gs = sres(x, gbuf)

        # Inner Newton corrections
        for inner in 1:10
            MOI.eval_constraint_jacobian(evaluator, Jval, x)
            vals = Jval[keep]
            J = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
            colnrm = [norm(J[:, j]) for j in 1:N]
            colnrm_s = max.(colnrm, 1e-3)
            Jsc = J * spdiagm(1.0 ./ colnrm_s)
            f0_inf = norm(gs, Inf)
            xt = copy(x); tbuf2 = zeros(NC); accepted = false
            for μ in (0.0, 1e-4, 1e-2)
                Jeq = μ > 0 ? (Jsc + μ * I) : Jsc
                local d_try
                try
                    d_eq = try -(lu(Jeq) \ gs) catch; -(qr(Jeq) \ gs); end
                    d_try = d_eq ./ colnrm_s
                    @inbounds for i in 1:N
                        colnrm[i] < 1e-3 && (d_try[i] = 0.0)
                    end
                catch
                    continue
                end
                any(!isfinite, d_try) && continue
                @inbounds for i in 1:N
                    maxstep = 0.10 * max(abs(x[freecols[i]]), 1e-3)
                    d_try[i] > maxstep && (d_try[i] = maxstep)
                    d_try[i] < -maxstep && (d_try[i] = -maxstep)
                end
                a = 1.0
                for i in 1:N
                    if isfinite(lb[i]) && d_try[i] < -1e-10
                        gap = x[freecols[i]] - lb[i]
                        if gap > 1e-8
                            cap = 0.99 * gap / (-d_try[i])
                            cap > 0 && cap < a && (a = cap)
                        else
                            d_try[i] = 0.0
                        end
                    end
                end
                best_a = 0.0; best_g = f0_inf; best_gt = gs; best_d = d_try
                for _ in 1:20
                    @inbounds for i in 1:N; xt[freecols[i]] = x[freecols[i]] + a * d_try[i]; end
                    gt = sres(xt, tbuf2)
                    ginf = norm(gt, Inf)
                    if isfinite(ginf) && ginf < best_g
                        best_g = ginf; best_a = a; best_gt = copy(gt); best_d = copy(d_try)
                    end
                    a *= 0.5
                    a < 1e-6 && break
                end
                if best_a > 0 && best_g < f0_inf * 0.999
                    @inbounds for i in 1:N
                        x[freecols[i]] += best_a * best_d[i]
                    end
                    gs = best_gt; accepted = true; break
                end
            end
            accepted || break
            iters_total += 1
            norm(gs, Inf) < inner_tol && break
        end

        if s == 1 || s == nsteps || s % max(1, div(nsteps, 10)) == 0
            log("  Euler step $s/$nsteps: ||F|| = $(norm(gs, Inf))")
        end
    end

    # ═══ Final Newton cleanup ═══
    log("  Newton cleanup...")
    for it in 1:15
        iters_total += 1
        MOI.eval_constraint_jacobian(evaluator, Jval, x)
        vals = Jval[keep]
        J = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
        colnrm = [norm(J[:, j]) for j in 1:N]
        colnrm_s = max.(colnrm, 1e-3)
        Jsc = J * spdiagm(1.0 ./ colnrm_s)
        f0_inf = norm(gs, Inf)
        xt = copy(x); tbuf3 = zeros(NC); accepted = false; alpha = 0.0
        for μ in (0.0, 1e-4)
            Jeq = μ > 0 ? (Jsc + μ * I) : Jsc
            local d_try
            try
                d_eq = try -(lu(Jeq) \ gs) catch; -(qr(Jeq) \ gs); end
                d_try = d_eq ./ colnrm_s
                @inbounds for i in 1:N; colnrm[i] < 1e-3 && (d_try[i] = 0.0); end
            catch
                continue
            end
            any(!isfinite, d_try) && continue
            @inbounds for i in 1:N
                maxstep = 0.25 * abs(x[freecols[i]]) + 1e-3
                d_try[i] > maxstep && (d_try[i] = maxstep)
                d_try[i] < -maxstep && (d_try[i] = -maxstep)
            end
            a = 1.0
            for i in 1:N
                if isfinite(lb[i]) && d_try[i] < -1e-10
                    gap = x[freecols[i]] - lb[i]
                    if gap > 1e-8
                        cap = 0.99 * gap / (-d_try[i])
                        cap > 0 && cap < a && (a = cap)
                    else
                        d_try[i] = 0.0
                    end
                end
            end
            for _ in 1:30
                @inbounds for i in 1:N; xt[freecols[i]] = x[freecols[i]] + a * d_try[i]; end
                gt = sres(xt, tbuf3)
                ginf = norm(gt, Inf)
                if isfinite(ginf) && ginf < f0_inf * (1 - 1e-8)
                    copyto!(x, xt); gs = copy(gt); accepted = true; alpha = a; break
                end
                a *= 0.5
                a < 1e-14 && break
            end
            accepted && break
        end
        log("  cleanup it $it: ||F||_inf=$(norm(gs, Inf)) alpha=$alpha")
        if norm(gs, Inf) < tol
            log("  Converged after $(it) cleanup iterations"); break
        elseif !accepted && norm(gs, Inf) < 10tol
            log("  Near-converged, accepting"); break
        elseif !accepted
            log("  No progress, stopping cleanup"); break
        end
    end

    # Extract values
    gbuf2 = zeros(NC); residual!(gbuf2, x, evaluator, rhs)
    values = Dict{String,Any}()
    for (nm, v) in vars
        if v isa AbstractArray
            out = Array{Float64}(undef, size(v))
            for idx in eachindex(v); out[idx] = x[pos[JuMP.index(v[idx]).value]]; end
            values[nm] = out
        else
            values[nm] = x[pos[JuMP.index(v).value]]
        end
    end
    for (nm, v) in vars, (k, vr) in enumerate(v isa AbstractArray ? vec(v) : [v])
        JuMP.is_fixed(vr) || JuMP.set_start_value(vr, (values[nm] isa AbstractArray ?
                                                       vec(values[nm])[k] : values[nm]))
    end

    status = norm(gs, Inf) < tol ? :converged : :maxit
    return NewtonResult(status, status == :converged, iters_total,
                        norm(gs, Inf), norm(gbuf2, Inf), values, N, NC)
end

"""
    residual!(g, xv, evaluator, rhs)

Evaluate the constraint residuals `F(xv)` into `g` using the given evaluator.
"""
function residual!(g::Vector{Float64}, xv::Vector{Float64},
                   evaluator::MOI.Nonlinear.Evaluator, rhs::Vector{Float64})
    MOI.eval_constraint(evaluator, g, xv)
    @inbounds for i in eachindex(g); g[i] -= rhs[i]; end
    return g
end
