function step1f_trade(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        tem = 0.0
        @inbounds for r in 1:NR; tem += TRADE[c,s,r,d]; end
        sc = tem > 0 ? BASIC_U[c,s,d] / tem : 1.0
        for r in 1:NR
            TRADE[c,s,r,d] *= sc
            for m in MAR_idx_full
                TRADMAR[c,s,m,r,d] *= sc
            end
        end
    end
end

NC=185; NS=2; NM=9; NR=34
TRADE = rand(NC, NS, NR, NR)
TRADMAR = rand(NC, NS, NM, NR, NR)
BASIC_U = rand(NC, NS, NR)
MAR_idx_full = collect(1:NM)

# warmup
step1f_trade(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)

GC.gc()
t0 = time()
step1f_trade(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
println("step1f time: $(time()-t0)s")

# with @inbounds update
function step1f_trade_inbounds(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        tem = 0.0
        @inbounds for r in 1:NR; tem += TRADE[c,s,r,d]; end
        sc = tem > 0 ? BASIC_U[c,s,d] / tem : 1.0
        @inbounds for r in 1:NR
            TRADE[c,s,r,d] *= sc
            for m in MAR_idx_full
                TRADMAR[c,s,m,r,d] *= sc
            end
        end
    end
end

step1f_trade_inbounds(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
GC.gc()
t0 = time()
step1f_trade_inbounds(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
println("step1f inbounds time: $(time()-t0)s")

# with @inbounds update and type stable
function step1f_typed(TRADE::Array{Float64,4}, TRADMAR::Array{Float64,5}, BASIC_U::Array{Float64,3}, MAR_idx_full::Vector{Int}, NC::Int, NS::Int, NR::Int)
    for c in 1:NC, s in 1:NS, d in 1:NR
        tem = 0.0
        @inbounds for r in 1:NR; tem += TRADE[c,s,r,d]; end
        sc = tem > 0 ? BASIC_U[c,s,d] / tem : 1.0
        @inbounds for r in 1:NR
            TRADE[c,s,r,d] *= sc
            for m in MAR_idx_full
                TRADMAR[c,s,m,r,d] *= sc
            end
        end
    end
end

step1f_typed(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
GC.gc()
t0 = time()
step1f_typed(TRADE, TRADMAR, BASIC_U, MAR_idx_full, NC, NS, NR)
println("step1f typed time: $(time()-t0)s")
