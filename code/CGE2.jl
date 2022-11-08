#  A Simple CGE Model with three sectors and intermediate input

using JuMP, Complementarity, DataFrames

sec = ["agri", "manu", "serv"]
sc = [1, 2, 3]
sam = [
    260     320     150     missing missing 635
    345     390     390     missing missing 600
    400     365     320     missing missing 385
    200     250     400     missing missing missing
    160     400     210     missing missing missing
    missing missing missing 850     770     missing 
]
samList = ["agri", "manu", "serv", "lab", "cap", "hh"]

qint0 = sam[sc, sc]
k0 = sam[length(sc) + 2, sc]
l0 = sam[length(sc) + 1, sc]
ks = sum(k0)
ls = sum(l0)
va0 = k0 + l0
ak = k0 ./ va0
al = l0 ./ va0

tint0 = sum(qint0, dims=1)[1, :]
aint = qint0 ./ (sum(qint0, dims=1))

q0 = tint0 + va0
atint = tint0 ./ q0
ava = va0 ./ q0

ks = sum(k0)
ls = sum(l0)
y0 = ks + ls
qh0 = sam[sc, 6]
alphah = qh0 ./ y0

rho = 0.6

# 3. Generate CGE Model
function solve_cge()
    m = MCPModel()
    wl = 1
    @variables m begin
        p[i = sc], (start = 1)
        wk, (start = 1) 
        # wl, (start = 1)
        pint[i = sc], (start = 1)
        pva[i = sc], (start = 1)
        q[i = sc], (start = q0[i])
        tint[i = sc], (start = tint0[i])
        va[i = sc], (start = va0[i])
        k[i = sc], (start = k0[i])
        l[i = sc], (start = k0[i])
        qint[j = sc, i = sc], (start = qint0[j, i])
        qh[i = sc], (start=qh0[i])
        y, (start=y0)
        walras, (start=0)
    end
    # tint + va = q
    @mapping(m, eqtint[i in sc], tint[i] - atint[i] * q[i])
    @complementarity(m, eqtint, tint)

    @mapping(m, eqva[i in sc], va[i] - ava[i] * q[i])
    @complementarity(m, eqva, va)

    @mapping(m, eqp[i in sc], atint[i] * pint[i] + ava[i] * pva[i] - p[i])
    @complementarity(m, eqp, p)

    # k + l = va
    @mapping(m, eqk[i in sc], k[i] * wk ^ rho - ak[i] * va[i] * pva[i] ^ rho)
    @complementarity(m, eqk, k)

    @mapping(m, eql[i in sc], l[i] * wl ^ rho - al[i] * va[i] * pva[i] ^ rho)
    @complementarity(m, eql, l)

    @mapping(m, eqpva[i in sc], al[i] * wl ^ (1 - rho) + ak[i] * wk ^ (1 - rho) - pva[i] ^ (1 - rho))
    @complementarity(m, eqpva, pva)

    # xigma qint = tint
    @mapping(m, eqint[j in sc, i in sc], qint[j, i] - aint[j, i] * tint[i])
    @complementarity(m, eqint, qint)

    @mapping(m, eqpint[i in sc], sum(aint[j, i] * p[j] for j in sc) - pint[i])
    @complementarity(m, eqpint, pint)

    # income
    @mapping(m, eqy, y - (wl * ls + wk * ks))
    @complementarity(m, eqy, y)

    # demand
    @mapping(m, eqqh[i in sc], alphah[i] * y - p[i] * qh[i])
    @complementarity(m, eqqh, qh)

    # demand - supply
    @mapping(m, eqq[i in sc], q[i] - (qh[i] + sum(qint[i, j] for j in sc)))
    @complementarity(m, eqq, q)

    # factor supply
    @mapping(m, eqwk, ks - sum(k[i] for i in sc))
    @complementarity(m, eqwk, wk)

    @mapping(m, eqwl, ls + walras - sum(l[i] for i in sc))
    @complementarity(m, eqwl, walras)

    # Model Solver
    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", ITERATION_LIMIT=10000)
    @show result_value.(walras)
end

solve_cge()

# second experiment by adjusting the supply of labor and capital
ls = sum(l0) * 2
ks = sum(k0) * 0.8
solve_cge()
