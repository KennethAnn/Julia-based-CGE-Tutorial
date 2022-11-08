#  A Simple CGE Model with government
using JuMP, Complementarity, DataFrames

sec = ["agri", "manu", "serv"]
sc = [1, 2, 3]
sam = [
    260     320     150     missing missing 630     10
    345     390     390     missing missing 590     15
    400     365     320     missing missing 385     5
    200     250     400     missing missing missing missing
    160     400     210     missing missing missing missing
    missing missing missing 850     770     missing 10
    5       5       5       missing missing 25      missing                             
]
samList = ["agri", "manu", "serv", "lab", "cap", "hh"]

# production block
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

q0 = tint0 + va0 + sam[7, sc]
atint = tint0 ./ q0 
ava = va0 ./ q0

# income block
transfr = sam[6,7]

ks = sum(k0)
ls = sum(l0)
yh0 = ks + ls + transfr
th = sam[7,6] / yh0

yd0 = yh0 * (1 - th)

ts = sam[7, sc] ./ q0
yg0 = th * yh0 + sum(sam[7,sc])

qg0 = sam[sc, 7]
shrg = qg0 ./ (yg0-transfr)
rho = 0.6

# LES system
qh0 = sam[sc, 6]
alphah = qh0 ./ yd0
# LESelas = [0.7, 1.0, 1.5]
# Frisch = -3
# LESbeta = (LESelas .* alphah) ./ sum(LESelas .* alphah)
# LESbetachk = sum(LESbeta)
# LESsub = qh0 + LESbeta * (yd0 / Frisch)


# 3. Generate CGE Model
function solve_cge()
    
    m = MCPModel()
    wl = 1
    @variables m begin
        p[i = sc], (start = 1)
        wk, (start = 1) 
        pint[i = sc], (start = 1)
        pva[i = sc], (start = 1)
        q[i = sc], (start = q0[i])
        tint[i = sc], (start = tint0[i])
        va[i = sc], (start = va0[i])
        k[i = sc], (start = k0[i])
        l[i = sc], (start = k0[i])
        qint[j = sc, i = sc], (start = qint0[j, i])
        qh[i = sc], (start=qh0[i])
        yh, (start=yh0)
        yd, (start=yd0)
        yg, (start=yg0)
        eg, (start=yg0)
        qg[i = sc], (start=qg0[i])
        walras, (start = 0)
        # gdp, (start=y0)
        # pgdp, (start=1)
    end
    # tint + va = q

    @mapping(m, eqtint[i in sc], tint[i] - atint[i] * q[i])
    @complementarity(m, eqtint, tint)

    @mapping(m, eqva[i in sc], va[i] - ava[i] * q[i])
    @complementarity(m, eqva, va)

    @mapping(m, eqp[i in sc], (atint[i] * pint[i] + ava[i] * pva[i]) / (1 - ts[i]) - p[i])
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
    @mapping(m, eqyh, (wl * ls + wk * ks + transfr) - yh)
    @complementarity(m, eqyh, yh)

    @mapping(m, eqyd, yh * (1 - th) - yd)
    @complementarity(m, eqyd, yd)

    @mapping(m, eqyg, yh * th + sum(q[i] *p[i] * ts[i] for i in sc) - yg)
    @complementarity(m, eqyg, yg)

    @mapping(m, eqeg, yg - eg)
    @complementarity(m, eqeg, eg)

    @mapping(m, eqqg[i in sc], qg[i] * p[i] - shrg[i] * (eg - transfr))
    @complementarity(m, eqqg, qg)

    # demand
    @mapping(m, eqqh[i in sc], p[i] * LESsub[i] + LESbeta[i] * (yd - sum(p[i] * LESsub[i] for i in sc))- p[i] * qh[i])
    # @mapping(m, eqqh[i in sc], alphah[i] * yd - p[i] * qh[i])
    @complementarity(m, eqqh, qh)

    # demand - supply
    @mapping(m, eqq[i in sc], q[i] - (qh[i] + qg[i] + sum(qint[i, j] for j in sc)))
    @complementarity(m, eqq, q)

    # factor supply
    @mapping(m, eqwk, ks - sum(k[i] for i in sc))
    @complementarity(m, eqwk, wk)

    @mapping(m, eqwl, ls + walras - sum(l[i] for i in sc))
    @complementarity(m, eqwl, walras)

    # @mapping(m, eqgdp, gdp - sum(qh[i] for i in sc))
    # @complementarity(m, eqgdp, gdp)

    # @mapping(m, eqpdgp, pgdp * gdp - sum(p[i] * qh[i] for i in sc))
    # @complementarity(m, eqpdgp, pgdp)

    # Model Solver
    print(m)
    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", ITERATION_LIMIT=10000)
    # @show result_value.(pgdp)
    # @show result_value.(gdp)
    # println("Walras is:")
    @show result_value.(walras)
    # @show result_value.(wk)
end

solve_cge()


ls = sum(l0) * 2
ks = sum(k0) * 3
solve_cge()
