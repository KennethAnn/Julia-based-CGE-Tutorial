#  A Simple CGE Model with two sectors and no intermediate input
using JuMP, Complementarity, DataFrames

sec = ["sec1", "sec2"]
sc = [1, 2]
sam = [
    missing missing missing missing 12      12
    missing missing missing missing 21      21
    9       7       missing missing missing 16
    3       14      missing missing missing 17
    missing missing 16      17      missing 33
    12      21      16      17      33      missing
]
samList = ["sec1", "sec2", "lab", "cap", "hh", "total"]

Q0 = sam[6, 1:2]
P0 = [1, 1]
LD0 = sam[3, 1:2]
KD0 = sam[4, 1:2]
LS = sum(LD0)[1]
KS = sum(KD0)[1]
WL0 = 1
WK0 = 1
QH0 = sam[1:2, 5]
Y0 = WL0 * LS + WK0 * KS
rho = 0.75

# 校准生产函数参数
al = LD0 ./ (LD0 + KD0)
ak = KD0 ./ (LD0 + KD0)
alphah = P0 .* QH0 / Y0

# 3. Generate CGE Model
function solve_cge()
    m = MCPModel()
    @variables m begin
        P[i = sc], (start = P0[i])
        WK, (start = WK0) 
        WL, (start = WL0)
        Q[i = sc], (start = Q0[i])
        LD[i = sc], (start = LD0[i])
        KD[i = sc], (start = KD0[i])
        Y, (start = Y0)
        QH[i = sc], (start = QH0[i])
    end

    fix(WK, 1.0)
    # production function - CES
    @mapping(m, eqL[i in sc], LD[i] * WL ^ rho - al[i] * Q[i] * P[i] ^ rho)
    @complementarity(m, eqL, LD)

    @mapping(m, eqK[i in sc], KD[i] * WK ^ rho - ak[i] * Q[i] * P[i] ^ rho)
    @complementarity(m, eqK, KD)

    @mapping(m, eqP[i in sc], al[i] * WL ^ (1 - rho) + ak[i] * WK ^ (1 - rho) - P[i] ^ (1 - rho))
    @complementarity(m, eqP, P)

    # income distribution
    @mapping(m, eqY, Y - (WL * LS + WK * KS))
    @complementarity(m, eqY, Y)

    # consumption structure
    @mapping(m, eqQH[i in sc], alphah[i] * Y - P[i] * QH[i])
    @complementarity(m, eqQH, QH)

    # market clearance: goods, production factors
    @mapping(m, eqQ[i in sc], Q[i] - QH[i])
    @complementarity(m, eqQ, Q)

    @mapping(m, eqWK, KS - sum(KD[i] for i in sc))
    @complementarity(m, eqWK, WK)

    @mapping(m, eqWL, LS - sum(LD[i] for i in sc))
    @complementarity(m, eqWL, WL)

    # Model Solver
    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=3000)
    @show result_value.(WL)
end

solve_cge()

# second experiment by adjusting the supply of labor and capital
LS = sum(LD0) * 2
solve_cge()