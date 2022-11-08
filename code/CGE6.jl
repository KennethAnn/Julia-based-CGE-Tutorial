#  A Simple CGE Model with one household and one government
using JuMP, Complementarity, DataFrames, CSV
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

samList = ["agri", "manu", "serv", "lab", "cap", "tax", "hh", "govt"]
sector = ["agri", "manu", "serv"]
sectors = collect(1:1:length(sector))
household = ["hh"]
households = collect(1:1:length(household))
va = ["lab","cap","tax"]
NumSector = length(sectors)
NumHouseholds = length(households)

sam = [
    260     320     150     missing missing missing 630     10
    345     390     390     missing missing missing 590     15
    400     365     320     missing missing missing 385     5
    200     250     400     missing missing missing missing missing
    160     400     210     missing missing missing missing missing
    5       5       5       missing missing missing missing missing
    missing missing missing 850     770     missing missing 10
    missing missing missing missing missing 15      25      missing                         
]
# production
qint0 = sam[sectors, sectors]
l0 = sam[length(sectors) + 1, sectors]
k0 = sam[length(sectors) + 2, sectors]
ls = sum(l0)
ks = sum(k0)
kl0 = k0 + l0
ak = k0 ./ kl0
al = l0 ./ kl0

int0 = sum(qint0, dims=1)[1, :]
aqint = qint0 ./ (sum(qint0, dims=1))

prod0 = int0 + kl0 + sam[NumSector + 3, sectors]
pprod0 = 1 .- sam[NumSector + 3,sectors] ./ prod0

aint = int0 ./ prod0 
akl = kl0 ./ prod0

# income block
transfr = sam[7, 8]
inch0 = ks + ls + transfr
taxh = sam[8, 7] / inch0

dinch0 = inch0 * (1 - taxh)

prodtr = sam[NumSector + 3, sectors]./ prod0
# pprod0 = 1 ./ (1 .+ prodtr)

incg0 = taxh * inch0 + sum(sam[6,sectors])
consg0 = sam[sectors, 8]
thetag = consg0 ./ (incg0-transfr)

rho = 0.6

# LES system
consh0 = sam[sectors, 7]
thetah = consh0 ./ dinch0
# LESelas = [0.7, 1.0, 1.5]
# Frisectorsh = -3
# LESbeta = (LESelas .* alphah) ./ sum(LESelas .* alphah)
# LESbetachk = sum(LESbeta)
# LESsub = qh0 + LESbeta * (yd0 / Frisectorsh)


elasdir = joinpath(@__DIR__, "..", "data", "Elasticities.csv")
elas = CSV.read(elasdir, DataFrames.DataFrame, header=1)

# 3. Calibrate the model
# 3.0 Elasticities
sigmap = elas[sectors, "sigmap"]
sigmakel = elas[sectors, "sigmakel"]
# sigmakel = [0.6, 0.6, 0.6]
# sigmakl = elas[sectors, "sigmakl"]
sigmakl = [0.6, 0.6, 0.6]

sigmaene = elas[sectors, "sigmaene"]
sigmafe = elas[sectors, "sigmafe"]
sigmav = elas[54, "sigmav"]
sigmam = elas[sectors, "sigmam"]
sigmae = elas[sectors, "sigmae"]
sigmaff = elas[sectors, "sigmaff"]
sigmafes = elas[sectors, "sigmafes"]
sigmanr = elas[sectors, "sigmanr"]
eta = elas[sectors, "eta"]

# 4. Generate CGE Model
function solve_cge()  
    m = MCPModel()
    pl = 1 # exogenous setting of labor price as one
    @variables m begin
        pprod[i = sectors], (start = pprod0[i]) # production cost of each sector
        pprodt[i = sectors], (start = 1)  # production cost by adding tax
        pk, (start = 1) # price of capital
        pint[i = sectors], (start = 1) # price of aggregated intermediate goods
        pkl[i = sectors], (start = 1) # price of capital-labor bundle
        prod[i = sectors], (start = prod0[i]) # production lelve of each sector
        int[i = sectors], (start = int0[i]) # amount of aggregated intermediate goods
        kl[i = sectors], (start = kl0[i]) # capital-labor bundle
        k[i = sectors], (start = k0[i]) # capital bundle
        l[i = sectors], (start = l0[i]) # labor bundle
        qint[j = sectors, i = sectors], (start = qint0[j, i]) # segmented intermediate goods
        consh[i = sectors], (start=consh0[i]) # consumption of households
        inch, (start=inch0) # income of households 
        dinch, (start=dinch0) # disposable income of households
        incg, (start=incg0) # income of governments
        expg, (start=incg0) # expenditure of governments
        consg[i = sectors], (start=consg0[i]) # consumption of government
        walras, (start = 0) # walras variable for testing the balance
        # gdp, (start=y0)
        # pgdp, (start=1)
    end
    # int + kl = prod
    @mapping(m, eq_int[i in sectors], int[i] * pint[i] ^ sigmap[i] - aint[i] * prod[i] * pprod[i] ^ sigmap[i])
    @complementarity(m, eq_int, int)

    @mapping(m, eq_kl[i in sectors], kl[i] * pkl[i] ^ sigmap[i] - akl[i] * prod[i] * pprod[i] ^ sigmap[i])
    @complementarity(m, eq_kl, kl)

    @mapping(m, eq_pprod[i in sectors], aint[i] * pint[i] ^ (1- sigmap[i]) + akl[i] * pkl[i] ^ (1 - sigmap[i]) - pprod[i] ^ (1 - sigmap[i]))
    @complementarity(m, eq_pprod, pprod)

    # @mapping(m, eq_int[i in sectors], int[i] - aint[i] * prod[i])
    # @complementarity(m, eq_int, int)

    # @mapping(m, eq_kl[i in sectors], kl[i]- akl[i] * prod[i])
    # @complementarity(m, eq_kl, kl)

    # @mapping(m, eq_pprodt[i in sectors], (aint[i] * pint[i]  + akl[i] * pkl[i]) / (1-prodtr[i]) - pprodt[i])
    # @complementarity(m, eq_pprodt, pprodt)

    # @mapping(m, eq_pprod[i in sectors], pint[i] * int[i] + pkl[i] * kl[i] - pprod[i] * prod[i])
    # @complementarity(m, eq_pprod, pprod)

    @mapping(m, eq_pprodt[i in sectors], pprod[i]/(1 - prodtr[i]) - pprodt[i])
    @complementarity(m, eq_pprodt, pprodt)

    # @mapping(m, eq_int[i in sc], int[i] - aint[i] * prod[i])
    # @complementarity(m, eq_int, int)

    # @mapping(m, eq_kl[i in sc], kl[i] - akl[i] * prod[i])
    # @complementarity(m, eq_kl, kl)

    # @mapping(m, eq_pprodt[i in sc], (aint[i] * pint[i] + akl[i] * pkl[i]) / (1 - prodtr[i]) - pprodt[i])
    # @complementarity(m, eq_pprodt, pprodt)

    # k + l = kl
    @mapping(m, eq_k[i in sectors], k[i] * pk ^ sigmakl[i] - ak[i] * kl[i] * pkl[i] ^ sigmakl[i])
    @complementarity(m, eq_k, k)

    @mapping(m, eq_l[i in sectors], l[i] * pl ^ sigmakl[i] - al[i] * kl[i] * pkl[i] ^ sigmakl[i])
    @complementarity(m, eq_l, l)

    @mapping(m, eq_pkl[i in sectors], ak[i] * pk ^ (1- sigmakl[i]) + al[i] * pl ^ (1 - sigmakl[i]) - pkl[i] ^ (1 - sigmakl[i]))
    @complementarity(m, eq_pkl, pkl)

    # qint ++ = int
    @mapping(m, eq_qint[j in sectors, i in sectors], qint[j, i] - aqint[j, i] * int[i])
    @complementarity(m, eq_qint, qint)

    @mapping(m, eq_pint[i in sectors], sum(aqint[j, i] * pprodt[j] for j in sectors) - pint[i])
    @complementarity(m, eq_pint, pint)

    # income
    @mapping(m, eq_inch, (pl * ls + pk * ks + transfr) - inch)
    @complementarity(m, eq_inch, inch)

    @mapping(m, eq_dinch, inch * (1 - taxh) - dinch)
    @complementarity(m, eq_dinch, dinch)

    @mapping(m, eq_incg, inch * taxh + sum(prod[i] *pprodt[i] * prodtr[i] for i in sectors) - incg)
    @complementarity(m, eq_incg, incg)

    @mapping(m, eq_expg, incg - expg)
    @complementarity(m, eq_expg, expg)

    @mapping(m, eq_consg[i in sectors], consg[i] * pprodt[i] - thetag[i] * (expg - transfr))
    @complementarity(m, eq_consg, consg)

    # demand
    @mapping(m, eq_consh[i in sectors], thetah[i] * dinch - consh[i] * pprodt[i])
    @complementarity(m, eq_consh, consh)

    # demand - supply
    @mapping(m, eq_prod[i in sectors], prod[i] - (consh[i] + consg[i] + sum(qint[i, j] for j in sectors)))
    @complementarity(m, eq_prod, prod)

    # factor supply
    @mapping(m, eq_pk, ks - sum(k[i] for i in sectors))
    @complementarity(m, eq_pk, pk)

    @mapping(m, eq_pl, ls + walras - sum(l[i] for i in sectors))
    @complementarity(m, eq_pl, walras)

    # @mapping(m, eqgdp, gdp - sum(qh[i] for i in sectors))
    # @complementarity(m, eqgdp, gdp)

    # @mapping(m, eqpdgp, pgdp * gdp - sum(p[i] * qh[i] for i in sectors))
    # @complementarity(m, eqpdgp, pgdp)

    # Model Solver
    # print(ComplementarityType[1])
    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=600)
    # @show result_value.(pgdp)
    # @show result_value.(gdp)
    # println("Walras is:")
    @show result_value.(pk)
    @show result_value.(walras)
end

solve_cge()

ls = sum(l0) * 6
ks = sum(k0) * 10
solve_cge()
