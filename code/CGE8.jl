#  A Large CGE Model with 42 sectors, 2 households and one government
using JuMP, Complementarity, DataFrames, CSV
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

samList = ["AgFFF", "Coal", "OilGas", "MetMin", "NMetMin", "FoodPr", "Textile", "Apparel",
            "WoodPr", "PaperPr", "RefPet", "Chemical", "NMetPr", "Metals", "MetalPr",
            "GenEqp", "SpecEqp", "TransEqp", "ElecEqp", "ICTEqp", "PrecInst", "OthMfg",
            "Waste", "MachRep", "ElecDist","GasDist", "WatDist", "Constr", "WhRetTr", "TranspSrv",
            "HotRest", "ICTServ", "Finance", "RealEst", "BusServe", "ResTech", "EnvServ",
            "ResServ", "Education", "Health", "RecEnt", "PubAdm", "lab", "cap", "tax", "hh", "govt"]

sector = ["AgFFF", "Coal", "OilGas", "MetMin", "NMetMin", "FoodPr", "Textile", "Apparel",
            "WoodPr", "PaperPr", "RefPet", "Chemical", "NMetPr", "Metals", "MetalPr",
            "GenEqp", "SpecEqp", "TransEqp", "ElecEqp", "ICTEqp", "PrecInst", "OthMfg",
            "Waste", "MachRep", "ElecDist","GasDist", "WatDist", "Constr", "WhRetTr", "TranspSrv",
            "HotRest", "ICTServ", "Finance", "RealEst", "BusServe", "ResTech", "EnvServ",
            "ResServ", "Education", "Health", "RecEnt", "PubAdm"]

sectors = collect(1:1:length(sector))
household = ["RuralHH", "UrbanHH"]
households = collect(1:1:length(household))
va = ["L","K","T"]
NumSector = length(sectors)
NumHouseholds = length(households)

samdir = joinpath(@__DIR__, "..", "data", "SAM42HH2.csv")
sam = CSV.read(samdir, DataFrames.DataFrame, header=1)
sam = Matrix(sam)[1:48, 2:49] # SAM table

# production
qint0 = sam[sectors, sectors]
l0 = sam[NumSector + 1, sectors]
k0 = sam[NumSector + 2, sectors]
ld0 = sum(l0)
kd0 = sum(k0)
kl0 = k0 + l0
ak = k0 ./ kl0
al = l0 ./ kl0

int0 = sum(qint0, dims=1)[1, :]
aqint = qint0 ./ (sum(qint0, dims=1))

prod0 = int0 + kl0 + sam[NumSector + 3, sectors]
aint = int0 ./ prod0 
akl = kl0 ./ prod0
prodtr = sam[NumSector + 3, sectors]./ prod0
pprod0 =(1 .- prodtr)

# income block
ls = sam[households .+ (NumSector + 3), NumSector + 1]
khs = sam[households .+ (NumSector + 3), NumSector + 2]
ks = kd0
akh = khs ./ kd0

transfr = sam[(NumSector + 3) .+ households, NumSector + 3 + NumHouseholds + 1]
inch0 = khs + ls + transfr
taxh = sam[NumSector + 3 + NumHouseholds + 1, (NumSector + 3) .+ households] ./ inch0
dinch0 = inch0 .* (1 .- taxh)

incg0 = sum(taxh[h] .* inch0[h] for h in households) + sum(sam[NumSector+3, sectors])
consg0 = sam[sectors, NumSector + 3 + NumHouseholds + 1]
thetag = consg0 ./ sum(consg0)

rho = 0.6

# LES system
consh0 = sam[sectors, households .+ (NumSector + 3)]
thetah = consh0 ./ sum(consh0, dims=1)
# LESelas = [0.7, 1.0, 1.5]
# Frisectorsh = -3
# LESbeta = (LESelas .* thetah) ./ sum(LESelas .* thetah)
# LESbetachk = sum(LESbeta)
# LESsub = consh0 + LESbeta * (dinch0 / Frisectorsh)


elasdir = joinpath(@__DIR__, "..", "data", "Elasticities.csv")
elas = CSV.read(elasdir, DataFrames.DataFrame, header=1)
# 3. Calibrate the model
# 3.0 Elasticities
sigmap = elas[sectors, "sigmap"]
sigmakel = elas[sectors, "sigmakel"]
sigmakl = elas[sectors, "sigmakl"]
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
    pl = 1
    @variables m begin
        pprod[i = sectors], (start = pprod0[i])
        pprodt[i = sectors], (start = 1)
        pk, (start = 1) 
        pint[i = sectors], (start = 1)
        pkl[i = sectors], (start = 1)
        prod[i = sectors], (start = prod0[i])
        int[i = sectors], (start = int0[i])
        kl[i = sectors], (start = kl0[i])
        k[i = sectors], (start = k0[i])
        l[i = sectors], (start = l0[i])
        qint[j = sectors, i = sectors], (start = qint0[j, i])
        consh[j = sectors, h = households], (start=consh0[j, h])
        inch[h = households], (start=inch0[h])
        dinch[h = households], (start=dinch0[h])
        incg, (start=incg0)
        expg, (start=incg0)
        consg[i = sectors], (start=consg0[i])
        walras, (start = 0)
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

    @mapping(m, eq_pprodt[i in sectors], pprod[i] - (1 - prodtr[i]) * pprodt[i])
    @complementarity(m, eq_pprodt, pprodt)

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
    @mapping(m, eq_inch[h in households], (pl * ls[h] + pk * ks * akh[h] + transfr[h]) - inch[h])
    @complementarity(m, eq_inch, inch)

    @mapping(m, eq_dinch[h in households], inch[h] * (1 - taxh[h]) - dinch[h])
    @complementarity(m, eq_dinch, dinch)

    @mapping(m, eq_incg, sum(inch[h] * taxh[h] for h in households) + sum(prod[i] * pprodt[i] * prodtr[i] for i in sectors) - incg)
    @complementarity(m, eq_incg, incg)

    @mapping(m, eq_expg, incg - expg)
    @complementarity(m, eq_expg, expg)

    @mapping(m, eq_consg[i in sectors], thetag[i] * (expg - sum(transfr[h] for h in households)) - consg[i] * pprodt[i])
    @complementarity(m, eq_consg, consg)

    # demand
    @mapping(m, eq_consh[i in sectors, h in households], thetah[i, h] * dinch[h] - consh[i, h] * pprodt[i])
    @complementarity(m, eq_consh, consh)

    # demand - supply
    @mapping(m, eq_prod[i in sectors], prod[i] - (sum(consh[i, h] for h in households) + consg[i] + sum(qint[i, j] for j in sectors)))
    @complementarity(m, eq_prod, prod)

    # factor supply
    @mapping(m, eq_pk, ks - sum(k[i] for i in sectors))
    @complementarity(m, eq_pk, pk)

    @mapping(m, eq_pl, sum(ls[h] for h in households) + walras - sum(l[i] for i in sectors))
    @complementarity(m, eq_pl, walras)

    # Model Solver
    # print(ComplementarityType[1])
    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=6000)
    # @show result_value.(pgdp)
    # @show result_value.(gdp)
    # println("Walras is:")
    @show result_value.(walras)
    @show result_value.(pk)
end

solve_cge()

ls = ls .* 2
ks = ks * 4
solve_cge()
