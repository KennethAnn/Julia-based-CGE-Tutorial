# A energy-economy-climate CGE model, adding the carbon accounting and carbon pricing module
using JuMP, Complementarity, DataFrames, CSV
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

samList = ["AgFFF", "Coal", "OilGas", "MetMin", "NMetMin", "FoodPr", "Textile", "Apparel",
            "WoodPr", "PaperPr", "RefPet", "Chemical", "NMetPr", "Metals", "MetalPr",
            "GenEqp", "SpecEqp", "TransEqp", "ElecEqp", "ICTEqp", "PrecInst", "OthMfg",
            "Waste", "MachRep", "ElecDist","GasDist", "WatDist", "Constr", "WhRetTr", "TranspSrv",
            "HotRest", "ICTServ", "Finance", "RealEst", "BusServe", "ResTech", "EnvServ",
            "ResServ", "Education", "Health", "RecEnt", "PubAdm", 
            "L", "K", "T", "RuralHH", "UrbanHH", "S", "E"]

sector = ["AgFFF", "Coal", "OilGas", "MetMin", "NMetMin", "FoodPr", "Textile", "Apparel",
            "WoodPr", "PaperPr", "RefPet", "Chemical", "NMetPr", "Metals", "MetalPr",
            "GenEqp", "SpecEqp", "TransEqp", "ElecEqp", "ICTEqp", "PrecInst", "OthMfg",
            "Waste", "MachRep", "ElecDist","GasDist", "WatDist", "Constr", "WhRetTr", "TranspSrv",
            "HotRest", "ICTServ", "Finance", "RealEst", "BusServe", "ResTech", "EnvServ",
            "ResServ", "Education", "Health", "RecEnt", "PubAdm"]

sectors = collect(1:1:length(sector))

energy = ["Coal", "RefPet", "ElecDist", "GasDist"]
ens = [2, 11, 25, 26] # energy sectors or sectors
fens = [2, 11, 26] # fossil fuel sectors or sectors
non_ens = setdiff(sectors, ens) # non-energy sectors
non_fens = setdiff(sectors, fens) # non-energy sectors

household = ["RuralHH", "UrbanHH"]
households = collect(1:1:length(household))
va = ["L","K","T"]
NumSector = length(sectors)
NumHouseholds = length(households)

samdir = joinpath(@__DIR__, "..", "data", "SAM2017.csv")
sam = CSV.read(samdir, DataFrames.DataFrame, header=1)
sam = Matrix(sam)[1:50, 2:51] # SAM table

# production
qint0 = sam[sectors, sectors]
l0 = sam[NumSector + 1, sectors]
k0 = sam[NumSector + 2, sectors]
ld0 = sum(l0)
kd0 = sum(k0)
ks = kd0
kl0 = k0 + l0
ak = k0 ./ kl0
al = l0 ./ kl0

## Energy
ene0 = sum(sam[ens, sectors], dims = 1)[1, :] # total energy input
aens = sam[sectors, sectors] ./ sum(sam[ens, sectors], dims = 1) # share of energy sources in sector i
fene0 = sum(sam[fens, sectors], dims = 1)[1, :] # total fossil energy input
elec0 = sam[25, sectors] # electricity input
afene = fene0 ./ ene0 # fossil fuel share
aelec = elec0 ./ ene0 # electricity share
afens = sam[sectors, sectors] ./ sum(sam[fens, sectors], dims = 1) # each fossil fuel share in the total fossil fuel inputs

## KEL
kel0 = kl0 + ene0 # KEL bundle
akl = kl0 ./ kel0 # share of KL in KEL bundle
aene = ene0 ./ kel0 # share of ene in KEL bundle

## Int
int0 = sum(sam[non_ens, sectors], dims=1)[1, :] # total non-energy intermediate input of sector in
aqint = sam[sectors, sectors] ./ sum(sam[non_ens, sectors], dims=1) # share of intermediate input

prod0 = int0 + kel0 + sam[NumSector + 3, sectors]
aint = int0 ./ prod0 
akel = kel0 ./ prod0
prodtr = sam[NumSector + 3, sectors]./ prod0
pprod0 =(1 .- prodtr)

# income block
ls = sam[households .+ (NumSector + 3), NumSector + 1]
akh = sam[households .+ (NumSector + 3), NumSector + 2] ./ kd0
aksav = sam[NumSector + 3 + NumHouseholds + 2, NumSector + 2] / kd0
akg = sam[NumSector + 3 + NumHouseholds + 1, NumSector + 2] / kd0

transfr = sam[(NumSector + 3) .+ households, NumSector + 3 + NumHouseholds + 1]
inch0 = akh * ks + ls + transfr
taxh = sam[NumSector + 3 + NumHouseholds + 1, (NumSector + 3) .+ households] ./ inch0
dinch0 = inch0 .* (1 .- taxh)
savh0 = sam[NumSector + 3 + NumHouseholds + 2, (NumSector + 3) .+ households]
asavh = savh0 ./ dinch0

# LES system
consh0 = sam[sectors, households .+ (NumSector + 3)]
thetah = consh0 ./ sum(consh0, dims=1)
# LESelas = [0.7, 1.0, 1.5]
# Frisectorsh = -3
# LESbeta = (LESelas .* thetah) ./ sum(LESelas .* thetah)
# LESbetachk = sum(LESbeta)
# LESsub = consh0 + LESbeta * (dinch0 / Frisectorsh)

incg0 = sum(taxh[h] .* inch0[h] for h in households) + sum(sam[NumSector+3, sectors]) + akg * ks
consg0 = sam[sectors, NumSector + 3 + NumHouseholds + 1]
thetag = consg0 ./ sum(consg0)
savg0 = sam[NumSector + 3 + NumHouseholds + 2, NumSector + 3 + NumHouseholds + 1]
asavg = savg0 / incg0

inv0 = sam[sectors, NumSector + 3 + NumHouseholds + 2]
tinv0 = sum(inv0)
ainv = inv0 ./ tinv0

# foreigner system
inve0 = (-1) * sam[NumSector + 3 + NumHouseholds + 2, NumSector + 3 + NumHouseholds + 3]
tdd0 = sum(qint0, dims=2)[:, 1] + sum(consh0, dims=2)[:, 1] + consg0 + inv0 # total domestic demand
imt0 = sam[NumSector + 3 + NumHouseholds + 3, sectors] # import from ROC
dmd0 = tdd0 - imt0

admd = dmd0 ./ tdd0
aimt = imt0 ./ tdd0

ext0 = sam[sectors, NumSector + 3 + NumHouseholds + 3] # export to ROC
adms = dmd0 ./ prod0
aext = ext0 ./ prod0

elasdir = joinpath(@__DIR__, "..", "data", "Elasticities.csv")
elas = CSV.read(elasdir, DataFrames.DataFrame, header=1)

# 2.2 Elasticities
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

# 2.3 Emissions
emisdir = joinpath(@__DIR__, "..", "data", "Emissions2017.csv") # emission factors
emission_factor =  CSV.read(emisdir, DataFrames.DataFrame, header=1)
emission_factor = Matrix(emission_factor)[1:44, 2:43] # emission_fac table

emit0 = zeros(Float32, NumSector)
emith0 = zeros(Float32, NumHouseholds)
for j in fens
    for i in sectors
        emit0[i] = sum(emission_factor[i, j] * qint0[j, i] for j in fens)
    end
    for h in households
        emith0[h] = sum(emission_factor[42 + h, j] * consh0[j, h] for j in fens)
    end
end
temit0 = sum(emit0[i] for i in sectors) + sum(emith0[h] for h in households)
emitcap = temit0 * 0.15

# 4. Generate CGE Model
function generate_CGE()  
    m = MCPModel()
    ptinv = 1
    @variables m begin
        pprod[i = sectors], (start = pprod0[i])
        pprodt[i = sectors], (start = 1)
        parm[i = sectors], (start = 1)
        pk, (start = 1) 
        pl, (start = 1)        
        pint[i = sectors], (start = 1)
        pkl[i = sectors], (start = 1)
        pene[i = sectors], (start = 1)
        pfene[i = sectors], (start = 1)
        pkel[i = sectors], (start = 1)
        pdm[i = sectors], (start = 1)
        pimt[i = sectors], (start = 1)
        pext[i = sectors], (start = 1)
        prod[i = sectors], (start = prod0[i])
        int[i = sectors], (start = int0[i])
        kl[i = sectors], (start = kl0[i])
        k[i = sectors], (start = k0[i])
        l[i = sectors], (start = l0[i])
        kel[i = sectors], (start = kel0[i])
        ene[i = sectors], (start = ene0[i])
        fene[i = sectors], (start = fene0[i])
        qint[j = sectors, i = sectors], (start = qint0[j, i])
        consh[j = sectors, h = households], (start=consh0[j, h])
        savh[h = households], (start=savh0[h])
        inch[h = households], (start=inch0[h])
        dinch[h = households], (start=dinch0[h])
        incg, (start=incg0)
        savg, (start=savg0)
        expg, (start=incg0)
        consg[i = sectors], (start=consg0[i])
        inv[j = sectors], (start=inv0[j])
        tinv, (start=tinv0)
        er, (start=1)
        # inve, (start=inve0)
        tdd[i=sectors], (start=tdd0[i])
        imt[i=sectors], (start=imt0[i])
        dmd[i=sectors], (start=dmd0[i])
        dms[i=sectors], (start=dmd0[i])
        ext[i=sectors], (start=ext0[i])
        emit[i=sectors], (start=emit0[i])
        emith[h=households], (start=emith0[h])
        temit, (start=temit0)
        pcarbon>=0, (start=0)
        etax[j=fens, i=sectors], (start=0)
        # etaxh[j=fens, i=households], (start=0.1)
        # rgdp, (start=rgdp0)
        # pgdp, (start=1)
    end
    @NLparameter(m, cap == temit0)
    # int + kel = prod
    @mapping(m, eq_int[i in sectors], int[i] * pint[i] ^ sigmap[i] - aint[i] * prod[i] * pprod[i] ^ sigmap[i])
    @complementarity(m, eq_int, int)

    @mapping(m, eq_kel[i in sectors], kel[i] * pkel[i] ^ sigmap[i] - akel[i] * prod[i] * pprod[i] ^ sigmap[i])
    @complementarity(m, eq_kel, kel)

    @mapping(m, eq_pprod[i in sectors], pprod[i] ^ (1 - sigmap[i])- aint[i] * pint[i] ^ (1- sigmap[i]) - akel[i] * pkel[i] ^ (1 - sigmap[i]))
    @complementarity(m, eq_pprod, pprod)

    @mapping(m, eq_pprodt[i in sectors], pprod[i] - (1 - prodtr[i]) * pprodt[i])
    @complementarity(m, eq_pprodt, pprodt)

    # kl + ene = kel
    @mapping(m, eq_kl[i in sectors], kl[i] * pkl[i] ^ sigmakel[i] - akl[i] * kel[i] * pkel[i] ^ sigmakel[i])
    @complementarity(m, eq_kl, kl)

    @mapping(m, eq_ene[i in sectors], ene[i] * pene[i] ^ sigmakel[i] - aene[i] * kel[i] * pkel[i] ^ sigmakel[i])
    @complementarity(m, eq_ene, ene)

    @mapping(m, eq_pkel[i in sectors], pkel[i] ^ (1 - sigmakel[i])- akl[i] * pkl[i] ^ (1- sigmakel[i]) - aene[i] * pene[i] ^ (1 - sigmakel[i]))
    @complementarity(m, eq_pkel, pkel)

    # k + l = kl
    @mapping(m, eq_k[i in sectors], k[i] * pk ^ sigmakl[i] - ak[i] * kl[i] * pkl[i] ^ sigmakl[i])
    @complementarity(m, eq_k, k)

    @mapping(m, eq_l[i in sectors], l[i] * pl ^ sigmakl[i] - al[i] * kl[i] * pkl[i] ^ sigmakl[i])
    @complementarity(m, eq_l, l)

    @mapping(m, eq_pkl[i in sectors], ak[i] * pk ^ (1- sigmakl[i]) + al[i] * pl ^ (1 - sigmakl[i]) - pkl[i] ^ (1 - sigmakl[i]))
    @complementarity(m, eq_pkl, pkl)

    # qint ++ = int
    @mapping(m, eq_qint[j in non_ens, i in sectors], qint[j, i] - aqint[j, i] * int[i])
    @complementarity(m, eq_qint, qint[non_ens, sectors])

    @mapping(m, eq_pint[i in sectors], sum(aqint[j, i] * parm[j] for j in non_ens) - pint[i])
    @complementarity(m, eq_pint, pint)

    # fene + elec = ene
    @mapping(m, eq_fene[i in sectors], fene[i] * pfene[i] ^ sigmaene[i] - afene[i] * ene[i] * pene[i] ^ sigmaene[i])
    @complementarity(m, eq_fene, fene)

    @mapping(m, eq_ielec[i in sectors], qint[25, i] * parm[25] ^ sigmaene[i] - aelec[i] * ene[i] * pene[i] ^ sigmaene[i])
    @complementarity(m, eq_ielec, qint[25, sectors])

    @mapping(m, eq_pene[i in sectors], pene[i] ^ (1 - sigmaene[i])- afene[i] * pfene[i] ^ (1- sigmaene[i]) - aelec[i] * parm[25] ^ (1 - sigmaene[i]))
    @complementarity(m, eq_pene, pene)

    # fens ++ = fene
    @mapping(m, eq_fens[j in fens, i in sectors], qint[j, i] * (parm[j] + etax[j, i]) ^ sigmafe[i] - afens[j, i] * fene[i] * pfene[i] ^ sigmafe[i])
    @complementarity(m, eq_fens, qint[fens, sectors])

    @mapping(m, eq_pfene[i in sectors], pfene[i] ^ (1 - sigmafe[i])- sum(afens[j, i] * (parm[j] + etax[j, i]) ^ (1- sigmafe[i]) for j in fens))
    @complementarity(m, eq_pfene, pfene)

    # income
    @mapping(m, eq_inch[h in households], (pl * ls[h] + pk * ks * akh[h] + transfr[h]) - inch[h])
    @complementarity(m, eq_inch, inch)

    @mapping(m, eq_dinch[h in households], inch[h] * (1 - taxh[h]) - dinch[h])
    @complementarity(m, eq_dinch, dinch)

    @mapping(m, eq_incg, sum(inch[h] * taxh[h] for h in households) + sum(prod[i] * pprodt[i] * prodtr[i] for i in sectors) + akg * pk * ks + pcarbon * temit - incg)
    @complementarity(m, eq_incg, incg)

    @mapping(m, eq_expg, incg - expg)
    @complementarity(m, eq_expg, expg)

    # households demand
    @mapping(m, eq_consh[i in sectors, h in households], thetah[i, h] * dinch[h] * (1 - asavh[h]) - consh[i, h] * parm[i])
    @complementarity(m, eq_consh, consh)
    
    @mapping(m, eq_savh[h in households], savh[h] - (dinch[h] - sum(consh[j, h] * parm[j] for j in sectors)))
    @complementarity(m, eq_savh, savh)

    # government demand
    @mapping(m, eq_consg[i in sectors], thetag[i] * (expg - sum(transfr[h] for h in households)) * (1 - asavg)- consg[i] * parm[i])
    @complementarity(m, eq_consg, consg)

    @mapping(m, eq_savg, savg - (incg - sum(consg[j] * parm[j] for j in sectors)))
    @complementarity(m, eq_savg, savg)

    # investment
    @mapping(m, eq_inv[j in sectors], inv[j] * parm[j] ^ sigmav - ainv[j] * tinv * ptinv ^ sigmav)
    @complementarity(m, eq_inv, inv)

    @mapping(m, eq_ptinv, ptinv ^ (1 - sigmav) - sum(ainv[j] * parm[j] ^ (1 - sigmav) for j in sectors))
    @complementarity(m, eq_ptinv, tinv)

    # import CES
    @mapping(m, eq_tdd[j in sectors], tdd[j] - (sum(qint[j, i] for i in sectors) + sum(consh[j, h] for h in households) + consg[j] + inv[j]))
    @complementarity(m, eq_tdd, tdd)

    @mapping(m, eq_imt[j in sectors], imt[j] * pimt[j] ^ sigmam[j] - aimt[j] * tdd[j] * parm[j] ^ sigmam[j])
    @complementarity(m, eq_imt, imt)

    @mapping(m, eq_dmd[j in sectors], dmd[j] * pdm[j] ^ sigmam[j] - admd[j] * tdd[j] * parm[j] ^ sigmam[j])
    @complementarity(m, eq_dmd, dmd)

    @mapping(m, eq_parm[j in sectors], aimt[j] * pimt[j] ^ (1 - sigmam[j]) + admd[j] * pdm[j] ^ (1 - sigmam[j]) - parm[j] ^ (1- sigmam[j]))
    @complementarity(m, eq_parm, parm)
    @mapping(m, eq_ext[j in sectors], ext[j] * pext[j] ^ sigmae[j] -  aext[j] * prod[j] * pprodt[j] ^ sigmae[j])
    @complementarity(m, eq_ext, ext)

    @mapping(m, eq_dms[j in sectors], dms[j] * pdm[j] ^ sigmae[j] - adms[j] *  prod[j] * pprodt[j] ^ sigmae[j])
    @complementarity(m, eq_dms, dms)

    @mapping(m, eq_prod[j in sectors], aext[j] * pext[j] ^ (1 - sigmae[j]) + adms[j] * pdm[j] ^ (1 - sigmae[j]) - pprodt[j] ^ (1 - sigmae[j]))
    @complementarity(m, eq_prod, prod)

    @mapping(m, eq_pdm[j in sectors], dms[j] - dmd[j])
    @complementarity(m, eq_pdm, pdm)

    @mapping(m, eq_pext[j in sectors], pext[j] - er)
    @complementarity(m, eq_pext, pext)

    @mapping(m, eq_pimt[j in sectors], pimt[j] - er)
    @complementarity(m, eq_pimt, pimt)
    
    # investment - saving
    # @mapping(m, eq_is, sum(savh[h] for h in households) + savg + aksav * pk * ks - tinv * ptinv - er * inve0)
    # @complementarity(m, eq_is, walras)

    @mapping(m, eq_inve, sum(ext[j] for j in sectors) - (sum(imt[j] for j in sectors) + inve0))
    @complementarity(m, eq_inve, er)

    # factor supply
    @mapping(m, eq_pk, ks - sum(k[i] for i in sectors))
    @complementarity(m, eq_pk, pk)

    @mapping(m, eq_pl, sum(ls[h] for h in households) - sum(l[i] for i in sectors))
    @complementarity(m, eq_pl, pl)

    # climate policy
    @mapping(m, eq_emit[i in sectors], emit[i] - sum(emission_factor[i, j] * qint[j, i] for j in fens))
    @complementarity(m, eq_emit, emit)

    @mapping(m, eq_emith[h in households], emith[h] - sum(emission_factor[42 + h, j] * consh[j, h] for j in fens))
    @complementarity(m, eq_emith, emith)

    @mapping(m, eq_temit, temit - (sum(emit[i] for i in sectors) + sum(emith[h] for h in households)))
    @complementarity(m, eq_temit, temit)
   
    @mapping(m, eq_pcarbon, cap - temit)
    @complementarity(m, eq_pcarbon, pcarbon)

    @mapping(m, eq_etax[j in fens, i in sectors], etax[j, i] - emission_factor[i, j] * pcarbon)
    @complementarity(m, eq_etax, etax)

    # @mapping(m, eq_etaxh[j in fens, h in households], etaxh[j, h] - emission_factor[h, j] * pcarbon)
    # @complementarity(m, eq_etaxh, etaxh)


    # for t in [1, 0.9, 0.8]
    #     # Model Solver
    #     set_value(cap, temit0 * t)
    #     status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=6000)
    #     # @show result_value.(pgdp)
    #     # @show result_value.(gdp)
    #     # println("Walras is:")
    #     @show result_value.(walras)
    #     @show result_value.(pcarbon)
    # end
    set_value(cap, temit0 * 10)
    status = solveMCP(m; convergence_tolerance=1e-8, output="yes", time_limit=6000)
    walras = result_value(savh[1]) + result_value(savh[2]) + result_value(savg) + aksav * result_value(pk) * ks - result_value(tinv) - result_value(er) * inve0
    println("walras value = ", walras)
    @show result_value.(pk)
    @show result_value.(temit)

end
solve_cge()

ls = ls * 2
ks = ks * 3
solve_cge()