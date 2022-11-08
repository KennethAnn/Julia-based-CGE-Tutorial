# Julia-based-CGE-Tutorial
A tutorial for Julia-based CGE modeling with lots of cases

### *1. Simple Model*

#### CGE1: *A simple model with two sectors and no intermediate input*

* **Notes: For the initial solution, the non-linear equations should be balanced**

#### CGE2: *A Simple CGE Model with three sectors and intermediate input*

#### CGE3: *A Simple CGE Model with three sectors and intermediate input, as well as Linear Expenditure System*

#### CGE4: *A Simple CGE Model with three sectors and intermediate input, as well as Linear Expenditure System, add the price base*

#### CGE5: *A Simple CGE Model with government*

#### CGE6: *A Simple CGE Model with one household and one government*

* ***Keynote: testing the balance of walras variable is the key for CGE model. It is the strong evidence that the model specification is wrong when the walras is not zero automatically. Therefore, checking the balance of walras variable is quite important at any time. 2022.09.07 17:45*** 

### *2. National Model*

#### CGE7: *A Large CGE Model with 42 sectors, one household and one government*

* **Key note: some variables may be negative. To do the comlementarity condition, please do not declare the variable as positive, i.e. no bound** (Since of the export and import, here the aggregation demand may be negative. Another case is about the investment, some of them are negative)
* ***Walras is right***

#### CGE8: *A Large CGE Model with 42 sectors, two households and one government*

#### CGE9:  *A Large CGE Model with 42 sectors, 2 households, 1 government and investment*

#### CGE10: *A open-economy CGE model with 42 sectors, 2 households, 1 government, 1 investment and 1 foreigner*

* ***change the expression of export CET function  (2022.9.10)***

* ```
  # export CET
  @mapping(m, eq_ext[j in sectors], ext[j] * pext[j] ^ sigmae[j] -  aext[j] * prod[j] * pprodt[j] ^ sigmae[j])
  @complementarity(m, eq_ext, ext)
  
  @mapping(m, eq_dms[j in sectors], dms[j] * pdm[j] ^ sigmae[j] - adms[j] *  prod[j] * pprodt[j] ^ sigmae[j])
  @complementarity(m, eq_dms, dms)
  
  @mapping(m, eq_prod[j in sectors], aext[j] * pext[j] ^ (1 - sigmae[j]) + adms[j] * pdm[j] ^ (1 - sigmae[j]) - pprodt[j] ^ (1 - sigmae[j]))
  @complementarity(m, eq_prod, prod)
  ```

* ***change the international capital flow  (2022.9.10)***

* ```
  @mapping(m, eq_pext[j in sectors], pext[j] - er)
  @complementarity(m, eq_pext, pext)
  
  @mapping(m, eq_pimt[j in sectors], pimt[j] - er)
  @complementarity(m, eq_pimt, pimt)
  
  # investment - saving
  @mapping(m, eq_is, sum(savh[h] for h in households) + savg + aksav * pk * ks - tinv * ptinv - er * inve0)
  @complementarity(m, eq_is, tinv)
  
  @mapping(m, eq_inve, sum(ext[j] for j in sectors) - (sum(imt[j] for j in sectors) + inve0))
  @complementarity(m, eq_inve, er)
  ```

#### CGE11: *A energy-economy CGE model, considering the energy-KL bundle*

#### CGE12: *A energy-economy-climate CGE model, adding the carbon accounting and carbon pricing module*

* **Key note: use the @NLparameter to adjust the parameters in expressions dynamically**
* **Key note: walras could be calculated after simulation, if walras value is almost 0, CGE is balanced and solved normally**

#### CGE13: *A recursive-dynamic energy-economy-climate CGE model, adding the carbon accounting and carbon pricing module*

#### CGE 14: Using labor efficiency and capital efficiency of each sector to reflect the tfp growth

```
# k + l = kl
@mapping(m, eq_k[i in sectors], lambdak[i] * k[i] * pk ^ sigmakl[i] - ak[i] * kl[i] * (lambdak[i] * pkl[i]) ^ sigmakl[i])
@complementarity(m, eq_k, k)

@mapping(m, eq_l[i in sectors], lambdal[i] * l[i] * pl ^ sigmakl[i] - al[i] * kl[i] * (lambdal[i] * pkl[i]) ^ sigmakl[i])
@complementarity(m, eq_l, l)

@mapping(m, eq_pkl[i in sectors], ak[i] * (pk / lambdak[i])^ (1- sigmakl[i]) + al[i] * (pl / lambdal[i]) ^ (1 - sigmakl[i]) - pkl[i] ^ (1 - sigmakl[i]))
@complementarity(m, eq_pkl, pkl)

# determine the tfp edogenously by setting rgdp
@mapping(m, eq_rgdp_exg, rgdp_exg - rgdp_cs)
@complementarity(m, eq_rgdp_exg, tfp)

# determine the tfp exogenously by setting tfp_exg
@mapping(m, eq_rgdp_exg, tfp - tfp_exg)
@complementarity(m, eq_rgdp_exg, tfp)
```
