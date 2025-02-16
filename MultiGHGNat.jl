# Adapted from WiNDC National Model
using MPSGE
using DataFrames, JLD2
using JuMP
using MPSGE.JuMP.Containers
using CSV, Plots
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
# P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
# P = MultiNatdata
# P = WplusSpAgCdata2020  
P = WplusSpAgCdata2022 # Re-aggregated/disaggregated from BEA detailed via WiNDC.jl, and my AggWiNDCdata.jl script

S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S
Sectors = CSV.read(joinpath(@__DIR__,"./Sectorsplus.csv"), DataFrame);
## New sectors:
# uel (utility: Electric power generation, transmission, and distribution), 
# ugs (utility: Natural gas distribution), 
# uwt (utility: Water, sewage and other systems), 
# coa (coal mining)
# min (other non-oil/gas mining [not coal])
### min => coa for CH4sectors, oil => oil & gas.
### emissions and abatement for gas split btw 'gas' and 'pip'?

I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
Ip = [[x for x in I if x∉[:uti]]; [:uel,:ugs, :uwt, :coa, :gas, :rnw]]
Jp = deepcopy(Ip) # Index for WiNDC BEA Sectors
CH4sectors = [:agr,:coa,:gas,:pip,:oil,:wst] #:uti? # subset index for relevant CH4 mitigation sectors (VA slack in benchmark)
VA = [va for va∈S[:va] if va!=:othtax] # Index Value Added (compen = returns to labour/wage, 'surplus' = returns to Kapital)
FD = S[:fd] # Index for final demand categories (pce, and investment types)
TS = S[:ts] #index for taxes/subsidies
YR = S[:yr] # Index for years for potential multi year runs
M = S[:m] # Index for margins (transport and trade)
M = [:trans, :trade]
a_0 = P[:a_0] #	    "Armington supply",
id_0 = P[:id_0] #	"Intermediate demand",
ys_0 = P[:ys_0]#	"Sectoral supply",
va_0 = P[:va_0] #	"Value added",
md_0 = P[:md_0] #	"Margin demand",
fd_0 = P[:fd_0] #	"Final demand",
pce_0 = deepcopy(P[:pce_0])
m_0 = P[:m_0] #	    "Imports",
ms_0 = P[:ms_0] #	"Margin supply",
bopdef_0 = P[:bopdef_0] #	"Balance of payments deficit",
x_0 = P[:x_0] #	    "Exports of goods and services",
y_0 = P[:y_0]  #	"Gross output",
y_0[:fbt,:value]=0; y_0[:mvt,:value]=0; ; y_0[:gmt,:value]=0

## Existing Taxes
ty_0 = P[:ty_0] #	"Output tax rate"
tm_0 = P[:tm_0] #	"Import tariff"; Initial, for price 
ta_0 = P[:ta_0] #	"Tax net subsidy rate on intermediate demand", benchmark data also for price level
ty_0DAA = DenseAxisArray([ty_0[i,:value] for i in Ip], Ip)
tm_0DAA = DenseAxisArray([tm_m0[i,:value] for i in Ip], Ip)
ta_0DAA = DenseAxisArray([ta_m0[i,:value] for i in Ip], Ip)

## Base Marginal Abatement Cost EPA data (2022). $/t and MMtCO2eq
MAC_CH4_data=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2022_data.csv"), DataFrame, header=2, limit=14)
MAC_CH4_totemiss=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2022_data.csv"), DataFrame, header=2, skipto=17)
# # MAC_CH4_data=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2020_data.csv"), DataFrame, header=2, limit=14)
# MAC_CH4_totemiss=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2020_data.csv"), DataFrame, header=2, skipto=17)

#####################
# SET SCENARIO
#####################
yr = Symbol(2022)

ReductTargetbase = 1117.47074039339 # 2022 imputed from Paris (= linear 2005 to 2005*(1-0.61) in 2035, gross emissions reduction assuming 2022 sink) diff to actual 2022
CH4toGHG2005 = .113943621 #fraction of CH4 to gross GHG in 2005, GHG Inv CO2eq=28
CH4oftarget = CH4toGHG2005*ReductTargetbase# % of CH4 in 2005/2035 extrapolation to adjust for different CH4 emissions scenarios

## Set tax rates 
# Paris Target reduction from 2022 = ReductTarget t  (= linear 2005 to 2005*(1-0.61) in 2035, gross emissions reduction assuming 2022 sink) diff to actual 2022
# CO2_taxrate = 45.271#<=re target w oilgas multi
# CH4_taxrate = 60.765#<=re target w oilgas multi
# CO2_taxrate = 52.1865#<=updated target w 20-yr GWP
# CH4_taxrate = 73.21#<=updated target w 20-yr GWP 
# CO2_taxrate = 49.9033#<=updated target w 20-yr GWP & oil/gas x 5 & oil gas only, Reduction targ=1845.07
# CH4_taxrate = 13.4084#<=updated target w 20-yr GWP & oil/gas x 5 & oil gas only, Reduction targ=1845.07 
CO2_taxrate = 45.361#<=target for all CH4
CH4_taxrate = 339.565#<=target for all CH4
# CH4_taxrate = 491.4 #<= target for CH$ oil/gas/pip
# CH4_taxrate = 339.55#<=re target w oilgas multi w/o abatement 
# CO2_taxrate = 200 * 1.130480652# SC CO2 EPA 2023 SCGHG report, 2022 year, 2020US$, central 2% near-term discount rate x BLS CPI adjustment from 2020$ # 45.361 re 1117 target                   64.046 (GHG for Paris, 2022)  2020->8.97 # 9.04 4.5# 68.9330(target just CO2)
# CH4_taxrate = 200 * 1.130480652 <= using SC CO2 because CH4 data is in MtCO2eq # 240.064 w abatement re gross and sink target #  339.565 re 1117 target    old->493.535 (GHG for Paris,2022) 39.725 #  39.99 12.69#

CH4abatement="yes" # Comment out CH4abatement="no" to allow CH4 abatment
# CH4abatement="no" # Until there's also CO2 abatemment, no CH4 abatement by default
Kmobile="yes" # Allow kapital & Labor to flow between sectors (original WiNDC)
# Kmobile="no" # Fix kapital in sectors, allow Labor to flow between
print("CO2tax: $CO2_taxrate, "); println("CH4tax: $CH4_taxrate")
print("$CH4abatement CH4 Abatement: "); print("$Kmobile mobile Kapital: ")

CH4_tax_switch_on_sectors = DenseAxisArray([0 for i in Ip], Ip)
CH4_tax_switch_on_sectors[CH4sectors] .=1 # Neutral: all sectors have tax (but only CH4sectors have non-zero CH4Intensity)
### Turn OFF Methane tax for coal, agriculture, and waste => Tax ONLY for oil, gas, and pip
# CH4_tax_switch_on_sectors[[:coa,:agr,:wst]] .=0 ; print("CH4 gas/oil/pip ONLY✓, ")# sectors here will have NO CH4 tax

CO2eqNonCO2_to_GHGInv_multiplier = 28/25 # Only FYI, not used: The non-CO2 MAC data (emissions and abatement) uses 25, the GHG Inv uses 28
GWP20_multiplier = (81.2)/25 # Emissions and abatement from non-MAC (100-yr, 25 CO2eq) to IPCC AR6 20-year, 81.2 CO2eq.
    GWPmulti = GWP20_multiplier / CO2eqNonCO2_to_GHGInv_multiplier; print("GWPmulti yes?->")
GWPmulti = 1; print("NO")# Comment out for 20-year GWP simulation 
    Undercount_oilgasch4_multiplier = 5; print(" undercount yes?->")#(Plant et al 2022   #Alvarez 2018 said 1.6)
Undercount_oilgasch4_multiplier = 1; print("NO ")# default (comment out for oilgas multi)
    multioil_CH4ratio = (sum(only(MAC_CH4_totemiss[:,["agr_livestock", "agr_rice", "COL", "wst_land", "wst_water"]]))+ only(MAC_CH4_totemiss[:,:GAS])*Undercount_oilgasch4_multiplier)/sum(only(MAC_CH4_totemiss[:,2:end]))
ReductTarget = ReductTargetbase - CH4oftarget + GWPmulti * multioil_CH4ratio *CH4oftarget;print(": Reduction targ=",ReductTarget,", ")
# GWPmulti = 1 ; ReductTarget=copy(ReductTargetbase); print("NO GWP multiplier, Reduction targ = ", ReductTarget)#Default is CH4 with GWP = 25, and CH4 emissions as reported


# MAC_CH4_totemiss[:,:GAS]/sum(only(MAC_CH4_totemiss[:,2:end]))
# Undercount_oilgasch4_multiplier*MAC_CH4_totemiss[:,:GAS]/(sum(only(MAC_CH4_totemiss[:,2:end]))-only(MAC_CH4_totemiss[:,:GAS])+only(MAC_CH4_totemiss[:,:GAS])*Undercount_oilgasch4_multiplier)


# % split for pipelines and oil from MAC for 'GAS' by proportion of combined economic output
gas_of_GAS = sum(ys_0[:gas,:])/(sum(ys_0[:pip,:])+sum(ys_0[:oil,:])+sum(ys_0[:gas,:]))
pip_of_GAS = sum(ys_0[:pip,:])/(sum(ys_0[:pip,:])+sum(ys_0[:oil,:])+sum(ys_0[:gas,:]))# 0.3693442274100097 from just pip and oil b4
oil_of_GAS = sum(ys_0[:oil,:])/(sum(ys_0[:pip,:])+sum(ys_0[:oil,:])+sum(ys_0[:gas,:]))# 0.6306557725899902 from just pip and oil b4
# Aggregate/disaggregate for WiNDC sectors, $/t and MMtCO2eq 
MAC_CH4_WiNDC=DataFrame([MAC_CH4_data[:,:cost_per_t], 
 (MAC_CH4_data[:,:agr_livestock]+MAC_CH4_data[:,:agr_rice]) * GWPmulti,
 MAC_CH4_data[:,:COL] * GWPmulti, 
 MAC_CH4_data[:,:GAS]*gas_of_GAS * GWPmulti * Undercount_oilgasch4_multiplier,
 MAC_CH4_data[:,:GAS]*pip_of_GAS * GWPmulti * Undercount_oilgasch4_multiplier,
 MAC_CH4_data[:,:GAS]*oil_of_GAS * GWPmulti * Undercount_oilgasch4_multiplier,
 (MAC_CH4_data[:,:wst_land]+MAC_CH4_data[:,:wst_water]) * GWPmulti ],
 [:cost_per_t; CH4sectors])
# Aggregate/disaggregate Total CH4 Emissions for WiNDC sectors, and include total VA for convenience, MMt and BillUS$
MAC_CH4_WiNDC_tot_MMt = DataFrame([MAC_CH4_totemiss[:,:cost_per_t],
 (MAC_CH4_totemiss[:,:agr_livestock]+MAC_CH4_totemiss[:,:agr_rice]) * GWPmulti,
 MAC_CH4_totemiss[:,:COL] * GWPmulti, 
 MAC_CH4_totemiss[:,:GAS]*gas_of_GAS * GWPmulti * Undercount_oilgasch4_multiplier,
 MAC_CH4_totemiss[:,:GAS]*pip_of_GAS * GWPmulti * Undercount_oilgasch4_multiplier,
 MAC_CH4_totemiss[:,:GAS]*oil_of_GAS * GWPmulti * Undercount_oilgasch4_multiplier,
 (MAC_CH4_totemiss[:,:wst_land]+MAC_CH4_totemiss[:,:wst_water]) * GWPmulti],
 [:cost_per_t; CH4sectors])
 MAC_CH4_WiNDC_tot = copy(MAC_CH4_WiNDC_tot_MMt); MAC_CH4_WiNDC_tot[:, 2:end] = MAC_CH4_WiNDC_tot_MMt[:,2:end].*10^-3
 push!(MAC_CH4_WiNDC_tot, ["Total va_cost"; [sum(va_0[VA,Symbol(sector)]) for sector in names(MAC_CH4_WiNDC[:,2:end])]])

# Initialise df with 0s to fill with calculations, Million U$ = $/t x MMt
CH4_cost_per_tier = copy(MAC_CH4_WiNDC[:,2:end]); CH4_cost_per_tier[:,:].=0
# For the first row only: Mulitply costs by tons of abatement potential
    for c in CH4sectors
        CH4_cost_per_tier[1,c] = MAC_CH4_WiNDC[1,:cost_per_t]*MAC_CH4_WiNDC[1,c]
    end
# Calculate Costs per tier of marginal abatement: subtract all emissions with abatement potential below tier, and mulitply by cost at that tier
    for c in CH4sectors
        for i in 2:length(CH4_cost_per_tier[:,1])
            CH4_cost_per_tier[i,c] = MAC_CH4_WiNDC[i,:cost_per_t]*(MAC_CH4_WiNDC[i,c]-MAC_CH4_WiNDC[i-1,c])
        end
    end
### Calculate the weighted cost for cumulative abatement at each level: Million U$ = $/t x MMt
    CH4_cumul_costAll = copy(MAC_CH4_WiNDC[:,2:end]); CH4_cumul_costAll.=0 # Intitialize df
    CH4_cumul_costAll[1,:]=copy(CH4_cost_per_tier[1,:])  # 1st Row copied as just tier costs
    # Rest of dataframe sums cost of each tier with all tiers below
    for c in CH4sectors
        for i in 2:length(CH4_cumul_costAll[:,1])
            CH4_cumul_costAll[i,c] = CH4_cost_per_tier[i,c]+CH4_cumul_costAll[i-1,c]
        end
    end
# Value-Added block names for DenseAxisArray indexing
VAMset = [:VAM5,:VAM10,:VAM15,:VAM20,:VAM30,:VAM40,:VAM50,:VAM100,:VAM500,:VAM1000]
# dfs to DenseAxisArrays, update units to BillUS$, and filter to cumulative Marginal costs which are positive for at least one sector i.e. 5th row, $5/t and up
    CH4_cumul_cost = DenseAxisArray(Matrix(CH4_cumul_costAll[5:end,:])*10^-3,VAMset, CH4sectors)
    CH4_EmMitigated = DenseAxisArray(Matrix(MAC_CH4_WiNDC[5:end,2:end])*10^-3, VAMset, CH4sectors)
## CH4 Emissions intensity of Standard value-added activities: total emissions/total value-added cost: MMt/BillUS$
    VASInt = DenseAxisArray(fill(0.,length(Jp)),Jp); for c in CH4sectors; VASInt[c] =  MAC_CH4_WiNDC_tot[1,c]/MAC_CH4_WiNDC_tot[2,c] end
## CH4 emission intensity (emissions/total va cost) for each sector at each abatement tier, subtracts cumulatively abated emissions and adds additional cost of abatement activities MMt/BillUS$
    VAM_CH4EmInt = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Initializing DenseAxisArray
    [VAM_CH4EmInt[v,c] = (MAC_CH4_WiNDC_tot[1,c]-CH4_EmMitigated[v,c])/(sum(va_0[VA,c])+CH4_cumul_cost[v,c]) for v in VAMset for c in CH4sectors]
## Relative cost of VA: (standard va cost + mitigation cost)/(standard cost) - used to multiply BOTH the va[:surplus] and va[:compen] equally in the blocks
    VAM_costover = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Initialise DenseAxisArray
    [VAM_costover[cost,c] = (sum(va_0[VA,c])+CH4_cumul_cost[cost,c])/sum(va_0[VA,c]) for cost in VAMset for c in CH4sectors]


CO2Int = DenseAxisArray(zeros(length(Jp)),Jp)
# 2024 GHG Inv Table 3-5: CO2 Emissions from Fossil Fuel Combustion by Fuel Type and Sector (MMT CO2 Eq.)
TotalCO2EmMMt_coal =  895.9 # <=2022 # 835.6 <=2020  EPA Inventory CO2 Stationary Combustion - Coal sum(Electricity, Industrial, Commercial, & Residential=0)
# Option with no transport emissions in the model: assumption that direct
# TotalCO2EmMMt_gas_oil = 2104.80
Natural_gasCO2 =  1740.70 #<=2022(incl 70.2 transport, Table 3-5:)  1615.7 #<=2020 (incl 58.7 transport, Table 3-5:)
PetroleumCO2 =  2115.40 #<=2022 (incl 1,681.10 transport, Table 3-5) # 1890.0 #<=2020 (incl 1,514.2 transport, Table 3-5)
# Option with all Transport CO2 emissions attributed to oil inputs: assumption that forms of oil fuel all transport that has direct CO2 emissions, and so taxing CO2 is total emissions from all oil as an input 
# TotalCO2EmMMt_gas_oil = Natural_gasCO2 + PetroleumCO2 # EPA inventory all CO2 transport + CO2 Stationary Combustion - both Oil & Natural Gas sum(Electricity, Industrial, Commercial, Residential & oil U.S. Territories)Sta
# (MMtCO2eq/$Bill x 10-^3 -> Gt/$B=t/$) EPA Inventory 2022 sum of CO2 MMt for coal, and for gas & oil per Billion of total benchmark intermediate input from sector
TotalCO2EmGt_coal = TotalCO2EmMMt_coal*10^-3
TotalCO2EmGt_gas = Natural_gasCO2*10^-3
TotalCO2EmGt_oil = PetroleumCO2*10^-3

CO2Int[:coa] =  TotalCO2EmGt_coal/sum(id_0[:coa,:]) 
CO2Int[:oil] =   TotalCO2EmGt_oil/sum(id_0[:oil,:])
CO2Int[:gas] =   TotalCO2EmGt_gas/sum(id_0[:gas,:])

TotCO2bnchmk =  TotalCO2EmGt_coal + TotalCO2EmGt_oil + TotalCO2EmGt_gas
# TotCH4bnchmk = sum(CH4emiss[:EPAemiss,:]) # from single step data version, same data, same value
TotCH4bnchmk = sum(MAC_CH4_WiNDC_tot[1,2:end]) # 0.6902703880400002
TotGHGbnchmk =  TotCO2bnchmk + TotCH4bnchmk # 5.0315703880400005

### NOT USED: 
# Crude oil barrels production 2020: 4144184x10^3
# Average price per barrel 2020: $36.86
# value_of_oil_2020 = 36.86 * 4144184*10^3
# # Natural gas markets production 2020: 36,520,826 x 10^6 ft^3  (total withdrawals: 40,729,927 x 10^6 ft^3)
# # Average natural gas spot price 2020 #3.32/thousand ft^3 $2.03 / million Btu
# value_of_gas_2020 = 3.32 * 36520826 * 10^3
# oil_fraction = value_of_oil_2020/(value_of_gas_2020+value_of_oil_2020)
# imports of crude oil 2020: 2150267 X10^3 barrel ; imports of oil products 727623 x 10^3 barrels
# imports of natural gas 2929: 2.55 trillion ft^3 ; 5.25 trillion ft^3
## End data preparations
###

# ## Upload table of elasticity parameter values drawn from SAGE 2.1.1 documentation and E3 book, with manual concordence
Elasdf=CSV.read(joinpath(@__DIR__,"./data/Elasticities_SAGE_E3.csv"), DataFrame, header=1)
Elas=DenseAxisArray(Matrix(Elasdf[:,2:end]),Symbol.(Elasdf[:,1]),Symbol.(names(Elasdf)[2:end]))
# Test to check significance of specific super high Armington elasticity for oil and gas by setting back to WiNDC 2 
#TODO # Need to get a better reference to determine which end is appropriate (or in the middle) 
# Elas[:oil,:SAGE_nf_Armington] = 2 ;println("ArmingtonOil=2")
# Elas[:gas,:SAGE_nf_Armington] = 2 ;println("ArmingtonGas=2")
# Match uel while uel/rnw split is just a 50/50 split all over to aviod Y[issue
Elas[:rnw,:E3_m_ID] = deepcopy(Elas[:uel,:E3_m_ID])
Elas[:rnw,:E3_e_E_El] = deepcopy(Elas[:uel,:E3_e_E_El])


# function Multiloop(CO2_taxrate,CH4_taxrate) 

MultiNat = MPSGEModel()

@parameters(MultiNat, begin
    ta[Jp], ta_0DAA[Jp]#ta_0[Jp]
    ty[Jp], ty_0DAA[Jp] #ty_0[Jp]
    tm[Jp], tm_0DAA[Jp] #tm_0[Jp]
    CH₄_tax, 0.
    CH4_secs[Jp], CH4_tax_switch_on_sectors[Jp]  #Boolean vector only CH4sectors matter, 1 = tax, 0 = no tax 
    CO₂_tax, 0.     
    # t_elas_y, 0.
    # elas_y, 0.05
    # t_elas_a, 2.            
    # elas_a, 0.
    # elas_dm, 2.
    # d_elas_ra, 1.
end)

@sectors(MultiNat,begin
    Y[Jp],      (description = "Sectoral Production",)
    A[Ip],      (description = "Armington Supply",)
    VAS[Jp],    (description = "Value Added, standard",)
    MS[M],     (description = "Margin Supply",)
    FDem,      (description = "Consumption goods final demand elasticities",)
end)
if (CH4abatement=="yes")
    @sectors(MultiNat,begin
    VAM5[CH4sectors],   (description = "Value Added, with mitigating activity up to \$5/t")
    VAM10[CH4sectors],  (description = "Value Added, with mitigating activity up to \$10/t")
    VAM15[CH4sectors],  (description = "Value Added, with mitigating activity up to \$15/t")
    VAM20[CH4sectors],  (description = "Value Added, with mitigating activity up to \$20/t")
    VAM30[CH4sectors],  (description = "Value Added, with mitigating activity up to \$30/t")
    VAM40[CH4sectors],  (description = "Value Added, with mitigating activity up to \$40/t")
    VAM50[CH4sectors],  (description = "Value Added, with mitigating activity up to \$50/t")
    VAM100[CH4sectors], (description = "Value Added, with mitigating activity up to \$100/t")
    VAM500[CH4sectors], (description = "Value Added, with mitigating activity up to \$500/t")
    VAM1000[CH4sectors],(description = "Value Added, with max mitigating activity, up to \$1000/t")
end)
end
if (Kmobile=="yes")
    @commodities(MultiNat,begin
    PVA[VA], (description = "Value-added Input to VA blocks",) # fully mobile, single price for all compen, and for surplus
end)
elseif (Kmobile=="no")
    @commodities(MultiNat, begin
    PVAK[Jp], (description = "Kapital Input to VA blocks",) # separate from L and fix within sector
    PVAL, (description = "Labour Input to VA blocks",) # separate from K, and mobile btw sectors (1 price)
    end)

end

@commodities(MultiNat,begin
    PA[Ip],   (description = "Armington Price")
    PY[Jp],   (description = "Supply",)
    PVAM[Jp], (description = "Value-added output - Input to Y",)
    PM[M],   (description = "Margin Price",)
    PFX,     (description = "Foreign Exachange",)
    PC,      (description = "Final Consumption price of demand goods",)
end)

# # Variables to track and report levels of CO2 emissions
@auxiliary(MultiNat, CO2em, index = [[:gas, :oil, :coa]])
@auxiliary(MultiNat, CO2TotEm, description = "Total CO2 emissions from fossil fuels")
# Variables to track and report levels of CH4 emissions
@auxiliary(MultiNat, CH4em, index = [[:agr,:coa,:gas,:pip,:oil,:wst]])
@auxiliary(MultiNat, CH4TotEm, description = "Total CH4 emissions")
@auxiliary(MultiNat, TotEm, description = "Total both emissions")


@consumer(MultiNat, RA, description = "Representative Agent")
"""
ID = Ip = [i for i∈Ip if i∉[:oil,:min, :uti]] # separate oilandgas and min, i without oil or min AND separate uti, so i without oil, min, or uti.
:wtr = uti (360) - ide_0 (290) ~ 70
:ele -> uti (360) - %ele = ide_0 (290)
:rnw -> :ele * .4
:elc -> coal for electrity : ele * .2
:elg -> gas for electricity : ele *.4
:oil -> oil - :elg : oil and gas for energy which is all oil and gas except gas for electrity

e is coming from id_0. %s. and then split to direct and elec. So e is a new DAA, with ide_0, input for each.
then split that data for elect, gas, and oil.
    then split elec for rnw, min (for now) and oil (for gas for now)
"""
# Domestic production for all sectors

# for j∈Jp
#         @production(MultiNat, Y[j], [t= 0, s=Elas[j, :SAGE_klem_Y], va=> s=Elas[j,:SAGE_kle_VAE],sm => s = Elas[j,:E3_m_ID]],begin # [t=0, s = 0, sv=> s = 0]
#         [@output(PY[i],ys_0[i,j], t, taxes = [Tax(RA,ty[j])]) for i∈Ip]... 
#         [@input(PA[i], id_0[i,j], sm, taxes = [Tax(RA,CO₂_tax * CO2Int[i])]) for i∈Ip]...
#         @input(PVAM[j], sum(va_0[VA,j]), va)
#         end)
# end

# Test of re-nesting within Y block without disaggregation
ID = [i for i ∈ Ip if i∉[:oil, :coa, :gas, :uel, :pet, :rnw] ] # Intermediate inputs EXCEPT oil and min
for j∈Jp
    @production(MultiNat, Y[j], [t= 0, s=Elas[j, :SAGE_klem_Y], vae=>s=Elas[j,:SAGE_kle_VAE], sm=>s= Elas[j,:E3_m_ID],
    va=>vae=0, En=>vae=Elas[j,:SAGE_ene], PrimENRG=>En=Elas[j,:SAGE_en], Elec=>En=Elas[j,:SAGE_en], # Which of these matter?
    oilgas=>PrimENRG=Elas[j,:E3_e_E_El], 
    inElec=>Elec=Elas[j,:SAGE_en]
    ],begin
    [@output(PY[i],ys_0[i,j], t, taxes = [Tax(RA,ty[j])]) for i∈Ip]...  # CO2 tax out here bc it will go in lower nests
    [@input(PA[i], id_0[i,j], sm) for i∈ID]... # elasticity btw inputs in vector 
     @input(PVAM[j], sum(va_0[VA,j]), va)
        @input(PA[:pet], id_0[:pet,j], PrimENRG) ##for f in [:oil, :gas]]...
          @input(PA[:oil], id_0[:oil,j],   oilgas, taxes=[Tax(RA,CO₂_tax * CO2Int[:oil])]) ##for f in [:oil, :gas]]...
          @input(PA[:gas], id_0[:gas,j]/2, oilgas, taxes=[Tax(RA,CO₂_tax * CO2Int[:gas])])## for f in [:oil, :gas]]...
        @input(PA[:uel], id_0[:uel,j] , Elec) #j==:uel ? 1 : 
          @input(PA[:rnw], id_0[:rnw,j],   inElec)
          @input(PA[:gas], id_0[:gas,j]/2, inElec, taxes=[Tax(RA,CO₂_tax * CO2Int[:gas])])
          @input(PA[:coa], id_0[:coa,j],   inElec, taxes=[Tax(RA,CO₂_tax * CO2Int[:coa])])
end)
end
# # Total value added cost as a function labor (compen) and kapital (surplus), standard (no mitigation)
# print("ElasVA = SAGEkl:: ")
if (Kmobile=="yes")
for j∈Jp
        @production(MultiNat, VAS[j], [t=0, s = 0, va => s = Elas[j,:SAGE_kl_VA]], begin # #     @production(MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin 
        [@output(PVAM[j],sum(va_0[VA,j]), t)]... 
        [@input(PVA[va], va_0[va,j], va, taxes = [Tax(RA,CH₄_tax* CH4_secs[j]* VASInt[j])]) for va∈VA]...
        end)
end
elseif (Kmobile=="no")
    for j∈Jp
        @production(MultiNat, VAS[j], [t=0, s = 0, va => s = Elas[j,:SAGE_kl_VA]], begin # #     @production(MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin 
        [@output(PVAM[j],sum(va_0[VA,j]), t)]... 
        @input(PVAK[j], va_0[:surplus,j], va, taxes = [Tax(RA,CH₄_tax* CH4_secs[j]* VASInt[j])])
        @input(PVAL, va_0[:compen,j], va, taxes = [Tax(RA,CH₄_tax*CH4_secs[j]* VASInt[j])])
        end)
end
end
## Loop over all the Marginal Abatement tiers as Value-Added production blocks
if CH4abatement=="yes"
    VAMcommodSet = [VAM5,VAM10,VAM15,VAM20,VAM30,VAM40,VAM50,VAM100,VAM500,VAM1000]
    if (Kmobile=="yes")
        for c∈CH4sectors
            for vam in VAMcommodSet
                if VAM_costover[vam.name,c]>1 # Some sectors are still cumulatively -negative costs at $5/t, so filtering those out.
                    @production(MultiNat, vam[c], [t=0, s = 0, va => s = Elas[c,:SAGE_kl_VA]], begin
                                                    # @production(MultiNat, vam[j], [t=0, s = 0, va => s = 1], begin
                        [@output(PVAM[c],sum(va_0[VA,c]), t)]... 
                        [@input(PVA[va], va_0[va,c]*VAM_costover[vam.name,c], va, taxes = [Tax(RA, CH₄_tax*CH4_secs[c]*VAM_CH4EmInt[vam.name,c])]) for va∈VA]...
                    end)
                end
            end
        end
    elseif (Kmobile=="no")
        for c∈CH4sectors
            for vam in VAMcommodSet
                if VAM_costover[vam.name,c]>1 # Some sectors are still cumulatively -negative costs at $5/t, so filtering those out.
                    @production(MultiNat, vam[c], [t=0, s = 0, va => s = Elas[c,:SAGE_kl_VA]], begin
                    # @production(MultiNat, vam[j], [t=0, s = 0, va => s = 1], begin
                        [@output(PVAM[c],sum(va_0[VA,c]), t)]... 
                        @input(PVAK[c],va_0[:surplus,c]*VAM_costover[vam.name,c], va, taxes = [Tax(RA, CH₄_tax*CH4_secs[c]*VAM_CH4EmInt[vam.name,c])])
                        @input(PVAL,    va_0[:compen,c]*VAM_costover[vam.name,c], va, taxes = [Tax(RA, CH₄_tax*CH4_secs[c]*VAM_CH4EmInt[vam.name,c])])
                    end)
                end
            end
        end
    end
end

for m∈M
    @production(MultiNat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[i,m] for i∈Ip), t)]...
        [@input(PY[i], ms_0[i,m], s) for i∈Ip]...
    end)
end
 
# for i∈Ip
    # @production(MultiNat, A[i], [t = t_elas_a, s = elas_a, dm => s = elas_dm], begin
    # @production(MultiNat, A[i], [t = 2, s = 0, dm => s = 2], begin
# println("elasdm=SAGEnf")
for i∈Ip
    @production(MultiNat, A[i], [t = 2, s = 0, dm => s = Elas[i,:SAGE_E3_Av_Armington]], begin
        [@output(PA[i], a_0[i,:value], t, taxes=[Tax(RA,ta[i])],reference_price=1-ta_0[i,:value])]... 
        [@output(PFX, x_0[i,:exports], t)]...
        [@input(PM[m], md_0[i,m], s) for m∈M]...
        @input(PY[i], y_0[i,:value], dm)
## Emissions tariff on CH4 goods because inputs
#         @input(PFX, m_0[i,:imports], dm, taxes = [Tax(RA,tm[i]),Tax(RA,CH₄_tax* VASInt[i])],reference_price=1+tm_0[i,:value])# No tariff on CO2 bc oil and gas are taxed as inputs to production, which includes these imports: Tax(RA,CO₂_tax * CO2Int[i]),
#     end)
# end;  println("Yes CH4 tariff: ")  
## Alternative, no tariff on CH4 goods        
@input(PFX, m_0[i,:imports], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[i,:value]) # without excise tariff CH4 goods (or oil, 'coal' i.e. 'min
    end)
end;
println(" No CH4 tariff:: ")

### Final Consumption with elasticity of Demand
    @production(MultiNat, FDem, [t=0, s= 1, p=>s=1], begin
        @output(PC, sum(pce_0),t) #for i∈Ip]    
        [@input(PA[i], pce_0[i,:pce],s) for i in Ip if i∉[:hos, :pet]]...
        @input(PA[:pet], pce_0[:pet,:pce],s)
        @input(PA[:hos], pce_0[:hos,:pce],s)
    end)

if (Kmobile=="yes")
    @demand(MultiNat, RA, begin
    # [@final_demand(PA[i], pce_0[i,:pce]) for i∈Ip]...
    @final_demand(PC, sum(pce_0))
    @endowment(PFX, only(bopdef_0))
    [@endowment(PA[i], -sum(fd_0[i,xfd] for xfd∈FD if xfd!=:pce)) for i∈Ip]...
    [@endowment(PVA[va], sum(va_0[va,j] for j∈Jp)) for va∈VA]...
    end, elasticity = 1)
    # end, elasticity = d_elas_ra)
elseif (Kmobile=="no")
    @demand(MultiNat, RA, begin
    [@final_demand(PA[i], pce_0[i,:pce]) for i∈Ip]...
    # [@final_demand(PA[i], fd_0[i,:pce]) for i∈Ip]...
    @endowment(PFX, only(bopdef_0))
    [@endowment(PA[i], -sum(fd_0[i,xfd] for xfd∈FD if xfd!=:pce)) for i∈Ip]...
    [@endowment(PVAK[j], sum(va_0[:surplus,j])) for j∈Jp]...
    [@endowment(PVAL, sum(va_0[:compen,j])) for j∈Jp]...
    end, elasticity = 1)
end

# ## CO2 emissions for fossil fuel sectors are the activity levels times the (base) total emissions intensity 
@aux_constraint(MultiNat, CO2em[:coa],  CO2em[:coa] - Y[:coa]*TotalCO2EmGt_coal)
@aux_constraint(MultiNat, CO2em[:oil],  CO2em[:oil] - Y[:oil]*TotalCO2EmGt_oil)
@aux_constraint(MultiNat, CO2em[:gas],  CO2em[:gas] - Y[:gas]*TotalCO2EmGt_gas)
## Total CO2 emissions are the sum of emissions from the 2 fossil fuel sectors (constraint expressed as equantion = 0)
@aux_constraint(MultiNat, CO2TotEm, CO2TotEm - (CO2em[:coa] + CO2em[:oil] + CO2em[:gas]))
## CH4 emissions for each CH4 emitting sector are the sum of (either): VA Standard activity levels x standard CH4 emissions intensity (benchmark = 1 x base emissions)
## +/or VA Mitigating activities at each tier, activity level x base emissions x mitigated emissions factor (mitigated intensity/baseline intensity)
for c in CH4sectors
    if CH4abatement=="yes"
        @aux_constraint(MultiNat, CH4em[c],  CH4em[c] - (VAS[c]*MAC_CH4_WiNDC_tot[1,c] +
        ifelse(VAM_costover[:VAM5,c]>1, VAM5[c]   *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM5,c]/VASInt[c] , 0) + # The $5/t tier includes sectors with negative costs, so have to filter those out here (AS WELL as in the production block)
                                        VAM10[c]  *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM10,c]/VASInt[c]+
                                        VAM15[c]  *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM15,c]/VASInt[c]+
                                        VAM20[c]  *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM20,c]/VASInt[c]+
                                        VAM30[c]  *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM30,c]/VASInt[c]+
                                        VAM40[c]  *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM40,c]/VASInt[c]+
                                        VAM50[c]  *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM50,c]/VASInt[c]+
                                        VAM100[c] *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM100,c]/VASInt[c]+
                                        VAM500[c] *MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM500,c]/VASInt[c]+
                                        VAM1000[c]*MAC_CH4_WiNDC_tot[1,c]*VAM_CH4EmInt[:VAM1000,c]/VASInt[c]
        ))
    else
        @aux_constraint(MultiNat, CH4em[c],  CH4em[c] - (VAS[c]*MAC_CH4_WiNDC_tot[1,c]  )) # +
    end
end
 
## Total CH4 Emissions are the sum of emissions from CH4 emitting sectors
@aux_constraint(MultiNat, CH4TotEm, CH4TotEm - (CH4em[:agr] + CH4em[:coa] + CH4em[:gas] + CH4em[:oil] + CH4em[:pip] + CH4em[:wst] ))
## Total GHG (CO2 & CH4) emissions in Mt CO2eq
@aux_constraint(MultiNat, TotEm, TotEm - (CH4TotEm + CO2TotEm))
set_silent(MultiNat)

# Benchmark 
# fix(PA[:oil],1)
# fix(PA[:rec], 1)
fix(RA, sum(pce_0[Ip,:pce])) # Numeraire, fixed at benchmark
# fix(RA, 14972.3)
### Note: Benchmark doesn't solve at 0 interation because of margins of slack activity. 
### Does balance with interactions or slack vars and production commented out.
# solve!(MultiNat , cumulative_iteration_limit = 0)
solve!(MultiNat)
fullvrbnch = generate_report(MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var])
# print(sort(fullvrbnch, :bmkmarg))

TaxRev = DataFrame(solve=Symbol[], Total_Revenue=Float64[], Change_in_Revenue=Float64[], Income=Float64[])
EqVar  = DataFrame(solve=Symbol[], utility=[], utility2=[], Mev=[], Mev2=[], Equiv_Variarion=Float64[], Equiv_Variarion2=Float64[], EV_pcnt= [], EV_pcnt2= [], Excess_Burden=Float64[])

totrevbnch  = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+ # taxes are negative in production, but positive for revenue
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
income      = totrevbnch + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:bnch totrevbnch totrevbnch-totrevbnch income])

totdem = -value(FDem)*value(compensated_demand(FDem,PC))
Mevbnch2 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
Mevbnch = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
EVbnch2  = Mevbnch2 - value(RA) 
EVbnch  = Mevbnch - value(RA) 
elasRA = MPSGE.elasticity(MultiNat.productions[FDem].input)
utilCES = (sum([pce_0[i,:pce]/sum(pce_0)^(1/elasRA)*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/(elasRA)) for i in Ip]))^(elasRA/(elasRA-1)) #
util2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])
util    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip])
push!(EqVar, [:bnch util util2 Mevbnch Mevbnch2 EVbnch EVbnch2 EVbnch/value(RA)*100 EVbnch2/value(RA)*100 -EVbnch])

# Initialize a Dataframe to save final demand results
FDemand = DataFrame(index=Vector{Symbol}(undef, length(Ip)),
desc=Vector{Symbol}(undef, length(Ip)), 
bnch=Vector{Float64}(undef, length(Ip)), 
WiNcntfact=Vector{Float64}(undef, length(Ip)), 
ch4tax=Vector{Float64}(undef, length(Ip)),
co2tax=Vector{Float64}(undef, length(Ip)),
both=Vector{Float64}(undef, length(Ip)),
ch4Qdelta=Vector{Float64}(undef, length(Ip)),
CO2Qdelta=Vector{Float64}(undef, length(Ip)),
bothQdelta=Vector{Float64}(undef, length(Ip)),
bncQpc=Vector{Float64}(undef, length(Ip)),
ch4Qpc=Vector{Float64}(undef, length(Ip)),
CO2Qpc=Vector{Float64}(undef, length(Ip)),
bothQpc=Vector{Float64}(undef, length(Ip))
)
for (n,i) in enumerate(Ip)
    FDemand[n,:index]= i
    FDemand[n,:desc] = Symbol(Sectors[Sectors.index.==string(i),2][1])
    FDemand[n,:bnch] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end
# unfix(RA)
## Check against WiNDC standard counterfactual
set_value!(ta,0)
set_value!(tm,0)

solve!(MultiNat)

totrevWiNcntfac = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
income          = totrevWiNcntfac +  sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:WiNcntfac totrevWiNcntfac totrevWiNcntfac - totrevbnch income])

totdemWiNcntfac = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
MevWiNcntfac2 = prod([(1/value(PA[i]))^((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) for i in Ip])* only(filter(x->x.solve==:WiNcntfac,TaxRev)[!,:Income]) #])
MevWiNcntfac    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:WiNcntfac,TaxRev)[!,:Income]) #]
EVWiNcntfac2  = MevWiNcntfac2 - value(RA) 
EVWiNcntfac  = MevWiNcntfac - value(RA) 
utilWiNcntfac2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) for i in Ip])
utilWiNcntfac    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
push!(EqVar, [:WiNcntfac utilWiNcntfac utilWiNcntfac2 MevWiNcntfac MevWiNcntfac2 EVWiNcntfac EVWiNcntfac2 EVWiNcntfac/value(RA)*100 EVWiNcntfac2/value(RA)*100 -EVWiNcntfac])

fullvrWiNcntfact = generate_report(MultiNat)
rename!(fullvrWiNcntfact, :value => :WiNcntfact, :margin => :WiNcntfactmarg)
fullvrWiNcntfact[!,:var] = Symbol.(fullvrWiNcntfact[:,:var])

for (n,i) in enumerate(Ip)
    FDemand[n,:WiNcntfact] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end

## DataFrame to hold the industry indices with descriptions
Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

## Re-set to benchmark taxes
set_value!(ta,ta_0DAA[Jp])
set_value!(tm,tm_0DAA[Jp])
solve!(MultiNat , cumulative_iteration_limit = 0)# Temp measure to address residual price changes
# sort(generate_report(MultiNat),:value)
# tax are at $/t of CH4(CO2eq)
## "EPA SC CH4 is $1600/t. 
set_value!(CH₄_tax, CH4_taxrate)
set_value!(CO₂_tax,0.) # Set CO2 tax to 0 for running separately.
solve!(MultiNat)

totrevch4   = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
income      = totrevch4 +  sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:ch4 totrevch4 totrevch4-totrevbnch income])

totdemch4 = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
Mevch42 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:ch4,TaxRev)[!,:Income]) #])
Mevch4    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:ch4,TaxRev)[!,:Income]) #]
EVch42  = Mevch42 - value(RA) 
EVch4  = Mevch4 - value(RA) 
utilch42    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemch4) for i in Ip])
utilch4    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
push!(EqVar, [:ch4 utilch4 utilch42 Mevch4 Mevch42 EVch4 EVch42 EVch4/value(RA)*100 EVch42/value(RA)*100 -EVch4])

for (n,i) in enumerate(Ip)
    FDemand[n,:ch4tax] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end

Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

fullvrch4 = generate_report(MultiNat)
rename!(fullvrch4, :value => :ch4tax, :margin => :ch4marg)
fullvrch4[!,:var] = Symbol.(fullvrch4[:,:var])

# print(sort(df, :margin, by= abs, rev=true))
# print(df)
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
set_value!(CO₂_tax, 0.0)
set_value!(ta,ta_0DAA[Jp]);
set_value!(tm,tm_0DAA[Jp]);
solve!(MultiNat)# Temp measure to address residual price changes 
set_value!(CO₂_tax, CO2_taxrate)
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
solve!(MultiNat, cumulative_iteration_limit=10000) #;

totrevco2   = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
income      = totrevco2 +  sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:co2 totrevco2 totrevco2-totrevbnch income])

totdemco2 = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
Mevco22 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:co2,TaxRev)[!,:Income]) #])
Mevco2    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:co2,TaxRev)[!,:Income]) #]
EVco22  = Mevco22 - value(RA) 
EVco2  = Mevco2 - value(RA) 
utilco22    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemco2) for i in Ip])
utilco2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
push!(EqVar, [:co2 utilco2 utilco22 Mevco2 Mevco22 EVco2 EVco22 EVco2/value(RA)*100 EVco22/value(RA)*100 -EVco2])


for (n,i) in enumerate(Ip)
    FDemand[n,:co2tax] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:4]
# Rs[:,2][68:71]

fullvrco2 = generate_report(MultiNat)
rename!(fullvrco2, :value => :CO2tax, :margin => :CO2marg)
fullvrco2[!,:var] = Symbol.(fullvrco2[:,:var])

## Re-set to benchmark taxes
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
set_value!(CO₂_tax, 0.0)
set_value!(ta,ta_0DAA[Jp])
set_value!(tm,tm_0DAA[Jp])
solve!(MultiNat)# Temp measure to address residual price changes
# # Counterfactual Fossil fuel extraction is ALSO taxed at emissions intensitiy of input x tax in $/ton
set_value!(CO₂_tax, CO2_taxrate)
set_value!(CH₄_tax, CH4_taxrate)
solve!(MultiNat, cumulative_iteration_limit=10000) #;

totrevboth  = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
income      = totrevboth + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:both totrevboth totrevboth-totrevbnch income])

totdemboth = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
Mevboth2 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:both,TaxRev)[!,:Income]) #])
Mevboth    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:both,TaxRev)[!,:Income]) #]
EVboth2  = Mevboth2 - value(RA) 
EVboth  = Mevboth - value(RA) 
utilboth2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemboth) for i in Ip])
utilboth    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
push!(EqVar, [:both utilboth utilboth2 Mevboth Mevboth2 EVboth EVboth2 EVboth/value(RA)*100 EVboth2/value(RA)*100 -EVboth])


for (n,i) in enumerate(Ip)
    FDemand[n,:both] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end
# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:5]
# Rs[:,2][68:71]
fullvrboth = generate_report(MultiNat)
rename!(fullvrboth, :value => :bothtaxes, :margin => :bothmarg)
fullvrboth[!,:var] = Symbol.(fullvrboth[:,:var])

#Generate Dataframe with all results (including names expressions)
FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrco2, fullvrboth, fullvrWiNcntfact, on = [:var], makeunique=true);
Compare = FullResults[1:end,[1,2,4,6,8,10]];
## Sum the difference of each tax applied individually
Compare.sum = Compare.ch4tax .- 1 + Compare.CO2tax .-1

## sum difference between the sum of the individual taxes to the combined taxes
Compare.diff = Vector{Union{Missing, Float64}}(undef, length(Compare[:,1]))
for n in 1:length(Compare[:,1])
    if Compare.ch4tax[n] == 0 || Compare.CO2tax[n] == 0
        Compare.diff[n] = missing
    elseif Compare.bothtaxes[n] >= 1 &&  Compare.sum[n] >= 0
        Compare.diff[n] = (Compare.bothtaxes[n] - 1) - Compare.sum[n]
    elseif Compare.bothtaxes[n] >= 1 &&  Compare.sum[n] <= 0
        Compare.diff[n] = 1 - Compare.bothtaxes[n] + Compare.sum[n]
    elseif Compare.bothtaxes[n] <= 1 &&  Compare.sum[n] <= 0
        Compare.diff[n] = 1 - Compare.bothtaxes[n] + Compare.sum[n]
    elseif Compare.bothtaxes[n] <= 1 &&  Compare.sum[n] >= 0
        Compare.diff[n] = 1 - Compare.bothtaxes[n] - Compare.sum[n]
    end
end

Compare[!,:var] = Symbol.(Compare[:,:var]);
# print(sort!(Compare, :var));
# println(sort!(Compare, :diff, by = abs, rev=true))#[1:25,:])
# CSV.write("C:\\Users\\Eli\\Box\\CGE\\MPSGE-JL\\First Mulit GHG taxes Paper\\MultiResults$(Dates.format(now(),"yyyy-mm-d_HhM")).csv", Compare, missingstring="missing", bom=true)
compCO2em = filter(:var => ==(:CO2TotEm), Compare);
compCH4em = filter(:var => ==(:CH4TotEm), Compare);
compTotem = filter(:var => ==(:TotEm), Compare);
Emissions = DataFrame(
["CO2" "Gt" TotCO2bnchmk ;;
only(compCO2em[:,:CO2tax]) ;;
only(compCO2em[:,:ch4tax]) ;;
only(compCO2em[:,:ch4tax]) + only(compCO2em[:,:CO2tax]) - TotCO2bnchmk;;
only(compCO2em[:,:bothtaxes]) ;;
only(compCO2em[:,:ch4tax]) + only(compCO2em[:,:CO2tax]) - TotCO2bnchmk - only(compCO2em[:,:bothtaxes]); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk ;;
only(compCH4em[:,:CO2tax]) ;;
only(compCH4em[:,:ch4tax]) ;;
only(compCH4em[:,:ch4tax]) + only(compCH4em[:,:CO2tax]) - TotCH4bnchmk ;;
only(compCH4em[:,:bothtaxes]) ;;
only(compCH4em[:,:ch4tax]) + only(compCH4em[:,:CO2tax]) - TotCH4bnchmk - only(compCH4em[:,:bothtaxes]); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk ;;
 only(compTotem[:,:CO2tax]) ;;
 only(compTotem[:,:ch4tax]) ;;
 only(compTotem[:,:ch4tax]) + only(compTotem[:,:CO2tax]) - TotGHGbnchmk ;;
 only(compTotem[:,:bothtaxes]) ;;
 only(compTotem[:,:ch4tax]) + only(compTotem[:,:CO2tax]) - TotGHGbnchmk - only(compTotem[:,:bothtaxes])],
["Emissions", "Unit", "Bnchmrk_Emissions", "CO2tax", "CH4tax","sum_of_taxes", "taxes_combined" ,"Interactions"])
Emissions_Mt = hcat(Emissions[:,1:2],Emissions[:,3:end].*10^3); Emissions_Mt[:,2] = ["Mt" ,"MtCO2eq" ,"MtCO2eq"];

EmissionReductionResults = DataFrame(
["CO2" "Gt" TotCO2bnchmk ;;
TotCO2bnchmk - only(compCO2em[:,:CO2tax]) ;;
TotCO2bnchmk - only(compCO2em[:,:ch4tax]) ;;
TotCO2bnchmk - only(compCO2em[:,:ch4tax]) + TotCO2bnchmk - only(compCO2em[:,:CO2tax]) ;;
TotCO2bnchmk - only(compCO2em[:,:bothtaxes]) ;;
TotCO2bnchmk - only(compCO2em[:,:ch4tax]) + TotCO2bnchmk - only(compCO2em[:,:CO2tax]) - (TotCO2bnchmk - only(compCO2em[:,:bothtaxes])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk ;;
TotCH4bnchmk - only(compCH4em[:,:CO2tax]) ;;
TotCH4bnchmk - only(compCH4em[:,:ch4tax]) ;;
TotCH4bnchmk - only(compCH4em[:,:ch4tax]) + TotCH4bnchmk - only(compCH4em[:,:CO2tax]) ;;
TotCH4bnchmk - only(compCH4em[:,:bothtaxes]) ;;
TotCH4bnchmk - only(compCH4em[:,:ch4tax]) + TotCH4bnchmk - only(compCH4em[:,:CO2tax]) - (TotCH4bnchmk - only(compCH4em[:,:bothtaxes])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk ;;
TotGHGbnchmk - only(compTotem[:,:CO2tax]) ;;
TotGHGbnchmk - only(compTotem[:,:ch4tax]) ;;
TotGHGbnchmk - only(compTotem[:,:ch4tax]) + TotGHGbnchmk - only(compTotem[:,:CO2tax]) ;;
TotGHGbnchmk - only(compTotem[:,:bothtaxes]) ;;
TotGHGbnchmk - only(compTotem[:,:ch4tax]) + TotGHGbnchmk - only(compTotem[:,:CO2tax]) - (TotGHGbnchmk - only(compTotem[:,:bothtaxes]))],
["Em_reductions", "Unit", "Bnchmrk_Emissions", "CO2tax_reduc", "CH4tax_reduc","Sum_of_each_tax", "both_taxes_combined" ,"Interactions"])
EmissionReductionResults_Mt = hcat(EmissionReductionResults[:,1:2],EmissionReductionResults[:,3:end].*10^3); EmissionReductionResults_Mt[:,2] = ["Mt" ,"MtCO2eq" ,"MtCO2eq"];
EmissionReductionResults_Mt.pc_diff .= -EmissionReductionResults_Mt.Interactions ./ EmissionReductionResults_Mt.Sum_of_each_tax .* 100 

function plottaxemisscurve(tax1, tax2, start, interval, finish ,vec, RAval, isfixed, cnst=1)
    """# runs a loop increasing each tax by \$1/t and then plotting Total GHG (CO2 & CH4) **incorporated** emissions 
    # Arguments are: which tax to change, other tax to either change simultaneously OR keep at 0, st=initial \$ tax value, fin= final \$ tax value,
    # and final (optional) argument can be set to 0 to remove other tax, by default """
    margemiss = DataFrame(tax1=Float64[], tax2=Float64[], Emissions=Float64[], utility=[],utility2=[], Mev=[], Mev2=[], Equiv_Variarion=Float64[], Equiv_Variarion2=Float64[], 
        EV_pcnt=[], EV_pcnt2=[], CH4Emissions=Float64[],CO2Emissions=Float64[], CH4perc_red=[], CO2perc_red=[])
    rename!(margemiss,:tax1=>tax1.name) , rename!(margemiss,:tax2=>tax2.name)
    Testvars = DataFrame(taxrt=Float64[], 
    # Yagr=Float64[],Ymin=Float64[],Ypip=Float64[],Yoil=Float64[],Apip=Float64[],Aoil=Float64[],Ywst=Float64[],
    # CompDYPApip=Float64[],CompDApipPApip=Float64[],DemRAPApip=Float64[],
    # PAagr=Float64[],PAmin=Float64[],PApip=Float64[],PAoil=Float64[],PAwst=Float64[],PAuti=Float64[],compDApipPAoil=Float64[],compdDAoilPAoil=Float64[],
    # VASagr=Float64[],VAMagr=Float64[],VASmin=Float64[],VAMmin=Float64[],VASpip=Float64[],VAMpip=Float64[],VASoil=Float64[],VAMoil=Float64[],VASwst=Float64[],VAMwst=Float64[],
    VASagr=Float64[],VAM10agr=Float64[],VAM100agr=Float64[],VAM500agr=Float64[],VAMkagr=Float64[],Apip=[],VASmin=Float64[],VAM10min=Float64[],VAM100min=Float64[],VAM500min=Float64[],VAMkmin=Float64[],VASpip=Float64[],VAM10pip=Float64[],VAM100pip=Float64[],VAM500pip=Float64[],VAMkpip=Float64[],VASoil=Float64[],VAM10oil=Float64[],VAM100oil=Float64[],VAM500oil=Float64[],VAMkoil=Float64[],VASwst=Float64[],VAM10wst=Float64[],VAM100wst=Float64[],VAM500wst=Float64[],VAMkwst=Float64[],
    TotEm=Float64[],CH4TotEm=Float64[],CO2TotEm=Float64[]
    # CH4emoil=Float64[],CH4empip=Float64[],CO2emin=Float64[],CO2emoil=Float64[]
    )
    ResultsTroubleshoot = DataFrame(var=[], value=Float64[], margin=Float64[], x1=Float64[]) 
    for (i,j) in zip(start:interval:finish,vec)
        print(i,", ",j,": ")
        set_value!(tax1, i)
        set_value!(tax2, cnst*j)
        solve!(MultiNat, output="no");
        Results = generate_report(MultiNat)
        Results[!,:var] = Symbol.(Results[:,:var]);
       
        totrevboth  = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
        sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
        sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
        sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
        income      = totrevboth + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
        
        totdem = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
        util    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip])
        util2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])
        Mev = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])*income
        Mev2 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])*income
        EV  = Mev - value(RA)
        EV2  = Mev2 - value(RA)

        ResultsTroubleshoot =vcat(ResultsTroubleshoot, [Results fill(i,length(Results[:,1]))])
        push!(margemiss, [i j only(filter(:var => ==(:TotEm), Results)[:, :value]) util util2 Mev Mev2 EV EV2 EV/value(RA)*100 EV2/value(RA)*100 only(filter(:var => ==(:CH4TotEm), Results)[:, :value]) only(filter(:var => ==(:CO2TotEm), Results)[:, :value])  100*(TotCH4bnchmk-only(filter(:var => ==(:CH4TotEm), Results)[:, :value]))/(TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:,:value]))  100*(TotCO2bnchmk-only(filter(:var => ==(:CO2TotEm), Results)[:, :value]))/(TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:,:value]))     ])
        push!(Testvars, [i,                
        # value(Y[:agr]),value(Y[:coa]),value(Y[:pip]),value(Y[:oil]),value(A[:pip]),value(A[:oil]),value(Y[:wst]),
        # value(compensated_demand(MultiNat[:Y][:pip],MultiNat[:PA][:pip])),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:pip])),value(demand(MultiNat[:RA],MultiNat[:PA][:pip])),
        # value(PA[:agr]),value(PA[:coa]),value(PA[:pip]),value(PA[:oil]),value(PA[:wst]),value(PA[:uti]),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:oil])),value(compensated_demand(MultiNat[:A][:oil],MultiNat[:PA][:oil])),
        # value(VAS[:agr]),value(VAM[:agr]),value(VAS[:coa]),value(VAM[:coa]),value(VAS[:pip]),value(VAM[:pip]),value(VAS[:oil]),value(VAM[:oil]),value(VAS[:wst]),value(VAM[:wst]),
        value(VAS[:agr]),value(VAM10[:agr]),value(VAM100[:agr]),value(VAM500[:agr]),value(VAM1000[:agr]),value(A[:pip]),value(VAS[:coa]),value(VAM10[:coa]),value(VAM100[:coa]),value(VAM500[:coa]),value(VAM1000[:coa]),value(VAS[:pip]),value(VAM10[:pip]),value(VAM100[:pip]),value(VAM500[:pip]),value(VAM1000[:pip]),value(VAS[:oil]),value(VAM10[:oil]),value(VAM100[:oil]),value(VAM500[:oil]),value(VAM1000[:oil]),value(VAS[:wst]),value(VAM10[:wst]),value(VAM100[:wst]),value(VAM500[:wst]),value(VAM1000[:wst]),
        value(TotEm),value(CH4TotEm),value(CO2TotEm)
        # value(CH4em[:oil]),value(CH4em[:pip]),value(CO2em[:coa]),value(CO2em[:oil])
        ] 
            )
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    plt = plot(margemiss[!,tax1.name], margemiss[!,:Emissions].*10^3, title= "RA:\$$RAval fxd:$isfixed", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
    plt = plot!([CO2_taxrate], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3))
    # plch4 = plot(margemiss[!,tax1.name], margemiss[!,:CH4Emissions].*10^3, label=false, title="CH4 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
    # plch4 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
    # plco2 = plot(margemiss[!,tax1.name], margemiss[!,:CO2Emissions].*10^3, label=false, title="CO2 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
    # plco2 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
    # plt3 = plot(margemiss[!,tax1.name], margemiss[!,:Emissions].*10^3, title= "Emissions with $tax1 $tax2in", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="\$/t", 
    # xlims=(0,finish), color=:black, linewidth=2, ylab="MMt CO2eq",legend=:left, linestyle=:dashdot)
    # plt3 = plot!(margemiss[!,tax1.name], margemiss[!,:CH4Emissions].*10^3, label="CH4 Emissions", color=:green, linestyle=:dot)
    # plt3 = plot!(margemiss[!,tax1.name], margemiss[!,:CO2Emissions].*10^3, label="CO2 Emissions", color=:blue, linestyle=:dash)
    # # plt3 = plot!([CO2_taxrate], seriestype=:vline, label="SCCO2", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
    # plt3 = plot!(margemiss[!,tax1.name],repeat([TotGHGbnchmk*10^3-ReductTarget],length(margemiss[:,tax1.name])), color=:yellow, label="target",linewidth=1.5)
    # plt3 = plot!(twinx(),margemiss[!,tax1.name],margemiss[!,:EV_pcnt], label="EV%", linewidth=2., xlim=(0,finish), ylim=(0,1),legend=:right)

    ### Next line for CH4 spillovers, only valid for CO2-only policy
    # plt3 = plot!(twinx(),margemiss[2:end,tax1.name],margemiss[2:end,:CH4perc_red], label="Ch4 spillover %\nright axis", color=:orange, yticks=([4.99,5,5.01,5.02,5.03],["4.99%","5%","5.01%","5.02%","5.03%"]), legend=:bottomright)
    ### Next line for CO2 spillovers, only valid for CH4-only policy
    # plt3 = plot!(twinx(),margemiss[2:end,tax1.name],margemiss[2:end,:CO2perc_red], label="CO2 spillover %\nright axis", color=:orange, yticks=([0,10,20,30,40,50,60,70],["0%","10%","20%","30%","40%","50%","60%","70%"]), legend=:topright)
    # plcdyp = plot(margemiss[!,tax1.name],Testvars[!,:CompDYPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(Y:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDYPApip]),maximum(Testvars[!,:CompDYPApip])), xlabel="$tax1 $tax2in \$/t")
    # plcdap = plot(margemiss[!,tax1.name],Testvars[!,:CompDApipPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(A:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDApipPApip]),maximum(Testvars[!,:CompDApipPApip])), xlabel="$tax1 $tax2in \$/t")
    # plfdrap = plot(margemiss[!,tax1.name],Testvars[!,:DemRAPApip], title= "RA:\$$RAval fxd:$isfixed", label="final_dem(RA,PA:pip)", ylim=(minimum(Testvars[!,:DemRAPApip]),maximum(Testvars[!,:DemRAPApip])), xlabel="$tax1 $tax2in \$/t")
    # pla = plot(margemiss[!,tax1.name],Testvars[!,:Yagr], label=false, title="Y:agr", ylim=(minimum(Testvars[!,:Yagr]),maximum(Testvars[!,:Yagr])), xlabel="$tax1 $tax2in \$/t")
    # plm = plot(margemiss[!,tax1.name],Testvars[!,:Ymin], title= "RA:\$$RAval fxd:$isfixed", label="Y:coa", ylim=(minimum(Testvars[!,:Ymin]),maximum(Testvars[!,:Ymin])), xlabel="$tax1 $tax2in \$/t")
    # plp = plot(margemiss[!,tax1.name],Testvars[!,:Ypip], title= "RA:\$$RAval fxd:$isfixed", label="Y:pip", ylim=(minimum(Testvars[!,:Ypip]),maximum(Testvars[!,:Ypip])), xlabel="$tax1 $tax2in \$/t")
    # plo = plot(margemiss[!,tax1.name],Testvars[!,:Yoil], title= "RA:\$$RAval fxd:$isfixed", label="Y:oil", ylim=(minimum(Testvars[!,:Yoil]),maximum(Testvars[!,:Yoil])), xlabel="$tax1 $tax2in \$/t")
    # plw = plot(margemiss[!,tax1.name],Testvars[!,:Ywst], title= "RA:\$$RAval fxd:$isfixed", label="Y:wst", ylim=(minimum(Testvars[!,:Ywst]),maximum(Testvars[!,:Ywst])), xlabel="$tax1 $tax2in \$/t")
    # plpa = plot(margemiss[!,tax1.name],Testvars[!,:PAagr], title= "RA:\$$RAval fxd:$isfixed", label="PA:agr", ylim=(minimum(Testvars[!,:PAagr]),maximum(Testvars[!,:PAagr])), xlabel="$tax1 $tax2in \$/t")
    # plpm = plot(margemiss[!,tax1.name],Testvars[!,:PAmin], title= "RA:\$$RAval fxd:$isfixed", label="PA:coa", ylim=(minimum(Testvars[!,:PAmin]),maximum(Testvars[!,:PAmin])), xlabel="$tax1 $tax2in \$/t")
    # plpp = plot(margemiss[!,tax1.name],Testvars[!,:PApip], title= "RA:\$$RAval fxd:$isfixed", label="PA:pip", ylim=(minimum(Testvars[!,:PApip]),maximum(Testvars[!,:PApip])), xlabel="$tax1 $tax2in \$/t")
    # plpo = plot(margemiss[!,tax1.name],Testvars[!,:PAoil], title= "RA:\$$RAval fxd:$isfixed", label="PA:oil", ylim=(minimum(Testvars[!,:PAoil]),maximum(Testvars[!,:PAoil])), xlabel="$tax1 $tax2in \$/t")
    # plpw = plot(margemiss[!,tax1.name],Testvars[!,:PAwst], title= "RA:\$$RAval fxd:$isfixed", label="PA:wst", ylim=(minimum(Testvars[!,:PAwst]),maximum(Testvars[!,:PAwst])), xlabel="$tax1 $tax2in \$/t")
    # plAp = plot(margemiss[!,tax1.name],Testvars[!,:Apip], title= "RA:\$$RAval fxd:$isfixed", label="A:pip", ylim=(minimum(Testvars[!,:Apip]),maximum(Testvars[!,:Apip])), xlabel="$tax1 $tax2in \$/t")
    # plAo = plot(margemiss[!,tax1.name],Testvars[!,:Aoil], title= "RA:\$$RAval fxd:$isfixed", label="A:oil", ylim=(minimum(Testvars[!,:Aoil]),maximum(Testvars[!,:Aoil])), xlabel="$tax1 $tax2in \$/t")     # Or label=false, title="price of oil commodity"
    return margemiss, plt, Testvars, ResultsTroubleshoot#, plch4, plco2, plt3 #,pla, plm, plp, plo, plw, plpa, plpm, plpp, plpo, plpw, 
    # plcdyp, plcdap, plfdrap, plAp, plAo
end

# set_upper_bound(MultiNat[:A][:pip], 8)
# set_upper_bound(MultiNat[:A][:pet], 4) # 30 is fine past $1000/t
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pet]))
# set_upper_bound(MultiNat[:Y][:sle], 12)
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(Y[:sle]))
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pip]))
# set_silent(MultiNat)
# co2vec = [45.38,41.3,40.9,40.75,39.52,39.4,39.21,39.05,38.9,38.75,38.6,38.45,37.9,37.75,37.6,37.5,37.35,37.19,36.75,36.65,36.5,36.35,36.2,36.05,35.3,35.17,35,34.9,34.75,34.6,34.45,34.35,34.2,34.05,33.91,33.78,32.73,32.6,32.45,32.35,32.2,32.1,31.95,31.85,31.7,31.57,31.45,31.3,30.55,30.45,30.3,30.2,30.1,29.94,29.85,29.7,29.6,29.45,29.35,29.25,28.8,28.68,28.56,28.44,28.33,28.21,28.09,27.97,27.85,27.74,27.62,27.51,27.39,27.28,27.16,27.05,26.93,26.82,26.7,26.58,26.46,26.35,26.23,26.11,25.99,25.87,25.75,25.64,25.52,25.4,25.29,25.18,25.06,24.95,24.84,24.73,24.61,24.5,24.38,24.27,24.15,24.04,23.92,23.81,23.69,23.58,23.46,23.35,23.23,23.12,23,22.89,22.78,22.66,22.55,22.44,22.33,22.21,22.1,22,20.65,20.55,20.45,20.32,20.25,20.15,20.03,19.9,19.82,19.72,19.6,19.48,19.37,19.27,19.16,19.06,18.95,18.85,18.75,18.64,18.54,18.43,18.34,18.25,18.12,18.02,17.92,17.82,17.73,17.63,17.5,17.4,17.3,17.2,17.1,17,16.9,16.8,16.7,16.6,16.5,16.39,16.28,16.18,16.09,16,15.88,15.79,15.69,15.6,15.5,15.39,15.28,15.18,15.07,14.98,14.9,14.79,14.68,14.59,14.5,14.4,14.3,14.2,14.1,14,13.9,13.8,13.7,13.6,13.5,13.4,13.3,13.2,13.1,13,12.9,12.8,12.73,12.64,12.54,12.42,12.32,12.22,12.13,12.04,11.95,11.86,11.76,11.67,11.58,11.47,11.36,11.27,11.18,11.09,11,10.91,10.82,10.71,10.62,10.53,10.44,10.34,10.25,10.15,10.05,9.95,9.85,9.75,9.66,9.58,9.5,9.4,9.31,9.22,9.12,9.03,8.94,8.84,8.75,8.66,8.56,8.47,8.37,8.28,8.19,8.09,8,7.9,7.81,7.72,7.63,7.54,7.46,7.37,7.28,7.19,7.1,7.01,6.92,6.83,6.74,6.65,6.56,6.46,6.37,6.28,6.19,6.09,6,5.91,5.83,5.74,5.65,5.57,5.48,5.39,5.3,5.22,5.13,5.04,4.94,4.85,4.77,4.68,4.6,4.51,4.42,4.33,4.24,4.15,4.06,3.97,3.89,3.8,3.71,3.62,3.53,3.44,3.35,3.27,3.18,3.1,3.01,2.93,2.85,2.76,2.68,2.6,2.51,2.43,2.34,2.26,2.18,2.09,2.01,1.92,1.84,1.75,1.67,1.58,1.5,1.41,1.32,1.23,1.14,1.05,0.99,0.91,0.83,0.74,0.66,0.57,0.49,0.41,0.33,0.24,0.16,0.08,0]
                                # Nothing co2vec = [45.361; 41.325; 40.8; collect(40.4:-((40.4-30.1)/50):30.1) ; collect(30:-((30-23.5)/50):23.5); collect(23:-((23-14.5)/76):14.5) ;collect(14.1:-((14.1-7.4)/76):7.4) ; collect(7:-((7-.36)/76):.36);.26;.175;.1;.03]  # Smoother, for all CH4
# # #4329.824705730001
# co2vec = [45.361; 42.01; 41.8; collect(41.5:-((41.5-26)/150):26) ; collect(26:-((26-11.55)/170):11.55) ; collect(11.5:-((11.5-.2)/163):.2);.18;.15;.1;.01]  # For oil, gas, pip only (up to 491) #491.43 with 0, but not a step
# GWP 81.2, oil gas x 5, oil/gas CH4 tax only
# co2vec = [52.1865; 51; 50; collect(49:-((49-26)/23):26) ; collect(26:-((26-11.55)/23):11.55) ; collect(11.5:-((11.5-.2)/23):.2);0;0;0;.0]  # For oil, gas, pip only (up to 491) #491.43 with 0, but not a step
# checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 1, 74, 
# co2vec
# , round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# # # # # EVdf2 = checkch4CO2[1]
# # # # # EVdf_slice2= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget) && x.Emissions >(TotGHGbnchmk*10^3-ReductTarget-1),EVdf2[:,[1,2,3,11]])# 4329.824705730001
# # # # # print(EVdf_slice2)
# # # # # # ##### checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 10, 400, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# # # checkch4CO2[7]
# # # # print(EVdf2)
# # # checkch4CO2[1]
# # # ##### png(checkch4CO2[7], "./Results/Bothtax-Allemiss")
# # ##### checkch4CO2[2]
# # # EVdf = checkch4CO2[1]
# # # ##### print("checkch4CO2",checkch4CO2[3])
# # # EVdf_slice= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget) && x.Emissions >(TotGHGbnchmk*10^3-ReductTarget-1),EVdf)# 4329.824705730001
# # # EVdf_cap= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget),EVdf)

# resultdf = copy(checkch4CO2[1])
# tax1 = names(resultdf)[1]; tax2 = names(resultdf)[2]

# # plt = plot(resultdf[!,tax1], resultdf[!,:Emissions].*10^3,  label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) & $(replace(tax2,"_"=>" "))  \$/t", xlims=(0,resultdf[end,:CH₄_tax]))#title= "RA:\$$(value(RA)) fxd:$isfixed",)
# # plt = plot!([CO2_taxrate], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3))
# # plch4 = plot(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, label=false, title="CH₄ Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(names(resultdf)[1]) $(names(resultdf)[2]) \$/t", xlims=(0,resultdf[end,:CH₄_tax]))
# # plch4 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
# # plco2 = plot(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label=false, title="CO₂ Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(names(resultdf)[1]) $(names(resultdf)[2])  \$/t", xlims=(0,resultdf[end,:CH₄_tax]))
# # plco2 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
# pltEV = plot(resultdf[!,tax1], resultdf[!,:Emissions].*10^3,  label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", #title= "Emissions with $(names(resultdf)[1]) $(names(resultdf)[2])",
# xlims=(0,resultdf[end,:CH₄_tax]), color=:black, linewidth=2, ylab="MMt CO₂eq",legend=:left, linestyle=:dashdot,yguidefontsize=9, legendfont=font(9,"Palatino Roman"),guidefont=font(10,"Palatino Roman"))
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, label="CH₄ Emissions", color=:green, linestyle=:dot)
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label="CO₂ Emissions", color=:blue, linestyle=:dash)
# pltEV = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
# pltEV = plot!(resultdf[!,tax1],repeat([TotGHGbnchmk*10^3-ReductTarget],length(resultdf[:,tax1])), color=:yellow, label="Reduction target",linewidth=1.5)
# pltEV = plot!(twinx(),resultdf[!,tax1],resultdf[!,:EV_pcnt2], label="Equivalent Variation:\nRight axis", ylabel= "Percentage change", guidefont=font(9,"Palatino Roman"), linewidth=2., xlim=(0,resultdf[end,tax1]), ylim=(0,1),legend=:right, yticks=([0.0,0.2,0.4,0.6,0.8,1.0],["0.0%","0.2%","0.4%","0.6%","0.8%","1.0%"]))
# pltEV = plot!(twiny(),resultdf[!,tax2],resultdf[!,:Emissions].*10^3,  xflip=true, xlim=(resultdf[end,tax2],resultdf[1,tax2]),xlabel="$(replace(tax2,"_"=>" ")) \$/t", linewidth=0, legend=false, guidefont=font(10,"Palatino Roman"))#, ylim=(0,1))
# 1
# # # filter(:EV_pcnt2 => ==(minimum(EVdf_slice.EV_pcnt2)),EVdf_slice)
# # # print(sort(EVdf,:EV_pcnt2)[1:80,:])
# # # print(sort(EVdf_cap,:EV_pcnt2)[1:40,:])
# png(pltEV, joinpath(@__DIR__,"./Results/EqVarTarget2"))
# savefig(pltEV, joinpath(@__DIR__,"./Results/EqVarTarget-oilgaspiponly.svg"))
# 1
# filter(:EV_pcnt => ==(minimum(EVdf_slice.EV_pcnt)),EVdf_slice)
# print(sort(EVdf,:EV_pcnt)[1:120,:])
# print(sort(EVdf_cap,:EV_pcnt)[1:40,:])

##### checkch4CO2[7]
###### Add Quant diff and % of consumption columns for Final Demand report dataframe
# FDemand[:,:ch4Qdelta]=FDemand[:,:ch4tax].-FDemand[:,:bnch]
# FDemand[:,:CO2Qdelta]=FDemand[:,:co2tax].-FDemand[:,:bnch]
# FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
# FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
# FDemand[:,:ch4Qpc]=FDemand[:,:ch4tax]./sum(FDemand[:,:ch4tax])*100
# FDemand[:,:CO2Qpc]=FDemand[:,:co2tax]./sum(FDemand[:,:co2tax])*100
# FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

# fix(RA, sum(fd_0[Ip,:pce]))


# checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, 400, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2[7]
# checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2[2] # Total emissions
# # # println("PApip up lim =", upper_bound(MultiNat[:A][:pip]))
# png(checkCO2[7], joinpath(@__DIR__,"./Results/CO2_spills"))

# checkCO2[5] # CH4 Emissions, was Y agr
# checkCO2[6] # CO2 emissions, was Y min
# checkCO2[7]
# png(checkCO2[7], "./Results/CO2tax-Allemiss")

# checkCO2[7] # Y pip
# checkCO2[8] # Y oil
# checkCO2[9] # Y wst
# checkCO2[10] # PA agr
# checkCO2[11] # PA min
# checkCO2[12] # PA pip
# checkCO2[13] # PA oil
# checkCO2[14] # PA wst
# checkCO2[15] # ch4 emissions
# checkCO2[16] # co2 emissions
# checkCO2[17] # comp demand Y pip
# checkCO2[18] # comp demand A pip
# checkCO2[19] # Final demand pip
# checkCO2[20] # A pip
# checkCO2[21] # A oil

# # print("checkCO2",checkCO2[3])
# # # checkCO2[4]
# fix(RA,16426.2) # RA value at $190/t
# checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 20, 1200, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# checkch4CO2[7]

# checkch4 = plottaxemisscurve(CH₄_tax,CO₂_tax, 0, 1, 500, zeros(501),round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# resultdfSpill = copy(checkch4[1])
# # # # checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, 400, zeros(401), round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# # # # resultdfSpill = copy(checkCO2[1])
1
# tax1 = names(resultdfSpill)[1]; tax2 = names(resultdfSpill)[2]
# pltSpill = plot(resultdfSpill[1:end,tax1], resultdfSpill[1:end,:Emissions].*10^3,  label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", #title= "Emissions with $(names(resultdfSpill)[1]) $(names(resultdfSpill)[2])",
# xlims=(0,resultdfSpill[end,tax1]), color=:black, linewidth=2, ylab="MMt CO₂eq",legend=:bottom, yguidefontsize=9, guidefont=font(10,"Palatino Roman"), legendfont=font(9,"Palatino Roman"))
# pltSpill = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
# pltSpill = plot!(resultdfSpill[1:end,tax1], resultdfSpill[1:end,:CO2Emissions].*10^3, label="CO₂ Emissions", color=:blue)
# pltSpill = plot!(resultdfSpill[1:end,tax1], resultdfSpill[1:end,:CH4Emissions].*10^3, label="CH₄ Emissions", color=:green)
# ### Next line for CH4 spillovers, only valid for CO2-only policy
# # pltSpill = plot!(twinx(),resultdfSpill[:,tax1],resultdfSpill[:,:CH4perc_red], ylim = (4.96,5.04),label="Ch4 spillover %\nright axis", color=:orange, xlims=(0,resultdfSpill[end,tax1]),legend=:right,guidefont=font(10,"Palatino Roman"),yticks=([4.96,4.98,5,5.02,5.04,],["4.96%","4.98%","5%","5.02%","5.04%"]))
# ### Next line for CO2 spillovers, only valid for CH4-only policy
# pltSpill = plot!(twinx(),resultdfSpill[2:end,tax1],resultdfSpill[2:end,:CO2perc_red], label="CO2 spillover %\nright axis", legendfont=font(9,"Palatino Roman"),color=:orange, xlims=(0,resultdfSpill[end,:CH₄_tax]),ylim=(0,100),yticks=([0,10,20,30,40,50,60,70,80,90,100],["0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"]), legend=:right,guidefont=font(10,"Palatino Roman"))
# 1
# savefig(pltSpill, joinpath(@__DIR__,"./Results/CH4-oilgas-Spillovers.svg"))

# checkch4[2]
# checkch4[5]
# checkch4[6]
# checkch4[7]
# png(checkch4[7], joinpath(@__DIR__,"./Results/CH4_spills"))
# png(checkch4[7], "./Results/CH4tax-Allemiss")

# checkch4[8]
# checkch4[9]
# checkch4[10]
# checkch4[11]
# checkch4[12]
# checkch4[13]
# checkch4[14]
# checkch4[15] #ch3 emissions
# checkch4[16] # Co2 emissions
# print("checkch4",checkch4[3])

# png(checkch4[2], "./Results/CH4to500-Totemiss")
# png(checkch4[15], "./Results/CH4to500-ch4emiss")
# png(checkch4[16], "./Results/CH4to500-co2emiss")
# png(checkch4[8], "./Results/CH4to500-oilactivity")
# png(checkch4[13], "./Results/CH4to500-oilprice")

# print("checkch4",checkch4[3])
# print(check[1])
# png(checkCO2[2], "./Results/CO2to1600RAUnfxd")
# return Emissions_Mt, 
println(
sort(filter(x -> x.var in [Symbol("Y[rec]"), Symbol("PVAM[uel]"),Symbol("PVAM[rnw]"),
# Symbol("Y[oil]"), Symbol("Y[gas]"), Symbol("Y[pet]"), Symbol("Y[coa]"), Symbol("Y[agr]"), 
Symbol("Y[uel]"), Symbol("Y[rnw]"), Symbol("Y[oil]"), Symbol("Y[gas]") 
# Symbol("Y[uwt]"), Symbol("Y[wst]"), Symbol("PVA[compen]"), Symbol("PVA[surplus]"), Symbol("RA")
], Compare))
,
)
println(Emissions_Mt)
println(EmissionReductionResults_Mt)# ;println("RA check, it's the 'both' bc that's the last simulation:",
# sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip]))
println(outerjoin(TaxRev,EqVar, on=:solve))
#Check Benchmark
# fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var]);
# print(sort!(fullvrbnch, [:bmkmarg]))
## End for Multiloop function version of model
# end
