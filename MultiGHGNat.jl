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
# P = WplusAgCSpdata2022 # Re-aggregated/disaggregated from BEA detailed via WiNDC.jl, and my AggWiNDCdata.jl script

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
tm_0DAA = DenseAxisArray([tm_0[i,:value] for i in Ip], Ip)
ta_0DAA = DenseAxisArray([ta_0[i,:value] for i in Ip], Ip)

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
CO2_taxrate = 45.361#<=target for all CH4
CH4_taxrate = 339.565#<=target for all CH4# CO2_taxrate = 45.271#<=re target w oilgas multi
    # CO2_taxrate = 22 #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement
    # CH4_taxrate = 119 #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement
    # CO2_taxrate = 10.39 #~ value for optimum combimation, all CH4, YES GWP, 1x oil/gas, w abatement
    # CH4_taxrate =47.75 #~ value for optimum combimation, all CH4, YES GWP, 1x oil/gas, w abatement
    # CO2_taxrate = 0.2095 #~ value for optimum combimation, all CH4, YES GWP, 1x oil/gas, w abatement
    # CH4_taxrate =9.75 #~ value for optimum combimation, all CH4, YES GWP, 1x oil/gas, w abatement
    # CO2_taxrate = 2 #~ value for optimum combimation, all CH4, 5 x oil/gas, w abatement
    # CH4_taxrate = 47.5 #~ value for optimum combimation, all CH4, 5 x oil/gas, w abatement
    # CO2_taxrate = 8.07#<=target for all CH4 combination with no abatement
    # CH4_taxrate = 188#<=target for all CH4# combination with no abatement
#     CO2_taxrate = 200 * 1.130480652# SC CO2 EPA 2023 SCGHG report, 2022 year, 2020US$, central 2% near-term discount rate x BLS CPI adjustment from 2020$
# CH4_taxrate = 200 * 1.130480652 #<= using SC CO2 because CH4 data is in MtCO2eq # 240.064 w abatement re gross and sink target #  
# CH4_taxrate = 240.062#<=target for all CH4, w/o CH4 abatement        
# CH4_taxrate = 60.765#<=re target w oilgas multi
# CO2_taxrate = 52.1865#<=updated target w 20-yr GWP
# CH4_taxrate = 73.2105#<=updated target w 20-yr GWP 
# CH4_taxrate = 163.925#<=updated target w 20-yr GWP, oil/gas tax 
# CH4_taxrate = 79.4405#<=updated target w 20-yr GWP, no abatement 
# CO2_taxrate = 45.2715#<=updated target x 5 oil/gas
# CH4_taxrate = 49.9205#<=updated target x 5 oil/gas 

# CO2_taxrate = 49.9033#<=updated target w 20-yr GWP & oil/gas x 5, Reduction targ=1845.07
# CH4_taxrate = 9.9224#<=updated target w 20-yr GWP & oil/gas x 5, Reduction targ=1845.07 
# CH4_taxrate = 13.4085#<=updated target w 20-yr GWP & oil/gas x 5 & oil gas taxed only, Reduction targ=1845.07 
# CH4_taxrate = 20.5398#<=updated target w 20-yr GWP & oil/gas x 5 & oil gas taxed only, No abatement, Reduction targ=1845.07 

# CO2_taxrate = 48.54#<=target for CO2 emissions only to meet target
# CH4_taxrate = 491.43 #<= target for CH$ oil/gas/pip
# CH4_taxrate = 318.64#<=re target w oilgas multi w/o abatement 

CH4abatement="yes" # Comment out CH4abatement="no" to allow CH4 abatment
# CH4abatement="no" 
Kmobile="yes" # Allow kapital & Labor to flow between sectors (original WiNDC)
# Kmobile="no" # Fix kapital in sectors, allow Labor to flow between
print("CO2tax: $CO2_taxrate, "); println("CH4tax: $CH4_taxrate")
print("$CH4abatement CH4 Abatement: "); print("$Kmobile mobile Kapital: ")

CH4_tax_switch_on_sectors = DenseAxisArray([0 for i in Ip], Ip)
CH4_tax_switch_on_sectors[CH4sectors] .=1 # Neutral: all sectors have tax (but only CH4sectors have non-zero CH4Intensity)
CH4_tax_switch_on_sectors[[:coa,:agr,:wst]] .=0 ; print("CH4 gas/oil/pip ONLY Yes?,or->")# sectors here will have NO CH4 tax
#########> Comment out below to turn OFF Methane tax for coal, agriculture, and waste => Tax ONLY for oil, gas, and pip
CH4_tax_switch_on_sectors[CH4sectors] .=1; print("NO,") # Neutral: all sectors have tax (but only CH4sectors have non-zero CH4Intensity)

EPACO2eqNonCO2_to_GHGInv_multiplier = 28/25 # The non-CO2 MAC data (emissions and abatement) uses 25, the GHG Inv uses 28
GWP20_multiplier = (81.2)/25 # Emissions and abatement from non-MAC (100-yr, 25 CO2eq) to IPCC AR6 20-year, 81.2 CO2eq.
GWPmulti = GWP20_multiplier / EPACO2eqNonCO2_to_GHGInv_multiplier; print(" GWPmulti yes?, or->")
######> Comment out for 20-year GWP simulation 
GWPmulti = 1; print("NO ") # Comment out for 20-year GWP simulation 
 
Undercount_oilgasch4_multiplier = 5; print(" Undercount CH₄ yes?, or->")#(Plant et al 2022   #Alvarez 2018 said 1.6)
Undercount_oilgasch4_multiplier = 1; print("NO ")# default (comment out for oilgas multi)

multioil_CH4ratio = (sum(only(MAC_CH4_totemiss[:,["agr_livestock", "agr_rice", "COL", "wst_land", "wst_water"]]))+ only(MAC_CH4_totemiss[:,:GAS])*Undercount_oilgasch4_multiplier)/sum(only(MAC_CH4_totemiss[:,2:end]))
ReductTarget = ReductTargetbase - CH4oftarget + GWPmulti * multioil_CH4ratio *CH4oftarget;print(": Reduction targ=",ReductTarget,", ")

#########################
### End scenario, start data
#########################
# % split for pipelines and oil from MAC for 'GAS' by proportion of combined economic output
gas_of_GAS = sum(ys_0[:,:gas])/(sum(ys_0[:,:pip])+sum(ys_0[:,:oil])+sum(ys_0[:,:gas]))
pip_of_GAS = sum(ys_0[:,:pip])/(sum(ys_0[:,:pip])+sum(ys_0[:,:oil])+sum(ys_0[:,:gas]))
oil_of_GAS = sum(ys_0[:,:oil])/(sum(ys_0[:,:pip])+sum(ys_0[:,:oil])+sum(ys_0[:,:gas]))
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
 MAC_CH4_WiNDC_tot = deepcopy(MAC_CH4_WiNDC_tot_MMt); MAC_CH4_WiNDC_tot[:, 2:end] = MAC_CH4_WiNDC_tot_MMt[:,2:end].*10^-3
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
    @production(MultiNat, Y[j], [t= 0, # Data re-allocation of oil and gas fixed this! re .05, # Can't be zero with slack activities. Need values from lit, but in absence, 0.05 is lowest that solves ~1 for Y[gas]&Y[oil] in the benchmark
    s=Elas[j, :SAGE_klem_Y], vae=>s=Elas[j,:SAGE_kle_VAE], sm=>s= Elas[j,:E3_m_ID],
    va=>vae=0, En=>vae=Elas[j,:SAGE_ene], PrimENRG=>En=Elas[j,:SAGE_en], Elec=>En=Elas[j,:SAGE_en],
    oilgas=>PrimENRG=Elas[j,:E3_e_E_El], inElec=>Elec=Elas[:uel,:SAGE_en]### Elec SHOULDN'T CHANGE WITH J!!!! 
    ],begin
    [@output(PY[i],ys_0[i,j], t, taxes = [Tax(RA,ty[j])]) for i∈Ip]...
    [@input(PA[i], id_0[i,j], sm) for i∈ID]... 
     @input(PVAM[j], sum(va_0[VA,j]), va)
        @input(PA[:pet], id_0[:pet,j], PrimENRG) 
          @input(PA[:oil], id_0[:oil,j],   oilgas, taxes=[Tax(RA,CO₂_tax * CO2Int[:oil])]) 
          @input(PA[:gas], id_0[:gas,j]/2, oilgas, taxes=[Tax(RA,CO₂_tax * CO2Int[:gas])])
        @input(PA[:uel], id_0[:uel,j] , Elec)
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

## Final Consumption with elasticity of Demand
@production(MultiNat, FDem, [t=0, s= .9999999, p=>s=.9999999], begin
    @output(PC, sum(pce_0),t) #for i∈Ip]    
    [@input(PA[i], pce_0[i,:pce],s) for i in Ip if i∉[:hos, :pet]]...
    @input(PA[:pet], pce_0[:pet,:pce],p)
    @input(PA[:hos], pce_0[:hos,:pce],p)
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
    @final_demand(PC, sum(pce_0))
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
### Note: Benchmark doesn't solve to appropriate residual at 0 interation because of margins of slack activity. 
### Does balance with interactions or slack vars and production commented out.
# solve!(MultiNat , cumulative_iteration_limit = 0)
solve!(MultiNat)#, output_minor_iterations_frequency=1)
fullvrbnch = generate_report(MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var])
# print(sort(fullvrbnch, :var))#:bmkmarg))#:bnchmrk))#

TaxRev = DataFrame(solve=Symbol[], Total_Revenue=Float64[], Change_in_Rev=Float64[], Income=Float64[])
EqVar  = DataFrame(solve=Symbol[], utility=Float64[], utilCES=Float64[], Mev=Float64[], OldEq_Var=Float64[], OldEV_pc=Float64[], Ut_perc=Float64[])

totrevbnch  = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+ # taxes are negative in production, but positive for revenue
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
income      = totrevbnch + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:bnch totrevbnch totrevbnch-totrevbnch income])

totdem = -value(FDem)*value(compensated_demand(FDem,PC))
# Mevbnch2 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
Mevbnch = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
MevbnchCES = prod([(1/value(PA[i])) for i in Ip if value(PA[i])>0])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
# EVbnch2  = Mevbnch2 - value(RA) 
EVbnch  = Mevbnch - value(RA) 
elasRA = MPSGE.elasticity(MultiNat.productions[FDem].input)
# priceIndCES = (sum([((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) * value(PA[i])^(1-elasRA) for i in Ip]))^(1/(1-elasRA))
PrIndexbnchCES = (sum([pce_0[i,:pce]/sum(pce_0)*value(PA[i])^(1-elasRA) for i in Ip]))^(1/(1-elasRA))
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
utilCES = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
# utilCES = sum([(pce_0[i,:pce]/sum(pce_0))^(1/elasRA)*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
# utilCD    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])
utilCD    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip])

push!(EqVar, [:bnch utilCD utilCES Mevbnch EVbnch EVbnch/value(RA)*100 -(utilCES-utilCES)/utilCES*100] )
sum([value(FDem)*value(compensated_demand(FDem,PA[i]))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])

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

LaspeyresPrIndex = sum([value(PA[i])*pce_0[i,:pce] for i in Ip])/value(RA)
PaaschePrIndex = sum([value(PA[i])*value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])/sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
Wmaybe = (value(RA)/LaspeyresPrIndex - value(RA))/value(RA) 
totdemWiNcntfac = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
# EVWiNcntfac2  = MevWiNcntfac2 - value(RA) 
utilWiNcntfac    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
utilCESWiNcntfac = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
# utilWiNcntfac2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) for i in Ip])
# MevWiNcntfac2 = prod([(1/value(PA[i]))^((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) for i in Ip])* only(filter(x->x.solve==:WiNcntfac,TaxRev)[!,:Income]) #])
MevWiNcntfac    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:WiNcntfac,TaxRev)[!,:Income]) #]
MevWiNcntfacCES = prod([(1/value(PA[i])) for i in Ip if value(PA[i])>0])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
EVWiNcntfac  = MevWiNcntfac - value(RA) 
EVWiNcntfacCES  = MevWiNcntfacCES - value(RA) 
EVWiNcntfacCES/value(RA) * 100
push!(EqVar, [:WiNcntfac utilWiNcntfac utilCESWiNcntfac MevWiNcntfac EVWiNcntfac EVWiNcntfac/value(RA)*100 (utilCESWiNcntfac-utilCES)/utilCES*100])

fullvrWiNcntfact = generate_report(MultiNat)
rename!(fullvrWiNcntfact, :value => :WiNcntfact, :margin => :WiNcntfactmarg)
fullvrWiNcntfact[!,:var] = Symbol.(fullvrWiNcntfact[:,:var])
# print(sort(fullvrWiNcntfact, :var))#:bmkmarg))#:bnchmrk))#

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
# Mevch42 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:ch4,TaxRev)[!,:Income]) #])
Mevch4    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:ch4,TaxRev)[!,:Income]) #]
# EVch42  = Mevch42 - value(RA) 
EVch4  = Mevch4 - value(RA) 
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
utilCESch4 = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
Mevch4CES = prod([(1/value(PA[i])) for i in Ip if value(PA[i])>0])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
EVch4CES  = Mevch4CES - value(RA) 
EVch4CES/value(RA) * 100
PrIndexCESch4 = (sum([pce_0[i,:pce]/sum(pce_0)*value(PA[i])^(1-elasRA) for i in Ip]))^(1/(1-elasRA))
PrIndexCESch4 * utilCESch4 - value(RA)
# utilch42    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemch4) for i in Ip])
utilch4    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])

push!(EqVar, [:ch4 utilch4 utilCESch4 Mevch4 EVch4 EVch4/value(RA)*100 (utilCESch4-utilCES)/utilCES*100])

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
# print(sort(fullvrch4, :ch4tax, by= abs, rev=true))#:ch4tax))#:var):ch4marg

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
# Mevco22 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:co2,TaxRev)[!,:Income]) #])
Mevco2    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:co2,TaxRev)[!,:Income]) #]
# EVco22  = Mevco22 - value(RA) 
EVco2  = Mevco2 - value(RA) 
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
utilCESco2 = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
# utilco22    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemco2) for i in Ip])
utilco2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
push!(EqVar, [:co2 utilco2 utilCESco2 Mevco2 EVco2 EVco2/value(RA)*100 (utilCESco2-utilCES)/utilCES*100])


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
# Mevboth2 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:both,TaxRev)[!,:Income]) #])
Mevboth    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:both,TaxRev)[!,:Income]) #]
# EVboth2  = Mevboth2 - value(RA) 
EVboth  = Mevboth - value(RA) 
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
utilCESboth = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
# utilboth2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemboth) for i in Ip])
utilboth    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
push!(EqVar, [:both utilboth utilCESboth Mevboth EVboth EVboth/value(RA)*100 (utilCESboth-utilCES)/utilCES*100])


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

function plottaxemisscurve(tax1, tax2, start, interval, finish ,vec, cnst=1)
    """# runs a loop increasing each tax by \$1/t and then plotting Total GHG (CO2 & CH4) **incorporated** emissions 
    # Arguments are: which tax to change, other tax to either change simultaneously OR keep at 0, st=initial \$ tax value, fin= final \$ tax value,
    # and final (optional) argument can be set to 0 to remove other tax, by default """
    margemiss = DataFrame(tax1=Float64[], tax2=Float64[], Emissions=Float64[], MevCES=[], EVCES2=Float64[], #EV_pcnt=[], EV_pcntCES=[], 
    # margemiss = DataFrame(tax1=Float64[], tax2=Float64[], Emissions=Float64[], utilityCES=[], Mev=[], MevCES=[], Equiv_Variarion=Float64[], EVCES2=Float64[], #EV_pcnt=[], EV_pcntCES=[], 
        CH4Emissions=Float64[],CO2Emissions=Float64[], CH4perc_red=[], CO2perc_red=[])
    rename!(margemiss,:tax1=>tax1.name) , rename!(margemiss,:tax2=>tax2.name)
    Testvars = DataFrame(taxrt=Float64[], 
    # Yagr=Float64[],Ycoa=Float64[], 
    Ypip=Float64[], Yoil=Float64[],Ygas=Float64[],
    # Ypet=Float64[], Apip=Float64[],Aoil=Float64[],Ywst=Float64[], CompDYPApip=Float64[],CompDApipPApip=Float64[],DemRAPApip=Float64[],
    # PAagr=Float64[],PAcoa=Float64[],PApip=Float64[],PAoil=Float64[],PAwst=Float64[],PAuel=Float64[],
    # compDApipPAoil=Float64[],compdDAoilPAoil=Float64[],
    # VAMagr=Float64[],VAScoa=Float64[],VAMcoa=Float64[],VASpip=Float64[],VAMpip=Float64[], VASwst=Float64[],VAMwst=Float64[],
    # VASagr=Float64[],VAM10agr=Float64[],VAM100agr=Float64[],VAM500agr=Float64[],VAMkagr=Float64[],Apip=[],
    # VAScoa=Float64[],VAM10coa=Float64[],VAM100coa=Float64[],VAM500coa=Float64[],VAMkcoa=Float64[],
    # VASpip=Float64[],VAM10pip=Float64[],
    VAM100pip=Float64[],VAM500pip=Float64[],
    # VAMkpip=Float64[],
    VASoil=Float64[],VAM10oil=Float64[],VAM20oil=Float64[],VAM50oil=Float64[],VAM100oil=Float64[],VAM500oil=Float64[],VAMkoil=Float64[],
    VASgas=Float64[],VAM10gas=Float64[],VAM20gas=Float64[],VAM50gas=Float64[],VAM100gas=Float64[],VAM500gas=Float64[],VAMkgas=Float64[],
    # VASwst=Float64[],VAM10wst=Float64[],VAM100wst=Float64[],VAM500wst=Float64[],VAMkwst=Float64[],
    TotEm=Float64[],CH4TotEm=Float64[],CO2TotEm=Float64[]
    # CH4emoil=Float64[],CH4empip=Float64[],CO2ecoa=Float64[],CO2emoil=Float64[]
    )
    ResultsTroubleshoot = DataFrame(var=[], value=Float64[], margin=Float64[], x1=Float64[]) 
    for (i,j) in zip(start:interval:finish,vec)
        print(i,", ",j,": ")
        set_value!(tax1, i)
        set_value!(tax2, cnst*j)
        solve!(MultiNat, output="no");
        Results = generate_report(MultiNat)
        Results[!,:var] = Symbol.(Results[:,:var]);
       
        push!(Testvars, [i,                
        # value(Y[:agr]), value(Y[:coa]),
        value(Y[:pip]),value(Y[:oil]),value(Y[:gas]),
        # value(Y[:pet]), value(A[:pip]),value(A[:oil]), value(Y[:wst]), value(compensated_demand(MultiNat[:Y][:pip],MultiNat[:PA][:pip])),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:pip])),value(demand(MultiNat[:RA],MultiNat[:PA][:pip])),
        # value(PA[:agr]),value(PA[:coa]),value(PA[:pip]),value(PA[:oil]),value(PA[:wst]),value(PA[:uel]),
        # value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:oil])),value(compensated_demand(MultiNat[:A][:oil],MultiNat[:PA][:oil])),
        # value(VAM[:agr]),value(VAS[:coa]),value(VAM[:coa]),value(VAS[:pip]),value(VAM[:pip]),
        # value(VAS[:wst]),value(VAM[:wst]), value(VAS[:agr]),value(VAM10[:agr]),value(VAM100[:agr]),value(VAM500[:agr]),value(VAM1000[:agr]),value(A[:pip]),
        # value(VAS[:coa]),value(VAM10[:coa]),value(VAM100[:coa]),value(VAM500[:coa]),value(VAM1000[:coa]),
        # value(VAS[:pip]),value(VAM10[:pip]),
        value(VAM100[:pip]),value(VAM500[:pip]),
        # value(VAM1000[:pip]),
        value(VAS[:oil]),value(VAM10[:oil]),value(VAM20[:oil]),value(VAM50[:oil]),value(VAM100[:oil]),value(VAM500[:oil]),value(VAM1000[:oil]),
        value(VAS[:gas]),value(VAM10[:gas]),value(VAM20[:gas]),value(VAM50[:gas]),value(VAM100[:gas]),value(VAM500[:gas]),value(VAM1000[:gas]),# value(VAS[:wst]),value(VAM10[:wst]),value(VAM100[:wst]),value(VAM500[:wst]),value(VAM1000[:wst]),
        value(TotEm),value(CH4TotEm),value(CO2TotEm)
        # value(CH4em[:oil]),value(CH4em[:pip]),value(CO2em[:coa]),value(CO2em[:oil])
        ]  )
        totrevboth  = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
        sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
        sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
        sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[vam.name,c]>1]) for vam in VAMcommodSet]))
        income      = totrevboth + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
                elasRA = MPSGE.elasticity(MultiNat.productions[FDem].input)
        totdem = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
        # util    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip])
        utilCESf    = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
        # utilCESf    = sum([(pce_0[i,:pce]/sum(pce_0))^(1/elasRA)*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
        # Mev = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])*income
        MevCES = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])*income
        # EV  = Mev - value(RA)
        # EVCES  =  value(RA) - MevCES
        EVCES2 = -((utilCESf-utilCES)/utilCES)*100
        ResultsTroubleshoot =vcat(ResultsTroubleshoot, [Results fill(i,length(Results[:,1]))])
        push!(margemiss, [i j only(filter(:var => ==(:TotEm), Results)[:, :value])  MevCES  EVCES2  ;;
               only(filter(:var => ==(:CH4TotEm), Results)[:, :value])    only(filter(:var => ==(:CO2TotEm), Results)[:, :value])    ;;
              100*(TotCH4bnchmk-only(filter(:var => ==(:CH4TotEm), Results)[:, :value]))/(TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:,:value])) ;;
              100*(TotCO2bnchmk-only(filter(:var => ==(:CO2TotEm), Results)[:, :value]))/(TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:,:value]))     ])
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    return margemiss, Testvars, ResultsTroubleshoot#, plch4, plco2, plt3 #,pla, plm, plp, plo, plw, plpa, plpm, plpp, plpo, plpw, 
    # plcdyp, plcdap, plfdrap, plAp, plAo
end

# set_upper_bound(MultiNat[:A][:pip], 8)
# set_upper_bound(MultiNat[:A][:pet], 4) # 30 is fine past $1000/t
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pet]))
# set_upper_bound(MultiNat[:Y][:sle], 12)
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(Y[:sle]))
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pip]))
set_silent(MultiNat)
#
###############################
## This Section for manaully tweaked vectors of carbon taxes that fit with increasing CH4 taxes to meet reduction targets in each setting
###############################
#  co2vec = [45.38,41.3,40.9,40.75,39.52,39.4,39.21,39.05,38.9,38.75,38.6,38.45,37.9,37.75,37.6,37.5,37.35,37.19,36.75,36.65,36.5,36.35,36.2,36.05,35.3,35.17,35,34.9,34.75,34.6,34.45,34.35,34.2,34.05,33.91,33.78,32.73,32.6,32.45,32.35,32.2,32.1,31.95,31.85,31.7,31.57,31.45,31.3,30.55,30.45,30.3,30.2,30.1,29.94,29.85,29.7,29.6,29.45,29.35,29.25,28.8,28.68,28.56,28.44,28.33,28.21,28.09,27.97,27.85,27.74,27.62,27.51,27.39,27.28,27.16,27.05,26.93,26.82,26.7,26.58,26.46,26.35,26.23,26.11,25.99,25.87,25.75,25.64,25.52,25.4,25.29,25.18,25.06,24.95,24.84,24.73,24.61,24.5,24.38,24.27,24.15,24.04,23.92,23.81,23.69,23.58,23.46,23.35,23.23,23.12,23,22.89,22.78,22.66,22.55,22.44,22.33,22.21,22.1,22,20.65,20.55,20.45,20.32,20.25,20.15,20.03,19.9,19.82,19.72,19.6,19.48,19.37,19.27,19.16,19.06,18.95,18.85,18.75,18.64,18.54,18.43,18.34,18.25,18.12,18.02,17.92,17.82,17.73,17.63,17.5,17.4,17.3,17.2,17.1,17,16.9,16.8,16.7,16.6,16.5,16.39,16.28,16.18,16.09,16,15.88,15.79,15.69,15.6,15.5,15.39,15.28,15.18,15.07,14.98,14.9,14.79,14.68,14.59,14.5,14.4,14.3,14.2,14.1,14,13.9,13.8,13.7,13.6,13.5,13.4,13.3,13.2,13.1,13,12.9,12.8,12.73,12.64,12.54,12.42,12.32,12.22,12.13,12.04,11.95,11.86,11.76,11.67,11.58,11.47,11.36,11.27,11.18,11.09,11,10.91,10.82,10.71,10.62,10.53,10.44,10.34,10.25,10.15,10.05,9.95,9.85,9.75,9.66,9.58,9.5,9.4,9.31,9.22,9.12,9.03,8.94,8.84,8.75,8.66,8.56,8.47,8.37,8.28,8.19,8.09,8,7.9,7.81,7.72,7.63,7.54,7.46,7.37,7.28,7.19,7.1,7.01,6.92,6.83,6.74,6.65,6.56,6.46,6.37,6.28,6.19,6.09,6,5.91,5.83,5.74,5.65,5.57,5.48,5.39,5.3,5.22,5.13,5.04,4.94,4.85,4.77,4.68,4.6,4.51,4.42,4.33,4.24,4.15,4.06,3.97,3.89,3.8,3.71,3.62,3.53,3.44,3.35,3.27,3.18,3.1,3.01,2.93,2.85,2.76,2.68,2.6,2.51,2.43,2.34,2.26,2.18,2.09,2.01,1.92,1.84,1.75,1.67,1.58,1.5,1.41,1.32,1.23,1.14,1.05,0.99,0.91,0.83,0.74,0.66,0.57,0.49,0.41,0.33,0.24,0.16,0.08,0
# ] # all ch4, 0:1:400
 # attempt at all ch4 smooth, but not as good actually co2vec = [45.361; 41.325; 40.8; collect(40.4:-((40.4-30.1)/50):30.1) ; collect(30:-((30-23.5)/50):23.5); collect(23:-((23-14.5)/76):14.5) ;collect(14.1:-((14.1-7.4)/76):7.4) ; collect(7:-((7-.36)/76):.36);.26;.175;.1;.03]
# co2vec = [45.361; 43; 42; collect(41.5:-((41.5-26)/150):26) ; collect(26:-((26-11.55)/169):11.55) ; collect(11.5:-((11.5-.2)/163):.2);.18;.15;.1;.01
# ]  # For oil, gas, pip only (up to 491) #491.43 with 0, but not a step
# co2vec = [45.361; 42.01; 41.8; collect(41.5:-((41.5-26)/150):26) ; collect(26:-((26-11.55)/169):11.55) ; collect(11.5:-((11.5-.2)/163):.2);.18;.15;.1;.01
# ] 492, not sure what for
# co2vec = [45.38,45.2,45.1,44.95,44.62,44.5,44.31,44.15,44,43.75,43.6,43.35,43.1,42.95,42.8,42.7,42.55,42.39,42.15,41.85,41.7,41.55,41.4,41.25,41.1,40.97,40.8,40.5,40.35,40.2,40.05,39.75,39.6,39.45,39.31,39.18,39.23,39.1,38.95,38.85,38.7,38.1,37.95,37.85,37.7,37.37,37.25,37.1,36.95,36.85,36.7,36.4,36.6,36.44,36.35,35.7,35.6,35.45,35.35,35.05,35.3,35.18,34.56,34.44,34.33,34.21,33.89,33.77,33.65,33.54,33.42,33.11,32.99,32.88,32.76,32.45,32.33,32.22,32.1,31.78,31.66,31.55,31.43,31.31,31.19,30.97,30.85,30.64,30.52,30.3,30.19,30.08,29.96,29.85,29.64,29.53,29.31,29.2,29.08,28.87,28.75,28.54,28.32,28.21,28.09,27.98,27.76,27.65,27.53,27.32,27.2,26.99,26.88,26.66,26.55,26.34,26.13,26.01,25.9,25.7,25.55,25.45,25.35,25.22,25.05,24.95,24.83,24.6,24.52,24.32,24.1,23.98,23.87,23.67,23.56,23.36,23.35,23.15,23.05,22.84,22.64,22.53,22.34,22.25,22.12,21.82,21.72,21.62,21.43,21.23,21.2,21,20.9,20.7,20.6,20.5,20.3,20.2,20,19.9,19.7,19.59,19.48,19.28,19.19,19,18.88,18.79,18.59,18.5,18.3,18.19,17.98,17.88,17.77,17.68,17.5,17.39,17.28,17.09,17,16.8,16.7,16.5,16.4,16.2,16.1,16,15.9,15.7,15.6,15.5,15.3,15.2,15.1,14.9,14.8,14.7,14.53,14.44,14.24,14.12,14.02,13.82,13.73,13.64,13.45,13.26,13.16,13.07,12.98,12.77,12.66,12.57,12.38,12.29,12.2,12.01,11.92,11.81,11.62,11.53,11.44,11.24,11.15,11.05,10.85,10.75,10.55,10.55,10.36,10.28,10.1,10,9.91,9.72,9.62,9.53,9.34,9.24,9.15,9.06,8.86,8.77,8.67,8.48,8.39,8.29,8.1,8,7.91,7.82,7.63,7.54,7.46,7.27,7.18,7.09,6.9,6.81,6.72,6.53,6.44,6.25,6.26,6.06,5.97,5.88,5.79,5.59,5.5,5.41,5.23,5.14,5.05,4.97,4.78,4.69,4.6,4.42,4.33,4.14,4.04,3.95,3.87,3.78,3.7,3.51,3.42,3.23,3.24,3.05,2.96,2.87,2.79,2.6,2.51,2.32,2.23,2.14,2.05,1.97,1.88,1.8,1.61,1.43,1.45,1.26,1.18,1.1,1.01,0.83,0.74,0.56,0.48,0.39,0.31,0.22,0
# ] # For oil/gas/pip only, No abatement. (up to 319) 
# co2vec = [49.8 49 19 17.7 17 16.5 16 15.3 14.9 14.4 14 13.5 13 12.5 12.1 11.6 9.45 9 8.5 8 7.7 7.2 6.8 6.4 6 5.5 5.2 4.7 4.4 3.8 3.5 3.1 2.8 2.5 2.1 1.7 1.3 .9 .6 .21 0
# ]# 41 steps 0:0.25:10 GWP=81.2, CH4 oil x 5, econ-wide CH4 tax, w Abatement
# co2vec = [49.9,49.1,19,18,17.8,17.4,16.9,16.6,16.1,15.7,15.5,15,14.6,14.3,13.8,13.3,12.9,12.6,12.2,11.8,11.5,11.1,10.25,10,9.5,9,8.7,8.5,8.2,7.8,7.6,7.2,6.8,6.3,6.2,6,5.5,5.2,4.8,4.6,4.3,3.9,3.5,3.3,3,2.8,2.6,2.2,1.8,1.6,1.1,0.8,0.5,0.21,0
# ]# 0:0.25:13.5 GWP=81.2, CH4 oil x 5, oil/gas CH4 tax, w Abatement
# co2vec = [49.9,49.1,48.4,47.4,46.8,45.8,45.1913,44.5826,43.7739,43.1652,42.3565,41.7478,41.1391,40.3304,39.7217,39.013,38.2043,37.4957,36.687,36.0783,35.4696,34.7609,34.1522,33.5435,32.8348,32.2261,31.4174,30.8087,30.2,29.5913,28.9826,28.3739,27.5652,26.9565,26.3478,25.7391,25.1304,24.5217,23.913,23.3043,22.6957,22.087,21.4783,20.8696,20.2609,19.6522,19.0435,18.4348,17.8261,17.2174,16.6087,16.2,15.5913,14.9826,14.3739,13.7652,13.3565,12.7478,12.2391,11.6304,11.0217,10.513,9.9043,9.49565,8.98696,8.57826,8.06957,7.46087,7.05217,6.44348,6.03478,5.42609,5.01739,4.4087,3.8,3.4,3,2.6,2.1,1.6,1.1,0.61,0.1
# ]# 0:0.25:13.5 GWP=81.2, CH4 oil x 5, oil/gas CH4 tax, No Abatement
# co2vec = [45.2715,31.3,30.45,30.13,29.63,29.3,29,28.7,27.47,27.2,26.89,26.5871,26.2943,25.9914,25.6986,25.3957,25.1229,24.8,24.5071,24.2343,23.9314,23.6586,23.3357,23.0629,22.44,22.1671,21.8943,21.6214,21.3286,21.0257,20.7529,20.48,20.2071,19.9043,19.6314,19.3586,18.7957,18.5229,18.24,17.97,17.6982,17.4164,17.1446,16.8929,16.6111,16.3393,16.0675,15.8357,14.5339,14.2821,14.0304,13.7786,13.5268,13.275,13.0232,12.7714,12.5196,12.2679,12.0161,11.7743,11.5125,11.2607,11.0089,10.7771,10.5354,10.2836,10.0318,9.8,9.54,9.304,9.068,8.782,7.006,6.79,6.564,6.348,6.112,5.896,5.67,5.454,5.218,5.012,4.796,4.56,4.354,4.138,3.912,3.716,3.48,3.28,3.058,2.852,2.636,2.42,2.22,2,0.777,0.56,0.36,	0.18,	0
# ] # x 5, & standard

# co2vec = [45.361; 32; 31.5; collect(30.8:-((30.8-20.5)/38):20.5) ; collect(20.1:-((20.1-9.1)/38):9.1) ; collect(8.8:-((8.8-.45)/38):.45);.24;.17;.1;.0
# ] # x 5 oil/gas tax only

# co2vec = [45.2715; 31.9; 31; 30.8; 30; 29.7; 29.3; 28.8; 27.9; 27.35; collect(26.96:-((26.96-17.92)/28):17.92) ; collect(17.65:-((17.65-9.2)/28):9.2) ; collect(8.9:-((8.9-1.75)/25):1.75);1.5;1.3;.747;.36;.26;0.1 ;0
# ]  # Just Set up  GWP, 0:1:75
# co2vec = [52.1865,51.9,43.3,41.1,40.9,40.7,40.5,39.8,39.7,39.4,39.5,39.1,38.8,38.8222,38.6444,38.4667,36.0889,35.5111,35.5333,35.3556,35.1778,35.1,34.9222,34.8144,34.6367,34.4689,34.3011,34.1333,33.9556,33.7778,33.7,33.6222,33.4444,33.2667,33.2889,33.1111,32.9333,32.7556,32.6778,32.5,32.4222,32.2444,32.0667,31.9889,31.9111,31.7333,31.6556,31.4778,30.38,30.2122,30.0444,29.8667,29.69,29.7111,29.5333,29.3556,29.3778,29.2,29.0222,28.9344,28.7667,28.6889,28.5111,28.4333,28.3556,28.1778,28,27.9222,27.7444,27.6667,27.4889,27.4111,26.6233,26.5356,26.3778,26.2,26.1222,26.0444,25.8667,25.7889,25.6111,25.5333,25.3556,25.2778,25.2,25.0222,24.9444,24.7667,24.6889,24.6111,24.5333,24.3556,24.2678,24.2,24.0222,23.8444,22.2467,22.0889,22.0011,21.9133,21.7556,21.6678,21.6,21.4,21.2744,21.1489,21.1133,20.9978,20.8722,20.7467,20.6211,20.4956,20.46,20.3444,20.2189,20.0933,19.9678,19.9322,19.8167,19.6911,19.5656,19.44,19.3144,19.2789,19.1633,18.9378,19.0122,18.8867,18.7611,18.6356,18.51,18.3844,18.2589,18.1333,18.0078,17.8822,17.9567,17.8311,17.7056,17.58,17.4544,17.3289,17.2933,17.1778,14.7522,14.8267,14.7011,14.6756,14.55,14.4244,14.3689,14.2433,14.1178,14.0022,13.8867,13.7711,13.6456,13.52,13.3944,13.2689,13.3433,13.2178,13.0922,12.9667,12.9311,12.8156,12.77,12.6544,12.5389,12.4133,12.3778,12.2622,12.2167,12.1011,11.9756,11.86,11.7344,11.6989,11.5833,11.5378,11.4322,11.3067,11.1811,11.1456,11.03,10.9044,10.8689,10.8333,10.7178,10.5922,10.5567,10.3511,9.4256,8.9,8.69,8.59667,8.49333,8.39,8.28667,8.27333,8.17,8.07667,8.05333,7.95,7.85667,7.75333,7.66,7.55667,7.45333,7.35,7.33667,7.23333,7.14,7.11667,7.01333,6.92,6.82667,6.72333,6.62,6.51667,6.41333,6.4,6.30667,6.27333,6.1,5.99667,5.98333,5.89,5.85667,5.68333,5.58,5.56667,5.47333,5.44,5.26667,5.16333,5.15,5.05667,5.02333,4.25,4.14667,4.04333,4.03,3.93667,3.90333,3.81,3.71667,3.62333,3.52,3.41667,3.31333,3.3,3.20667,3.10333,3,2.98667,2.89333,2.87,2.68667,2.67333,2.58,2.54667,2.37333,2.36,2.26667,2.23333,2.06,2.04667,1.95333,1.92,1.74667,1.73333,1.64,1.55,1.43333,1.33,1.22667,1.32333,1.22,1.11667,1.01333,0.91,0.80667,0.703333,0.77,0.6,0.6,0.4,0.4,0.4,0.2,0.2,0.1,0,0,0,0
# ]  # GWP, 0:.25:75
# 5408.
# co2vec = [52.1865;48;47;collect(46:-(46-2.1)/150:2.1); 2;1.9;1.4;;2;1.4;1;.6;.3;.2;0

# ]
# co2vec = [52.1865,51.97,51.76,51.52,51.3,51.11,50.89,50.68,50.458,50.25,50.047,49.835,49.59,49.3965,49.2031,48.9596,48.7661,48.5227,48.3292,48.1357,47.9223,47.7288,47.4853,47.2919,47.0984,46.855,46.6615,46.458,46.2646,46.0711,45.8776,45.6342,45.4407,45.2472,45.0538,44.8603,44.6268,44.4334,44.2199,43.9964,43.803,43.6095,43.406,43.2126,43.0191,42.8256,42.6322,42.4387,42.2452,42.0318,41.8383,41.6249,41.4314,41.2379,41.0445,40.851,40.6575,40.4641,40.2706,40.0771,39.8837,39.6902,39.4967,39.3033,39.1098,38.9163,38.7229,38.5294,38.3359,38.1425,37.949,37.7755,37.5821,37.3886,37.1951,37.0017,36.8282,36.6348,36.4413,36.2478,36.0544,35.8709,35.6774,35.484,35.3405,35.147,34.9536,34.7801,34.5866,34.4132,34.2197,34.0262,33.8428,33.6493,33.4958,33.3024,33.1089,32.9254,32.782,32.5885,32.415,32.2016,32.0581,31.8647,31.6912,31.4977,31.3343,31.1608,30.9673,30.7739,30.6304,30.4369,30.2735,30.1,29.9,29.7436,29.5871,29.3807,29.2243,29.0678,28.8614,28.705,28.5485,28.3421,28.1856,28.0292,27.8528,27.6663,27.5099,27.3335,27.177,26.9706,26.8042,26.6477,26.4813,26.2949,26.1284,25.972,25.7955,25.6391,25.4827,25.3062,25.1498,24.9734,24.7869,24.6205,24.4941,24.3176,24.1312,23.9748,23.8183,23.6619,23.5054,23.329,23.1426,22.9861,22.8297,22.6733,22.5168,22.3604,22.204,22.0475,21.8711,21.7147,21.5282,21.3718,21.2153,21.0589,20.9025,20.746,20.5896,20.4332,20.2767,20.1203,19.9639,19.8074,19.651,19.4946,19.3381,19.1817,19.0252,18.8688,18.7124,18.5559,18.3995,18.2431,18.0866,17.9302,17.7738,17.6173,17.4909,17.3345,17.198,17.0416,16.8851,16.7287,16.5723,16.4158,16.2594,16.103,15.9465,15.8201,15.6837,15.5272,15.3708,15.2144,15.0579,14.9315,14.795,14.6386,14.4822,14.3257,14.1693,14.0429,13.8864,13.75,13.6,13.4463,13.3127,13.179,12.9953,12.8816,12.728,12.5943,12.4406,12.2969,12.1633,12.0296,11.8759,11.7122,11.5786,11.4449,11.2912,11.1576,10.9939,10.8602,10.7265,10.5929,10.4592,10.3255,10.1718,10.0382,9.8745,9.7408,9.6071,9.47347,9.3398,9.20612,9.07245,8.93878,8.7851,8.62143,8.48776,8.35408,8.22041,8.08673,7.95306,7.81939,7.68571,7.55204,7.41837,7.28469,7.15102,7.01735,6.88367,6.75,6.61633,6.48265,6.34898,6.21531,6.08163,5.94796,5.81429,5.68061,5.54694,5.41327,5.27959,5.14592,5.01224,4.87857,4.7449,4.61122,4.47755,4.37388,4.2602,4.12653,3.99286,3.85918,3.72551,3.59184,3.45816,3.35449,3.22082,3.08714,2.95347,2.8198,2.68612,2.57245,2.43878,2.3051,2.20143,2.08776,1.95408,1.82041,1.68673,1.58306,1.44939,1.33571,1.20204,1.06837,0.93469,0.83102,0.697347,0.583673,0.45,0.35,0.23,0.09,0
# ] # GWP, No Abatement
# co2vec = [45.38,45.1,44.9,44.7,44.4,44.2,43.9,43.7,43.4,43.2,43,42.7,42.5,42.3,42.1,42.1,41.6,41.4125,41.2375,41.05,40.8625,40.675,40.4875,40,39.6889,39.4778,39.3667,39.0056,38.8944,38.5833,38.4722,38.1611,37.95,37.72,37.49,37.36,37.13,36.9,36.67,36.44,36.21,35.98,35.65,35.4786,35.3071,35.1357,34.8643,34.5929,34.4214,34.25,34.0182,33.8364,33.5545,33.3727,33.0909,32.9091,32.7273,32.5455,32.3636,32.0318,31.85,31.6208,31.3917,31.2625,31.0333,30.7542,30.625,30.3958,30.2667,30.0375,29.8083,29.5792,29.35,29.1208,28.9917,28.7625,28.6333,28.3042,28.175,27.9458,27.7167,27.4875,27.3583,27.1292,26.9,26.7708,26.5417,26.3125,26.0833,25.9542,25.725,25.5958,25.3667,25.1375,24.9083,24.7792,24.55,24.4208,24.1917,23.9625,23.7333,23.6042,23.325,23.1958,23.0667,22.8375,22.6083,22.3792,22.35,22.05,21.85,21.65,21.45,21.3,21.1,20.9,20.7,20.55,20.6,20.5,20,19.8,19.6,19.4,19.2,19,18.8,18.6,18.5,18.3111,18.0722,17.9333,17.9944,17.5556,17.4167,17.1778,17.0389,17.1,16.9611,16.7222,16.5833,16.1444,16.2056,15.7667,15.5278,15.3889,15.25,15.0625,14.875,14.6875,14.6,14.4125,14.225,14.0375,13.85,13.6625,13.475,13.2875,13.15,12.9625,12.775,12.5875,12.4,12.2722,12.1444,12.0167,11.6889,11.5611,11.3833,11.2556,11.1278,10.9,10.6722,10.5444,10.4167,10.1889,10.0611,9.8333,9.80556,9.57778,9.35,9.22059,8.99118,8.86176,8.73235,8.50294,8.37353,8.24412,8.01471,7.88529,7.70588,7.57647,7.44706,7.21765,7.08824,6.95882,6.72941,6.6,6.38095,6.3619,6.14286,5.92381,5.80476,5.63571,5.51667,5.39762,5.17857,4.95952,4.94048,4.72143,4.50238,4.38333,4.26429,4.04524,3.92619,3.70714,3.5881,3.41905,3.4,3.1,3.09,2.9,2.69,2.55,2.40143,2.25,2.09,1.96714,1.80571,1.65,1.48286,1.36,1.21,1.05,0.9,0.765714,0.604286,0.47,0.32,0.16,.0,0
# ]  # Standard but No CH4 abatement, 0:1:241


# checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 1, 1600, 
# zeros(1601)) 
# # # collect(0:1:400))
# co2vec)

# checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 10, 480, zeros(51), 0)
# checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 600, 1, 800, ze ros(201), 0)

# # # # # # # # EVdf2 = checkch4CO2[1]
# # # # # # EVdf_slice2= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget) && x.Emissions >(TotGHGbnchmk*10^3-ReductTarget-1),EVdf2[:,[1,2,3,11]])# 4329.824705730001
# # # # # # print(EVdf_slice2)
# # # # # # # ##### checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 10, 400, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# # # # checkch4CO2[7]
# # # # # print(EVdf2)
# # # # checkch4CO2[1]
# # # # ##### png(checkch4CO2[7], "./Results/Bothtax-Allemiss")
# # # ##### checkch4CO2[2]
# # # # EVdf = checkch4CO2[1]
# # # # ##### print("checkch4CO2",checkch4CO2[3])
# # # # EVdf_slice= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget) && x.Emissions >(TotGHGbnchmk*10^3-ReductTarget-1),EVdf)# 4329.824705730001
# # # # EVdf_cap= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget),EVdf)

# resultdf = copy(checkCO2[1])#[1:1300,:]
# # resultdf = copy(checkch4CO2[1])
# tax1 = names(resultdf)[1]; tax2 = names(resultdf)[2]

# # plt = plot(resultdf[!,tax1], resultdf[!,:Emissions].*10^3,  label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) & $(replace(tax2,"_"=>" "))  \$/t", xlims=(0,resultdf[end,:CH₄_tax]))#title= "RA:\$$(value(RA)) fxd:$isfixed",)
# # plt = plot!([CO2_taxrate], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3))
# # plch4 = plot(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, label=false, title="CH₄ Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(names(resultdf)[1]) $(names(resultdf)[2]) \$/t", xlims=(0,resultdf[end,:CH₄_tax]))
# # plch4 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
# # plco2 = plot(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label=false, title="CO₂ Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(names(resultdf)[1]) $(names(resultdf)[2])  \$/t", xlims=(0,resultdf[end,:CH₄_tax]))
# # plco2 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
# pltEV = plot(resultdf[!,tax1], resultdf[!,:Emissions].*10^3, legend=:left, label="Total GHG Emissions", ylim=(0,maximum(resultdf[:,:Emissions])*10^3+500), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", #title= "Emissions with $(names(resultdf)[1]) $(names(resultdf)[2])",
# xlims=(0,resultdf[end,tax1]), color=:black, linewidth=1, ylab="MMt CO₂eq", linestyle=:dashdot,
# yguidefontsize=9, legendfont=font(9,"Palatino Roman"),guidefont=font(10,"Palatino Roman"))
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, label="CH₄ Emissions", color=:darkgreen, linestyle=:dash)
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label="CO₂ Emissions", color=:blue, linestyle=:dot)
# pltEV = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(0,maximum(resultdf[:,:Emissions])*10^3+500), color=:red, linewidth=0.2)
# # pltEV = plot!(resultdf[!,tax1],repeat([TotGHGbnchmk*10^3-ReductTarget],length(resultdf[:,tax1])), color=:yellow, label="Reduction target",linewidth=1.5)
# # pltEV = plot!(twinx(),resultdf[!,tax1],resultdf[!,:EVCES2], legend=:right, label="Excess Burden:\nRight axis", ylabel= "Percentage change", guidefont=font(9,"Palatino Roman"), linewidth=2., xlim=(0,resultdf[end,tax1]), ylim=(0,1), yticks=([0.0,0.2,0.4,0.6,0.8,1.0],["0.0%","0.2%","0.4%","0.6%","0.8%","1.0%"]), legendfont=font(9,"Palatino Roman"))
# # pltEV = plot!(twiny(),resultdf[!,tax2],resultdf[!,:Emissions].*10^3,  xflip=true, xlim=(resultdf[end,tax2],resultdf[1,tax2]),xlabel="$(replace(tax2,"_"=>" ")) \$/t", linewidth=0, legend=false, guidefont=font(10,"Palatino Roman"),xticks=([0,45],["0","\$45"]))
# #Not working
# # pltEV = Plots.scatter!(twiny(),[resultdf[120,tax2]],[resultdf[120,:CO2Emissions].*10^3],  xflip=true, xlim=(resultdf[end,tax2],resultdf[1,tax2]),xticks=false,legend=false, guidefont=font(10,"Palatino Roman"))

# # # # ### Next line for CH4 spillovers, only valid for CO2-only policy
# pltEV = plot!(twinx(),resultdf[2:end,tax1],resultdf[2:end,:CH4perc_red],legend=:right, label="CH₄ spillover %\nright axis", color=:orange, linewidth=2,legendfont=font(9,"Palatino Roman"), 
# ylim=(minimum(resultdf[2:end,:CH4perc_red])-.1,maximum(resultdf[2:end,:CH4perc_red])+.1), yticks=([4.5,4.6,4.7,4.8,4.9],["4.5%","4.6%","4.7%","4.8%","4.9%"]))
# pltEV = plot!(twinx(),resultdf[2:end,tax1],(resultdf[1,:CH4Emissions] .-resultdf[2:end,:CH4Emissions])./(resultdf[1,:Emissions] .-resultdf[2:end,:Emissions]), label="CO2 spillover %\nright axis", color=:orange, legend=:left)#, ylim=(0,1),xticks=([0,.10,.20,.30,.40,.50,.60,.70],["0%","10%","20%","30%","40%","50%","60%","70%"])
# ## Next line for CO2 spillovers, only valid for CH4-only policy
# pltEV = plot!(twinx(),resultdf[2:end,tax1],(resultdf[1,:CO2Emissions] .-resultdf[2:end,:CO2Emissions])./(resultdf[1,:Emissions] .-resultdf[2:end,:Emissions]), legend=:topright,label="CO2 spillover %\nright axis", color=:orange, ylim=(0,1),xticks=([0,.10,.20,.30,.40,.50,.60,.70],["0%","10%","20%","30%","40%","50%","60%","70%"]))
# pltEV = plot!(twinx(),resultdf[2:end,tax1],resultdf[2:end,:CO2perc_red],legend=:right, label="CO₂ spillover %\nright axis", legendfont=font(9,"Palatino Roman"), color=:orange, xlim=(0,maximum(resultdf[:,tax1])), ylim=(0,100),yticks=([0,10,20,30,40,50,60,70,80,90,100],["0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"]))

# 1
# # TestVars = copy(checkch4CO2[2])
# # TestVars = copy(checkCO2[2])
# print(TestVars)
# # # filter(:EV_pcnt2 => ==(minimum(EVdf_slice.EV_pcnt2)),EVdf_slice)
# # # print(sort(EVdf,:EV_pcnt2)[1:80,:])
# # # print(sort(EVdf_cap,:EV_pcnt2)[1:40,:])

# png(pltEV, joinpath(@__DIR__,"./Results/EVTarget"))
# savefig(pltEV, joinpath(@__DIR__,"./Results/EVTarget.svg"))
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

# checkCO2[7]
1
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
# 1
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

##############################################################
## OTHER plots from function
# plt = plot(margemiss[!,tax1.name], margemiss[!,:Emissions].*10^3, title= "RA:\$$RAval fxd:$isfixed", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
# plt = plot!([CO2_taxrate], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3))
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
##############################################

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
