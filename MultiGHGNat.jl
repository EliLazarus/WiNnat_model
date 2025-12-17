# Adapted from WiNDC National Model
using MPSGE
using DataFrames, JLD2
using JuMP
using MPSGE.JuMP.Containers
using CSV

### Run the economic data aggregation/disaggregation script if it hasn't been
if !@isdefined(WplusSpAgCdata2022)
    include("./data/AggWiNDCdata.jl")
end

# WiNDC v0.1.1 `https://github.com/uw-windc/WiNDC.jl.git#aggregation`   ...C:\Users\Eli\.julia\packages\WiNDC\htB5C
P = WplusSpAgCdata2022 # Re-aggregated/disaggregated from BEA detailed via WiNDC.jl, and my AggWiNDCdata.jl script

S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S
Sectors = CSV.read(joinpath(@__DIR__,"./Sectorsplus.csv"), DataFrame);
## New sectors:
# uel (utility: Electric power generation, transmission, and distribution), ugs (utility: Natural gas distribution), uwt (utility: Water, sewage and other systems), 
# coa (coal mining), min (other non-oil/gas mining [not coal])
### oil => oil & gas emissions and abatement for gas split btw 'gas' and 'pip'?

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
tx_0DAA = DenseAxisArray([0 for i in Ip], Ip)
## Base Marginal Abatement Cost EPA data (2022). $/t and MMtCO2eq
MAC_CH4_data=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2022_data.csv"), DataFrame, header=2, limit=14)
MAC_CH4_totemiss=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2022_data.csv"), DataFrame, header=2, skipto=17)

#####################
# SET SCENARIO
#####################
yr = Symbol(2022)

ReductTargetbase = 1117.47074039339 # 2022 imputed from Paris (= linear 2005 to 2005*(1-0.61) in 2035, gross emissions reduction assuming 2022 sink) diff to actual 2022
CH4toGHG2005 = .113943621 #fraction of CH4 to gross GHG in 2005, GHG Inv CO2eq=28
CH4oftarget = CH4toGHG2005*ReductTargetbase# % of CH4 in 2005/2035 extrapolation to adjust for different CH4 emissions scenarios

SCGHG = 200 * 1.130480652

## Set tax rates 
# Paris Target reduction from 2022 = ReductTarget t  (= linear 2005 to 2005*(1-0.61) in 2035, gross emissions reduction assuming 2022 sink) diff to actual 2022
CO2_taxrate = 48.198#<=target for all CH4 (Nested CES consumption utility, elasticities from Marc)
CH4_taxrate = 292.954#<=target for all CH (Nested CES consumption utility, elasticities from Marc)
### Paris Target lowest cost combination (Nested CES consumption utility, elasticities from Marc)
CO2_taxrate = 21.1719    #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement (post Consumption utility elasticities from Marc)
CH4_taxrate = 119 #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement (post Consumption utility elasticities from Marc)
                #### Redundant, temp for reference
                #### CO2_taxrate = 29.2774    #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement (from Cobb-Douglas Utilities)
                #### CH4_taxrate = 60 #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement (from Cobb-Douglas Utilities)
                #### CO2_taxrate = 42.68#<=target for all CH4
                #### CH4_taxrate = 282.885#<=target for all CH
                #### CO2_taxrate = 5.6134#<= Paris target with SCCH4 for all CH4
                #### CO2_taxrate = 45.271#<=re target w oilgas multi
### Optimal combinations for Paris
    # CO2_taxrate = 18.44176473    #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement
    # CH4_taxrate = 119 #~ value for optimum combimation, all CH4, no GWP, 1x oil/gas, w abatement
    # CO2_taxrate = 10.39 #~ value for optimum combimation, all CH4, YES GWP, 1x oil/gas, w abatement
    # CH4_taxrate = 39.92 #~ value for optimum combimation, all CH4, YES GWP, 1x oil/gas, w abatement 362.55
    # CO2_taxrate = 0.2095 #~ value for optimum combimation, all CH4, YES GWP, 5x oil/gas, w abatement
    # CH4_taxrate =9.75 #~ value for optimum combimation, all CH4, YES GWP, 5x oil/gas, w abatement
    # CO2_taxrate = 2 #~ value for optimum combimation, all CH4, 5 x oil/gas, w abatement
    # CH4_taxrate = 38.271 #~ value for optimum combimation, all CH4, 5 x oil/gas, w abatement # 42.3 what is this?
    # CO2_taxrate = 35.19435  + 200 * 1.130480652 # CO2 tax to match CH4 reductions at SCC
    # CH4_taxrate = .08 ## Test combo to reach 3053.83 SCCO2 reductions with combo
    # CO2_taxrate = .19393 ## Test combo to reach 3053.83 SCCO2 reductions with combo

    ### SCCCO2
### Version ΔCH4/ΔCO2 x \Delta SCCO2
# Optimum CO2 tax ($237/t) incorporating the additional cost of the CO2 associated with that amount of CH4
# CO2_taxrate =  SCGHG  #* 2# SC CO2 EPA 2023 SCGHG report, 2022 year, 2020US$, central 2% near-term discount rate x BLS CPI adjustment from 2020$
# CO2_taxrate = SCGHG + 127.11002457204812/2719.3414026721402 * SCGHG # -# Δ: 10.56839890075912 = 236.6645293007591-(200*1.130480652) [0.046742944614138617]
#### (127.39569955056795-127.11002457204812)/(2724.65463406823-2719.3414026721402) Ratio at the margin (SCCO2 + 1) -->  0.05376671129551907
# ((EmissionReductionResults_Mt[1,:CO2tax_reduc]/EmissionReductionResults_Mt[1,:Sum_of_each_tax])*
        # (EmissionReductionResults_Mt[1,:Interactions]/EmissionReductionResults_Mt[3,:Sum_of_each_tax])
        # +
        # (EmissionReductionResults_Mt[2,:CO2tax_reduc]/EmissionReductionResults_Mt[2,:Sum_of_each_tax])*
        # (EmissionReductionResults_Mt[2,:Interactions]/EmissionReductionResults_Mt[3,:Sum_of_each_tax]))*     
        # SCGHG
# CO2_taxrate = 207.52891160027443
        # CO2_taxrate = 200 * 1.130480652 + .1 # CH4/CO2: [0.04675664135837693] 127.39569955056795/2724.65463406823 * 200 * 1.130480652
            # CO2_taxrate =  200 * 1.130480652 - 1 # CH4/CO2: [0.04672922057023451] (126.82327535984061/2714.003653650159* 200 * 1.130480652 +127.39569955056795/2724.65463406823 * 200 * 1.130480652)/2
            # CO2_taxrate = 238.25257576825084 # (127.39569955056795-127.11002457204812)/(2724.65463406823-2719.3414026721402) * 200 * 1.130480652 + 200 * 1.130480652
            ## Oltimal combination of taxes:: marginal benefit (tons reduced x SCC) = marginal cost (Δ MEV)
            # CH4_taxrate = 270 #(from 269.9) with CO2_Tax 21.835  SCCO2 & SCCH4 with CO2 EV_ton, both taxes =    -458.73443613181297
            # CO2_taxrate = 17.776 #(from 17.676) * 1.130480652  #   -102.65914688792492

            ## Optimal price, the cost (with CH4 in CO2 eq = reductions) of the direct CH4 emissions reductions from the SCC, but incorporating the CO2 reductions
            ## So it's the price that reduces total emissions by the amount of CH4 emission reductions (301.27768273240076) from the price at the SCC
            # CO2_taxrate = 204.560589 # (2719.3414026721402 MMt)
### SCCH4
    # CH4_taxrate = SCGHG # 200 * 1.130480652 #* 2# SC CH4 EPA 2023 SCGHG report, 2022 year, 2020US$, central 2% near-term discount rate x BLS CPI adjustment from 2020$
### Optimal CH4 tax:: Version ΔCO2/ΔCH4 x \Delta SCCO2
    # CH4_taxrate = SCGHG
    # CH4_taxrate = CH4_taxrate + 637.9479632557219/301.2776827323572 * SCGHG # -# 704.8490356286785 (478.7529052286785 SCC of the CO2 portion [2.1174750066783035 x SCC])

    ###    (640.4027205455742-637.9479632557219)/(301.5650649871174-301.2776827323572) Ratio at the margin (SCHGH + 1) --> 8.541784501971808
    #     (EmissionReductionResults_Mt[1,:Interactions]/EmissionReductionResults_Mt[3,:Sum_of_each_tax])
    #     +
    #     (EmissionReductionResults_Mt[2,:CH4tax_reduc]/EmissionReductionResults_Mt[2,:Sum_of_each_tax])*
    #     (EmissionReductionResults_Mt[2,:Interactions]/EmissionReductionResults_Mt[3,:Sum_of_each_tax]))*     
    #     SCGHG
    #     CH4_taxrate = 685.2873207942273
### Quantity of GHG reductions from CO2 tax to adjust from overlap
# CO2_taxrate = 142.78404 # GHG reductions -> 2267.5078529824714 = 2904.364030775726 (total GHG from SCCO2*) - 615.7451507979843  - 21.111026995270763
# Adjusted CO2 
#    EmissionReductionResults_Mt[3,:CO2tax_reduc]-
#                 (EmissionReductionResults_Mt[1,:CO2tax_reduc]/EmissionReductionResults_Mt[1,:Sum_of_each_tax]*EmissionReductionResults_Mt[1,:Interactions]+
#                 EmissionReductionResults_Mt[2,:CO2tax_reduc]/EmissionReductionResults_Mt[2,:Sum_of_each_tax]*EmissionReductionResults_Mt[2,:Interactions])
### Quantity of GHG reductions from CH4 tax to adjust from overlap
# CH4_taxrate = 499.365 # GHG reductions -> 1610.13502093482 = 2037.7215944758027 - 360.4366412345451 -67.1499323064377

# (EmissionReductionResults_Mt[2,:CH4tax_reduc]/EmissionReductionResults_Mt[2,:Sum_of_each_tax])*EmissionReductionResults_Mt[2,:Interactions]
# EmissionReductionResults_Mt[3,:CH4tax_reduc]-
#                 (EmissionReductionResults_Mt[1,:CH4tax_reduc]/EmissionReductionResults_Mt[1,:Sum_of_each_tax]*EmissionReductionResults_Mt[1,:Interactions]+
#                 EmissionReductionResults_Mt[2,:CH4tax_reduc]/EmissionReductionResults_Mt[2,:Sum_of_each_tax]*EmissionReductionResults_Mt[2,:Interactions])


                    #     CH4_taxrate = 200 * 1.130480652  #CO2/CH4: [2.1235971765261725] 
                # CH4_taxrate = 200 * 1.130480652 +.1 #CO2/CH4: [2.1113381073050443]
                # CH4_taxrate = 200 * 1.130480652 #CO2/CH4: [2.1113381073050443]
                # 2.1174750066783035 @ SCCH4
                # 2.117467641915608 = (2.1235971765261725+2.1113381073050443)/2 * 200 * 1.130480652 +200 * 1.130480652
                # 478.7529050900957 +200 * 1.130480652
###############################################
# Marginal cost = marginal benefit, Each individual tax on its own
# CO2_taxrate = 31.67; CH4_taxrate = 0 # (from 31.57)
# CH4_taxrate = 351.48 #;CO2_taxrate = 0 #(from 351.38)
###############################################
### Work to find point where marginal benefit (Δ $/t taxes) = marginal cost (ΔEV/t)
# EV_ton, both taxes =    -461.3295040164493 SCCO2 & SCCH4
# CH4_taxrate = 200 * 1.130480652 +1 
# EV_ton, both taxes  SCCO2 & SCCH4+1 461.3655192985324-461.3295040164493 (0.03601528208309901)
# CH4_taxrate =478.75290523009437 # SCCO2 & SCCH4 with CO2 -475.9501996593062
# CH4_taxrate =478.75290523009437+1 # SCCO2 & SCCH4+1 with CO2 476.022838515209-475.9501996593062 (0.0726388559028237)

        #### 600, 44.65 balance benefits and costs to 4 digits   601 & 44.5555 : CH4 up 1, CO2 down .0944999   44.65-44.5555 OR CH4 down 1, CO2 up 
        #### 800, 21.985 balance benefits and costs to 4 digits 801, 21.835: CH4 up 1, CO2 down 0.14999 
# CH4_taxrate = 267.435 #(from 267.43) with CO2_Tax 18.4052  SCCO2 & SCCH4 with CO2 EV_ton, both taxes =    -458.73443613181297
# CO2_taxrate = 18.4052 #(from 18.4051) * 1.130480652  #   -102.65914688792492
## 267.44 from 267.43      18.4052 from 18.4051

    # CO2_taxrate = 18.50515 &
    # CO2 tax: 18.40515 -> 1.1434253651706994e9 CO2, CO2 tax: 18.50515 -> 1.1453564678965917e9 CO2
    #  1.1434253651706994e9 - 1.1453564678965917e9
    # CO2 tax: 18.40515 ->  3.215047073057263e8 CH4,   CO2 tax: 18.50515 -> 3.2154707660363823e8 CH4 from the optimum, increase just the CH4 by $0.1/t
    #  3.215047073057263e8 - 3.2154707660363823e8
    ## Optimal price, the cost (with CH4 in CO2 eq = reductions) of the direct CH4 emissions reductions from the SCC, but incorporating the CO2 reductions
                            ## So it's the price that reduces total emissions by the amount of CH4 emission reductions (301.27768273240076) from the price at the SCC
                        ### Right UNDER the line of a hyper sensitive jump from 298.864 to 309.826 
                            # CH4_taxrate =  35.813478 # (298.864 ≈ 301.27768273240076 Total GHG reductions)
                ### Starting points where each individual tax acheives 3020.619085404541 MMt reductions
                # CO2_taxrate = 259.486; CH4_taxrate = 1373.61135 # (# 3020.619085404541 total GHG reductions, min EV: CO2 from SCCO2 + CH4 from SCCH4)
                ### Optimal (lowest EV) combination that achieves 3020.619085404541 total GHG reductions (min EV: CO2 from SCCO2 + CH4 from SCCH4)
                # CO2_taxrate = 33.292; CH4_taxrate = 1030
 
                ### Optimum CH4 tax ($216/t) subtract the spilllover benefit of the CH4 reduction associated with that tax rate.
                ####### CH4iter=2
                ####### CH4_taxrate = CH4price[CH4iter-1,:CH4tax] - 200 * 1.130480652 * 127.11002457204812 / (2719.3414026721402)  # 10.568398900759114-> 215.52773149924087
                ####### CH4iter=3
                ####### CH4_taxrate = CH4price[CH4iter-1,:CH4tax] - (127.11002457204812 - 124.02382054981165) / 2661.6547040891874 * 200 * 1.130480652 # 0.2621597707548472-> 215.26557172848604
                ####### CH4iter=4
                ####### CH4_taxrate = CH4price[CH4iter-1,:CH4tax] - (124.02382054981165 - 123.94566491601378) / 2660.1870457751643 * 200 * 1.130480652 # 0.006642648079472584-> 215.25892908040657
                ####### CH4iter=5
                ####### CH4_taxrate = CH4_taxrate - (123.94566491601378 - 123.94368356521501) / 2660.1498343131825 * 200 * 1.130480652 # 0.00016840244966237751-> 215.2587606779569

        ### Optimum CH4 tax ($468/t) incorporating the additional cost of the CO2 associated with that amount of CH4
        ### Version over total difference and x SCC
            # CH4_taxrate = CH4_taxrate + (993.739-637.948)/1334.35 * 200 * 1.130480652 # 153.570678618798 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1122.603-993.739)/1476.464 * 200 * 1.130480652 # 19.73339800216301 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1163.562-1122.603)/1521.543 * 200 * 1.130480652 # 6.086368512131147 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1176.077-1163.562)/1535.307 * 200 * 1.130480652 # 1.8430145058649652 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1179.856-1176.077)/1539.463 * 200 * 1.130480652 # 0.5550099461835712 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1180.993-1179.856)/1540.713 * 200 * 1.130480652 # 0.16685216537069994 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.335-1180.993)/1541.089 * 200 * 1.130480652 # 0.05017547759851781 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.437-1181.335)/1541.202 * 200 * 1.130480652 # 0.014963518929231045 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.468-1181.437)/1541.235 * 200 * 1.130480652 # 0.004547638771790086 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.477-1181.468)/1541.246 * 200 * 1.130480652 # 0.0013202728010994287 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.48-1181.477)/1541.249 * 200 * 1.130480652 # 0.00044009007706344696 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.481-1181.48)/1541.25 * 200 * 1.130480652 # 0.00014669659717414665 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (1181.481-1181.481)/1541.25 * 200 * 1.130480652 # 0.00014669659717414665 SCC of the CO2 portion

            ### Combined optimum of both (CO2: 234.9, CH4: 277.6)
            # CO2_taxrate =  200 * 1.130480652 #* 2# SC CO2 EPA 2023 SCGHG report, 2022 year, 2020US$, central 2% near-term discount rate x BLS CPI adjustment from 2020$
            # CH4_taxrate = 200 * 1.130480652# -2 #* 2#<= using SC CO2 because CH4 data is in MtCO2eq #
            # CO2_taxrate = CO2_taxrate + 127.11/3330.8 * 200 * 1.130480652 # 8.628281234281253
            # CH4_taxrate = CH4_taxrate + 637.948/3330.8 * 200 * 1.130480652 # 43.30418343833889 SCC of the CO2 portion
            # CO2_taxrate = CO2_taxrate + (129.54-127.11)/3423.223 * 200 * 1.130480652 # 0.16049599949287507
            # CH4_taxrate = CH4_taxrate + (742.524-637.948)/3423.223 * 200 * 1.130480652 # 6.907008083525497 SCC of the CO2 portion
            # CO2_taxrate = CO2_taxrate + (129.585-129.54)/3432.093 * 200 * 1.130480652 # 0.0029644668335046856
            # CH4_taxrate = CH4_taxrate + (758.884-742.524)/3432.093 * 200 * 1.130480652 # 1.0777483865804345 SCC of the CO2 portion
            # CO2_taxrate = CO2_taxrate + (129.586-129.585)/3433.375 * 200 * 1.130480652 # 6.585244268426245e-5
            # CH4_taxrate = CH4_taxrate + (761.429-758.884)/3433.375 * 200 * 1.130480652 # 0.16759446663064498 SCC of the CO2 portion
            # CO2_taxrate = CO2_taxrate + (129.586-129.586)/3433.572 * 200 * 1.130480652 # 6.585244268426245e-5
            # CH4_taxrate = CH4_taxrate + (761.824-761.429)/3433.572 * 200 * 1.130480652 # 0.16759446663064498 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (761.886-761.824)/3433.603 * 200 * 1.130480652 # 0.004082580334652163 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (761.895-761.886)/3433.608 * 200 * 1.130480652 # 0.0005926317662363584 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (761.897-761.895)/3433.608 * 200 * 1.130480652 # 0.000131695948056683 SCC of the CO2 portion
            # CH4_taxrate = CH4_taxrate + (761.897-761.897)/3433.609 * 200 * 1.130480652 # 0.0 SCC of the CO2 portion

#######################################
# Average Costs = SCGHG = 200 * 1.130480652
### Costs of rates that equal average benefit   
# CO2_taxrate = 72.62654854 # marginal tax rate = EV/t (cost): $226.0961/t 
# CH4_taxrate = 880.1551194 # tax rate (benefit) = EV/t (cost): $226.0961304/t 
### Combined rates that equal average benefit
# CH4_taxrate = 733;CO2_taxrate = 30.7826009;# 2579.01 # tax rate (benefit) = EV/t (cost): $226.0961304/t 
    # CH4_taxrate = 735;CO2_taxrate = 30.54468;# 2579.00 # other combinations that also ~= EV/t (cost): $226.0961304/t
    # CH4_taxrate = 731;CO2_taxrate = 31.019373;# 2579.0 # other combinations that also ~= EV/t (cost): $226.0961304/t

        # CH4_taxrate = 1356.16 # target = 
            # CH4_taxrate = 1356.16 + 200 * 1.130480652# Rate that gets same GHG reduction as CO2 at SCC 
### Alt specs
# CH4_taxrate = 60.765#<=re target w 5 x oil/gas, oil/gas tax
# CO2_taxrate = 49.8915#<=updated target w 20-yr GWP
# # CH4_taxrate = 62.539#<=updated target w 20-yr GWP 
# CH4_taxrate = 48.4625#42.3005#<=updated target x 5 oil/gas , oil/gas tax
# CO2_taxrate = 50.479#<=updated target w 20-yr GWP x 5 oil/gas, Reduction targ=1845.07
# CH4_taxrate = 8.219#<=updated target w 20-yr GWP x 5 oil/gas, Reduction targ=1845.07 
# CH4_taxrate = 10.9021#<=updated target w 20-yr GWP & x 5 oil/gas & oil gas taxed only, Reduction targ=1845.07 
# CO2_taxrate = 48.54#<=target for CO2 emissions only to meet target
# CH4_taxrate = 318.64#<=re target w oilgas multi w/o abatement 

CH4abatement="yes" # Comment out CH4abatement="no" to allow CH4 abatment
# CH4abatement="no" 
# CO2_taxrate = 48.198#<=target for all CH4 combination with no abatement (post Consumption utility elasticities from Marc)
# CH4_taxrate = 209.325#<=target for all CH4# combination with no abatement (post Consumption utility elasticities from Marc)
    ###Optimal
# CO2_taxrate = 5.532#<=target for all CH4 (post Consumption utility elasticities from Marc)
# CH4_taxrate = 180#<=target for all CH (post Consumption utility elasticities from Marc)
    	

# CO2_taxrate = 8.07#<=target for all CH4 combination with no abatement
# CH4_taxrate = 202.34#<=target for all CH4# combination with no abatement
# CH4_taxrate = 79.4405#<=updated target w 20-yr GWP, no abatement 

Kmobile="yes" # Allow kapital & Labor to flow between sectors (original WiNDC)
# Kmobile="no" # Fix kapital in sectors, allow Labor to flow between

only_oilgaspip = true
only_oilgaspip = false ####> Comment out to implement tax ONLY on oil, gas, pip =>turn **OFF** Methane tax for coal, agriculture, and waste
CH4_tax_switch_on_sectors = DenseAxisArray([0 for i in Ip], Ip)
if !only_oilgaspip;  ; print("[All CH4 tax] "); end
CH4_tax_switch_on_sectors[CH4sectors] .=1 # All sectors switched on by default, but only CH4sectors in the loop for the taxes and have non-zero CH4Intensity
if only_oilgaspip; CH4_tax_switch_on_sectors[[:coa,:agr,:wst]] .=0 ; print("[CH4 gas/oil/pip ONLY] "); end# sectors here will have NO CH4 tax

# CH4_taxrate = 406.364 #<= target for CH$ oil/gas/pip (was $393)

### If I was to represent methane in non CO2 eq units:
### Multiply emissions and the emission part of abatement cost by 25
### Then the SCCH4 rate would be "EPA 2022 SCCH4/SCCO2" (1799/200) = 8.995 higher. 
### That translate to approx 0.3598. 1799 * 1.130480652
GWP20year = true
GWP20year = false###> Comment out to implement 20-year GWP simulation 
EPACO2eqNonCO2_to_GHGInv_multiplier = 28/25 # The non-CO2 MAC data (emissions and abatement) uses 25, as a step to convert to 20yr from GHG Inventory value which has 28/81.2
GWP20_multiplier = (81.2)/25 # Emissions and abatement from non-MAC (100-yr, 25 CO2eq) to IPCC AR6 20-year, 81.2 CO2eq.
if GWP20year; GWPmulti = GWP20_multiplier / EPACO2eqNonCO2_to_GHGInv_multiplier; print("{20 year GWP} "); end
if !GWP20year; GWPmulti = 1; print("{NO GWP multiple} "); end # Comment out for 20-year GWP simulation 
# CO2_taxrate = 56.6298#<= 20-yr GWP (post Consumption utility elasticities from Marc)
# CH4_taxrate = 64.8137#<= w 20-yr GWP (post Consumption utility elasticities from Marc) 
    ###Optimal 
    # CO2_taxrate = 9.122#<= 20-yr GWP (post Consumption utility elasticities from Marc)
    # CH4_taxrate = 47#<= w 20-yr GWP (post Consumption utility elasticities from Marc) 
### Oil/Gas tax, 20-yr GWP
# CO2_taxrate =  56.6298 #<=oil/gas tax, 20=yr GWP (post Consumption utility elasticities from Marc)
# CH4_taxrate =  134.89 #<=oil/gas tax, 20=yr GWP (post Consumption utility elasticities from Marc)
    ###Optimal oil/gas tax, 20=yr GWP (post Consumption utility elasticities from Marc)
    # CO2_taxrate = 4.345 #<=oil/gas tax, 20=yr GWP (post Consumption utility elasticities from Marc)
    # CH4_taxrate = 119 #<=oil/gas tax, 20=yr GWP (post Consumption utility elasticities from Marc)

CH4x5 = true
CH4x5 = false ###> Comment out to implement x 5 multiple CH4 emissions
if CH4x5; Undercount_oilgasch4_multiplier = 5; print("::5 x CH₄::"); end#(Plant et al 2022   #Alvarez 2018 said 1.6)
if !CH4x5; Undercount_oilgasch4_multiplier = 1; print("::1 x CH₄::");end# default (comment out for oilgas multi)
# CO2_taxrate = 49.7935#<=updated target x 5 oil/gas (post Consumption utility elasticities from Marc)
# CH4_taxrate = 43.803#79999#<=updated target x 5 oil/gas (post Consumption utility elasticities from Marc)
    ###Optimal x 5 combination
    # CO2_taxrate = 0#<=updated target x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CH4_taxrate = 43.803#79999#<=updated target x 5 oil/gas (post Consumption utility elasticities from Marc)

# CO2_taxrate = 58.117 #<=Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
# CH4_taxrate = 8.504 #<=Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
    ###Optimal Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CH4_taxrate = 8.504 #<=Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CO2_taxrate = 0 #<=Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)

# CO2_taxrate = 49.7931 #<=oil/gas tax, x 5 oil/gas (post Consumption utility elasticities from Marc)
# CH4_taxrate = 50.2133 #<=oil/gas tax, x 5 oil/gas (post Consumption utility elasticities from Marc)
    # ###Optimal oil/gas tax, x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CO2_taxrate = 2.309 #<=oil/gas tax, x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CH4_taxrate = 47 #<=oil/gas tax, x 5 oil/gas (post Consumption utility elasticities from Marc)

# CO2_taxrate = 58.117 #<=oil/gas tax, Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
# CH4_taxrate = 11.2758 #<=oil/gas tax, Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
    # ###Optimal oil/gas tax, Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CO2_taxrate = 0 #<=oil/gas tax, Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)
    # CH4_taxrate = 11.2758 #<=oil/gas tax, Both 20=yr GWP & x 5 oil/gas (post Consumption utility elasticities from Marc)

multioil_CH4ratio = (sum(only(MAC_CH4_totemiss[:,["agr_livestock", "agr_rice", "COL", "wst_land", "wst_water"]]))+ only(MAC_CH4_totemiss[:,:GAS])*Undercount_oilgasch4_multiplier)/sum(only(MAC_CH4_totemiss[:,2:end]))
ReductTarget = ReductTargetbase - CH4oftarget + GWPmulti * multioil_CH4ratio *CH4oftarget;println(": Reduction target=",ReductTarget)
print("CO2tax: $CO2_taxrate, "); println("CH4tax: $CH4_taxrate")
print("$CH4abatement CH4 Abatement: "); print("$Kmobile mobile Kapital: ")
#########################
### End scenario, => Start data set up
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

CO2Int = DenseAxisArray(zeros(length(Jp)),Jp) # Initialise DAA
# 2024 GHG Inv Table 3-5: CO2 Emissions from Fossil Fuel Combustion by Fuel Type and Sector (MMT CO2 Eq.)
TotalCO2EmMMt_coal =  895.9 # <=2022 # 835.6 <=2020  EPA Inventory CO2 Stationary Combustion - Coal sum(Electricity, Industrial, Commercial, & Residential=0)
### (Option with no transport emissions in the model: assumption that direct)
### TotalCO2EmMMt_gas_oil = 2104.80
TotalCO2EmMMt_Natural_gas =  1740.70 #<=2022(incl 70.2 transport, Table 3-5:)  1615.7 #<=2020 (incl 58.7 transport, Table 3-5:)
TotalCO2EmMMt_Petroleum =  2115.40 #<=2022 (incl 1,681.10 transport, Table 3-5) # 1890.0 #<=2020 (incl 1,514.2 transport, Table 3-5)
# Option with all Transport CO2 emissions attributed to oil inputs: assumption that forms of oil fuel all transport that has direct CO2 emissions, and so taxing CO2 is total emissions from all oil as an input 
# TotalCO2EmMMt_gas_oil = TotalCO2EmMMt_Natural_gas + TotalCO2EmMMt_Petroleum # EPA inventory all CO2 transport + CO2 Stationary Combustion - both Oil & Natural Gas sum(Electricity, Industrial, Commercial, Residential & oil U.S. Territories)Sta
# (MMtCO2eq/$Bill x 10-^3 -> Gt/$B=t/$) EPA Inventory 2022 sum of CO2 MMt for coal, and for gas & oil per Billion of total benchmark intermediate input from sector
TotalCO2EmGt_coal = TotalCO2EmMMt_coal*10^-3
TotalCO2EmGt_gas  = TotalCO2EmMMt_Natural_gas*10^-3
TotalCO2EmGt_oil  = TotalCO2EmMMt_Petroleum*10^-3

CO2Int[:coa] =  TotalCO2EmGt_coal/sum(id_0[:coa,:]) # Sum of value of coal as an input into all sectors
CO2Int[:gas] =   TotalCO2EmGt_gas/sum(id_0[:gas,:])
CO2Int[:oil] =   TotalCO2EmGt_oil/sum(id_0[:oil,:])

TotCO2bnchmk =  TotalCO2EmGt_coal + TotalCO2EmGt_oil + TotalCO2EmGt_gas
TotCH4bnchmk = sum(MAC_CH4_WiNDC_tot[1,2:end]) # 2020 was 0.6902703880400002
TotGHGbnchmk =  TotCO2bnchmk + TotCH4bnchmk # 2020 was 5.0315703880400005

realCO2taxpc = []
for ff in [:coa,:oil,:gas]
push!(realCO2taxpc,CO2Int[ff]*CO2_taxrate)
end
CO2_realtaxpc = DenseAxisArray(realCO2taxpc,([CO2_taxrate])) # The index is the $/t rate, the values at tax % on inputs

### This returns the total % tax tax on each CH4 sector
function get_real_ch4tax(solve::Symbol)
    df = DataFrame(Solve=Symbol[],Sector=String[],Index=Symbol[],tax_rate=Float64[])
    for vam in [[VAS] ;VAMcommodSet]
        for c in CH4sectors
            if (vam in VAMcommodSet && VAM_costover[MPSGE.name(vam),c]>1 && value(MPSGE.tax_revenue(vam[c],RA))<0) || (vam == VAS && value(MPSGE.tax_revenue(vam[c],RA))<0)
                push!(df, [solve split(string(vam),"[")[1] c value(MPSGE.taxes(MPSGE.netputs(vam[c],PVA[:compen])[1],RA))])#value(MPSGE.tax_revenue(vam[c], RA))])
            end
        end
    end
    return df
end

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

### Proportion of residential electricity consumption for vehicle charging
veh_chrg_pc = 5252/1509233 # https://www.eia.gov/totalenergy/data/monthly/pdf/sec7_19.pdf
# Alternative, 2020, similar, veh_chrg_pc = 3.2/824.3 # electric vehicle charging/total https://www.eia.gov/consumption/residential/data/2020/c&e/pdf/ce5.1a.pdf 
# function Multiloop(CO2_taxrate,CH4_taxrate) 

MultiNat = MPSGEModel()

@parameters(MultiNat, begin
    ta[j=Jp], ta_0DAA[j]
    ty[j=Jp], ty_0DAA[j]
    tm[j=Jp], tm_0DAA[j]
    tx[j=Jp], tx_0DAA[j] # all 0s, set up for export tax to represent foreign tariffs
    CH₄_tax, 0.
    CH4_secs[j=Jp], CH4_tax_switch_on_sectors[j]  #Boolean vector only CH4sectors matter, 1 = tax, 0 = no tax 
    CO₂_tax, 0.
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
# @auxiliary(MultiNat, CO2em, index = [[:gas, :oil, :coa]])
@auxiliary(MultiNat, CO2em[[:gas, :oil, :coa]])
@auxiliary(MultiNat, CO2TotEm, description = "Total CO2 emissions from fossil fuels")
# Variables to track and report levels of CH4 emissions
# @auxiliary(MultiNat, CH4em, index = [[:agr,:coa,:gas,:pip,:oil,:wst]])
@auxiliary(MultiNat, CH4em[[:agr,:coa,:gas,:pip,:oil,:wst]])
@auxiliary(MultiNat, CH4TotEm, description = "Total CH4 emissions")
@auxiliary(MultiNat, TotEm, description = "Total both emissions")


@consumer(MultiNat, RA, description = "Representative Agent")
@consumer(MultiNat, ROW, description = "Sink for export tax")
# Domestic production for all sectors
# for j∈Jp
#         @production(MultiNat, Y[j], [t=0, s=Elas[j, :SAGE_klem_Y], va=> s=Elas[j,:SAGE_kle_VAE],sm => s=Elas[j,:E3_m_ID]],begin # [t=0, s=0, sv=> s=0]
#         [@output(PY[i],ys_0[i,j], t, taxes = [Tax(RA,ty[j])]) for i∈Ip]... 
#         [@input(PA[i], id_0[i,j], sm, taxes = [Tax(RA,CO₂_tax * CO2Int[i])]) for i∈Ip]...
#         @input(PVAM[j], sum(va_0[VA,j]), va)
#         end)
# end

ID = [i for i ∈ Ip if i∉[:oil, :coa, :gas, :uel, :pet, :rnw] ] # Intermediate inputs EXCEPT oil and min
for j∈Jp
    @production(MultiNat, Y[j], [t=0,
    s=Elas[j, :SAGE_klem_Y], vae=>s=Elas[j,:SAGE_kle_VAE], sm=>s=Elas[j,:E3_m_ID],
    va=>vae=0, En=>vae=Elas[j,:SAGE_ene], PrimENRG=>En=Elas[j,:SAGE_en], Elec=>En=Elas[j,:SAGE_en],
    oilgas=>PrimENRG=Elas[j,:E3_e_E_El], inElec=>Elec=Elas[:uel,:SAGE_en]### Elec preference between coa/gas/rnw doesn't change between sectors - common to all sectors 
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
if (Kmobile=="yes")
for j∈Jp
        @production(MultiNat, VAS[j], [t=0, s=0, va => s=Elas[j,:SAGE_kl_VA]], begin # #     @production(MultiNat, VAS[j], [t=0, s=0, va => s= 1], begin 
        [@output(PVAM[j],sum(va_0[VA,j]), t)]... 
        [@input(PVA[va], va_0[va,j], va, taxes = [Tax(RA,CH₄_tax* CH4_secs[j]* VASInt[j])]) for va∈VA]...
        end)
end
elseif (Kmobile=="no")
    for j∈Jp
        @production(MultiNat, VAS[j], [t=0, s=0, va => s=Elas[j,:SAGE_kl_VA]], begin # #     @production(MultiNat, VAS[j], [t=0, s=0, va => s= 1], begin 
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
        for c∈CH4sectors # Here we filter for sector-specific methane tax cases (CH$sectors is boolean)
            for vam in VAMcommodSet
                if VAM_costover[MPSGE.name(vam),c]>1 # Some sectors are still cumulatively -negative costs at $5/t, so filtering those out.
                    @production(MultiNat, vam[c], [t=0, s=0, va => s=Elas[c,:SAGE_kl_VA]], begin
                        [@output(PVAM[c],sum(va_0[VA,c]), t)]... 
                        [@input(PVA[va], va_0[va,c]*VAM_costover[MPSGE.name(vam),c], va, taxes = [Tax(RA, CH₄_tax*CH4_secs[c]*VAM_CH4EmInt[MPSGE.name(vam),c])]) for va∈VA]... #CH4_secs is a boolean to be able to turn of any sector from the CH4 tax
                    end)
                end
            end
        end
    elseif (Kmobile=="no")
        for c∈CH4sectors # Here we filter for sector-specific methane tax cases (CH$sectors is boolean)
            for vam in VAMcommodSet
                if VAM_costover[MPSGE.name(vam),c]>1 # Some sectors are still cumulatively -negative costs at $5/t, so filtering those out.
                    @production(MultiNat, vam[c], [t=0, s=0, va => s=Elas[c,:SAGE_kl_VA]], begin
                       [@output(PVAM[c],sum(va_0[VA,c]), t)]... 
                        @input(PVAK[c],va_0[:surplus,c]*VAM_costover[MPSGE.name(vam),c], va, taxes = [Tax(RA, CH₄_tax*CH4_secs[c]*VAM_CH4EmInt[MPSGE.name(vam),c])])
                        @input(PVAL,    va_0[:compen,c]*VAM_costover[MPSGE.name(vam),c], va, taxes = [Tax(RA, CH₄_tax*CH4_secs[c]*VAM_CH4EmInt[MPSGE.name(vam),c])])
                    end)
                end
            end
        end
    end
end

for m∈M
    @production(MultiNat, MS[m], [t=0, s=0], begin
        [@output(PM[m], sum(ms_0[i,m] for i∈Ip), t)]...
        [@input(PY[i], ms_0[i,m], s) for i∈Ip]...
    end)
end
 
for i∈Ip
    @production(MultiNat, A[i], [t= 2, s=0, dm => s=Elas[i,:SAGE_E3_Av_Armington]], begin
        [@output(PA[i], a_0[i,:value], t, taxes=[Tax(RA,ta[i])],reference_price=1-ta_0[i,:value])]... 
        [@output(PFX, x_0[i,:exports], t, taxes=[Tax(ROW,tx[i])])]...
        [@input(PM[m], md_0[i,m], s) for m∈M]...
        @input(PY[i], y_0[i,:value], dm)
        @input(PFX, m_0[i,:imports], dm, taxes = [Tax(RA,tm[i]),
## Emissions tariffs ## Does not work, does not make sense.
        # Tax(RA,CH₄_tax* VASInt[i]), # Actually INCREASES CO2 emissions (effect on CH4 seems reasonable) 
        # Tax(RA, CO₂_tax * CO2Int[i]) # massively REDUCES reductions
        ],
        reference_price=1+tm_0[i,:value])# No tariff on CO2 bc oil and gas are taxed as inputs to production, which includes these imports: Tax(RA,CO₂_tax * CO2Int[i]),
    end)
end;

## Final Consumption with CES elasticity of Demand (Consumption utility)
# @production(MultiNat, FDem, [t=0, s=.999999], begin
@production(MultiNat, FDem, [t=0,           s=0.5, # Big difference in utility change comes from here. 
transp=>s=0.35,                             
non_pers_transp=>transp=0.25,                  pers_transp=>transp=0.4, 
veh_exp=>pers_transp=0.25,                     fuels=>pers_transp=0.25,
veh_elect=>fuels=Elas[:uel,:SAGE_en],  petr=>fuels=0, # veh_fuel leaf, elasticity value v minor impact, but nesting definition needed
# non-transport expenditures  
non_transp=>s=0.5,
goods=>non_transp=.87,                        housing_exp=>non_transp=0.3,  # goods value is average of E3 intermediate demand 'materials' substitution elasticities
own_occ_exp=>housing_exp=0.25,                 home_nrg_expd=>housing_exp=0.75,
elect=>home_nrg_expd=Elas[:uel,:SAGE_en],   homefuels=>home_nrg_expd=0.5,
], begin
    @output(PC, sum(pce_0),t)
    # [@input(PA[i], pce_0[i,:pce],s) for i in Ip]...
    # transport expenditures
    [@input(PA[i], pce_0[i,:pce],non_pers_transp)        for i in [:air :trn :wtt :trk :grd :otr]]...  # non_pers_transp 
    [@input(PA[i], pce_0[i,:pce],veh_exp)                for i in [:mot :ote :mvt]]... # Vehicle and Service Expenditures pers_transp  0.4
    [@input(PA[i], pce_0[i,:pce]*veh_chrg_pc, veh_elect)     for i in [:uel, :rnw]]... # Electricity 0.25 small fraction, for vehicle charging
    [@input(PA[i], pce_0[i,:pce],petr)               for i in [:pet, :oil]]...   #Motor Vehicle Fuels 0.25 # pet
    # non-transport expenditures
    [@input(PA[i], pce_0[i,:pce],goods)                  for i in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa :oil ]),Ip)]... # Everything else, filter out all specified
    [@input(PA[i], pce_0[i,:pce],own_occ_exp)            for i in [:hou,:ore]]...    #Owner-occupied and Rental Expenditures 
    [@input(PA[i], pce_0[i,:pce]*(1-veh_chrg_pc),elect)  for i in [:uel, :rnw]]...    # 'Electricity', home_nrg_exp, 0.75
    [@input(PA[i], pce_0[i,:pce],homefuels)              for i in [:ugs, :coa]]... # natural gas, homefuels, 0.5   # pip and gas are 0
    
    # [@input(PA[:ugs], pce_0[:ugs,:pce],ugs)] # natural gas, homefuels, 0.5   # pip and gas are 0
    # [@input(PA[:coa], pce_0[:coa,:pce],coa)]  # 'fuel oil'  TODO this should really have PART of pet...but not much. As is, there is no oil, and coal extraction is very small.
end)

if (Kmobile=="yes")
    @demand(MultiNat, RA, begin
    # [@final_demand(PA[i], pce_0[i,:pce]) for i∈Ip]...
    @final_demand(PC, sum(pce_0))
    @endowment(PFX, only(bopdef_0))
    [@endowment(PA[i], -sum(fd_0[i,xfd] for xfd∈FD if xfd!=:pce)) for i∈Ip]...
    [@endowment(PVA[va], sum(va_0[va,j] for j∈Jp)) for va∈VA]...
    end, elasticity = .999999)
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

@demand(MultiNat, ROW, begin # Add Rest of World to distribute 'export tax' revenue.
    @final_demand(PC, 0)
    @endowment(PFX, 0)
    end, elasticity = .999999) # elasticity doesn't matter, only one demand good

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

### print for emissions excise taxes
if contains(string(MultiNat.productions[Symbol("A[oil]")]),"CH₄"); println("Yes CH4 tariff!!!"); end
if contains(string(MultiNat.productions[Symbol("A[oil]")]),"CO₂"); println("Yes CO2 tariff!!!"); end 
# Benchmark 
fix(RA, sum(pce_0[Ip,:pce])) # Numeraire, fixed at benchmark
### Note: Benchmark doesn't solve to appropriate residual at 0 interation because of margins of slack activity. 
### Does balance with interactions or slack vars and production commented out.
# solve!(MultiNat , cumulative_iteration_limit = 0)
solve!(MultiNat)#, output_minor_iterations_frequency=1)
# fullvrbnch = generate_report(MultiNat); rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg); fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var])
# print(sort(fullvrbnch, :var))#:bmkmarg))#:bnchmrk))#
if CH4abatement=="yes" ; realCH4taxpcc_bnch = get_real_ch4tax(:bnch) ;end

if !@isdefined(CESutility_multinat)
    include("CESUtility.jl")
end

TaxRev = DataFrame(solve=Symbol[], Total_Revenue=Float64[], Change_in_Rev=Float64[], Income=Float64[])
EqVar  = DataFrame(solve=Symbol[], Mev=Float64[], MEVdelta=Float64[],  utilCES=Float64[], Ut_perc=Float64[])#EV=Float64[],
function totalrevenue()
-(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+ # taxes are negative in production, but positive for revenue
sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
if CH4abatement=="yes"
    sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[MPSGE.name(vam),c]>1]) for vam in VAMcommodSet])
else 0
end)
end
totrevbnch = totalrevenue()
income      = totrevbnch + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:bnch totrevbnch totrevbnch-totrevbnch income])

totdem = -value(FDem)*value(compensated_demand(FDem,PC))
# Mevbnch2 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])

# Mevbnch = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
# MevbnchCES = prod([(1/value(PA[i])) for i in Ip if value(PA[i])>0])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
# EVbnch2  = Mevbnch2 - value(RA) 
# elasRA = MPSGE.elasticity(MultiNat.productions[:FDem].input)
# priceIndCES = (sum([((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) * value(PA[i])^(1-elasRA) for i in Ip]))^(1/(1-elasRA))
# PrIndexbnchCES = (sum([pce_0[i,:pce]/sum(pce_0)*value(PA[i])^(1-elasRA) for i in Ip]))^(1/(1-elasRA))
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.

# function CESutil(sector::ScalarSector,I::Vector, nest::Symbol)
#     sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasticity(sector,nest)-1)/elasticity(sector,nest)) for i in I if value(compensated_demand(FDem,PA[i]))>0])^(elasticity(sector,nest)/(elasticity(sector,nest)-1))
# end
# function utilCESfn()
#     sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
# end
# utilCESbnch = utilCESfn()
# utilCD    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip])

utilCESbnch = CESutility_multinat(pce_0)
Ut_pc = (utilCESbnch-utilCESbnch)/utilCESbnch
MEVdelt  = value(RA)*Ut_pc
Mevbnch =  value(RA)+MEVdelt

# push!(EqVar, [:bnch utilCESbnch Mevbnch EVbnch value(RA)*(utilCESbnch-utilCESbnch)/utilCESbnch  (utilCESbnch-utilCESbnch)/utilCESbnch*100] )
push!(EqVar, [:bnch  Mevbnch MEVdelt utilCESbnch Ut_pc*100] )
# sum([value(FDem)*value(compensated_demand(FDem,PA[i]))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])

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
totrevWiNcntfac = totalrevenue()

if CH4abatement=="yes"; realCH4taxpc_WiNcntfac = get_real_ch4tax(:WiNcntfac); end
income          = totrevWiNcntfac +  sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:WiNcntfac totrevWiNcntfac totrevWiNcntfac - totrevbnch income])

# LaspeyresPrIndex = sum([value(PA[i])*pce_0[i,:pce] for i in Ip])/value(RA)
# PaaschePrIndex = sum([value(PA[i])*value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])/sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
# Wmaybe = (value(RA)/LaspeyresPrIndex - value(RA))/value(RA) 
# totdemWiNcntfac = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
# EVWiNcntfac2  = MevWiNcntfac2 - value(RA) 
# utilWiNcntfac    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
# utilCESWiNcntfac = utilCESfn()
utilCESWiNcntfac = CESutility_multinat(pce_0)
Ut_pcWiNcntfac = (utilCESWiNcntfac-utilCESbnch)/utilCESbnch
MEVdeltWiNcntfac  = value(RA)*Ut_pcWiNcntfac
MevWiNcntfac =  value(RA)+MEVdeltWiNcntfac
push!(EqVar, [:WiNcntfac MevWiNcntfac MEVdeltWiNcntfac utilCESWiNcntfac Ut_pcWiNcntfac*100] )
# push!(EqVar, [:bnch utilCESbnch Mevbnch EVbnch value(RA)*(utilCESbnch-utilCESbnch)/utilCESbnch  (utilCESbnch-utilCESbnch)/utilCESbnch*100] )
# utilWiNcntfac2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) for i in Ip])
# MevWiNcntfac2 = prod([(1/value(PA[i]))^((value(FDem)*value(compensated_demand(FDem,PA[i])))/totdem) for i in Ip])* only(filter(x->x.solve==:WiNcntfac,TaxRev)[!,:Income]) #])
# MevWiNcntfac    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:WiNcntfac,TaxRev)[!,:Income]) #]
# MevWiNcntfacCES = prod([(1/value(PA[i])) for i in Ip if value(PA[i])>0])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
# EVWiNcntfac  = MevWiNcntfac - value(RA) 
# EVWiNcntfacCES  = MevWiNcntfacCES - value(RA) 
# EVWiNcntfacCES/value(RA) * 100
# push!(EqVar, [:WiNcntfac utilCESWiNcntfac MevWiNcntfac EVWiNcntfac value(RA)*(utilCESWiNcntfac-utilCESbnch)/utilCESbnch (utilCESWiNcntfac-utilCESbnch)/utilCESbnch*100])

# fullvrWiNcntfact = generate_report(MultiNat)
# rename!(fullvrWiNcntfact, :value => :WiNcntfact, :margin => :WiNcntfactmarg)
# fullvrWiNcntfact[!,:var] = Symbol.(fullvrWiNcntfact[:,:var])
# print(sort(fullvrWiNcntfact, :var))#:bmkmarg))#:bnchmrk))#

for (n,i) in enumerate(Ip)
    FDemand[n,:WiNcntfact] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end

## DataFrame to hold the industry indices with descriptions
Rs = DataFrame([Y value.(Y) first.(last.(string.(Y),4),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
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
if CH4abatement=="yes"; realCH4taxpc_ch4 = get_real_ch4tax(:CH₄); end

totrevch4 = totalrevenue()
income      = totrevch4 +  sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:ch4 totrevch4 totrevch4-totrevbnch income])

# totdemch4 = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
# # Mevch42 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:ch4,TaxRev)[!,:Income]) #])
# Mevch4    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:ch4,TaxRev)[!,:Income]) #]
# # EVch42  = Mevch42 - value(RA) 
# EVch4  = Mevch4 - value(RA) 
# No exponent ^(1/elasRA) on share. Results are the same with either consistently applied. Not sure which is most 'correct'.
# utilCESch4  = utilCESfn()
utilCESch4  = CESutility_multinat(pce_0)
Ut_pcch4 = (utilCESch4-utilCESbnch)/utilCESbnch
MEVdeltch4  = value(RA)*Ut_pcch4
Mevch4 =  value(RA)+MEVdeltch4
push!(EqVar, [:ch4 Mevch4 MEVdeltch4 utilCESch4 Ut_pcch4*100] )

# utilCESch4  = CESutility_multinatAlt()
# Mevch4CES = prod([(1/value(PA[i])) for i in Ip if value(PA[i])>0])* only(filter(x->x.solve==:bnch,TaxRev)[!,:Income])
# EVch4CES  = Mevch4CES - value(RA) 
# EVch4CES/value(RA) * 100
# PrIndexCESch4 = (sum([pce_0[i,:pce]/sum(pce_0)*value(PA[i])^(1-elasRA) for i in Ip]))^(1/(1-elasRA))
# PrIndexCESch4 * utilCESch4 - value(RA)
# utilch42    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemch4) for i in Ip])
# utilch4    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])

# push!(EqVar, [:ch4 utilCESch4 Mevch4 EVch4 value(RA)*(utilCESch4-utilCESbnch)/utilCESbnch (utilCESch4-utilCESbnch)/utilCESbnch*100])

for (n,i) in enumerate(Ip)
    FDemand[n,:ch4tax] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end
Rs = DataFrame([Y value.(Y) first.(last.(string.(Y),4),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

CO2_ch4tax = value(CO2TotEm)
CH4_ch4tax = value(CH4TotEm)
GHG_ch4tax = value(TotEm)
# fullvrch4 = generate_report(MultiNat)
# rename!(fullvrch4, :value => :ch4tax, :margin => :ch4marg)
# fullvrch4[!,:var] = Symbol.(fullvrch4[:,:var])
# print(sort(fullvrch4, :ch4tax, by= abs, rev=true))#:ch4tax))#:var):ch4marg

set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
set_value!(CO₂_tax, 0.0)
set_value!(ta,ta_0DAA[Jp]);
set_value!(tm,tm_0DAA[Jp]);
solve!(MultiNat)# Temp measure to address residual price changes 
set_value!(CO₂_tax, CO2_taxrate)
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
solve!(MultiNat, cumulative_iteration_limit=10000) #;

if CH4abatement=="yes"; realCH4taxpc_co2 = get_real_ch4tax(:CO₂); end
totrevco2 = totalrevenue()
income      = totrevco2 +  sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:co2 totrevco2 totrevco2-totrevbnch income])

# totdemco2 = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
# # Mevco22 = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])* only(filter(x->x.solve==:co2,TaxRev)[!,:Income]) #])
# Mevco2    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:co2,TaxRev)[!,:Income]) #]
# # EVco22  = Mevco22 - value(RA) 
# EVco2  = Mevco2 - value(RA) 
# utilCESco2  = utilCESfn()
utilCESco2  = CESutility_multinat(pce_0)
Ut_pcco2 = (utilCESco2-utilCESbnch)/utilCESbnch
MEVdeltco2  = value(RA)*Ut_pcco2
Mevco2 =  value(RA)+MEVdeltco2
push!(EqVar, [:co2 Mevco2 MEVdeltco2 utilCESco2 Ut_pcco2*100] )

# utilCESco2  = CESutility_multinatAlt()
# utilco22    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemco2) for i in Ip])
# utilco2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
# push!(EqVar, [:co2 utilCESco2 Mevco2 EVco2 value(RA)*(utilCESco2-utilCESbnch)/utilCESbnch (utilCESco2-utilCESbnch)/utilCESbnch*100])


for (n,i) in enumerate(Ip)
    FDemand[n,:co2tax] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:4]
# Rs[:,2][68:71]
CO2_co2tax = value(CO2TotEm)
CH4_co2tax = value(CH4TotEm)
GHG_co2tax = value(TotEm)
# fullvrco2 = generate_report(MultiNat)
# rename!(fullvrco2, :value => :CO2tax, :margin => :CO2marg)
# fullvrco2[!,:var] = Symbol.(fullvrco2[:,:var])

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

if CH4abatement=="yes"; realCH4taxpc_both = get_real_ch4tax(:both); end
totrevboth = totalrevenue()
income      = totrevboth + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
push!(TaxRev, [:both totrevboth totrevboth-totrevbnch income])

# totdemboth = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
# Mevboth    = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])* only(filter(x->x.solve==:both,TaxRev)[!,:Income]) #]
# EVboth  = Mevboth - value(RA) 
# utilCESboth  = utilCESfn()
utilCESboth  = CESutility_multinat(pce_0)
Ut_pcboth = (utilCESboth-utilCESbnch)/utilCESbnch
MEVdeltboth  = value(RA)*Ut_pcboth
Mevboth =  value(RA)+MEVdeltboth
push!(EqVar, [:both Mevboth MEVdeltboth utilCESboth Ut_pcboth*100] )
# utilCESboth  = CESutility_multinatAlt()
# utilboth2    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdemboth) for i in Ip])
# utilboth    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip ])
# push!(EqVar, [:both utilCESboth Mevboth EVboth value(RA)*(utilCESboth-utilCESbnch)/utilCESbnch (utilCESboth-utilCESbnch)/utilCESbnch*100])

for (n,i) in enumerate(Ip)
    FDemand[n,:both] = value(FDem)*value(compensated_demand(FDem,PA[i]))
end
# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:5]
# Rs[:,2][68:71]
CO2_bothtax = value(CO2TotEm)
CH4_bothtax = value(CH4TotEm)
GHG_bothtax = value(TotEm)

# fullvrboth = generate_report(MultiNat)
# rename!(fullvrboth, :value => :bothtaxes, :margin => :bothmarg)
# fullvrboth[!,:var] = Symbol.(fullvrboth[:,:var])

###### Add Quant diff and % of consumption columns for Final Demand report dataframe
# FDemand[:,:ch4Qdelta]=FDemand[:,:ch4tax].-FDemand[:,:bnch]
# FDemand[:,:CO2Qdelta]=FDemand[:,:co2tax].-FDemand[:,:bnch]
# FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
# FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
# FDemand[:,:ch4Qpc]=FDemand[:,:ch4tax]./sum(FDemand[:,:ch4tax])*100
# FDemand[:,:CO2Qpc]=FDemand[:,:co2tax]./sum(FDemand[:,:co2tax])*100
# FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

#Generate Dataframe with all results (including names expressions)
# FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrco2, fullvrboth, fullvrWiNcntfact, on = [:var], makeunique=true);
# Compare = FullResults[1:end,[1,2,4,6,8,10]];
# ## Sum the difference of each tax applied individually
# Compare.sum = Compare.ch4tax .- 1 + Compare.CO2tax .-1

## sum difference between the sum of the individual taxes to the combined taxes
# Compare.diff = Vector{Union{Missing, Float64}}(undef, length(Compare[:,1]))
# for n in 1:length(Compare[:,1])
#     if Compare.ch4tax[n] == 0 || Compare.CO2tax[n] == 0
#         Compare.diff[n] = missing
#     elseif Compare.bothtaxes[n] >= 1 &&  Compare.sum[n] >= 0
#         Compare.diff[n] = (Compare.bothtaxes[n] - 1) - Compare.sum[n]
#     elseif Compare.bothtaxes[n] >= 1 &&  Compare.sum[n] <= 0
#         Compare.diff[n] = 1 - Compare.bothtaxes[n] + Compare.sum[n]
#     elseif Compare.bothtaxes[n] <= 1 &&  Compare.sum[n] <= 0
#         Compare.diff[n] = 1 - Compare.bothtaxes[n] + Compare.sum[n]
#     elseif Compare.bothtaxes[n] <= 1 &&  Compare.sum[n] >= 0
#         Compare.diff[n] = 1 - Compare.bothtaxes[n] - Compare.sum[n]
#     end
# end

# Compare[!,:var] = Symbol.(Compare[:,:var]);
# print(sort!(Compare, :var));
# println(sort!(Compare, :diff, by = abs, rev=true))#[1:25,:])
# CSV.write("C:\\Users\\Eli\\Box\\CGE\\MPSGE-JL\\First Mulit GHG taxes Paper\\MultiResults$(Dates.format(now(),"yyyy-mm-d_HhM")).csv", Compare, missingstring="missing", bom=true)
# compCO2em = filter(:var => ==(:CO2TotEm), Compare);
# compCH4em = filter(:var => ==(:CH4TotEm), Compare);
# compTotem = filter(:var => ==(:TotEm), Compare);

Emissions = DataFrame(
["CO2" "Gt" TotCO2bnchmk ;;
### Need a CO2 emissions variable from each solve
# Instead of compCO2em[:,:CO2tax], there's a CO2 emissions from CO2 taxes variable, and then a CO2 emissions from ch4 taxes variable.
CO2_co2tax;; #only(compCO2em[:,:CO2tax]) ;;
CO2_ch4tax ;; #only(compCO2em[:,:ch4tax]) ;;
CO2_ch4tax + CO2_co2tax - TotCO2bnchmk;;
CO2_bothtax;; #only(compCO2em[:,:bothtaxes]) ;;
CO2_ch4tax + CO2_co2tax - TotCO2bnchmk - CO2_bothtax; # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk ;;
CH4_co2tax ;; #only(compCH4em[:,:CO2tax]) ;;
CH4_ch4tax ;; #only(compCH4em[:,:ch4tax]) ;;
CH4_ch4tax + CH4_co2tax - TotCH4bnchmk ;;
CH4_bothtax ;; #only(compCH4em[:,:bothtaxes]) ;;
CH4_ch4tax + CH4_co2tax - TotCH4bnchmk - CH4_bothtax; # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk ;;
 GHG_co2tax ;; #only(compTotem[:,:CO2tax]) ;;
 GHG_ch4tax ;; # only(compTotem[:,:ch4tax]) ;;
 GHG_ch4tax + GHG_co2tax - TotGHGbnchmk ;;
 GHG_bothtax ;; #only(compTotem[:,:bothtaxes]) ;;
 GHG_ch4tax + GHG_co2tax - TotGHGbnchmk - GHG_bothtax],
["Emissions", "Unit", "Bnchmrk_Emissions", "CO2tax", "CH4tax","sum_of_taxes", "taxes_combined" ,"Interactions"])
Emissions_Mt = hcat(Emissions[:,1:2],Emissions[:,3:end].*10^3); Emissions_Mt[:,2] = ["Mt" ,"MtCO2eq" ,"MtCO2eq"];

EmissionReductionResults = DataFrame(
["CO2" "Gt" TotCO2bnchmk ;;
TotCO2bnchmk - CO2_co2tax ;;
TotCO2bnchmk - CO2_ch4tax ;;
TotCO2bnchmk - CO2_ch4tax + TotCO2bnchmk - CO2_co2tax ;;
TotCO2bnchmk - CO2_bothtax ;;
TotCO2bnchmk - CO2_ch4tax + TotCO2bnchmk - CO2_co2tax - (TotCO2bnchmk - CO2_bothtax); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk ;;
TotCH4bnchmk - CH4_co2tax ;;
TotCH4bnchmk - CH4_ch4tax ;;
TotCH4bnchmk - CH4_ch4tax + TotCH4bnchmk - CH4_co2tax ;;
TotCH4bnchmk - CH4_bothtax ;;
TotCH4bnchmk - CH4_ch4tax + TotCH4bnchmk - CH4_co2tax - (TotCH4bnchmk - CH4_bothtax); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk ;;
TotGHGbnchmk - GHG_co2tax ;;
TotGHGbnchmk - GHG_ch4tax ;;
TotGHGbnchmk - GHG_ch4tax + TotGHGbnchmk - GHG_co2tax ;;
TotGHGbnchmk - GHG_bothtax ;;
TotGHGbnchmk - GHG_ch4tax + TotGHGbnchmk - GHG_co2tax - (TotGHGbnchmk - GHG_bothtax)],
["Em_reductions", "Unit", "Bnchmrk_Emissions", "CO2tax_reduc", "CH4tax_reduc","Sum_of_each_tax", "both_taxes_combined" ,"Interactions"])
EmissionReductionResults_Mt = hcat(EmissionReductionResults[:,1:2],EmissionReductionResults[:,3:end].*10^3); EmissionReductionResults_Mt[:,2] = ["Mt" ,"MtCO2eq" ,"MtCO2eq"];

### TODO Bring back later : not worth the dataframe generation right now, or figuing out the whole change to Compare, but maybe later
# EmissionReductionResults_Mt.pc_diff .= -EmissionReductionResults_Mt.Interactions ./ EmissionReductionResults_Mt.Sum_of_each_tax .* 100 

# println(
# sort(filter(x -> x.var in [Symbol("Y[rec]"), Symbol("PVAM[uel]"),Symbol("PVAM[rnw]"),
# # Symbol("Y[oil]"), Symbol("Y[gas]"), Symbol("Y[pet]"), Symbol("Y[coa]"), Symbol("Y[agr]"), 
# Symbol("Y[uel]"), Symbol("Y[rnw]"), Symbol("Y[oil]"), Symbol("Y[gas]") 
# # Symbol("Y[uwt]"), Symbol("Y[wst]"), Symbol("PVA[compen]"), Symbol("PVA[surplus]"), Symbol("RA")
# ], Compare))
# ,
# )
println(Emissions_Mt)
println(EmissionReductionResults_Mt)# ;println("RA check, it's the 'both' bc that's the last simulation:",
# sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip]))
println(outerjoin(TaxRev,EqVar, on=:solve))
#Check Benchmark
# fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var]);
# print(sort!(fullvrbnch, [:bmkmarg]))
## End for Multiloop function version of model
# end
# Billions of dollars/Gigatons = $/t (RA in billions, EV in %)/Reductions in MMt/10^3
EVpertonGHG_CH4 = value(RA)*(EqVar[3,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[3,:CH4tax_reduc]*10^-3) # Total EV change from CH4 tax /per ton GHG reduced
EmissionReductionResults_Mt[3,:CH4tax_reduc]
EVpertonCO2_CH4 = value(RA)*(EqVar[3,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[1,:CH4tax_reduc]*10^-3) # Total EV change from CH4 tax /per ton CO2 reduced
EmissionReductionResults_Mt[1,:CH4tax_reduc]
EVpertonCH4_CH4 = value(RA)*(EqVar[3,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[2,:CH4tax_reduc]*10^-3)# Total EV change from CH4 tax /per ton CH4 reduced
EmissionReductionResults_Mt[2,:CH4tax_reduc]
# println("EV_ton, CH4 tax =\t",EVperton_CH4)
EVpertonGHG_CO2 = value(RA)*(EqVar[4,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[3,:CO2tax_reduc]*10^-3) # Total EV change from CO2 tax /per ton GHG reduced
EmissionReductionResults_Mt[3,:CO2tax_reduc]
EVpertonCO2_CO2 = value(RA)*(EqVar[4,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[1,:CO2tax_reduc]*10^-3) # Total EV change from CO2 tax /per ton CO2 reduced
EmissionReductionResults_Mt[1,:CO2tax_reduc]
EVpertonCH4_CO2 = value(RA)*(EqVar[4,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[2,:CO2tax_reduc]*10^-3) # Total EV change from CO2 tax /per ton CH4 reduced
EmissionReductionResults_Mt[2,:CO2tax_reduc]
# println("EV_ton, CO2 tax =\t",EVperton_CO2)
EVpertonGHG_both = value(RA)*(EqVar[5,:Ut_perc]*10^-2)/(EmissionReductionResults_Mt[3,:both_taxes_combined]*10^-3)
println("EV/ton, CH4 tax =\t",EVpertonGHG_CH4)
println("EV/ton, CO2 taxes =\t",EVpertonGHG_CO2)
println("EV_ton, both taxes =\t",EVpertonGHG_both)
println("Benefits_ton, both taxes =\$$(200 * 1.130480652)\t")

####################################
### This section is for marginal changes between simulations. The marg_em variable calculations and prints must be commented out for the initial run 
####################################
# print("CO2tax: $CO2_taxrate, "); println("CH4tax: $CH4_taxrate")
# println(EqVar[5,:MEVdelta]*10^9,"\t= \$Cost: EV")
# println(" ",EmissionReductionResults_Mt[3,:both_taxes_combined]*10^6*200 * 1.130480652,"\t= \$Benefit: Total tons reduced x SCC")
        marg_em_reduct = prev_em - EmissionReductionResults_Mt[3,:both_taxes_combined]*10^6
        # marg_em_reductCO2 = prev_emCO2 - EmissionReductionResults_Mt[1,:both_taxes_combined]*10^6
        # marg_em_reductCH4 = prev_emCH4 - EmissionReductionResults_Mt[2,:both_taxes_combined]*10^6
# println(marg_em_reduct, " Δ GHG")
# println(marg_em_reductCO2, " Δ CO2")
# println(marg_em_reductCH4, " Δ CH4")
# println(marg_em_reductCO2/marg_em_reductCH4 * 200 * 1.130480652 + 200 * 1.130480652,": CH4 tax w CO2")
# println(marg_em_reductCO2/marg_em_reductCH4 * 200 * 1.130480652,": Δ on CH4 tax")
# println(marg_em_reductCO2/marg_em_reductCH4, " Δ fraction on margin CH4 SCC tax")
# println(marg_em_reductCO2/(TotCO2bnchmk*10^6), " marginal reduction as % of bnch CO2")
# println(marg_em_reductCH4/(TotCH4bnchmk*10^6), " marginal reduction as % of bnch CH4")
# println(marg_em_reductCH4/marg_em_reductCO2 * 200 * 1.130480652 + 200 * 1.130480652,": CO2 tax w CH4")
# println(marg_em_reductCH4/marg_em_reductCO2 * 200 * 1.130480652,": Δ on CO2 tax")
# println(marg_em_reductCH4/marg_em_reductCO2, " Δ fraction on margin CO2 SCC tax")
# println(marg_em_reduct, "= marginal reduction")
        marg_MEVdelta = prev_MEVdelta - EqVar[5,:MEVdelta]*10^9
        println(" ",marg_MEVdelta, "= marginal cost")
        println(marg_em_reduct * 200 * 1.130480652, "= marginal benefit")
        prev_em = EmissionReductionResults_Mt[3,:both_taxes_combined]*10^6;
        # prev_emCO2 = EmissionReductionResults_Mt[1,:both_taxes_combined]*10^6;
#          prev_emCH4 = EmissionReductionResults_Mt[2,:both_taxes_combined]*10^6;
        prev_MEVdelta = EqVar[5,:MEVdelta]*10^9;
# println("CO2 reductions: ",EmissionReductionResults_Mt[1,:CO2tax_reduc])
# println("CH4 reductions: ",EmissionReductionResults_Mt[2,:CO2tax_reduc]) 


############################## Redundant, pretty sure _
# if CO2_taxrate ==  200 * 1.130480652
#     println("benefit as a % of direct: ",(EmissionReductionResults_Mt[2,:CO2tax_reduc])/(EmissionReductionResults_Mt[2,:CO2tax_reduc]))
#     println(CO2_taxrate* (EmissionReductionResults_Mt[2,:CO2tax_reduc])/(EmissionReductionResults_Mt[2,:CO2tax_reduc])," difference in CO2 price")
#     println("CO2_taxrate = CO2_taxrate - ($(EmissionReductionResults_Mt[2,:CO2tax_reduc])) / ($(EmissionReductionResults_Mt[2,:CO2tax_reduc])) * 200 * 1.130480652 # ",
#         200 * 1.130480652 * EmissionReductionResults_Mt[2,:CO2tax_reduc]/EmissionReductionResults_Mt[2,:CO2tax_reduc],"-> ",
#         CO2_taxrate - (200 * 1.130480652*(EmissionReductionResults_Mt[2,:CO2tax_reduc])/(EmissionReductionResults_Mt[2,:CO2tax_reduc])))
#     else
#     println("benefit as a % of direct: ",(CO2price[CO2iter-1,:CH4reduc]-EmissionReductionResults_Mt[2,:CO2tax_reduc])/(EmissionReductionResults_Mt[2,:CO2tax_reduc]))    
#     println("CO2_taxrate = CO2_taxrate - ($(CO2price[CO2iter-1,:CH4reduc]) - $(EmissionReductionResults_Mt[1,:CO2tax_reduc])) / $(EmissionReductionResults_Mt[2,:CO2tax_reduc]) * 200 * 1.130480652 # ",
#         200 * 1.130480652 * (CO2price[CO2iter-1,:CH4reduc]-EmissionReductionResults_Mt[1,:CO2tax_reduc])/EmissionReductionResults_Mt[2,:CO2tax_reduc],"-> ",
#         CO2_taxrate - (200 * 1.130480652 *(CO2price[CO2iter-1,:CH4reduc]-EmissionReductionResults_Mt[1,:CO2tax_reduc])/EmissionReductionResults_Mt[2,:CO2tax_reduc]))
# end
# push!(CO2price,[CO2iter CO2_taxrate EmissionReductionResults_Mt[1,:CO2tax_reduc] EmissionReductionResults_Mt[2,:CO2tax_reduc]])

# if CH4_taxrate ==  200 * 1.130480652
#     println("benefit as a % of direct: ",(EmissionReductionResults_Mt[1,:CH4tax_reduc])/(EmissionReductionResults_Mt[2,:CH4tax_reduc]))
#     println(CH4_taxrate* (EmissionReductionResults_Mt[1,:CH4tax_reduc])/(EmissionReductionResults_Mt[2,:CH4tax_reduc])," difference in CO2 price")
#     println("CH4_taxrate = CH4_taxrate - ($(EmissionReductionResults_Mt[1,:CH4tax_reduc])) / ($(EmissionReductionResults_Mt[2,:CH4tax_reduc])) * 200 * 1.130480652 # ",
#         200 * 1.130480652 * EmissionReductionResults_Mt[1,:CH4tax_reduc]/EmissionReductionResults_Mt[2,:CH4tax_reduc],"-> ",
#         CH4_taxrate - (200 * 1.130480652*(EmissionReductionResults_Mt[1,:CH4tax_reduc])/(EmissionReductionResults_Mt[2,:CH4tax_reduc])))
#     else
#     println("benefit as a % of direct: ",(CH4price[CO2iter-1,:CH4reduc]-EmissionReductionResults_Mt[1,:CH4tax_reduc])/(EmissionReductionResults_Mt[2,:CH4tax_reduc]))    
#     println("CH4_taxrate = CH4_taxrate - ($(CH4price[CO2iter-1,:CH4reduc]) - $(EmissionReductionResults_Mt[1,:CH4tax_reduc])) / $(EmissionReductionResults_Mt[2,:CH4tax_reduc]) * 200 * 1.130480652 # ",
#         200 * 1.130480652 * (CH4price[CO2iter-1,:CH4reduc]-EmissionReductionResults_Mt[1,:CH4tax_reduc])/EmissionReductionResults_Mt[2,:CH4tax_reduc],"-> ",
#         CH4_taxrate - (200 * 1.130480652 *(CH4price[CO2iter-1,:CH4reduc]-EmissionReductionResults_Mt[1,:CH4tax_reduc])/EmissionReductionResults_Mt[2,:CH4tax_reduc]))
# end
# push!(CH4price,[CH4iter CH4_taxrate EmissionReductionResults_Mt[1,:CH4tax_reduc] EmissionReductionResults_Mt[2,:CH4tax_reduc]])
#    200 * 1.130480652 * EmissionReductionResults_Mt[2,:CH4tax_reduc]/(EmissionReductionResults_Mt[1,:CH4tax_reduc]+EmissionReductionResults_Mt[2,:CH4tax_reduc])
# (EmissionReductionResults_Mt[1,:CO2tax_reduc]+EmissionReductionResults_Mt[2,:CO2tax_reduc])
# EmissionReductionResults_Mt[2,:CH4tax_reduc]+EmissionReductionResults_Mt[1,:CO2tax_reduc]
# println(EmissionReductionResults_Mt[3,:CO2tax_reduc])
# println(EmissionReductionResults_Mt[3,:CH4tax_reduc])
############ <---- ####################################

###################################
#### for Paris Reduction etc tables for paper
###################################
# println(Int64(round(TotGHGbnchmk*10^3, digits=0)))
# println(Int64(round(TotCO2bnchmk*10^3, digits=0)))
# println(Int64(round(TotCH4bnchmk*10^3, digits=0)))
# println(Int64(round(ReductTarget, digits=0)))
# println("(",Int64(round(ReductTarget/(TotGHGbnchmk*10^1))),"%)")

# #### Full table of results from CO2 tax
# println("\$",Int64(round(value(CO2_taxrate),digits=0)),"/t")
# println("-")
# # println(Int64(round(ReductTarget, digits=0)))
# # For SCC
# println(Int64(round(EmissionReductionResults_Mt[3,:CO2tax_reduc],digits=0)))
# # println("(",Int64(round(EmissionReductionResults_Mt[3,:CO2tax_reduc]/Emissions_Mt[3,:Bnchmrk_Emissions]*100,digits=0)),"%)")# % CO2 of CO2 tax
# println(Int64(round(EmissionReductionResults_Mt[1,:CO2tax_reduc],digits=0)))
# println("(",Int64(round(EmissionReductionResults_Mt[1,:CO2tax_reduc]/Emissions_Mt[1,:Bnchmrk_Emissions]*100,digits=0)),"%)")# % CO2 of CO2 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:CO2tax_reduc],digits=0)))
# println("(",Int64(round(EmissionReductionResults_Mt[2,:CO2tax_reduc]/Emissions_Mt[2,:Bnchmrk_Emissions]*100,digits=0)),"%)") # % CH4 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[3,:CO2tax_reduc],digits=0)))
# # println(Int64(round(ReductTarget, digits=0)))
# println(Int64(round(EmissionReductionResults_Mt[1,:CO2tax_reduc]/EmissionReductionResults_Mt[3,:CO2tax_reduc]*100, digits=0)),"%") # % CO2 of GHG from CO2 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:CO2tax_reduc]/EmissionReductionResults_Mt[3,:CO2tax_reduc], digits=2)*100),"%") # % CH4 of GHG from CO2 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:CO2tax_reduc]/EmissionReductionResults_Mt[3,:CO2tax_reduc], digits=2)*100),"%") # Spillovers = above
# println("\$",Int64(round(TaxRev[4,:Change_in_Rev])))
# println("\$",-1*Int64(round(EqVar[4,:MEVdelta])))#    value(RA)*(EqVar[4,:Ut_perc]*10^-2))))
# println("\$",-1*Int64(round(EVperton_CO2)),"/t")
# println(round(EqVar[4,:Ut_perc], digits=2),"%")
# # println("\$",-1*Int64(round(EVperton_CO2)),"/t")

# # # #### Full table of results from CH4 tax
# println("-")
# println("\$",Int64(round(value(CH4_taxrate),digits=0)),"/t")
# # println(Int64(round(ReductTarget, digits=0)))
# # For SCC
# # println(Int64(round(EmissionReductionResults_Mt[3,:CH4tax_reduc],digits=0)))
# # println("(",Int64(round(EmissionReductionResults_Mt[3,:CH4tax_reduc]/Emissions_Mt[3,:Bnchmrk_Emissions]*100,digits=0)),"%)") # % CO2 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[1,:CH4tax_reduc],digits=0)))
# println("(",Int64(round(EmissionReductionResults_Mt[1,:CH4tax_reduc]/Emissions_Mt[1,:Bnchmrk_Emissions]*100,digits=0)),"%)") # % CO2 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:CH4tax_reduc],digits=0)))
# println("(",Int64(round(EmissionReductionResults_Mt[2,:CH4tax_reduc]/Emissions_Mt[2,:Bnchmrk_Emissions],digits=2)*100),"%)") # % CH4 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[3,:CH4tax_reduc],digits=0)))
# # println(Int64(round(ReductTarget, digits=0)))
# println(Int64(round(EmissionReductionResults_Mt[1,:CH4tax_reduc]/EmissionReductionResults_Mt[3,:CH4tax_reduc], digits=2)*100),"%") # % CO2 of GHG from CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:CH4tax_reduc]/EmissionReductionResults_Mt[3,:CH4tax_reduc]*100, digits=0)),"%") # % CH4 of GHG from CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[1,:CH4tax_reduc]/EmissionReductionResults_Mt[3,:CH4tax_reduc], digits=2)*100),"%") # % Spillovers = CH4 from above
# println("\$",Int64(round(TaxRev[3,:Change_in_Rev])))
# println("\$",-1*Int64(round(value(RA)*(EqVar[3,:Ut_perc]*10^-2))))
# println("\$",-1*Int64(round(EVperton_CH4)),"/t")
# println(round(EqVar[3,:Ut_perc], digits=2),"%")
# # println("\$",-1*Int64(round(EVperton_CH4)),"/t")

# # #### Full table of results from both/bocmbined
# println("\$",Int64(round(value(CO2_taxrate),digits=0)),"/t")
# println("\$",Int64(round(value(CH4_taxrate),digits=0)),"/t")
# # For SCC
# # println(Int64(round(EmissionReductionResults_Mt[3,:both_taxes_combined],digits=0)))
# # println("(",Int64(round(EmissionReductionResults_Mt[3,:both_taxes_combined]/Emissions_Mt[3,:Bnchmrk_Emissions]*100,digits=0)),"%)") # % CO2 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[1,:both_taxes_combined],digits=0)))
# println("(",Int64(round(EmissionReductionResults_Mt[1,:both_taxes_combined]/Emissions_Mt[1,:Bnchmrk_Emissions]*100,digits=0)),"%)") # % CO2 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:both_taxes_combined],digits=0)))
# println("(",Int64(round(EmissionReductionResults_Mt[2,:both_taxes_combined]/Emissions_Mt[2,:Bnchmrk_Emissions],digits=2)*100),"%)") # % CH4 of CH4 tax
# println(Int64(round(ReductTarget, digits=0)))
# println(Int64(round(EmissionReductionResults_Mt[1,:both_taxes_combined]/EmissionReductionResults_Mt[3,:both_taxes_combined], digits=2)*100),"%") # % CH4 of CH4 tax
# println(Int64(round(EmissionReductionResults_Mt[2,:both_taxes_combined]/EmissionReductionResults_Mt[3,:both_taxes_combined]*100, digits=0)),"%") # % CH4 of CH4 tax
# println("-")#println(Int64(round(EmissionReductionResults_Mt[1,:both_taxes_combined]/EmissionReductionResults_Mt[3,:both_taxes_combined], digits=2)*100),"%") # % CH4 of CH4 tax
# println("\$",Int64(round(TaxRev[5,:Change_in_Rev])))
# println("\$",-1*Int64(round(value(RA)*(EqVar[5,:Ut_perc]*10^-2))))
# println("\$",-1*Int64(round(EVperton_both)),"/t")
# println(round(EqVar[5,:Ut_perc], digits=2),"%")
# println("\$",-1*Int64(round(EVperton_both)),"/t")

# For alternative specification table
# println(Int64(round(TotGHGbnchmk*10^3, digits=0)))
# println(Int64(round(TotCH4bnchmk/TotGHGbnchmk, digits=2)*100),"%")
# println("\$",Int64(round(value(CO2_taxrate),digits=0)),"/t")
# println("\$",Int64(round(value(CH4_taxrate),digits=0)),"/t")
# println("\$",Int64(round(TaxRev[5,:Change_in_Rev])))
# println(round(EqVar[5,:Ut_perc], digits=2),"%")

# realCH4taxpc_WiNcntfac
# realCH4taxpc_ch4
# realCH4taxpc_co2
# realCH4taxpc_both