# Adapted from WiNDC National Model
using MPSGE
using DataFrames, JLD2
using JuMP
using MPSGE.JuMP.Containers
using CSV, Plots
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S
Sectors = CSV.read(joinpath(@__DIR__,"./Sectors.csv"), DataFrame);
## New sectors:
# uel (utility: Electric power generation, transmission, and distribution), 
# ugs (utility: Natural gas distribution), 
# uwt (utility: Water, sewage and other systems), 
# coa (coal mining)
# min (other non-oil/gas mining [not coal])
### min => coa for CH4sectors, oil => oil & gas.
### emissions and abatement for gas split btw 'gas' and 'pip'?

I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
J = [i for i∈S[:j] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
CH4sectors = [:agr,:min,:pip,:oil,:wst] #:uti? # subset index for relevant CH4 mitigation sectors (VA slack in benchmark)
VA = [va for va∈S[:va] if va!=:othtax] # Index Value Added (compen = returns to labour/wage, 'surplus' = returns to Kapital)
FD = S[:fd] # Index for final demand categories (pce, and investment types)
TS = S[:ts] #index for taxes/subsidies
YR = S[:yr] # Index for years for potential multi year runs
M = S[:m] # Index for margins (transport and trade)
a_0 = P[:a_0] #	    "Armington supply",
id_0 = P[:id_0] #	"Intermediate demand",

ys_0 = P[:ys_0]#	"Sectoral supply",
va_0 = P[:va_0] #	"Value added",
md_0 = P[:md_0] #	"Margin demand",
fd_0 = P[:fd_0] #	"Final demand",
m_0 = P[:m_0] #	    "Imports",
ms_0 = P[:ms_0] #	"Margin supply",
bopdef_0 = P[:bopdef_0] #	"Balance of payments deficit",
x_0 = P[:x_0] #	    "Exports of goods and services",
y_0 = P[:y_0]  #	"Gross output",
## Existing Taxes
ty_0 = P[:ty_0] #	"Output tax rate"
tm_0 = P[:tm_0] #	"Import tariff"; Initial, for price 
ta_0 = P[:ta_0] #	"Tax net subsidy rate on intermediate demand", benchmark data also for price level

yr = Symbol(2020)

## Base Marginal Abatement Cost EPA data (2020). $/t and MMtCO2eq
MAC_CH4_data=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2020_data.csv"), DataFrame, header=2, limit=14)
MAC_CH4_totemiss=CSV.read(joinpath(@__DIR__,"./data/EPA_CH4_MAC_2020_data.csv"), DataFrame, header=2, skipto=17)
# % split for pipelines and oil from MAC for 'GAS' by proportion of combined economic output
pip_of_GAS = sum(ys_0[yr,:pip,:])/((sum(ys_0[yr,:pip,:])+sum(ys_0[yr,:oil,:])))
oil_of_GAS = sum(ys_0[yr,:oil,:])/((sum(ys_0[yr,:pip,:])+sum(ys_0[yr,:oil,:])))
# Aggregate/disaggregate for WiNDC sectors, $/t and MMtCO2eq 
MAC_CH4_WiNDC=DataFrame([MAC_CH4_data[:,:cost_per_t], MAC_CH4_data[:,:agr_livestock]+MAC_CH4_data[:,:agr_rice],
 MAC_CH4_data[:,:min], 
 MAC_CH4_data[:,:GAS]*pip_of_GAS,
 MAC_CH4_data[:,:GAS]*oil_of_GAS,
 MAC_CH4_data[:,:wst_land]+MAC_CH4_data[:,:wst_water]],
 [:cost_per_t; CH4sectors])
# Aggregate/disaggregate Total CH4 Emissions for WiNDC sectors, and include total VA for convenience, MMt and BillUS$
MAC_CH4_WiNDC_tot_MMt = DataFrame([MAC_CH4_totemiss[:,:cost_per_t], MAC_CH4_totemiss[:,:agr_livestock]+MAC_CH4_totemiss[:,:agr_rice],
 MAC_CH4_totemiss[:,:min], 
 MAC_CH4_totemiss[:,:GAS]*pip_of_GAS,
 MAC_CH4_totemiss[:,:GAS]*oil_of_GAS,
 MAC_CH4_totemiss[:,:wst_land]+MAC_CH4_totemiss[:,:wst_water]],
 [:cost_per_t; CH4sectors])
 MAC_CH4_WiNDC_tot = copy(MAC_CH4_WiNDC_tot_MMt); MAC_CH4_WiNDC_tot[:, 2:end] = MAC_CH4_WiNDC_tot_MMt[:,2:end].*10^-3
 push!(MAC_CH4_WiNDC_tot, ["Total va_cost"; [sum(va_0[yr,:,Symbol(sector)]) for sector in names(MAC_CH4_WiNDC[:,2:end])]])

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
    VASInt = DenseAxisArray(fill(0.,length(J)),J); for c in CH4sectors; VASInt[c] =  MAC_CH4_WiNDC_tot[1,c]/MAC_CH4_WiNDC_tot[2,c] end
## CH4 emission intensity (emissions/total va cost) for each sector at each abatement tier, subtracts cumulatively abated emissions and adds additional cost of abatement activities MMt/BillUS$
    VAM_CH4EmInt = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Initializing DenseAxisArray
    [VAM_CH4EmInt[v,c] = (MAC_CH4_WiNDC_tot[1,c]-CH4_EmMitigated[v,c])/(sum(va_0[yr,:,c])+CH4_cumul_cost[v,c]) for v in VAMset for c in CH4sectors]
## Relative cost of VA: (standard va cost + mitigation cost)/(standard cost) - used to multiply BOTH the va[:surplus] and va[:compen] equally in the blocks
    VAM_costover = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Initialise DenseAxisArray
    [VAM_costover[cost,c] = (sum(va_0[yr,:,c])+CH4_cumul_cost[cost,c])/sum(va_0[yr,:,c]) for cost in VAMset for c in CH4sectors]


CO2Int = DenseAxisArray(zeros(length(J)),J)
# 2024 GHG Inv Table 3-5: CO2 Emissions from Fossil Fuel Combustion by Fuel Type and Sector (MMT CO2 Eq.)
TotalCO2EmMMt_coal = 835.6 # 2020 2022=>895.9 # EPA Inventory CO2 Stationary Combustion - Coal sum(Electricity, Industrial, Commercial, & Residential=0)
# Option with no transport emissions in the model: assumption that direct
# TotalCO2EmMMt_gas_oil = 2104.80
Natural_gasCO2 = 1615.7 #2020 (incl 58.7 transport)
PetroleumCO2 = 1890.0 # 2020 (incl 1,514.2 transport)
# Option with all Transport CO2 emissions attributed to oil inputs: assumption that forms of oil fuel all transport that has direct CO2 emissions, and so taxing CO2 is total emissions from all oil as an input 
TotalCO2EmMMt_gas_oil = Natural_gasCO2 + PetroleumCO2 # EPA inventory all CO2 transport + CO2 Stationary Combustion - both Oil & Natural Gas sum(Electricity, Industrial, Commercial, Residential & oil U.S. Territories)Sta
# (MMtCO2eq/$Bill x 10-^3 -> Gt/$B=t/$) EPA Inventory 2022 sum of CO2 MMt for coal, and for gas & oil per Billion of total benchmark intermediate input from sector
TotalCO2EmGt_coal = TotalCO2EmMMt_coal*10^-3
TotalCO2EmGt_gas_oil = TotalCO2EmMMt_gas_oil*10^-3
CO2Int[:min] =  TotalCO2EmGt_coal/sum(id_0[yr,:min,:]) 
CO2Int[:oil] =   TotalCO2EmGt_gas_oil/sum(id_0[yr,:oil,:])  
TotCO2bnchmk =  TotalCO2EmGt_coal + TotalCO2EmGt_gas_oil
# TotCH4bnchmk = sum(CH4emiss[:EPAemiss,:]) # from single step data version, same data, same value
TotCH4bnchmk = sum(MAC_CH4_WiNDC_tot[1,2:end])
TotGHGbnchmk =  TotCO2bnchmk + TotCH4bnchmk


# Crude oil barrels production 2020: 4144184x10^3
# Average price per barrel 2020: $36.86
value_of_oil_2020 = 36.86 * 4144184*10^3
# Natural gas markets production 2020: 36,520,826 x 10^6 ft^3  (total withdrawals: 40,729,927 x 10^6 ft^3)
# Average natural gas spot price 2020 #3.32/thousand ft^3 $2.03 / million Btu
value_of_gas_2020 = 3.32 * 36520826 * 10^3
oil_fraction = value_of_oil_2020/(value_of_gas_2020+value_of_oil_2020)
# imports of crude oil 2020: 2150267 X10^3 barrel ; imports of oil products 727623 x 10^3 barrels
# imports of natural gas 2929: 2.55 trillion ft^3 ; 5.25 trillion ft^3
## End data preparations

## Set tax rates
CO2_taxrate = 190 # SC CO2 EPA 2023 SCGHG report, 2020 year, 2020US$, central 2% discount rate
CH4_taxrate = 190 # using SC CO2 because CH4 data is in MtCO2eq
# CH4abatement="yes" # Comment out CH4abatement="no" to allow CH4 abatment
CH4abatement="no" # Until there's also CO2 abatemment, no CH4 abatement by default
# Kmobile="yes" # Allow kapital & Labor to flow between sectors (original WiNDC)
Kmobile="no" # Fix kapital in sectors, allow Labor to flow between
print("CO2tax: $CO2_taxrate, "); println("CH4tax: $CH4_taxrate")
print("$CH4abatement CH4 Abatement: "); print("$Kmobile mobile Kapital:: ")


## Upload table of elasticity parameter values drawn from SAGE 2.1.1 documentation and E3 book, with manual concordence
Elasdf=CSV.read(joinpath(@__DIR__,"./data/Elasticities_SAGE_E3.csv"), DataFrame, header=1)
Elas=DenseAxisArray(Matrix(Elasdf[:,2:end]),Symbol.(Elasdf[:,1]),Symbol.(names(Elasdf)[2:end]))
# Test to check significance of specific super high Armington elasticity for oil and gas by setting back to WiNDC 2 
#TODO # Need to get a better reference to determine which end is appropriate (or in the middle) 
# Elas[:oil,:SAGE_nf_Armington] = 2 ;println("ArmingtonOil=2")
# Elas[:gas,:SAGE_nf_Armington] # Gas is not yet disaggregated

# function Multiloop(CO2_taxrate,CH4_taxrate) 

MultiNat = MPSGEModel()

@parameters(MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    CH4_tax, 0.
    CO2_tax, 0.
    t_elas_y, 0.
    elas_y, 0.05
    t_elas_a, 2.            
    elas_a, 0.
    elas_dm, 2.
    d_elas_ra, 1.
end)

@sectors(MultiNat,begin
    Y[J],      (description = "Sectoral Production",)
    A[I],      (description = "Armington Supply",)
    VAS[J],    (description = "Value Added, standard")
    MS[M],     (description = "Margin Supply",)
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
    PVAK[J], (description = "Kapital Input to VA blocks",) # separate from L and fix within sector
    PVAL, (description = "Labour Input to VA blocks",) # separate from K, and mobile btw sectors (1 price)
    end)

end

@commodities(MultiNat,begin
    PA[I],   (description = "Armington Price")
    PY[J],   (description = "Supply",)
    PVAM[J], (description = "Value-added output - Input to Y",)
    PM[M],   (description = "Margin Price",)
    PFX,     (description = "Foreign Exachange",)
end)

# Variables to track and report levels of CO2 emissions
@auxiliary(MultiNat, CO2em, index = [[:min, :oil]])
@auxiliary(MultiNat, CO2TotEm, description = "Total CO2 emissions from fossil fuels")
# Variables to track and report levels of CH4 emissions
@auxiliary(MultiNat, CH4em, index = [[:agr,:min,:oil,:pip,:wst]])
@auxiliary(MultiNat, CH4TotEm, description = "Total CH4 emissions")
@auxiliary(MultiNat, TotEm, description = "Total both emissions")


@consumer(MultiNat, RA, description = "Representative Agent")
"""
ID = I = [i for i∈I if i∉[:oil,:min, :uti]] # separate oilandgas and min, i without oil or min AND separate uti, so i without oil, min, or uti.
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
# # println("ElasY=E3m")
# for j∈J
#         @production(MultiNat, Y[j], [t= 0, s=Elas[j, :SAGE_klem_Y], va=> s=Elas[j,:SAGE_kle_VAE],sm => s = Elas[j,:E3_m_ID]],begin # [t=0, s = 0, sv=> s = 0]
#         [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
#         [@input(PA[i], id_0[yr,i,j], sm, taxes = [Tax(RA,CO2_tax * CO2Int[i])]) for i∈I]...
#         @input(PVAM[j], sum(va_0[yr,VA,j]), va)
#         end)
# end

# Test of re-nesting within Y block without disaggregation
ID = [i for i ∈ I if i∉[:oil, :min] ] # Intermediate inputs EXCEPT oil and min
for j∈J
    @production(MultiNat, Y[j], [t= 0, s=Elas[j, :SAGE_klem_Y], vae=>s=Elas[j,:SAGE_kle_VAE],sm => s = Elas[j,:E3_m_ID],
    va=>vae=0, oilmin=> vae=1
    ],begin
    [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
    # take tax out here bc it will go in lower nests
    [@input(PA[i], id_0[yr,i,j], sm) for i∈ID]... # elasticity btw inputs in vector 
    @input(PVAM[j], sum(va_0[yr,VA,j]), va)
    [@input(PA[f], id_0[yr,f,j], oilmin, taxes = [Tax(RA,CO2_tax * CO2Int[f])]) for f in [:oil, :min]]...
end)
end

## Draft set up for full nesting after disaggregation
# ID = [i for i ∈ I if i∉[:pet, :oil, :min, :ele, :rnw, :gas, :col]
# for j∈J
#     # @production(MultiNat, Y[j], [t= 0, s=Elas[j, :SAGE_klem_Y], vae=> s=Elas[j,:SAGE_kle_VAE],sm => s = Elas[j,:E3_m_ID],
#     # va=>vae=   , e=> vae = [j], el => e = [j, value irrelevant just linking nests], oilgas=> e=(btw oil and gas energy),f => e = (part gas & oil), fe=>el= (coal & gas)
#     # ],begin
#     [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
#     [@input(PA[i], id_0[yr,i,j], sm, taxes = [Tax(RA,CO2_tax * CO2Int[i])]) for i∈I]...
    
#     # take tax out here bc it will go in lower nests
#     # [@input(PA[i], id_0[yr,i,j], sm) for i∈ID]... # elasticity btw inputs in vector 
#     @input(PVAM[j], sum(va_0[yr,VA,j]), va)
#     # @input(PA[:uti], sum(electricity), el)
#     # [@input(PA[f], sum(oil and part gas), oilgas) for f in [:oil, :gas]]...
#     # @input(PA[:rnw], sum(benchmark renewables), rel)
#     # [@input(PA[ff], sum(benchmark renewables), rel) for ff in [:gas, :col]]...
# end)
# end


# # One price of electricity, based on how all firms get energy. Taking electricity out of id_o and ys_0. And gas, coal, and oil out. 
# begin
#     @output(PEL, electriticy for all firms (sum elect for all firms ))for i in I...
#     @input(PY[r], renewable energy from Y) # somehow I need to get some kind of estimate of renewable generation separated out from uti->elec (it's under it in NAICS, but no data in BEA: 22111 Electric Power Generation -> 
#     # {221111	Hydroelectric Power Generation, 221112	Fossil Fuel Electric Power Generation, 221113	Nuclear Electric Power Generation, 221114	Solar Electric Power Generation, 221115	Wind Electric Power Generation, 221116	Geothermal Electric Power Generation, 221117	Biomass Electric Power Generation, 221118	Other Electric Power Generation}
#      @input(PY[gc], gas and coal input from Y, for electricity, with tax )
#     @input(PVAM) # pull out uti from VA too.
#     [@input(PA[i] id_0[el]) for i in I]... #material/services for electricity
# end


# Total value added cost as a function labor (compen) and kapital (surplus), standard (no mitigation)
print("ElasVA = SAGEkl:: ")
if (Kmobile=="yes")
for j∈J
        @production(MultiNat, VAS[j], [t=0, s = 0, va => s = Elas[j,:SAGE_kl_VA]], begin # #     @production(MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin 
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], va_0[yr,va,j], va, taxes = [Tax(RA,CH4_tax* VASInt[j])]) for va∈VA]...
        end)
end
elseif (Kmobile=="no")
    for j∈J
        @production(MultiNat, VAS[j], [t=0, s = 0, va => s = Elas[j,:SAGE_kl_VA]], begin # #     @production(MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin 
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        @input(PVAK[j], va_0[yr,:surplus,j], va, taxes = [Tax(RA,CH4_tax* VASInt[j])])
        @input(PVAL, va_0[yr,:compen,j], va, taxes = [Tax(RA,CH4_tax* VASInt[j])])
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
                        [@output(PVAM[c],sum(va_0[yr,:,c]), t)]... 
                        [@input(PVA[va], va_0[yr,va,c]*VAM_costover[vam.name,c], va, taxes = [Tax(RA, CH4_tax*VAM_CH4EmInt[vam.name,c])]) for va∈VA]...
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
                        [@output(PVAM[c],sum(va_0[yr,:,c]), t)]... 
                        @input(PVAK[c],va_0[yr,:surplus,c]*VAM_costover[vam.name,c], va, taxes = [Tax(RA, CH4_tax*VAM_CH4EmInt[vam.name,c])])
                        @input(PVAL,    va_0[yr,:compen,c]*VAM_costover[vam.name,c], va, taxes = [Tax(RA, CH4_tax*VAM_CH4EmInt[vam.name,c])])
                    end)
                end
            end
        end
    end
end

for m∈M
    @production(MultiNat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end
 
for i∈I
    @production(MultiNat, A[i], [t = t_elas_a, s = elas_a, dm => s = elas_dm], begin
    # @production(MultiNat, A[i], [t = 2, s = 0, dm => s = 2], begin
# println("elasdm=SAGEnf")
# for i∈I
#     @production(MultiNat, A[i], [t = 2, s = 0, dm => s = Elas[i,:SAGE_nf_Armington]], begin
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i])],reference_price=1-ta_0[yr,i])]... 
        [@output(PFX, x_0[yr,i], t)]...
        [@input(PM[m], md_0[yr,m,i], s) for m∈M]...
        @input(PY[i], y_0[yr,i], dm)
## Emissions tariff on CH4 goods because inputs
#         @input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i]),Tax(RA,CH4_tax* VASInt[i])],reference_price=1+tm_0[yr,i])# No tariff on CO2 bc oil and gas are taxed as inputs to production, which includes these imports: Tax(RA,CO2_tax * CO2Int[i]),
#     end)
# end; # println("Yes CH4 tariff: ")  
## Alternative, no tariff on CH4 goods        
@input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[yr,i]) # without excise tariff CH4 goods (or oil, 'coal' i.e. 'min
    end)
end; println("No CH4 tariff:: ")

if (Kmobile=="yes")
    @demand(MultiNat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVA[va], sum(va_0[yr,va,j] for j∈J)) for va∈VA]...
    end, elasticity = 1)
    # end, elasticity = d_elas_ra)
elseif (Kmobile=="no")
    @demand(MultiNat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVAK[j], sum(va_0[yr,:surplus,j])) for j∈J]...
    [@endowment(PVAL, sum(va_0[yr,:compen,j])) for j∈J]...
    end, elasticity = 1)
    # end, elasticity = d_elas_ra)
end

## CO2 emissions for fossil fuel sectors are the activity levels times the (base) total emissions intensity 
@aux_constraint(MultiNat, CO2em[:min],  CO2em[:min] - Y[:min]*TotalCO2EmGt_coal)
@aux_constraint(MultiNat, CO2em[:oil],  CO2em[:oil] - Y[:oil]*TotalCO2EmGt_gas_oil)
## Total CO2 emissions are the sum of emissions from the 2 fossil fuel sectors (constraint expressed as equantion = 0)
@aux_constraint(MultiNat, CO2TotEm, CO2TotEm - (CO2em[:min] + CO2em[:oil]))
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
@aux_constraint(MultiNat, CH4TotEm, CH4TotEm - (CH4em[:agr] + CH4em[:min] + CH4em[:oil] + CH4em[:pip] + CH4em[:wst] ))
## Total GHG (CO2 & CH4) emissions in Mt CO2eq
@aux_constraint(MultiNat, TotEm, TotEm - (CH4TotEm + CO2TotEm))
set_silent(MultiNat)

# Benchmark 
# fix(PA[:oil],1)
# fix(PA[:rec], 1)
# fix(RA, sum(fd_0[yr,I,:pce])) # Numeraire, fixed at benchmark
## Note: Benchmark doesn't solve at 0 interation because of margins of slack activity. Does balance with interactions or slack vars and production commented out.
solve!(MultiNat , cumulative_iteration_limit = 0)

fullvrbnch = generate_report(MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
# print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))
# fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var])
# fullvrch4[!,:var] = Symbol.(fullvrch4[:,:var])
# filter(x -> x.var in [Symbol("Y[oil]"), Symbol("Y[pet]"), Symbol("Y[min]"), Symbol("Y[agr]"),
#  Symbol("Y[wst]"), Symbol("PVA[compen]"), Symbol("PVA[surplus]"),
#  Symbol("PVAM[agr]"), Symbol("PVAM[oil]"), Symbol("PVAM[ppd]")], fullvrbnch)

#  filter(x -> x.var in [Symbol("Y[oil]"), Symbol("Y[pet]"), Symbol("Y[min]"), Symbol("Y[agr]"),
#  Symbol("Y[wst]"), Symbol("PVA[compen]"), Symbol("PVA[surplus]"),
#  Symbol("PVAM[agr]"), Symbol("PVAM[oil]"), Symbol("PVAM[ppd]")], fullvrch4)#fullvrbnch)

# Initialize a Dataframe to save final demand results
FDemand = DataFrame(index=Vector{Symbol}(undef, length(I)),
desc=Vector{Symbol}(undef, length(I)), 
bnch=Vector{Float64}(undef, length(I)), 
WiNcntfact=Vector{Float64}(undef, length(I)), 
ch4tax=Vector{Float64}(undef, length(I)),
co2tax=Vector{Float64}(undef, length(I)),
both=Vector{Float64}(undef, length(I)),
ch4Qdelta=Vector{Float64}(undef, length(I)),
CO2Qdelta=Vector{Float64}(undef, length(I)),
bothQdelta=Vector{Float64}(undef, length(I)),
bncQpc=Vector{Float64}(undef, length(I)),
ch4Qpc=Vector{Float64}(undef, length(I)),
CO2Qpc=Vector{Float64}(undef, length(I)),
bothQpc=Vector{Float64}(undef, length(I))
)
for (n,i) in enumerate(I)
    FDemand[n,:index]= i
    FDemand[n,:desc] = Symbol(Sectors[Sectors.index.==string(i),2][1])
    FDemand[n,:bnch] = value(demand(RA,PA[i]))
end
# unfix(RA)
## Check against WiNDC standard counterfactual
set_value!(ta,0)
set_value!(tm,0)

solve!(MultiNat)

fullvrWiNcntfact = generate_report(MultiNat)
rename!(fullvrWiNcntfact, :value => :WiNcntfact, :margin => :WiNcntfactmarg)

for (n,i) in enumerate(I)
    FDemand[n,:WiNcntfact] = value(demand(RA,PA[i]))
end

## DataFrame to hold the industry indices with descriptions
Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

## Re-set to benchmark taxes
set_value!(ta,ta_0[yr,J])
set_value!(tm,tm_0[yr,J])
solve!(MultiNat , cumulative_iteration_limit = 0)# Temp measure to address residual price changes

# tax are at $/t of CH4(CO2eq)
## "EPA SC CH4 is $1600/t. 
set_value!(CH4_tax, CH4_taxrate)
set_value!(CO2_tax,0.) # Set CO2 tax to 0 for running separately.
solve!(MultiNat)

for (n,i) in enumerate(I)
    FDemand[n,:ch4tax] = value(demand(RA,PA[i]))
end

Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

fullvrch4 = generate_report(MultiNat)
rename!(fullvrch4, :value => :ch4tax, :margin => :ch4marg)

# print(sort(df, :margin, by= abs, rev=true))
# print(df)
solve!(MultiNat , cumulative_iteration_limit = 0) # Temp measure to address residual price changes 
set_value!(CO2_tax, CO2_taxrate)
set_value!(CH4_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
solve!(MultiNat, cumulative_iteration_limit=10000) #;

for (n,i) in enumerate(I)
    FDemand[n,:co2tax] = value(demand(RA,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:4]
# Rs[:,2][68:71]

fullvrco2 = generate_report(MultiNat)
rename!(fullvrco2, :value => :CO2tax, :margin => :CO2marg)

# # Counterfactual Fossil fuel extraction is ALSO taxed at emissions intensitiy of input x tax in $/ton
set_value!(CO2_tax, CO2_taxrate)
set_value!(CH4_tax, CH4_taxrate)
solve!(MultiNat, cumulative_iteration_limit=10000) #;

for (n,i) in enumerate(I)
    FDemand[n,:both] = value(demand(RA,PA[i]))
end
# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:5]
# Rs[:,2][68:71]
fullvrboth = generate_report(MultiNat)
rename!(fullvrboth, :value => :bothtaxes, :margin => :bothmarg)

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
EmissionReductionResults = DataFrame(
["CO2" "Gt" TotCO2bnchmk TotCO2bnchmk - only(compCO2em[:,:CO2tax]) TotCO2bnchmk - only(compCO2em[:,:ch4tax]) TotCO2bnchmk - only(compCO2em[:,:ch4tax]) + TotCO2bnchmk - only(compCO2em[:,:CO2tax]) TotCO2bnchmk - only(compCO2em[:,:bothtaxes])  TotCO2bnchmk - only(compCO2em[:,:ch4tax]) + TotCO2bnchmk - only(compCO2em[:,:CO2tax]) - (TotCO2bnchmk - only(compCO2em[:,:bothtaxes])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk TotCH4bnchmk - only(compCH4em[:,:CO2tax]) TotCH4bnchmk - only(compCH4em[:,:ch4tax]) TotCH4bnchmk - only(compCH4em[:,:ch4tax]) + TotCH4bnchmk - only(compCH4em[:,:CO2tax]) TotCH4bnchmk - only(compCH4em[:,:bothtaxes]) TotCH4bnchmk - only(compCH4em[:,:ch4tax]) + TotCH4bnchmk - only(compCH4em[:,:CO2tax]) - (TotCH4bnchmk - only(compCH4em[:,:bothtaxes])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk TotGHGbnchmk - only(compTotem[:,:CO2tax]) TotGHGbnchmk - only(compTotem[:,:ch4tax]) TotGHGbnchmk - only(compTotem[:,:ch4tax]) + TotGHGbnchmk - only(compTotem[:,:CO2tax]) TotGHGbnchmk - only(compTotem[:,:bothtaxes]) TotGHGbnchmk - only(compTotem[:,:ch4tax]) + TotGHGbnchmk - only(compTotem[:,:CO2tax]) - (TotGHGbnchmk - only(compTotem[:,:bothtaxes]))], ["Emissions", "Unit", "Bnchmrk_Emissions", "CO2tax_reduc", "CH4tax_reduc","Sum_of_both_taxes", "taxes_combined" ,"Interactions"])

## Generate subset DataFrame with just the Value-Added activity for the emitting sectors, show those results
#filter(row -> row.var ∈ [Symbol("VAM[agr]"),Symbol("VAM[min]"),Symbol("VAM[pip]"),Symbol("VAS[oil]"),Symbol("VAM[oil]"),Symbol("VAS[min]"),Symbol("VAS[pip]"),Symbol("VAS[agr]"),Symbol("VAS[wst]"),Symbol("VAM[wst]"),], Compare)

function plottaxemisscurve(tax1, tax2, start, interval, finish, RAval, isfixed, cnst=1)
    """# runs a loop increasing each tax by \$1/t and then plotting Total GHG (CO2 & CH4) **incorporated** emissions 
    # Arguments are: which tax to change, other tax to either change simultaneously OR keep at 0, st=initial \$ tax value, fin= final \$ tax value,
    # and final (optional) argument can be set to 0 to remove other tax, by default """
    margemiss = DataFrame(tax=Float64[], Emissions=Float64[], CH4Emissions=Float64[],CO2Emissions=Float64[])
    Testvars = DataFrame(taxrt=Float64[], 
    # Yagr=Float64[],Ymin=Float64[],Ypip=Float64[],Yoil=Float64[],Apip=Float64[],Aoil=Float64[],Ywst=Float64[],
    # CompDYPApip=Float64[],CompDApipPApip=Float64[],DemRAPApip=Float64[],
    # PAagr=Float64[],PAmin=Float64[],PApip=Float64[],PAoil=Float64[],PAwst=Float64[],PAuti=Float64[],compDApipPAoil=Float64[],compdDAoilPAoil=Float64[],
    # VASagr=Float64[],VAMagr=Float64[],VASmin=Float64[],VAMmin=Float64[],VASpip=Float64[],VAMpip=Float64[],VASoil=Float64[],VAMoil=Float64[],VASwst=Float64[],VAMwst=Float64[],
    VASagr=Float64[],VAM10agr=Float64[],VAM100agr=Float64[],VAM500agr=Float64[],VAMkagr=Float64[],VASmin=Float64[],VAM10min=Float64[],VAM100min=Float64[],VAM500min=Float64[],VAMkmin=Float64[],VASpip=Float64[],VAM10pip=Float64[],VAM100pip=Float64[],VAM500pip=Float64[],VAMkpip=Float64[],VASoil=Float64[],VAM10oil=Float64[],VAM100oil=Float64[],VAM500oil=Float64[],VAMkoil=Float64[],VASwst=Float64[],VAM10wst=Float64[],VAM100wst=Float64[],VAM500wst=Float64[],VAMkwst=Float64[],
    TotEm=Float64[],CH4TotEm=Float64[],CO2TotEm=Float64[]
    # CH4emoil=Float64[],CH4empip=Float64[],CO2emin=Float64[],CO2emoil=Float64[]
    )
    ResultsTroubleshoot = DataFrame(var=[], value=Float64[], margin=Float64[], x1=Float64[]) 
    for i in start:interval:finish
        set_value!(tax1, i)
        set_value!(tax2, cnst*i)
        solve!(MultiNat, output="no");
        Results = generate_report(MultiNat)
        Results[!,:var] = Symbol.(Results[:,:var]);
        ResultsTroubleshoot =vcat(ResultsTroubleshoot, [Results fill(i,length(Results[:,1]))])
        push!(margemiss, [i only(filter(:var => ==(:TotEm), Results)[:, :value]) only(filter(:var => ==(:CH4TotEm), Results)[:, :value]) only(filter(:var => ==(:CO2TotEm), Results)[:, :value])])
        push!(Testvars, [i,                
        # value(Y[:agr]),value(Y[:min]),value(Y[:pip]),value(Y[:oil]),value(A[:pip]),value(A[:oil]),value(Y[:wst]),
        # value(compensated_demand(MultiNat[:Y][:pip],MultiNat[:PA][:pip])),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:pip])),value(demand(MultiNat[:RA],MultiNat[:PA][:pip])),
        # value(PA[:agr]),value(PA[:min]),value(PA[:pip]),value(PA[:oil]),value(PA[:wst]),value(PA[:uti]),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:oil])),value(compensated_demand(MultiNat[:A][:oil],MultiNat[:PA][:oil])),
        # value(VAS[:agr]),value(VAM[:agr]),value(VAS[:min]),value(VAM[:min]),value(VAS[:pip]),value(VAM[:pip]),value(VAS[:oil]),value(VAM[:oil]),value(VAS[:wst]),value(VAM[:wst]),
        value(VAS[:agr]),value(VAM10[:agr]),value(VAM100[:agr]),value(VAM500[:agr]),value(VAM1000[:agr]),value(VAS[:min]),value(VAM10[:min]),value(VAM100[:min]),value(VAM500[:min]),value(VAM1000[:min]),value(VAS[:pip]),value(VAM10[:pip]),value(VAM100[:pip]),value(VAM500[:pip]),value(VAM1000[:pip]),value(VAS[:oil]),value(VAM10[:oil]),value(VAM100[:oil]),value(VAM500[:oil]),value(VAM1000[:oil]),value(VAS[:wst]),value(VAM10[:wst]),value(VAM100[:wst]),value(VAM500[:wst]),value(VAM1000[:wst]),
        value(TotEm),value(CH4TotEm),value(CO2TotEm)
        # value(CH4em[:oil]),value(CH4em[:pip]),value(CO2em[:min]),value(CO2em[:oil])
        ] 
            )
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    plt = plot(margemiss[!,:tax], margemiss[!,:Emissions].*10^3, title= "RA:\$$RAval fxd:$isfixed", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,1600))
    plt = plot!([190], seriestype=:vline, label="SCCO2", ylim=(0,TotGHGbnchmk*10^3))
    plch4 = plot(margemiss[!,:tax], margemiss[!,:CH4Emissions].*10^3, label=false, title="CH4 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,1600))
    plch4 = plot!([190], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
    plco2 = plot(margemiss[!,:tax], margemiss[!,:CO2Emissions].*10^3, label=false, title="CO2 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,1600))
    plco2 = plot!([190], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
    plt3 = plot(margemiss[!,:tax], margemiss[!,:Emissions].*10^3, title= "Emissions with $tax1 $tax2in", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="\$/t", xlims=(0,1600), color=:black, linewidth=3)
    plt3 = plot!(margemiss[!,:tax], margemiss[!,:CH4Emissions].*10^3, label="CH4 Emissions", color=:green)
    plt3 = plot!(margemiss[!,:tax], margemiss[!,:CO2Emissions].*10^3, label="CO2 Emissions", color=:blue)
    plt3 = plot!([190], seriestype=:vline, label="SCCO2", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
    # plcdyp = plot(margemiss[!,:tax],Testvars[!,:CompDYPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(Y:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDYPApip]),maximum(Testvars[!,:CompDYPApip])), xlabel="$tax1 $tax2in \$/t")
    # plcdap = plot(margemiss[!,:tax],Testvars[!,:CompDApipPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(A:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDApipPApip]),maximum(Testvars[!,:CompDApipPApip])), xlabel="$tax1 $tax2in \$/t")
    # plfdrap = plot(margemiss[!,:tax],Testvars[!,:DemRAPApip], title= "RA:\$$RAval fxd:$isfixed", label="final_dem(RA,PA:pip)", ylim=(minimum(Testvars[!,:DemRAPApip]),maximum(Testvars[!,:DemRAPApip])), xlabel="$tax1 $tax2in \$/t")
    # pla = plot(margemiss[!,:tax],Testvars[!,:Yagr], label=false, title="Y:agr", ylim=(minimum(Testvars[!,:Yagr]),maximum(Testvars[!,:Yagr])), xlabel="$tax1 $tax2in \$/t")
    # plm = plot(margemiss[!,:tax],Testvars[!,:Ymin], title= "RA:\$$RAval fxd:$isfixed", label="Y:min", ylim=(minimum(Testvars[!,:Ymin]),maximum(Testvars[!,:Ymin])), xlabel="$tax1 $tax2in \$/t")
    # plp = plot(margemiss[!,:tax],Testvars[!,:Ypip], title= "RA:\$$RAval fxd:$isfixed", label="Y:pip", ylim=(minimum(Testvars[!,:Ypip]),maximum(Testvars[!,:Ypip])), xlabel="$tax1 $tax2in \$/t")
    # plo = plot(margemiss[!,:tax],Testvars[!,:Yoil], title= "RA:\$$RAval fxd:$isfixed", label="Y:oil", ylim=(minimum(Testvars[!,:Yoil]),maximum(Testvars[!,:Yoil])), xlabel="$tax1 $tax2in \$/t")
    # plw = plot(margemiss[!,:tax],Testvars[!,:Ywst], title= "RA:\$$RAval fxd:$isfixed", label="Y:wst", ylim=(minimum(Testvars[!,:Ywst]),maximum(Testvars[!,:Ywst])), xlabel="$tax1 $tax2in \$/t")
    # plpa = plot(margemiss[!,:tax],Testvars[!,:PAagr], title= "RA:\$$RAval fxd:$isfixed", label="PA:agr", ylim=(minimum(Testvars[!,:PAagr]),maximum(Testvars[!,:PAagr])), xlabel="$tax1 $tax2in \$/t")
    # plpm = plot(margemiss[!,:tax],Testvars[!,:PAmin], title= "RA:\$$RAval fxd:$isfixed", label="PA:min", ylim=(minimum(Testvars[!,:PAmin]),maximum(Testvars[!,:PAmin])), xlabel="$tax1 $tax2in \$/t")
    # plpp = plot(margemiss[!,:tax],Testvars[!,:PApip], title= "RA:\$$RAval fxd:$isfixed", label="PA:pip", ylim=(minimum(Testvars[!,:PApip]),maximum(Testvars[!,:PApip])), xlabel="$tax1 $tax2in \$/t")
    # plpo = plot(margemiss[!,:tax],Testvars[!,:PAoil], title= "RA:\$$RAval fxd:$isfixed", label="PA:oil", ylim=(minimum(Testvars[!,:PAoil]),maximum(Testvars[!,:PAoil])), xlabel="$tax1 $tax2in \$/t")
    # plpw = plot(margemiss[!,:tax],Testvars[!,:PAwst], title= "RA:\$$RAval fxd:$isfixed", label="PA:wst", ylim=(minimum(Testvars[!,:PAwst]),maximum(Testvars[!,:PAwst])), xlabel="$tax1 $tax2in \$/t")
    # plAp = plot(margemiss[!,:tax],Testvars[!,:Apip], title= "RA:\$$RAval fxd:$isfixed", label="A:pip", ylim=(minimum(Testvars[!,:Apip]),maximum(Testvars[!,:Apip])), xlabel="$tax1 $tax2in \$/t")
    # plAo = plot(margemiss[!,:tax],Testvars[!,:Aoil], title= "RA:\$$RAval fxd:$isfixed", label="A:oil", ylim=(minimum(Testvars[!,:Aoil]),maximum(Testvars[!,:Aoil])), xlabel="$tax1 $tax2in \$/t")     # Or label=false, title="price of oil commodity"
    return margemiss, plt, Testvars, ResultsTroubleshoot, plch4, plco2, plt3 #,pla, plm, plp, plo, plw, plpa, plpm, plpp, plpo, plpw, 
    # plcdyp, plcdap, plfdrap, plAp, plAo
end

## Add Quant diff and % of consumption columns for Final Demand report dataframe
FDemand[:,:ch4Qdelta]=FDemand[:,:ch4tax].-FDemand[:,:bnch]
FDemand[:,:CO2Qdelta]=FDemand[:,:co2tax].-FDemand[:,:bnch]
FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
FDemand[:,:ch4Qpc]=FDemand[:,:ch4tax]./sum(FDemand[:,:ch4tax])*100
FDemand[:,:CO2Qpc]=FDemand[:,:co2tax]./sum(FDemand[:,:co2tax])*100
FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

set_silent(MultiNat)
# checkch4CO2 = plottaxemisscurve(CH4_tax, CO2_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# checkch4CO2[7]
# png(checkch4CO2[7], "./Results/Bothtax-Allemiss")
# checkch4CO2[2]
# print("checkch4CO2",checkch4CO2[3])

# fix(RA, sum(fd_0[yr,I,:pce]))
set_upper_bound(MultiNat[:A][:pip], 10)
# # MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pip]))
# # fix(RA, 14008.668551652801)
# # unfix(RA)
# fix(RA,16030.7) # RA value at $190/t
# # checkCO2 = plottaxemisscurve(CO2_tax, CH4_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2 = plottaxemisscurve(CO2_tax, CH4_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2[2] # Total emissions
# # # println("PApip up lim =", upper_bound(MultiNat[:A][:pip]))

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
# checkch4CO2 = plottaxemisscurve(CH4_tax, CO2_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# # checkch4CO2[2]

# fix(RA,fix(RA,16030.7) # RA value at $190/t
# checkch4 = plottaxemisscurve(CH4_tax,CO2_tax, 0, 10, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkch4[2]
# checkch4[5]
# checkch4[6]
# checkch4[7]
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

println("Emissions (remaining, i.e. not mitigated)")
EmissUnits_mt = DataFrame();
EmissUnits_mt.Unit=["Mt"; "MtCO2eq"; "MtCO2eq"];
EmissionReductionResults_Mt =[EmissionReductionResults[:,1:1] EmissUnits_mt EmissionReductionResults[:,3:end].*10^3] 
Emissions_Mt = DataFrame(names(EmissionReductionResults_Mt)[1] => EmissionReductionResults_Mt[:,1],
names(EmissionReductionResults_Mt)[2] => EmissionReductionResults_Mt[:,2],
names(EmissionReductionResults_Mt)[3] => EmissionReductionResults_Mt[:,3],
chop(names(EmissionReductionResults_Mt)[4], tail=6) => EmissionReductionResults_Mt[:,3]-EmissionReductionResults_Mt[:,4],
chop(names(EmissionReductionResults_Mt)[5], tail=6) => EmissionReductionResults_Mt[:,3]-EmissionReductionResults_Mt[:,5],
names(EmissionReductionResults_Mt)[6] => EmissionReductionResults_Mt[:,3]-EmissionReductionResults_Mt[:,6],
names(EmissionReductionResults_Mt)[7] => EmissionReductionResults_Mt[:,3]-EmissionReductionResults_Mt[:,7])

# return Emissions_Mt, 
println(Emissions_Mt)
println(
filter(x -> x.var in [Symbol("Y[oil]"), Symbol("Y[pet]"), Symbol("Y[min]"), Symbol("Y[agr]")
, Symbol("Y[wst]"), Symbol("PVA[compen]"), Symbol("PVA[surplus]"), Symbol("RA")], Compare)
,
);println("RA check, it's the 'both' bc that's the last simulation:",
sum([value(demand(RA,PA[i]))*value(PA[i]) for i in I])

)
#Check Benchmark
# fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var]);
# print(sort!(fullvrbnch, [:bmkmarg]))
# end

# product (1/PA[i])^demand(RA,PA[i]/sum(demand(RA,PA[i]}* RA? - )))