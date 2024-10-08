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
Sectors = CSV.read("Sectors.csv", DataFrame);

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

yr = Symbol(2017)

## Base Marginal Abatement Cost EPA data (2020)
MAC_CH4_data=CSV.read("./data/EPA_CH4_MAC_2020_data.csv", DataFrame, header=2, limit=14)
MAC_CH4_totemiss=CSV.read("./data/EPA_CH4_MAC_2020_data.csv", DataFrame, header=2, skipto=17)
# % split for pipelines and oil from MAC for 'GAS' by proportion of combined economic output
pip_of_GAS = sum(ys_0[yr,:pip,:])/((sum(ys_0[yr,:pip,:])+sum(ys_0[yr,:oil,:])))
oil_of_GAS = sum(ys_0[yr,:oil,:])/((sum(ys_0[yr,:pip,:])+sum(ys_0[yr,:oil,:])))
# Aggregate/disaggregate for WiNDC sectors 
MAC_CH4_WiNDC=DataFrame([MAC_CH4_data[:,:cost_per_t], MAC_CH4_data[:,:agr_livestock]+MAC_CH4_data[:,:agr_rice],
 MAC_CH4_data[:,:min], 
 MAC_CH4_data[:,:GAS]*pip_of_GAS,
 MAC_CH4_data[:,:GAS]*oil_of_GAS,
 MAC_CH4_data[:,:wst_land]+MAC_CH4_data[:,:wst_water]],
 [:cost_per_t; CH4sectors])

MAC_CH4_WiNDC_tot=DataFrame([MAC_CH4_totemiss[:,:cost_per_t], MAC_CH4_totemiss[:,:agr_livestock]+MAC_CH4_totemiss[:,:agr_rice],
 MAC_CH4_totemiss[:,:min], 
 MAC_CH4_totemiss[:,:GAS]*pip_of_GAS,
 MAC_CH4_totemiss[:,:GAS]*oil_of_GAS,
 MAC_CH4_totemiss[:,:wst_land]+MAC_CH4_totemiss[:,:wst_water]],
 [:cost_per_t; CH4sectors])
 push!(MAC_CH4_WiNDC_tot, ["Total va_cost"; [sum(va_0[yr,:,Symbol(sector)]) for sector in names(MAC_CH4_WiNDC[:,2:end])]])

# Initialise df with 0s to fill with calculations
CH4_cost_per_tier = copy(MAC_CH4_WiNDC[:,2:end]); CH4_cost_per_tier[:,:].=0
for c in CH4sectors
    CH4_cost_per_tier[1,c] = MAC_CH4_WiNDC[1,:cost_per_t]*MAC_CH4_WiNDC[1,c]
end
# Calculate Costs per tier of marginal abatement
for c in CH4sectors
    for i in 2:length(CH4_cost_per_tier[:,1])
        CH4_cost_per_tier[i,c] = MAC_CH4_WiNDC[i,:cost_per_t]*(MAC_CH4_WiNDC[i,c]-MAC_CH4_WiNDC[i-1,c])
    end
end
# Calculate the cost for cumulative abatement at eavh level
CH4_cumul_costAll = copy(MAC_CH4_WiNDC[:,2:end]); CH4_cumul_costAll.=0
# Set up 1st Row as just tier costs
CH4_cumul_costAll[1,:]=copy(CH4_cost_per_tier[1,:])
# Rest of dataframe sums cost of each tier with all tiers below
for c in CH4sectors
    for i in 2:length(CH4_cumul_costAll[:,1])
        CH4_cumul_costAll[i,c] = CH4_cost_per_tier[i,c]+CH4_cumul_costAll[i-1,c]
    end
end
VAMset = [:VAM5,:VAM10,:VAM15,:VAM20,:VAM30,:VAM40,:VAM50,:VAM100,:VAM500,:VAM1000]
# Filter to cumulative Marginal costs are positive for at least one sector  
CH4_cumul_cost = DenseAxisArray(Matrix(CH4_cumul_costAll[5:end,:]),VAMset, CH4sectors)
## CH4 Emissions intensity of Standard value-added activities 
VASInt =  [MAC_CH4_WiNDC_tot[1,c]*10^-3/MAC_CH4_WiNDC_tot[2,c] for c in CH4sectors]
## Cumulative abated emissions at each level of abatement cost, filtered to +$5/t and up only                        
CH4_EmMitigated = DenseAxisArray(Matrix(MAC_CH4_WiNDC[5:end,2:end]), VAMset, CH4sectors)
VAM_CH4EmInt = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Just initializing DenseAxisArray
## CH4 emission intensity (emissions/total va cost) for each sector at each abatement tier, subtracts cumulatively abated emissions and adds additional cost of abatement activities
[VAM_CH4EmInt[v,c] = (MAC_CH4_WiNDC_tot[1,c]*10^-3-CH4_EmMitigated[v,c]*10^-3)/(sum(va_0[yr,:,c])+CH4_cumul_cost[v,c]*10^-3) for v in VAMset for c in CH4sectors]
## Initialise DenseAxisArray to fill
VAM_costover = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors)
## Relative cost of VA including max mitigation
[VAM_costover[cost,c] = (sum(va_0[yr,:,c])+CH4_cumul_cost[cost,c]*10^-3)/sum(va_0[yr,:,c]) for cost in VAMset for c in CH4sectors]

#####--------Single Mitigation step data set up-----------#####
## Base Data for reference
    # EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDs : EPAnonCO2-report-data-annex-9-30-19_0\NonCO2 MACs-Summary Tables.xlsx
    CH4emissdatadf = DataFrame(Wsector = [:agr,:agr,:min,:oil,:wst,:wst],
    ## EPA Total Emissions per sector, MMt CO2eq 
    EPAemiss =[260.483532*10^-3,13.70952225*10^-3,59.31302643*10^-3,224.8979059*10^-3,111.5049515*10^-3,20.36144996*10^-3],
    ## EPA maximum % of abatement per sector at <$1000/t
    MaxpercMit = [.304,.280,.645,.475,.050,.350],
    PAsector = ["AGRICULTURE, LIVE","AGRICULTURE, RICE","ENERGY, COL",
    "ENERGY, GAS","WASTE, LAN","WASTE, WWR"])
## Calculate CH4 Intensity factors from EPA data
# EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDs
# Emissions MMt, sum and weighted average with current aggregation (before potential disaggregation of BEA/WiNDC data at some point)
CH4emiss = DenseAxisArray([260.483532*10^-3+13.70952225*10^-3 59.31302643*10^-3 224.8979059*10^-3*pip_of_GAS 224.8979059*10^-3*oil_of_GAS 111.5049515*10^-3+20.36144996*10^-3
## Weighted average mitigated potential per sector
(260.483532*0.3039385+13.70952225*0.2803694)/(260.483532+13.70952225)  0.6452054 0.4749547 0.4749547 (111.5049515*0.049615808+20.36144996*0.122234272)/(111.5049515+20.36144996)
    # ## Typo! 300.85 in denominator (260.483532*.304+13.70952225*.280)/(300.8535461+260.483532+13.70952225)  .645 .475 .475 (111.5049515*.050+20.36144996*.350)/(111.5049515+20.36144996)
## Mitigation cost per sector in Million $US 2019, for the maximum % mititation potential
## EPA Non-CO2 MAC Sum of each $s/ton mit x tons mitigated at that wedge of abatement cost potential - calculated in Excel
## cost x sum(MMT for that sector) + cost x sum(MMT additional at that cost for that sector) + etc.
## :oil and :pip split EPA "Gas" by proportion of combined output from ys_0
    7090.604857 269.0188218 6862.194574*pip_of_GAS 6862.194574*oil_of_GAS 1795.700322],
# DenseAxisArray Indices, 2D
[:EPAemiss :MaxpercMit :MitCostTot],
[:agr,:min,:pip,:oil,:wst])

CH4calc = DenseAxisArray([
## CH4 Emissions (MMt), 2019 / $US Billion (2017) Value Added inputs (kapital and labor, i.e. productive actiity)
    [CH4emiss[:EPAemiss,c]/sum(va_0[yr,:,c]) for c in CH4sectors];;
## subtract (maximum) mitigated CH4, so CH4 of remaining emissions after maximum abatement (at <$1000/t)
    [CH4emiss[:EPAemiss,c]*(1-CH4emiss[:MaxpercMit,c])/(sum(va_0[yr,:,c])+CH4emiss[:MitCostTot,c]*10^-3) for c in CH4sectors];;
## Total Cost of Mitigation is the standard VA inputs + the additional cost of the mitigation in US$Bill
    [CH4emiss[:MitCostTot,c]*10^-3+sum(va_0[yr,:,c]) for c in CH4sectors]],
CH4sectors, # dimension 1 indexed by sector
[:CH4Intens :CH4MitIntens :TotCostwMit ]) # dimension 2 indexed by values
## Relative cost of VA including max mitigation
MitCostoverVA = DenseAxisArray([CH4calc[c,:TotCostwMit]/sum(va_0[yr,:,c]) for c in CH4sectors], CH4sectors)

vam_0 = deepcopy(va_0) #copy for slack mitigating activties
# # benchmark value added input levels for with additional % of costs for mitigation
for c in CH4sectors
    vam_0[:,:,c] = va_0[:,:,c].data .* MitCostoverVA[c]
end
## default of 0 for non-emitting and less significant emitting sectors
ch4VASInt = DenseAxisArray(zeros(length(J)),J)
## Set vector of standard CH4 intensities: CH4 (in CO2eq)/value-added factor inputs
for c in CH4sectors
    ch4VASInt[c] = CH4calc[c,:CH4Intens]
end
## Still default of 0 for non-emitting and less significant emitting sectors
ch4VAMInt = DenseAxisArray(zeros(length(J)),J)
## Set vector of CH4 *mitigated* intensities: Mitigated CH4 (in CO2eq)/value-added factor inputs
for c in CH4sectors
    ch4VAMInt[c] = CH4calc[c,:CH4MitIntens]
end

CO2Int = DenseAxisArray(zeros(length(J)),J)
# (MMtCO2eq/$Bill) EPA Inventory 2022 sum of CO2 MMt for coal, and for gas & oil per Billion of total benchmark intermediate input from sector
CO2Int[:min] =  895.9*10^-3/sum(id_0[yr,:min,:]) 
CO2Int[:oil] =  2104.80*10^-3/sum(id_0[yr,:oil,:]) 

TotCO2bnchmk =  895.9*10^-3 + 2104.8*10^-3
TotCH4bnchmk =  sum(CH4emiss[:EPAemiss,:])
TotGHGbnchmk =  TotCO2bnchmk + TotCH4bnchmk
## End data preparations

## Set tax rates
CO2_taxrate = 190
CH4_taxrate = 190

MultiNat = MPSGEModel()

@parameters(MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    ch4_tax, 0.
    tax_CO2, 0.
end)

@sectors(MultiNat,begin
    Y[J],      (description = "Sectoral Production",)
    A[I],      (description = "Armington Supply",)
    VAS[J],    (description = "Value Added, standard")
    # VAM[J],    (description = "Value Added, with additional max mitigating activity")
    MS[M],     (description = "Margin Supply",)
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

@commodities(MultiNat,begin
    PA[I],   (description = "Armington Price")
    PY[J],   (description = "Supply",)
    PVA[VA], (description = "Value-added Input to VA blocks",)
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

for j∈J
    @production(MultiNat, Y[j], [t=0, s = 0], begin
        [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
        [@input(PA[i], id_0[yr,i,j], s, taxes = [Tax(RA,tax_CO2 * CO2Int[i])]) for i∈I]...
         @input(PVAM[j], sum(va_0[yr,VA,j]), s)
    end)
end

for j∈J
    @production(MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], va_0[yr,va,j], va, taxes = [Tax(RA,ch4_tax* ch4VASInt[j])]) for va∈VA]...
    end)
end

# TODO set up to loop over all VAMs
for j∈CH4sectors
    if VAM_costover[:VAM5,j]>1 # Some sectors are still cumulatively -negative costs at $5/t, so filtering those out.
        @production(MultiNat, VAM5[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM5,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM5,j])]) for va∈VA]...
        end)
    end
end

for j∈CH4sectors
        @production(MultiNat, VAM10[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM10,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM10,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM15[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM15,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM15,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM20[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM20,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM20,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM30[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM30,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM30,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM40[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM40,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM40,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM50[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM50,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM50,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM100[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM100,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM100,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM500[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j]*VAM_costover[:VAM500,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM500,j])]) for va∈VA]...
        end)
end
for j∈CH4sectors
        @production(MultiNat, VAM1000[j], [t=0, s = 0, va => s = 1], begin
            [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
            [@input(PVA[va], va_0[yr,va,j].*VAM_costover[:VAM1000,j], va, taxes = [Tax(RA, ch4_tax*VAM_CH4EmInt[:VAM1000,j])]) for va∈VA]...
        end)
end

## First pass, all EPA mitigation potential up to $1,000/t
# Slack mitigating VA activities for main CH4 producing sectors: total weighted Marginal Abatment version
# for j∈CH4sectors
#     @production(MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
#         [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
#         [@input(PVA[va], vam_0[yr,va,j], va, taxes = [Tax(RA, ch4_tax* ch4VAMInt[j])]) for va∈VA]...
#     end)
# end


for m∈M
    @production(MultiNat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end

for i∈I
    @production(MultiNat, A[i], [t = 2, s = 0, dm => s = 2], begin
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i])],reference_price=1-ta_0[yr,i])]...
        [@output(PFX, x_0[yr,i], t)]...
        [@input(PM[m], md_0[yr,m,i], s) for m∈M]...
        @input(PY[i], y_0[yr,i], dm)
        @input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[yr,i])
    end)
end

@demand(MultiNat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    end,begin
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVA[va], sum(va_0[yr,va,j] for j∈J)) for va∈VA]...
end, elasticity = 1)

## CO2 emissions for fossil fuel sectors are the activity levels times the (base) total emissions intensity 
@aux_constraint(MultiNat, CO2em[:min],  CO2em[:min] - Y[:min]*895.9*10^-3)
@aux_constraint(MultiNat, CO2em[:oil],  CO2em[:oil] - Y[:oil]*2104.80*10^-3)
## Total CO2 emissions are the sum of emissions from the 2 fossil fuel sectors (constraint expressed as equantion = 0)
@aux_constraint(MultiNat, CO2TotEm, CO2TotEm - (CO2em[:min] + CO2em[:oil]))
## CH4 emissions for each CH4 emitting sector are the sum of (either): VA Standard activity levels x standard CH4 emissions intensity (benchmark = 1 x base emissions)
## +/or VA Mitigating activities at each tier, activity level x base emissions x mitigated emissions factor (mitigated intensity/baseline intensity)
for c in CH4sectors
    # VAM, the old one, for testing.
    # @aux_constraint(MultiNat, CH4em[c],  CH4em[c] - (VAS[c]*CH4emiss[:EPAemiss,c]+VAM[c]*CH4emiss[:EPAemiss,c]*ch4VAMInt[c]/ch4VASInt[c]))
    @aux_constraint(MultiNat, CH4em[c],  CH4em[c] - (VAS[c]*CH4emiss[:EPAemiss,c]+
    ifelse(VAM_costover[:VAM5,c]>1, VAM5[c]*CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM5,c]/ch4VASInt[c] , 0) + # The $5/t tier includes sectors with negative costs, so have to filter those out here (AS WELL as in the production block)
    VAM10[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM10,c]/ch4VASInt[c]+
    VAM15[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM15,c]/ch4VASInt[c]+
    VAM20[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM20,c]/ch4VASInt[c]+
    VAM30[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM30,c]/ch4VASInt[c]+
    VAM40[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM40,c]/ch4VASInt[c]+
    VAM50[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM50,c]/ch4VASInt[c]+
    VAM100[c] *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM100,c]/ch4VASInt[c]+
    VAM500[c] *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM500,c]/ch4VASInt[c]+
    VAM1000[c]*CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM1000,c]/ch4VASInt[c]
    ))
end
 
## Total CH4 Emissions are the sum of emissions from CH4 emitting sectors
@aux_constraint(MultiNat, CH4TotEm, CH4TotEm - (CH4em[:agr] + CH4em[:min] + CH4em[:oil] + CH4em[:pip] + CH4em[:wst] ))
## Total GHG (CO2 & CH4) emissions in Mt CO2eq
@aux_constraint(MultiNat, TotEm, TotEm - (CH4TotEm + CO2TotEm))
set_silent(MultiNat)

# Benchmark 
# fix(RA, sum(fd_0[yr,i,:pce] for i∈I))
## Note: Benchmark doesn't solve at 0 interation because of margins of slack activity. Does balance with interactions or slack vars and production commented out.
solve!(MultiNat)
#; cumulative_iteration_limit = 0)

fullvrbnch = generate_report(MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
# print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))
# fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var])

# Initialize a Dataframe to save final demand results
FDemand = DataFrame(index=Vector{Symbol}(undef, length(I)),
desc=Vector{Symbol}(undef, length(I)), 
bnch=Vector{Float64}(undef, length(I)), 
cntr=Vector{Float64}(undef, length(I)), 
ch4=Vector{Float64}(undef, length(I)),
cO2=Vector{Float64}(undef, length(I)),
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

fullvrCntr = generate_report(MultiNat)
rename!(fullvrCntr, :value => :cntr, :margin => :cntrmarg)

for (n,i) in enumerate(I)
    FDemand[n,:cntr] = value(demand(RA,PA[i]))
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

# tax are at $s per ton of CH4 (CO2eq)
## "OR 1600 * CO2 conversion rate back to per ton of 
## Or alternatively, re-work data replacing CH4 in actual tons"
set_value!(ch4_tax, CH4_taxrate)

solve!(MultiNat)

for (n,i) in enumerate(I)
    FDemand[n,:ch4] = value(demand(RA,PA[i]))
end

Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

fullvrch4 = generate_report(MultiNat)
rename!(fullvrch4, :value => :ch4, :margin => :ch4marg)

# print(sort(df, :margin, by= abs, rev=true))
# print(df)

# # Counterfactual Fossil fuel extraction is ALSO taxed at emissions intensitiy of input x tax in $/ton
set_value!(tax_CO2, CO2_taxrate)

solve!(MultiNat, cumulative_iteration_limit=10000) #;

for (n,i) in enumerate(I)
    FDemand[n,:both] = value(demand(RA,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:5]
# Rs[:,2][68:71]

fullvrboth = generate_report(MultiNat)
rename!(fullvrboth, :value => :both, :margin => :bothmarg)

## Then, set CH4 taxes back to 0 to generate CO2 tax ONLY
    set_value!(ch4_tax, 0.0)

solve!(MultiNat, cumulative_iteration_limit=10000) #;

for (n,i) in enumerate(I)
    FDemand[n,:cO2] = value(demand(RA,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:4]
# Rs[:,2][68:71]
fullvrCO2 = generate_report(MultiNat)
rename!(fullvrCO2, :value => :CO2, :margin => :CO2marg)

#Generate Dataframe with all results (including names expressions)
FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrCO2, fullvrboth, fullvrCntr, on = [:var], makeunique=true);
Compare = FullResults[1:end,[1,2,4,6,8,10]];
## Sum the difference of each tax applied individually
Compare.sum = Compare.ch4 .- 1 + Compare.CO2 .-1

## sum difference between the sum of the individual taxes to the combined taxes
Compare.diff = Vector{Union{Missing, Float64}}(undef, length(Compare[:,1]))
for n in 1:length(Compare[:,1])
    if Compare.ch4[n] == 0 || Compare.CO2[n] == 0
        Compare.diff[n] = missing
    elseif Compare.both[n] >= 1 &&  Compare.sum[n] >= 0
        Compare.diff[n] = (Compare.both[n] - 1) - Compare.sum[n]
    elseif Compare.both[n] >= 1 &&  Compare.sum[n] <= 0
        Compare.diff[n] = 1 - Compare.both[n] + Compare.sum[n]
    elseif Compare.both[n] <= 1 &&  Compare.sum[n] <= 0
        Compare.diff[n] = 1 - Compare.both[n] + Compare.sum[n]
    elseif Compare.both[n] <= 1 &&  Compare.sum[n] >= 0
        Compare.diff[n] = 1 - Compare.both[n] - Compare.sum[n]
    end
end

Compare[!,:var] = Symbol.(Compare[:,:var]);
# print(sort!(Compare, :var));
# println(sort!(Compare, :diff, by = abs, rev=true))#[1:25,:])
# CSV.write("C:\\Users\\Eli\\Box\\CGE\\MPSGE-JL\\First Mulit GHG taxes Paper\\MultiResults$(Dates.format(now(),"yyyy-mm-d_HhM")).csv", Compare, missingstring="missing", bom=true)

## Look at the sum of each reduction compared to the combined reduction
compCO2em = filter(:var => ==(:CO2TotEm), Compare);
println("CO2 Reduction Sum Mt: ", (TotCO2bnchmk - only(compCO2em[:,:ch4])+TotCO2bnchmk - only(compCO2em[:,:CO2]))*10^3) # Sum of the difference from benchmark CO2 emissions for both taxes applied separately
println("CO2 Reduction Combined: ", (TotCO2bnchmk - only(compCO2em[:,:both]))*10^3) # Total CO2 reduction with taxes combined
println("CO2 Reduction Interaction: ", (TotCO2bnchmk - only(compCO2em[:,:ch4])+TotCO2bnchmk - only(compCO2em[:,:CO2])- # Sum of the difference from benchmark CO2 emissions for both taxes applied separately
(TotCO2bnchmk - only(compCO2em[:,:both])))*10^3)# Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes.
compCH4em = filter(:var => ==(:CH4TotEm), Compare);
println("CH4 Reduction Sum: ", (TotCH4bnchmk - only(compCH4em[:,:ch4]) + TotCH4bnchmk - only(compCH4em[:,:CO2]))*10^3) # Sum of the difference from benchmark CH4 emissions for both taxes applied separately
println("CH4 Reduction Combined: ", (TotCH4bnchmk - only(compCH4em[:,:both]))*10^3) # Total CH4 reduction with taxes combined
println("CH4 Reduction Interaction: ", (TotCH4bnchmk - only(compCH4em[:,:ch4]) + TotCH4bnchmk - only(compCH4em[:,:CO2])- # Sum of the difference from benchmark CH4 emissions for both taxes applied separately
(TotCH4bnchmk - only(compCH4em[:,:both])))*10^3) # Total CH4 reduction with taxes combined
compTotem = filter(:var => ==(:TotEm), Compare);
println("Total GHG Emission Reduction Sum: ", (TotGHGbnchmk - only(compTotem[:,:ch4]) + TotGHGbnchmk - only(compTotem[:,:CO2]))*10^3) # Sum of the difference from benchmark GHG (CO2&CH4) emissions for both taxes applied separately
println("Total GHG Emission Reduction Combined: ", (TotGHGbnchmk - only(compTotem[:,:both]))*10^3) # Total GHG (CO2&CH4) reduction with taxes combined
println("Total GHG Emission Interaction: ", (TotGHGbnchmk - only(compTotem[:,:ch4]) + TotGHGbnchmk - only(compTotem[:,:CO2])- # Sum of the difference from benchmark GHG (CO2&CH4) emissions for both taxes applied separately
(TotGHGbnchmk - only(compTotem[:,:both])))*10^3) # Total GHG (CO2&CH4) reduction with taxes combined

EmissionReductionResults = DataFrame(
["CO2" "Gt" TotCO2bnchmk - only(compCO2em[:,:ch4]) + TotCO2bnchmk - only(compCO2em[:,:CO2]) TotCO2bnchmk - only(compCO2em[:,:both])  TotCO2bnchmk - only(compCO2em[:,:ch4]) + TotCO2bnchmk - only(compCO2em[:,:CO2]) - (TotCO2bnchmk - only(compCO2em[:,:both])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk - only(compCH4em[:,:ch4]) + TotCH4bnchmk - only(compCH4em[:,:CO2]) TotCH4bnchmk - only(compCH4em[:,:both]) TotCH4bnchmk - only(compCH4em[:,:ch4]) + TotCH4bnchmk - only(compCH4em[:,:CO2]) - (TotCH4bnchmk - only(compCH4em[:,:both])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk - only(compTotem[:,:ch4]) + TotGHGbnchmk - only(compTotem[:,:CO2]) TotGHGbnchmk - only(compTotem[:,:both]) TotGHGbnchmk - only(compTotem[:,:ch4]) + TotGHGbnchmk - only(compTotem[:,:CO2]) - (TotGHGbnchmk - only(compTotem[:,:both]))], ["Emissions", "Unit", "Sum of each tax", "taxes combined" ,"Interactions"])

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
        set_value!(tax1, i)#/10^3)
        set_value!(tax2, cnst*i)#/10^3)
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
        value(TotEm),value(CH4TotEm),value(CH4TotEm)
        # value(CH4em[:oil]),value(CH4em[:pip]),value(CO2em[:min]),value(CO2em[:oil])
        ] 
            )
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    plt = plot(margemiss[!,:tax], margemiss[!,:Emissions].*10^3, title= "RA:\$$RAval fxd:$isfixed", label="CH4 Emiss", ylim=(0,4000), xlabel="$tax1 $tax2in \$/t")
    plch4 = plot(margemiss[!,:tax], margemiss[!,:CH4Emissions].*10^3, label=false, title="CH4 Emissions", ylim=(0,4000), xlabel="$tax1 $tax2in \$/t")
    plco2 = plot(margemiss[!,:tax], margemiss[!,:CO2Emissions].*10^3, label=false, title="CO2 Emissions", ylim=(0,4000), xlabel="$tax1 $tax2in \$/t")
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
    return margemiss, plt, Testvars, ResultsTroubleshoot, plch4, plco2 #,pla, plm, plp, plo, plw, plpa, plpm, plpp, plpo, plpw, 
    # plcdyp, plcdap, plfdrap, plAp, plAo
end

## Add Quant diff and % of consumption columns for Final Demand report dataframe
# FDemand[:,:ch4Qdelta]=FDemand[:,:ch4].-FDemand[:,:bnch]
# FDemand[:,:CO2Qdelta]=FDemand[:,:cO2].-FDemand[:,:bnch]
# FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
# FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
# FDemand[:,:ch4Qpc]=FDemand[:,:ch4]./sum(FDemand[:,:ch4])*100
# FDemand[:,:CO2Qpc]=FDemand[:,:cO2]./sum(FDemand[:,:cO2])*100
# FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

# set_silent(MultiNat)
# checkch4CO2 = plottaxemisscurve(ch4_tax, tax_CO2, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# checkch4CO2[2]
# print("checkch4CO2",checkch4CO2[3])

fix(RA, sum(fd_0[yr,I,:pce]))
set_upper_bound(MultiNat[:A][:pip], 2)
# # MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pip]))
# # fix(RA, 14008.668551652801)
# # unfix(RA)
# # checkCO2 = plottaxemisscurve(tax_CO2, ch4_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2 = plottaxemisscurve(tax_CO2, ch4_tax, 0, 4, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2[2] # Total emissions
# # println("PApip up lim =", upper_bound(MultiNat[:A][:pip]))
# checkCO2[3]

# checkCO2[5] # Y agr
# checkCO2[6] # Y min
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
# # # checkch4CO2 = plottaxemisscurve(ch4_tax, tax_CO2, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# # checkch4CO2[2]

checkch4 = plottaxemisscurve(ch4_tax,tax_CO2, 0, 10, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
checkch4[2]
checkch4[5]
checkch4[6]
# checkch4[7]
# checkch4[8]
# checkch4[9]
# checkch4[10]
# checkch4[11]
# checkch4[12]
# checkch4[13]
# checkch4[14]
# checkch4[15] #ch3 emissions
# checkch4[16] # Co2 emissions
print("checkch4",checkch4[3])

# png(checkch4[2], "./Results/CH4to500-Totemiss")
# png(checkch4[15], "./Results/CH4to500-ch4emiss")
# png(checkch4[16], "./Results/CH4to500-co2emiss")
# png(checkch4[8], "./Results/CH4to500-oilactivity")
# png(checkch4[13], "./Results/CH4to500-oilprice")

# print("checkch4",checkch4[3])
# print(check[1])

# CH4sectors = [:agr,:min,:pip,:oil,:wst] #:uti? # subset index for relevant CH4 mitigation sectors (VA slack in benchmark)

# png(checkCO2[2], "./Results/CO2to1600RAUnfxd")

EmissUnits_mt = DataFrame();
EmissUnits_mt.Unit=["Mt"; "MtCO2eq"; "MtCO2eq"];
EmissionReductionResults_Mt =[EmissionReductionResults[:,1:1] EmissUnits_mt EmissionReductionResults[:,3:5].*10^3] 
