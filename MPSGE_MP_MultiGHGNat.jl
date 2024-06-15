# Adapted from WiNDC National Model
using MPSGE_MP
using DataFrames, JLD2
using JuMP
using MPSGE_MP.JuMP.Containers
using CSV
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S
Sectors = CSV.read("Sectors.csv", DataFrame);

I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
J = [i for i∈S[:j] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
#subset index for Slack CH4 mitigation production
C = [:agr,:oil,:pip,:wst,:min,] #:uti, :pip,
VA = [va for va∈S[:va] if va!=:othtax] # Index Value Added (compen = returns to labour/wage, 'surplus' = returns to Kapital)
FD = S[:fd]
TS = S[:ts]
YR = S[:yr] # Index for years for potential multi year runs
M = S[:m]

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

    ## Base Data for reference
    # EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDS
    CH4emissdatadf = DataFrame(Wsector = [:agr,:agr,:min,:oil,:wst,:wst],
    ## EPA Total Emissions per sector, MMt CO2eq 
    EPAemiss =[260.483532,13.70952225,59.31302643,224.8979059,111.5049515,20.36144996],
    ## EPA maximum % of abatement per sector at <$1000/t
    MaxpercMit = [.304,.280,.645,.475,.050,.350],
    PAsector = ["AGRICULTURE, LIVE","AGRICULTURE, RICE","ENERGY, COL",
    "ENERGY, GAS","WASTE, LAN","WASTE, WWR"])

## Calculate CH4 Intensity factors from EPA data
# EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDS
# Sum and weighted average with current aggregation, before disaggregation
CH4emiss = DenseAxisArray([260.483532+13.70952225 59.31302643 224.8979059/2 224.8979059/2 111.5049515+20.36144996
## Weighted average per sector
    (260.483532*.304+13.70952225*.280)/(300.8535461+260.483532+13.70952225)  .645 .475 .475 (111.5049515*.050+20.36144996*.350)/(111.5049515+20.36144996)
## Million $US 2019: EPA Non-CO2 MAC Sum of each $s/ton mit x tons mitigated at that wedge of abatement cost potential - calculated in Excel
## cost x sum(MMT for that sector) + cost x sum(MMT additional at that cost for that sector) + etc.
## :oil and :pip split EPA GAS in half for want of more specific MAC data or proportiona
    7090.604857 269.0188218 6862.194574/2 6862.194574/2 1795.700322],
[:EPAemiss :MaxpercMit :MitCostTot],
[:agr,:min,:pip,:oil,:wst])

# #Alternative, calculated with dataframe
# testdf = combine(groupby(CH4emissdatadf, :Wsector), [:EPAemiss, :Mitigated] .=> sum,  renamecols=false)
# testdf.MaxpercMit = testdf[!,:Mitigated]./testdf[!,:EPAemiss]
# #bHere add already calculated, but better to upload full table and calculate in Julia
# testdf.MitCostTot =  [8962.758884, 269.0188218, 6862.194574, 1795.700322]

CH4calc = DenseAxisArray([
## CH4 Emissions (MMt), 2019 / $US Billion (2017) Value Added inputs (kapital and labor, i.e. productive actiity)
    [CH4emiss[:EPAemiss,i]/sum(va_0[yr,:,i]) for i in axes(CH4emiss)[2]];;
## subtract (maximum) mitigated CH4, so CH4 of remaining emissions after maximum abatement (at <$1000/t)
    [(1-CH4emiss[:MaxpercMit,i])*CH4emiss[:EPAemiss,i]/sum(va_0[yr,:,i]) for i in axes(CH4emiss)[2]];;
## Total Cost of Mitigation is the standard VA inputs + the additional cost of the mitigation in US$Bill
    [CH4emiss[:MitCostTot,i]/10^3+sum(va_0[yr,:,i]) for i in axes(CH4emiss)[2]]],
[i for i in axes(CH4emiss)[2]],
[:CH4Intens :CH4MitIntens :TotCostwMit ])

## Relative cost of VA including max mitigation
MitCostoverVA = DenseAxisArray([CH4calc[i,:TotCostwMit]/sum(va_0[yr,:,i]) for i in axes(CH4emiss)[2]],
    [i for i in axes(CH4emiss)[2]])

vam_0 = deepcopy(va_0) #copy for slack mitigating activties
# benchmark value added input levels for with additional % of costs for mitigation
# vam_0[:,:,:] = va_0[:,:,:].data .*5 # Default for mitigation, back up to ensure slack (but redundant with filtered index for mitigating production)
vam_0[:,:,:agr] = va_0[:,:,:agr].data .* MitCostoverVA[:agr]
vam_0[:,:,:min] = va_0[:,:,:min].data .* MitCostoverVA[:min]   
vam_0[:,:,:oil] = va_0[:,:,:oil].data .* MitCostoverVA[:oil]   
## Use oil proportional mitigation cost for pipeline,no EPA MAC data
vam_0[:,:,:pip] = va_0[:,:,:pip].data .* MitCostoverVA[:pip]   
vam_0[:,:,:wst] = va_0[:,:,:wst].data .* MitCostoverVA[:wst]   

## Set vector of Values of CH4 intensities, default of 0 for less significant emitting sectors
ch4VASInt = DenseAxisArray(zeros(length(J)),J)
ch4VASInt[:agr] = CH4calc[:agr,:CH4Intens]
ch4VASInt[:min] = CH4calc[:min,:CH4Intens]
ch4VASInt[:oil] = CH4calc[:oil,:CH4Intens]
ch4VASInt[:pip] = CH4calc[:pip,:CH4Intens]
ch4VASInt[:wst] = CH4calc[:wst,:CH4Intens]

ch4VAMInt = DenseAxisArray(zeros(length(J)),J)
ch4VAMInt[:agr] = CH4calc[:agr,:CH4MitIntens]
ch4VAMInt[:min] = CH4calc[:min,:CH4MitIntens]
ch4VAMInt[:oil] = CH4calc[:oil,:CH4MitIntens]
ch4VAMInt[:pip] = CH4calc[:pip,:CH4MitIntens]
ch4VAMInt[:wst] = CH4calc[:wst,:CH4MitIntens]

CO2Int = DenseAxisArray(zeros(length(J)),J)
# (MMtCO2eq/$Bill) EPA Inventory 2022 sum of CO2 MMt for coal, and for gas & oil per Billion of total benchmark intermediate input from sector
CO2Int[:min] =   895.9/sum(id_0[yr,:min,:]) 
CO2Int[:oil] =   2104.80/sum(id_0[yr,:oil,:]) 

## End data preparations

MP_MultiNat = MPSGEModel()

@parameters(MP_MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    ch4_tax, 0.
    tax_co2, 0.
end)

@sectors(MP_MultiNat,begin
    Y[J],  (description = "Sectoral Production",)
    A[I],  (description = "Armington Supply",)
    VAS[J], (description = "Value Added, standard")
    VAM[J], (description = "Value Added, with additional mitigating activity")
    MS[M], (description = "Margin Supply",)

end)

@commodities(MP_MultiNat,begin
    PA[I],   (description = "Armington Price",)
    PY[J],   (description = "Supply",)
    PVA[VA], (description = "Value-added Input to VA blocks",)
    PVAM[J], (description = "Value-added output - Input to Y",)
    PM[M],   (description = "Margin Price",)
    PFX,     (description = "Foreign Exachange",)
end)

@auxiliary(MP_MultiNat, CO2em, index = [[:min, :oil]])
@auxiliary(MP_MultiNat, CO2TotEm, description = "Total CO2 emissions from fossil fuels")
@auxiliary(MP_MultiNat, CH4em, index = [[:agr,:min,:oil,:pip,:wst]])
@auxiliary(MP_MultiNat, CH4TotEm, description = "Total CH4 emissions")
@auxiliary(MP_MultiNat, TotEm, description = "Total both emissions")


@consumer(MP_MultiNat, RA, description = "Representative Agent")

for j∈J
    @production(MP_MultiNat, Y[j], [t=0, s = 0], begin
        [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
        [@input(PA[i], id_0[yr,i,j], s, taxes = [Tax(RA,tax_co2 * CO2Int[i])]) for i∈I]...
        @input(PVAM[j], sum(va_0[yr,VA,j]), s)
    end)
end

for j∈J
    @production(MP_MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], va_0[yr,va,j], va, taxes = [Tax(RA,ch4_tax* ch4VASInt[j])]) for va∈VA]...
    end)
end

# Slack mitigating VA activities for main CH4 producing sectors
for j∈C
    @production(MP_MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], vam_0[yr,va,j], va, taxes = [Tax(RA, ch4_tax* ch4VAMInt[j])]) for va∈VA]...
    end)
end

for m∈M
    @production(MP_MultiNat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end

for i∈I
    @production(MP_MultiNat, A[i], [t = 2, s = 0, dm => s = 2], begin
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i])],reference_price=1-(ta_0[yr,i]))]...
        [@output(PFX, x_0[yr,i], t)]...
        [@input(PM[m], md_0[yr,m,i], s) for m∈M]...
        @input(PY[i], y_0[yr,i], dm)
        @input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[yr,i])
    end)
end

@demand(MP_MultiNat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    end,begin
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVA[va], sum(va_0[yr,va,j] for j∈J)) for va∈VA]...
end, elasticity = 1)

@aux_constraint(MP_MultiNat, CO2em[:min],  CO2em[:min] - Y[:min]*895.9)
@aux_constraint(MP_MultiNat, CO2em[:oil],  CO2em[:oil] - Y[:oil]*2104.80)
@aux_constraint(MP_MultiNat, CO2TotEm, CO2TotEm - (CO2em[:min] + CO2em[:oil]))
@aux_constraint(MP_MultiNat, CH4em[:agr],  CH4em[:agr] - (VAS[:agr]*CH4emiss[:EPAemiss,:agr]+VAM[:agr]*CH4emiss[:EPAemiss,:agr]*ch4VAMInt[:agr]/ch4VASInt[:agr]))
@aux_constraint(MP_MultiNat, CH4em[:min],  CH4em[:min] - (VAS[:min]*CH4emiss[:EPAemiss,:min]+VAM[:min]*CH4emiss[:EPAemiss,:min]*ch4VAMInt[:min]/ch4VASInt[:min]))
@aux_constraint(MP_MultiNat, CH4em[:oil],  CH4em[:oil] - (VAS[:oil]*CH4emiss[:EPAemiss,:oil]+VAM[:oil]*CH4emiss[:EPAemiss,:oil]*ch4VAMInt[:oil]/ch4VASInt[:oil]))
@aux_constraint(MP_MultiNat, CH4em[:pip],  CH4em[:pip] - (VAS[:pip]*CH4emiss[:EPAemiss,:pip]+VAM[:pip]*CH4emiss[:EPAemiss,:pip]*ch4VAMInt[:pip]/ch4VASInt[:pip]))
@aux_constraint(MP_MultiNat, CH4em[:wst],  CH4em[:wst] - (VAS[:wst]*CH4emiss[:EPAemiss,:wst]+VAM[:wst]*CH4emiss[:EPAemiss,:wst]*ch4VAMInt[:wst]/ch4VASInt[:wst]))
@aux_constraint(MP_MultiNat, CH4TotEm, CH4TotEm - (CH4em[:agr] + CH4em[:min] + CH4em[:oil] + CH4em[:pip] + CH4em[:wst] ))
@aux_constraint(MP_MultiNat, TotEm, TotEm - (CH4TotEm + CO2TotEm))


# Benchmark 
# fix(RA, sum(fd_0[yr,i,:pce] for i∈I))
## Note: Benchmark doesn't solve at 0 interation because of margins of slack activity. Does balance with interactions or slack vars and production commented out.
solve!(MP_MultiNat)#; cumulative_iteration_limit = 0)

fullvrbnch = generate_report(MP_MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))

# tax are at $s per ton of CH4 (CO2eq)
taxrate_ch4 = 190. 
"OR 1600 * CO2 conversion rate back to per ton of 
Or alternatively, re-work data replcing CH4 in actual tons"
## Divided by 1,000 for $Bill/MMt
set_value!(ch4_tax, (taxrate_ch4 *10^-3))

unfix(RA)
solve!(MP_MultiNat)

Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
Rs[:,2][1:4]
Rs[:,2][67:71]

fullvrch4 = generate_report(MP_MultiNat)
rename!(fullvrch4, :value => :ch4, :margin => :ch4marg)

# print(sort(df, :margin, by= abs, rev=true))
# print(df)

# Base WinNat Counterfactual (remove import tariffs and net subsidies on intermediate inputs)
# unfix(RA)
# set_value!(tm,0)
# set_value!(ta,0)

# # Counterfactual Fossil fuel extraction is taxed at emissions intensitiy of input x tax in $/ton
CO2_taxrate = 190 
set_value!(tax_co2, CO2_taxrate * 10^-3)

solve!(MP_MultiNat, cumulative_iteration_limit=10000) #;

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:5]
# Rs[:,2][67:71]

fullvrboth = generate_report(MP_MultiNat)
rename!(fullvrboth, :value => :both, :margin => :bothmarg)

# #Then, set CH4 taxes back to 0 to generate CO2 tax only
    set_value!(ch4_tax, 0.0)

solve!(MP_MultiNat, cumulative_iteration_limit=10000) #;

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:4]
# Rs[:,2][67:71]
#Generate Dataframe with all results (including names expressions)
fullvrco2 = generate_report(MP_MultiNat)
rename!(fullvrco2, :value => :co2, :margin => :co2marg)

FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrco2, fullvrboth, on = [:var], makeunique=true)
Compare = FullResults[1:end,[1,2,4,6,8]]
## Sum the difference of each tax applied individually
Compare.sum = Compare.ch4 .- 1 + Compare.co2 .-1

## sum difference between the sum of the individual taxes to the combined taxes
Compare.diff = Vector{Union{Missing, Float64}}(undef, length(Compare[:,1]))
for n in 1:length(Compare[:,1])
    if Compare.ch4[n] == 0 || Compare.co2[n] == 0
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

Compare[!,:var] = Symbol.(Compare[:,:var])
print(sort!(Compare, :var))
println(sort!(Compare, :diff, by = abs, rev=true))#[1:25,:])
# CSV.write("C:\\Users\\Eli\\Box\\CGE\\MPSGE-JL\\First Mulit GHG taxes Paper\\MultiResults.csv", Compare, missingstring="missing", bom=true)

compCO2em = filter(:var => ==(:CO2TotEm), Compare);
println("CO2 Reduction Sum: ", only(compCO2em[:,:bnchmrk]-compCO2em[:,:ch4]+compCO2em[:,:bnchmrk]-compCO2em[:,:co2]))
println("CO2 Reduction Combined: ", only(compCO2em[:,:bnchmrk]-compCO2em[:,:both]))
compCH4em = filter(:var => ==(:CH4TotEm), Compare);
println("CO2 Reduction Sum: ", only(compCH4em[:,:bnchmrk]-compCH4em[:,:ch4]+compCH4em[:,:bnchmrk]-compCH4em[:,:co2]))
println("CO2 Reduction Combined: ", only(compCH4em[:,:bnchmrk]-compCH4em[:,:both]))
compTotem = filter(:var => ==(:TotEm), Compare);
println("CO2 Reduction Sum: ", only(compTotem[:,:bnchmrk]-compTotem[:,:ch4]+compTotem[:,:bnchmrk]-compTotem[:,:co2]))
println("CO2 Reduction Combined: ", only(compTotem[:,:bnchmrk]-compTotem[:,:both]))
