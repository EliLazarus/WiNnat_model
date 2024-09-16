# Adapted from WiNDC National Model
using MPSGE
using DataFrames, JLD2
using JuMP
using MPSGE.JuMP.Containers
using CSV
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
CH4emiss = DenseAxisArray([260.483532*10^-3+13.70952225*10^-3 59.31302643*10^-3 224.8979059*10^-3/2 224.8979059*10^-3/2 111.5049515*10^-3+20.36144996*10^-3
## Weighted average mitigated potential per sector
    (260.483532*.304+13.70952225*.280)/(300.8535461+260.483532+13.70952225)  .645 .475 .475 (111.5049515*.050+20.36144996*.350)/(111.5049515+20.36144996)
## Mitigation cost per sector in Million $US 2019, for the maximum % mititation potential
## EPA Non-CO2 MAC Sum of each $s/ton mit x tons mitigated at that wedge of abatement cost potential - calculated in Excel
## cost x sum(MMT for that sector) + cost x sum(MMT additional at that cost for that sector) + etc.
## :oil and :pip split EPA "Gas" in half for want of more specific MAC data or proportiona
    7090.604857 269.0188218 6862.194574/2 6862.194574/2 1795.700322],
[:EPAemiss :MaxpercMit :MitCostTot],
[:agr,:min,:pip,:oil,:wst])

CH4calc = DenseAxisArray([
## CH4 Emissions (MMt), 2019 / $US Billion (2017) Value Added inputs (kapital and labor, i.e. productive actiity)
    [CH4emiss[:EPAemiss,c]/sum(va_0[yr,:,c]) for c in CH4sectors];;
## subtract (maximum) mitigated CH4, so CH4 of remaining emissions after maximum abatement (at <$1000/t)
    [(1-CH4emiss[:MaxpercMit,c])*CH4emiss[:EPAemiss,c]/sum(va_0[yr,:,c]) for c in CH4sectors];;
## Total Cost of Mitigation is the standard VA inputs + the additional cost of the mitigation in US$Bill
    [CH4emiss[:MitCostTot,c]/10^3+sum(va_0[yr,:,c]) for c in CH4sectors]],
CH4sectors, # dimension 1 indexed by sector
[:CH4Intens :CH4MitIntens :TotCostwMit ]) # dimension 2 indexed by values

## Relative cost of VA including max mitigation
MitCostoverVA = DenseAxisArray([CH4calc[c,:TotCostwMit]/sum(va_0[yr,:,c]) for c in CH4sectors],
    CH4sectors)

vam_0 = deepcopy(va_0) #copy for slack mitigating activties
# benchmark value added input levels for with additional % of costs for mitigation
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

MultiNat = MPSGEModel()

@parameters(MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    ch4_tax, 0.
    tax_CO2, 0.
end)

@sectors(MultiNat,begin
    Y[J],  (description = "Sectoral Production",)
    A[I],  (description = "Armington Supply",)
    VAS[J], (description = "Value Added, standard")
    VAM[J], (description = "Value Added, with additional mitigating activity")
    MS[M], (description = "Margin Supply",)

end)

@commodities(MultiNat,begin
    PA[I],   (description = "Armington Price",)
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

# Slack mitigating VA activities for main CH4 producing sectors
for j∈CH4sectors
    @production(MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], vam_0[yr,va,j], va, taxes = [Tax(RA, ch4_tax* ch4VAMInt[j])]) for va∈VA]...
    end)
end

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
## CH4 emissions for each CH4 emitting sector are the sum of (either): VA Standard activity levels x standard CH4 emissions intensity 
## +/or VA Mitigating activity x 
@aux_constraint(MultiNat, CH4em[:agr],  CH4em[:agr] - (VAS[:agr]*CH4emiss[:EPAemiss,:agr]+VAM[:agr]*CH4emiss[:EPAemiss,:agr]*ch4VAMInt[:agr]/ch4VASInt[:agr]))
@aux_constraint(MultiNat, CH4em[:min],  CH4em[:min] - (VAS[:min]*CH4emiss[:EPAemiss,:min]+VAM[:min]*CH4emiss[:EPAemiss,:min]*ch4VAMInt[:min]/ch4VASInt[:min]))
@aux_constraint(MultiNat, CH4em[:oil],  CH4em[:oil] - (VAS[:oil]*CH4emiss[:EPAemiss,:oil]+VAM[:oil]*CH4emiss[:EPAemiss,:oil]*ch4VAMInt[:oil]/ch4VASInt[:oil]))
@aux_constraint(MultiNat, CH4em[:pip],  CH4em[:pip] - (VAS[:pip]*CH4emiss[:EPAemiss,:pip]+VAM[:pip]*CH4emiss[:EPAemiss,:pip]*ch4VAMInt[:pip]/ch4VASInt[:pip]))
@aux_constraint(MultiNat, CH4em[:wst],  CH4em[:wst] - (VAS[:wst]*CH4emiss[:EPAemiss,:wst]+VAM[:wst]*CH4emiss[:EPAemiss,:wst]*ch4VAMInt[:wst]/ch4VASInt[:wst]))
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
taxrate_ch4 = 190. 
## "OR 1600 * CO2 conversion rate back to per ton of 
## Or alternatively, re-work data replacing CH4 in actual tons"
set_value!(ch4_tax, taxrate_ch4)#*10^-3) # Divided by 1,000 for $Bill/MMt

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
CO2_taxrate = 190 
set_value!(tax_CO2, CO2_taxrate)# * 10^-3) # Divided by 1,000 for $Bill/MMt

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
    margemiss = DataFrame(tax=Float64[], Emissions=Float64[])
    Testvars = DataFrame(taxrt=Float64[], 
    
    Ypip=Float64[],Yoil=Float64[],CH4emoil=Float64[],CH4empip=Float64[],VASpip=Float64[],VAMpip=Float64[],CO2emin=Float64[],CO2emoil=Float64[]
    
    )
    for i in start:interval:finish
        set_value!(tax1, i)#/10^3)
        set_value!(tax2, cnst*i)#/10^3)
        solve!(MultiNat, output="no");
        Results = generate_report(MultiNat)
        Results[!,:var] = Symbol.(Results[:,:var]);
        push!(margemiss, [i only(filter(:var => ==(:TotEm), Results)[:, :value])])
        push!(Testvars, [i,                

    value(Y[:pip]),value(Y[:oil]),value(CH4em[:oil]),value(CH4em[:pip]),value(VAS[:pip]),value(VAM[:pip]),value(CO2em[:min]),value(CO2em[:oil])] 
        
        )
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    return margemiss,    plot(margemiss[!,:tax], margemiss[!,:Emissions].*10^3, title= "$tax1 $tax2in RA:\$$RAval fxd:$isfixed", ylim=(0,maximum(margemiss[!,:Emissions].*10^3))), Testvars
end

EmissUnits_mt = DataFrame();
EmissUnits_mt.Unit=["Mt"; "MtCO2eq"; "MtCO2eq"];
EmissionReductionResults_Mt =[EmissionReductionResults[:,1:1] EmissUnits_mt EmissionReductionResults[:,3:5].*10^3] 

FDemand[:,:ch4Qdelta]=FDemand[:,:ch4].-FDemand[:,:bnch]
FDemand[:,:CO2Qdelta]=FDemand[:,:cO2].-FDemand[:,:bnch]
FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
FDemand[:,:ch4Qpc]=FDemand[:,:ch4]./sum(FDemand[:,:ch4])*100
FDemand[:,:cO2Qpc]=FDemand[:,:cO2]./sum(FDemand[:,:cO2])*100
FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

# set_silent(MultiNat)
# checkch4CO2 = plottaxemisscurve(ch4_tax, tax_CO2, 1260, 1, 1280)
# checkch4CO2[2]
# print("checkch4CO2",checkch4CO2[3])
# checkch4CO2 = plottaxemisscurve(ch4_tax, tax_CO2, 0, 1, 1600)
# checkch4CO2[2]
# print("checkch4CO2",checkch4CO2[3])

# fix(RA, sum(fd_0[yr,i,:pce] for i∈I))
# fix(RA, 14008.668551652801)
# unfix(RA)
# checkCO2 = plottaxemisscurve(tax_CO2, ch4_tax, 1000, 1, 1600, round(value(RA),digits=2), is_fixed(RA), 0)
# checkCO2[2]
# print("checkCO2",checkCO2[3])
# checkch4CO2 = plottaxemisscurve(ch4_tax, tax_CO2, 1000, 1, 1600, round(value(RA),digits=2), is_fixed(RA))
# checkch4CO2[2]

# checkch4 = plottaxemisscurve(ch4_tax,tax_CO2, 900, 100, 1600, 0)
# checkch4[2]
# print("checkch4",checkch4[3])
# print(check[1])

# CH4sectors = [:agr,:min,:pip,:oil,:wst] #:uti? # subset index for relevant CH4 mitigation sectors (VA slack in benchmark)

# png(check[2], "CO2to1600")