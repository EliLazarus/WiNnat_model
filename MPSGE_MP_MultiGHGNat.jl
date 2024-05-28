# TODO 1) Done: (assuming the index was correct) Actualy CO2 intentisities (from EPA Inventory, TAble 1-4)
# TODO 2) Separate CO2 tax
# TODO 3) Re-check magnitudes and units for CH4 taxes
# TODO 4) Magniture for both - theoretically I want to be able to apply a tax at a reasonable 
# absolute magnitude, BUT taxes in MPSGE at % per input, or output

# Adapted from WiNDC National Model
using MPSGE_MP
using DataFrames, JLD2
using JuMP
using MPSGE_MP.JuMP.Containers
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S

I = [i for i∈S[:i] if i∉[:use,:oth]]
J = [i for i∈S[:j] if i∉[:use,:oth]]
#subset index for Slack CH4 mitigation production
C = [:agr,:oil,:wst,:min,] #:uti, :pip,
VA = [va for va∈S[:va] if va!=:othtax]
FD = S[:fd]
TS = S[:ts]
YR = S[:yr]
M = S[:m]

a_0 = P[:a_0]
id_0 = P[:id_0]
ys_0 = P[:ys_0]
tm_0 = P[:tm_0]
va_0 = P[:va_0]

yr = Symbol(2017)

## Base Data for reference
# EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDS
CH4emissdatadf = DataFrame(Wsector = [:agr,:agr,:agr,:min,:oil,:wst,:wst],
## EPA Total Emissions per sector, MMt CO2eq 
EPAemiss =[300.8535461,260.483532,13.70952225,59.31302643,224.8979059,111.5049515,20.36144996],
## EPA maximum % of abatement per sector at <$1000/t
MaxpercMit = [.038,.304,.280,.645,.475,.050,.350],
PAsector = ["AGRICULTURE, CROP","AGRICULTURE, LIVE","AGRICULTURE, RICE","ENERGY, COL","ENERGY, GAS","WASTE, LAN","WASTE, WWR"])

## Calculate CH4 Intensity factors from EPA data
# EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDS
# Sum and weighted average with current aggregation, before disaggregation
CH4emiss = DenseAxisArray([300.8535461+260.483532+13.70952225 59.31302643 224.8979059 111.5049515+20.36144996
## Weighted average per sector
(300.8535461*.038+260.483532*.304+13.70952225*.280)/(300.8535461+260.483532+13.70952225)  .645 .475 (111.5049515*.050+20.36144996*.350)/(111.5049515+20.36144996);
## Million $US 2019: EPA Non-CO2 MAC Sum of each $s/ton mit x tons mitigated at that wedge of abatement cost potential - calculated in Excel
## cost x sum(MMT for that sector) + cost x sum(MMT additional at that cost for that sector) + etc.
8962.758884 269.0188218 6944.389519 1795.700322],
[:EPAemiss :MaxpercMit :MitCostTot],
[:agr,:min,:oil,:wst])

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

vam_0 = deepcopy(va_0) #copy for slack mitigating activty
# inputs with additional % needed for mitigation
vam_0[:,:,:] = va_0[:,:,:].data .*5 # Default for mitigation, back up to ensure slack (but redundant with filtered index for mitigating production)
vam_0[:,:,:agr] = va_0[:,:,:agr].data .* MitCostoverVA[:agr]
vam_0[:,:,:min] = va_0[:,:,:min].data .* MitCostoverVA[:min]   
vam_0[:,:,:oil] = va_0[:,:,:oil].data .* MitCostoverVA[:oil]   
## Use oil proportional mitigation cost for pipeline,no EPA MAC data
vam_0[:,:,:pip] = va_0[:,:,:pip].data .* MitCostoverVA[:oil]   
vam_0[:,:,:wst] = va_0[:,:,:wst].data .* MitCostoverVA[:wst]   

## Set vector of Values of CH4 intensities, default of 0
ch4Int = DenseAxisArray(zeros(length(J)),J)
ch4Int[:agr] = CH4calc[:agr,:CH4Intens]
ch4Int[:min] = CH4calc[:min,:CH4Intens]
ch4Int[:oil] = CH4calc[:oil,:CH4Intens]
## Use oil proportional mitigation cost for pipeline,no EPA MAC data
ch4Int[:pip] = CH4calc[:oil,:CH4Intens]
ch4Int[:wst] = CH4calc[:wst,:CH4Intens]

ch4Intmit = DenseAxisArray(zeros(length(J)),J)
ch4Intmit[:agr] = CH4calc[:agr,:CH4MitIntens]
ch4Intmit[:min] = CH4calc[:min,:CH4MitIntens]
ch4Intmit[:oil] = CH4calc[:oil,:CH4MitIntens]
ch4Intmit[:pip] = CH4calc[:oil,:CH4MitIntens]
ch4Intmit[:wst] = CH4calc[:wst,:CH4MitIntens]

CO2Int = DenseAxisArray(zeros(length(J)),J)
# MMtCO2eq/$Bill
CO2Int[:min] =  8.121789674248 
CO2Int[:oil] =  5.113219404020

md_0 = P[:md_0]
fd_0 = P[:fd_0]
m_0 = P[:m_0]
ty_0 = P[:ty_0]
ms_0 = P[:ms_0]
bopdef_0 = P[:bopdef_0]
x_0 = P[:x_0]
ta_0 = P[:ta_0]
#s_0 = P[:s_0]
fs_0 = P[:fs_0]
y_0 = P[:y_0];

MP_MultiNat = MPSGEModel()

@parameters(MP_MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    ch4_tax[J], DenseAxisArray(zeros(length(J)),J)
    ch4_taxmit[C], DenseAxisArray(zeros(length(C)),C)
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

@consumer(MP_MultiNat, RA, description = "Representative Agent")

# @auxiliaries(MP_MultiNat, begin
#     CH4tax[J], (description = "CH4 Emissions, emissions Inensity x Inputs")
#     CH4mittax[C], (description = "CH4 Emissions, mitigated emissions Inensity x Inputs")

# end)


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
        [@input(PVA[va], va_0[yr,va,j], va, taxes = [Tax(RA,ch4_tax[j]* ch4Int[j])]) for va∈VA]...
    end)
end

# Slack mitigating VA activities for main CH4 producing sectors
for j∈C
    @production(MP_MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], vam_0[yr,va,j], va, taxes = [Tax(RA, ch4_taxmit[j]* ch4Intmit[j])]) for va∈VA]...
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
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i])],reference_price=1-(ta_0[yr,i]))]... # , Tax(RA,tax_co2*CO2Int[i]), Tax(RA,CH4tax[i]) # price +tax_co2*CO2Int[i]+CH4tax[i]
        [@output(PFX, x_0[yr,i], t)]...
        [@input(PM[m], md_0[yr,m,i], s) for m∈M]...
        @input(PY[i], y_0[yr,i], dm)
        @input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[yr,i])
    end)
end

@demand(MP_MultiNat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    end,begin
    [@endowment(PY[i], fs_0[yr,i]) for i∈I]...
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVA[va], sum(va_0[yr,va,j] for j∈J)) for va∈VA]...
end)

# for j ∈ J
#     @aux_constraint(MP_MultiNat, CH4tax[j],
#     (taxrate_ch4 * ch4Int[j]*sum(va_0[yr,:,j])*10^-3)/sum(va_0[yr,:,j])) #*VAS[j]
# end

# ## Tax rate (% of input value) is function of emissions, tax ($/t) activity level and input levels over the benchmark input value
# for j ∈ C
#     @aux_constraint(MP_MultiNat, CH4mittax[j],
#         (taxrate_ch4 * ch4Intmit[j]*sum(vam_0[yr,:,j])*10^-3)/sum(vam_0[yr,:,j])) #* VAM[j]
#     end
    
# Benchmark 
# fix(RA, sum(fd_0[yr,i,:pce] for i∈I))

solve!(MP_MultiNat)#; cumulative_iteration_limit = 0)

fullvrbnch = generate_report(MP_MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))

# tax per ton of CH4 (CO2eq??)
# set_value!(taxrate_ch4,  1.)
taxrate_ch4 = 30.
for j in J
set_value!(ch4_tax[j],      (taxrate_ch4 *10^-3))
end

unfix(RA)

solve!(MP_MultiNat)
fullvrch4 = generate_report(MP_MultiNat)
rename!(fullvrch4, :value => :ch4, :margin => :ch4marg)
# print(sort(df, :margin, by= abs, rev=true))
# print(df)

# WinNat Counterfactual
# unfix(RA)

# set_value!(ta,0)
# set_value!(tm,0)

# # Counterfactual Fossil fuel extraction is taxed at emissions intensitiy of input x tax in $/ton
CO2_taxrate = 80 
set_value!(tax_co2, CO2_taxrate * 10^-3)


solve!(MP_MultiNat, cumulative_iteration_limit=10000) #;
fullvrboth = generate_report(MP_MultiNat)
rename!(fullvrboth, :value => :both, :margin => :bothmarg)

# #Then, set CH4 taxes back to 0 to generate CO2 tax only
    set_value!(ch4_tax, 0.0)

solve!(MP_MultiNat, cumulative_iteration_limit=10000) #;
#Generate Dataframe with all results (including names expressions)
fullvrco2 = generate_report(MP_MultiNat)
rename!(fullvrco2, :value => :co2, :margin => :co2marg)

FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrco2, fullvrboth, on = [:var], makeunique=true)
CompareFullResults = FullResults[1:end,[1,4,6,8,2]]

CompareFullResults.check = Vector{Union{Missing, Float64}}(undef, length(CompareFullResults[:,1]))
for n in 1:length(CompareFullResults[:,1]) 
    if CompareFullResults[n,2] < 1 && CompareFullResults[n,3] > 1 #less than 1, more than 1
        CompareFullResults.check[n] = CompareFullResults[n,2]+(CompareFullResults[n,3]-1)
elseif CompareFullResults[n,2] < 1 && CompareFullResults[n,3] < 1 #less than 1, less than 1
    CompareFullResults.check[n] = CompareFullResults[n,2]-(1-CompareFullResults[n,3])
elseif CompareFullResults[n,2] > 1 && CompareFullResults[n,3] < 1 #less than 1, more than 1
    CompareFullResults.check[n] = CompareFullResults[n,2]-(1-CompareFullResults[n,3])
else CompareFullResults.check[n] =  CompareFullResults[n,2]+(CompareFullResults[n,3]-1)
end
end
CompareFullResults.diff = CompareFullResults.both - CompareFullResults.check

CompareFullResults[!,:var] = Symbol.(CompareFullResults[:,:var])
# print(CompareFullResults)#[359:1200,:])
# print(sort!(CompareFullResults, :bmkmarg, by = abs, rev =true))#[5800:6000,:])
print(sort!(CompareFullResults, :var))
print(sort!(CompareFullResults, :diff, by = abs, rev=true))#[1:15,:])