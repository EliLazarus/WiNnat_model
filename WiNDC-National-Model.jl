# WiNDC National Model
using MPSGE
using DataFrames, JLD2
using JuMP
using MPSGE.JuMP.Containers
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S

I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
J = [i for i∈S[:j] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
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
fs_0 = P[:fs_0] #	"Household supply", # All zeros
y_0 = P[:y_0];  #	"Gross output",

ty_0 = P[:ty_0] #	"Output tax rate"
tm_0 = P[:tm_0] #	"Import tariff"; Initial, for price 
ta_0 = P[:ta_0] #	"Tax net subsidy rate on intermediate demand", benchmark as data also for price level

yr = Symbol(2017)

WiNnat = MPSGEModel()

@parameters(WiNnat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    t_elas_y, 0            
    elas_y,   0            
    elas_va,  1          
    t_elas_m, 0            
    elas_m,   0         
    t_elas_a, 2           
    elas_a,   0         
    elas_dm,  2          
    d_elas_ra,1
end)

@sectors(WiNnat,begin
    Y[J],  (description = "Sectoral Production",)
    A[I],  (description = "Armington Supply",)
    MS[M], (description = "Margin Supply",)
end)

@commodities(WiNnat,begin
    PA[I],   (description = "Armington Price",)
    PY[J],   (description = "Supply",)
    PVA[VA], (description = "Value-added",)
    PM[M],   (description = "Margin Price",)
    PFX,     (description = "Foreign Exachange",)
end)

@consumer(WiNnat, RA, description = "Representative Agent")

for j∈J
    @production(WiNnat, Y[j], [t= t_elas_y, s = elas_y, va => s = elas_va], begin # 0, 0, 1
        [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
        [@input(PA[i], id_0[yr,i,j], s) for i∈I]...
        [@input(PVA[va], va_0[yr,va,j], va) for va∈VA]...
    end)
end



for m∈M
    @production(WiNnat, MS[m], [t = 0, s = 0], begin # 0, 0
    # @production(WiNnat, MS[m], [t = t_elas_m, s = elas_m], begin # 0, 0
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end

for i∈I
    # @production(WiNnat, A[i], [t = t_elas_a, s = elas_a, dm => s = elas_dm], begin # 2, 0, 2
    @production(WiNnat, A[i], [t = 2, s = 0, dm => s = 2], begin # 2, 0, 2
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i])],reference_price=1-ta_0[yr,i])]...
        [@output(PFX, x_0[yr,i], t)]...
        [@input(PM[m], md_0[yr,m,i], s) for m∈M]...
        @input(PY[i], y_0[yr,i], dm)
        @input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[yr,i])
    end)
end

@demand(WiNnat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    end,begin
    [@endowment(PY[i], fs_0[yr,i]) for i∈I]...
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVA[va], sum(va_0[yr,va,j] for j∈J)) for va∈VA]...
end, elasticity = d_elas_ra)

# Benchmark 
# fix(RA, sum(fd_0[yr,i,:pce] for i∈I))

solve!(WiNnat; cumulative_iteration_limit = 0)

df_benchmark = generate_report(WiNnat);

rename!(df_benchmark, :value => :bnchmrk, :margin => :bmkmarg)
df_benchmark[!,:var] = Symbol.(df_benchmark[:,:var]);

# print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))
Sectors = CSV.read("Sectors.csv", DataFrame);

# Initialize a Dataframe to save final demand results
FDemandWiNnat = DataFrame(index=Vector{Symbol}(undef, length(I)),
desc=Vector{Symbol}(undef, length(I)), 
bnch=Vector{Float64}(undef, length(I)), 
bnchPr=Vector{Float64}(undef, length(I)), 
cntr=Vector{Float64}(undef, length(I)),
cntrPr=Vector{Float64}(undef, length(I)),
cntrval=Vector{Float64}(undef, length(I)))
for (n,i) in enumerate(I)
    FDemandWiNnat[n,:index]= i
    FDemandWiNnat[n,:desc] = Symbol(Sectors[Sectors.index.==string(i),2][1])
    FDemandWiNnat[n,:bnch] = value(demand(WiNnat[:RA],WiNnat[:PA][i]))
    FDemandWiNnat[n,:bnchPr] = filter(:var => ==(Symbol("PA[$i]")),df_benchmark)[1,2]
end

# Counterfactual
# fix(RA,12453.896315446877)

# 12453.896315446877/sum(fd_0[yr,i,:pce] for i∈I)
# 13154.978277803244/sum(fd_0[yr,i,:pce] for i∈I)
# for i in I; fdW+=sum(fd_0[yr,i,FD[2:18]])*-value(PA[i])+sum(va_0[yr,VA[1],i])*value(PVA[VA[1]])+sum(va_0[yr,VA[2],i])*value(PVA[VA[2]])
# end

# for i in I; fdW+=sum(fd_0[yr,i,:pce])*-value(PA[i])
# end

set_value!(ta,0)
set_value!(tm,0)

solve!(WiNnat)
df_counter = generate_report(WiNnat)
rename!(df_counter, :value => :counter, :margin => :countermarg)
df_counter[!,:var] = Symbol.(df_counter[:,:var]);

for (n,i) in enumerate(I)
     FDemandWiNnat[n,:cntr] = value(demand(WiNnat[:RA],WiNnat[:PA][i]))
     FDemandWiNnat[n,:cntrPr] = filter(:var => ==(Symbol("PA[$i]")),df_counter)[1,2]
     FDemandWiNnat[n,:cntrval] = FDemandWiNnat[n,:cntr] * FDemandWiNnat[n,:cntrPr]

end

# [value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:pip]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:trn]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:trk]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:wtt]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:air]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:otr]))	,
# sum([value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:pip]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:trn]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:trk]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:wtt]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:air]))	,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PY][:otr]))]) ,
# value(compensated_demand(WiNnat[:MS][:trn],WiNnat[:PM][:trn]))  ,
# value(MS[:trn])]


# df = generate_report(WiNnat);
# df |>
#     x -> sort(x, :margin, rev=true)
# [Y value.(Y)][sortperm([Y value.(Y)][:,2], rev= true), :]