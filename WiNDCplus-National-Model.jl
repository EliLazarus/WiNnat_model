# WiNDC National Model
using MPSGE
using DataFrames, JLD2, CSV
using JuMP
using MPSGE.JuMP.Containers
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
P = MultiNatdata

S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S

I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
Ip = [[x for x in I if x∉[:uti]]; [:uel,:ugs, :uwt, :coa]]#, :gas]]
Jp = deepcopy(Ip) # Index for WiNDC BEA Sectors
VA = [va for va∈S[:va] if va!=:othtax] # Index Value Added (compen = returns to labour/wage, 'surplus' = returns to Kapital)
FD = S[:fd]
TS = S[:ts]
YR = S[:yr] # Index for years for potential multi year runs
M = S[:m]

a_m0 = P[:a_0] #	    "Armington supply",
id_m0 = P[:id_0] #	"Intermediate demand",
ys_m0 = P[:ys_0]#	"Sectoral supply",
va_m0 = P[:va_0] #	"Value added",
md_m0 = P[:md_0] #	"Margin demand",
fd_m0 = P[:fd_0] #	"Final demand",
m_m0 = P[:m_0] #	    "Imports",
ms_m0 = P[:ms_0] #	"Margin supply",
bopdef_m0 = P[:bopdef_0] #	"Balance of payments deficit",
x_m0 = P[:x_0] #	    "Exports of goods and services",
# fs_m0 = P[:fs_0] #	"Household supply", # All zeros
y_m0 = P[:y_0];  #	"Gross output",

ty_m0 = P[:ty_0] #	"Output tax rate"
tm_m0 = P[:tm_0] #	"Import tariff"; Initial, for price 
ta_m0 = P[:ta_0] #	"Tax net subsidy rate on intermediate demand", benchmark as data also for price level

yr = Symbol(2020)

WiNnat = MPSGEModel()

@parameters(WiNnat, begin
    ta_m[Jp], ta_m0[Jp]
    ty_m[Jp], ty_m0[Jp]
    tm_m[Jp], tm_m0[Jp]
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
    Y[Jp],  (description = "Sectoral Production",)
    A[Ip],  (description = "Armington Supply",)
    MS[M], (description = "Margin Supply",)
end)

@commodities(WiNnat,begin
    PA[Ip],   (description = "Armington Price",)
    PY[Jp],   (description = "Supply",)
    PVA[VA], (description = "Value-added",)
    PM[M],   (description = "Margin Price",)
    PFX,     (description = "Foreign Exachange",)
end)

@consumer(WiNnat, RA, description = "Representative Agent")

for j∈Jp
    @production(WiNnat, Y[j], [t= t_elas_y, s = elas_y, va => s = elas_va], begin # 0, 0, 1
        [@output(PY[i],ys_m0[j,i], t, taxes = [Tax(RA,ty_m[j])]) for i∈Ip]... 
        [@input(PA[i], id_m0[i,j], s) for i∈Ip]...
        [@input(PVA[va], va_m0[va,j], va) for va∈VA]...
    end)
end

for m∈M
    @production(WiNnat, MS[m], [t = 0, s = 0], begin # 0, 0
    # @production(WiNnat, MS[m], [t = t_elas_m, s = elas_m], begin # 0, 0
        [@output(PM[m], sum(ms_m0[i,m] for i∈Ip), t)]...
        [@input(PY[i], ms_m0[i,m], s) for i∈Ip]...
    end)
end

for i∈Ip
    # @production(WiNnat, A[i], [t = t_elas_a, s = elas_a, dm => s = elas_dm], begin # 2, 0, 2
    @production(WiNnat, A[i], [t = 2, s = 0, dm => s = 2], begin # 2, 0, 2
        [@output(PA[i], a_m0[i], t, taxes=[Tax(RA,ta_m[i])],reference_price=1-ta_m0[i])]...
        [@output(PFX, x_m0[i], t)]...
        [@input(PM[m], md_m0[m,i], s) for m∈M]...
        @input(PY[i], y_m0[i], dm)
        @input(PFX, m_m0[i], dm, taxes = [Tax(RA,tm_m[i])],reference_price=1+tm_m0[i])
    end)
end

@demand(WiNnat, RA, begin
    [@final_demand(PA[i], fd_m0[i,:pce]) for i∈Ip]...
    # [@endowment(PY[i], fs_m0[i]) for i∈Ip]...
    @endowment(PFX, only(bopdef_m0))
    [@endowment(PA[i], -sum(fd_m0[i,xfd] for xfd∈FD if xfd!=:pce)) for i∈Ip]...
    [@endowment(PVA[va], sum(va_m0[va,j] for j∈Jp)) for va∈VA]...
end, elasticity = d_elas_ra)

# Benchmark 
# fix(RA, sum(fd_m0[i,:pce] for i∈Ip))

solve!(WiNnat; cumulative_iteration_limit = 0)

df_benchmark = generate_report(WiNnat);
# sort!(df_benchmark, [:value])
rename!(df_benchmark, :value => :bnchmrk, :margin => :bmkmarg)
df_benchmark[!,:var] = Symbol.(df_benchmark[:,:var]);
print(sort(df_benchmark, :bmkmarg))
# print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))
# Sectors = CSV.read(joinpath(@__DIR__,"Sectorsplus.csv"), DataFrame);

# # Initialize a Dataframe to save final demand results
# FDemandWiNnat = DataFrame(index=Vector{Symbol}(undef, length(Ip)),
# desc=Vector{Symbol}(undef, length(Ip)), 
# bnch=Vector{Float64}(undef, length(Ip)), 
# bnchPr=Vector{Float64}(undef, length(Ip)), 
# cntr=Vector{Float64}(undef, length(Ip)),
# cntrPr=Vector{Float64}(undef, length(Ip)),
# cntrval=Vector{Float64}(undef, length(Ip)))
# for (n,i) in enumerate(Ip)
#     FDemandWiNnat[n,:index]= i
#     FDemandWiNnat[n,:desc] = Symbol(Sectors[Sectors.index.==string(i),2][1])
#     FDemandWiNnat[n,:bnch] = value(demand(WiNnat[:RA],WiNnat[:PA][i]))
#     FDemandWiNnat[n,:bnchPr] = filter(:var => ==(Symbol("PA[$i]")),df_benchmark)[1,2]
# end

# # Counterfactual
# # fix(RA,12453.896315446877)

# # 12453.896315446877/sum(fd_0[i,:pce] for i∈Ip)
# # 13154.978277803244/sum(fd_0[i,:pce] for i∈Ip)
# # for i in Ip; fdW+=sum(fd_0[i,FD[2:18]])*-value(PA[i])+sum(va_0[VA[1],i])*value(PVA[VA[1]])+sum(va_0[VA[2],i])*value(PVA[VA[2]])
# # end

# # for i in Ip; fdW+=sum(fd_0[i,:pce])*-value(PA[i])
# # end

# set_value!(ta_m,0)
# set_value!(tm_m,0)

# solve!(WiNnat)

# df = generate_report(WiNnat);
# df |>
#     x -> sort(x, :margin, rev=true)
# [Y value.(Y)][sortperm([Y value.(Y)][:,2], rev= true), :]