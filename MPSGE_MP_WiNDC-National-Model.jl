# WiNDC National Model
using MPSGE_MP
using DataFrames, JLD2
using JuMP
using MPSGE_MP.JuMP.Containers
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
    @production(WiNnat, Y[j], [t=0, s = 0, va => s = 1], begin
        [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
        [@input(PA[i], id_0[yr,i,j], s) for i∈I]...
        [@input(PVA[va], va_0[yr,va,j], va) for va∈VA]...
    end)
end



for m∈M
    @production(WiNnat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end

for i∈I
    @production(WiNnat, A[i], [t = 2, s = 0, dm => s = 2], begin
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
end)

1;
# Benchmark 
#fix(RA, sum(fd_0[yr,i,:pce] for i∈I))

solve!(WiNnat; cumulative_iteration_limit = 0)

df_benchmark = generate_report(WiNnat);
# Counterfactual
fix(RA,12453.896315446877)

set_value!(ta,0)
set_value!(tm,0)

solve!(WiNnat)

df = generate_report(WiNnat);
df |>
    x -> sort(x, :margin, rev=true)
value.(Y)
