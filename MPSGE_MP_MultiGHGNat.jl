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
C = [:agr,:oil,:uti,:wst,:pip,:min,]
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
vam_0 = deepcopy(va_0) #copy for slack mitigating activty
# inputs with additional % needed for mitigation
vam_0[:,:,:] = va_0[:,:,:].data .*5 # Default for mitigation, slack
vam_0[:,:,:agr] = va_0[:,:,:agr].data .*1.01213932526577 
vam_0[:,:,:min] = va_0[:,:,:min].data .*1.00217834542134 
vam_0[:,:,:oil] = va_0[:,:,:oil].data .*1.00389656859361 
vam_0[:,:,:pip] = va_0[:,:,:pip].data .*1.00389656859361 
vam_0[:,:,:wst] = va_0[:,:,:wst].data .*1.22546701239839 

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

yr = Symbol(2017)

#y_ = [j for j∈J if sum(ys_0[yr,j,i] for i∈I) !=0]
#a_ = [i_ for i_∈I if a_0[yr,i_]!=0]

MP_MultiNat = MPSGEModel()

@parameters(MP_MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    tch4[J], DenseAxisArray(zeros(length(J)),J)
	tco2[J], DenseAxisArray(zeros(length(J)),J)
    # tch4 = add!(MultiNat, Parameter(:tch4, indices = (y_,), value=zeros(length(y_))))
	# tco2 = add!(MultiNat, Parameter(:tco2, indices = (y_,), value=zeros(length(y_))))

end)

@sectors(MP_MultiNat,begin
    Y[J],  (description = "Sectoral Production",)
    A[I],  (description = "Armington Supply",)
    VALAD[J], (description = "Value Added, standard")
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

for j∈J
    @production(MP_MultiNat, Y[j], [t=0, s = 0], begin
        [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
        [@input(PA[i], id_0[yr,i,j], s, taxes = [Tax(RA,tco2[j])]) for i∈I]...
        @input(PVAM[j], sum(va_0[yr,VA,j]), s)
    end)
end

for j∈J
    @production(MP_MultiNat, VALAD[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], va_0[yr,va,j], va, taxes = [Tax(RA,tch4[j])]) for va∈VA]...
    end)
end

for j∈J
    @production(MP_MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], vam_0[yr,va,j], va) for va∈VA]...
    end)
end
# for j∈C
#     @production(MP_MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
#         [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
#         [@input(PVA[va], vam_0[yr,va,j], va) for va∈VA]...
#     end)
# end

for m∈M
    @production(MP_MultiNat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end

for i∈I
    @production(MP_MultiNat, A[i], [t = 2, s = 0, dm => s = 2], begin
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i]), Tax(RA,tch4[i]), Tax(RA,tco2[i])],reference_price=1-(ta_0[yr,i]+tch4[i]+tco2[j]))]... 
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

# Benchmark 
fix(RA, sum(fd_0[yr,i,:pce] for i∈I))

solve!(MP_MultiNat)#; cumulative_iteration_limit = 0)

fullvrbnch = generate_report(MP_MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))

set_value!(tch4[:agr], 0.4)
set_value!(tch4[:oil], 0.4)
set_value!(tch4[:uti], 0.1)
set_value!(tch4[:wst], 0.15)
set_value!(tch4[:pip], 0.3)
set_value!(tch4[:min], 0.2)
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


set_value!(tco2[:oil], 0.2) # nominal value of tax and carbon intensity combied
set_value!(tco2[:min], 0.4) # nominal value of tax and carbon intensity combied
solve!(MP_MultiNat, cumulative_iteration_limit=10000) #;
fullvrboth = generate_report(MP_MultiNat)
rename!(fullvrboth, :value => :both, :margin => :bothmarg)

# # Counterfactual Fossil fuel extraction is taxed at (VERY NOMINAL ) ~ carbon content of combustion
# #First, set CH4 taxes back to 0
for i in I
    set_value!(tch4[i], 0.0)
end

solve!(MP_MultiNat, cumulative_iteration_limit=10000) #;
#Generate Dataframe with all results (including names expressions)
fullvrco2 = generate_report(MP_MultiNat)
rename!(fullvrco2, :value => :co2, :margin => :co2marg)

FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrco2, fullvrboth, on = [:var], makeunique=true)
CompareFullResults = FullResults#[1:end,[1,2,4,6,8]]

CompareFullResults[!,:var] = Symbol.(CompareFullResults[:,:var])
# print(CompareFullResults)#[359:1200,:])
# print(sort!(CompareFullResults, :bmkmarg, by = abs, rev =true))#[5800:6000,:])
print(sort!(CompareFullResults, :var))

Testlimitedindex = outerjoin(CompareFullResults, CompareFullResults_j, on = [:var], makeunique=true)

# Testlimitedindex = outerjoin(CompareFullResults[1:end,[1,3,4,5]], CompareFullResults_j[1:end,[1,3,4,5]], on = [:var], makeunique=true)
Testlimitedindex.CH4diff = Testlimitedindex.ch4 .- Testlimitedindex.ch4_1
Testlimitedindex.co2diff = Testlimitedindex.co2 .- Testlimitedindex.co2_1
Testlimitedindex.bothdiff = Testlimitedindex.both .- Testlimitedindex.both_1
print(sort!(Testlimitedindex, by = abs, :bothdiff, rev=true))
