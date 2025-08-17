n=0
for i in keys(MultiNat.productions);n+=1;
       println(n,"; ",i)
       end

for (a,b) in enumerate(array/Dict)
       print("Numerical index $a, value $b")
end

for (a,b) in enumerate(MPSGE.production_sectors(MultiNat))
       println("$a, $b")
end

WiNnat.productions[collect(keys(WiNnat.productions))[11]]  #or
WiNnat.productions[Symbol("A[hos]")]  # This one!! with production name.
MPSGE.production_sectors(MultiNat)[1]


compensated_demand(MultiNat[:Y][:ppd], MultiNat[:PY][:air])
value(compensated_demand(WiNnat[:Y][:ppd], WiNnat[:PY][:air]))
demand(RA,PA[:agr])
value(demand(RA,PA[:agr]))
sum([value(demand(RA,PA[i])) for i in I])
compensated_demand(WiNnat[:A][:agr], WiNnat[:PFX], :dm) # dm here specifies the nest name, useful if the same commodity appears in different nests
## Model name prints Dict of MPSGEModelstruct

print([keys(ta_0[Symbol("2017"), :]) ta_0[Symbol("2017"),:]])

only(filter(:var => ==(:TotEm), fullvrboth)[:value])

print(filter(row -> row.var ∈ [Symbol("Y[alt]"),Symbol("Y[cep]"),Symbol("Y[ote]"),Symbol("Y[mmf]"),Symbol("Y[rec]"),Symbol("Y[min]"),Symbol("Y[pmt]"),Symbol("Y[sle]"),Symbol("Y[pet]"),Symbol("Y[oil]"),Symbol("Y[wst]"),Symbol("Y[agr]")], Compare))

print(filter(row -> row.var ∈ [Symbol("Y[alt]"),Symbol("Y[cep]"),Symbol("Y[ote]"),Symbol("Y[mmf]"),Symbol("Y[rec]"),Symbol("Y[min]"),Symbol("Y[pmt]"),Symbol("Y[sle]"),Symbol("Y[pet]"),Symbol("Y[oil]"),Symbol("Y[wst]"),Symbol("Y[agr]")], Compare))
print(sort(filter(row -> row.var ∈ Symbol.(MultiNat[:Y][CH4sectors]), fullvrch4),:ch4tax))

subset(fullvrbnch, :var => ByRow(x -> occursin(r"PVAM",string(x)) ))

subset(Compare, :var => ByRow(x -> occursin(r"PA",string(x)) ))

set_silent(MultiNat)
unset_silent(MultiNat)

collect(1:.2:10) #turn iterator into vector

# collects rows (filtered on index) of df into a new df
PAall=DataFrame([[],[],[],[],[],[],[],[]],[:var,      :bnchmrk,  :ch4,      :co2,      :both,     :cntr,      :sum,       :diff])
for i in Ip                                                                                                                       
       push!(PAall,Vector{Any}(filter(:var => ==(Symbol(MultiNat[:PA][i])),Compare)[1,:]))                                                       
       end
## get final demand, and benchmark fd etc (with the short description)
       for i in I
              println(WiNnat[:PA][i],": ",filter(:index => ==(string(i)), Sectors)[:,:Short_description], "\t: ",value(WiNnat[:PA][i]),"\t->",value(demand(WiNnat[:RA],WiNnat[:PA][i])),"\t;",fd_0[yr,i,:pce],"\t=>",value(demand(WiNnat[:RA],WiNnat[:PA][i])))
       end

       for i in I
              println(i,": ",round(fd_0[yr,i,:pce],digits=4),"\t",fd_0[yr,i,:pce]/sum(fd_0[yr,:,:pce]),"\t ",round(value(PA[i]),digits=4),"->", round(only(filter(:index => ==(i),FDemand)[:,:cO2tax]),digits=4),"\t PrxQ: ", value(PA[i])*only(filter(:index => ==(i), FDemand)[:,:cO2tax])/a)
       end

       CheckFDemand = DataFrame([Any[] for c in 1:6],[:sector, :fd_0, :fd_0proportion, :PAi, :FDmdCO2, :CO2proportion])
       for i in I
              push!(CheckFDemand,[i,round(fd_0[yr,i,:pce],digits=4),fd_0[yr,i,:pce]/sum(fd_0[yr,:,:pce]),round(value(PA[i]),digits=4), round(only(filter(:index => ==(i),FDemand)[:,:cO2tax]),digits=4), value(PA[i])*only(filter(:index => ==(i), FDemand)[:,:cO2tax])/a])
       end

       # Govt sectors, by guess
       GovtSectors = [:fnd, :fen, :slg, :sle, :soc, :fdd]

       println(i,": ",only(filter(y-> y.index ==(string(i)), Sectors)[!,:Short_description]),":\t",sum(ys_0[yr,i,:]))
       filter(y-> y.index ==(string(:ppd)), Sectors)[!,:Short_description]

    "   unstack([data_frame],
        [columns to use as row keys],
        [column storing column names],
        [column storing the values];
        combine=[function to apply to the values]) "

              for (n,i) in enumerate(sectors(MultiNat))                                                                                                                      
                     println(n,": ",i)                                                                                                                                              
                     end
              MultiNat.productions[sectors(MultiNat)[263]]

              sum([value(MPSGE.tax_revenue(A[i],RA)) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])
              sum([value(MPSGE.tax_revenue(Y[i],RA)) for i in [i for i in Ip]])
              # More specific to the model
              sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
              sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])

# Emissionsresults = deepcopy(Emissions_Mt)
# for z in [Emissionsresults, EmissionReductionResults_Mt, Incomeresults, TaxRevresults ]
# delete!(z,[i for i in 1:nrow(z)])
# end
println(filter(x->x.var in [Symbol("Y[$i]") for i in Ip],fullvrbnch))
filter(x->x.EVCES2==minimum(resultdf[:,:EVCES2]),resultdf)

subset(bnchM1, :var => ByRow(x -> occursin(r"RA", x))) # need vars as strings for regex

# Show entire DataFrame in repl
show(stdout,"text/plain",df)

# Copy to clipboard
clipboard(sprint(show, "text/tab-separated-values", resultdf))

@run fn() # goes straight into the debugger for the fn!
@code_warntype

### My hack to run compensated_demand with nests for MultiNat electricity for CES Utility without parent_name_chain
### Goes in build.jl
function compensated_demand(S::ScalarSector, C::ScalarCommodity, nest::Symbol; virtual = false)
if nest == :elect
N =   [netputs(S,C)[2]]#[:elect, :home_nrg_expd, :housing_exp, :non_transp, :s]
elseif nest == :veh_elect
N = [netputs(S,C)[1]] #[:veh_elect, :fuels, :pers_transp, :transp, :s]
else
    N = [n for n∈netputs(S,C) if nest∈parent_name_chain(n)] #Then it checks that the actual nest name is there.
end
    if isempty(N)
        error("This netput (input or output) doesn't appear in a nest of that name")
    end
    return sum(compensated_demand.(N, virtual = virtual); init = 0)
end


