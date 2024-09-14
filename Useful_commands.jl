n=0
for i in keys(WiNnat.productions);n+=1;
       println(n,"; ",i)
       end

for (a,b) in enumerate(array/Dict)
       print("Numerical index $a, value $b")
end

WiNnat.productions[collect(keys(WiNnat.productions))[11]]  #or
MPSGE.production_sectors(MultiNat)[1]

compensated_demand(Y[:ppd], PY[:air])
value(compensated_demand(Y[:ppd], PY[:air]))
demand(RA,PA[:agr])
value(demand(RA,PA[:agr]))
sum([value(demand(RA,PA[i])) for i in I])
compensated_demand(WiNnat[:A][:agr], WiNnat[:PFX], :dm) # dm here specifies the nest name, useful if the same commodity appears in different nests
## Model name prints Dict of MPSGEModelstruct

print([keys(ta_0[Symbol("2017"), :]) ta_0[Symbol("2017"),:]])

only(filter(:var => ==(:TotEm), fullvrboth)[:value])

print(filter(row -> row.var ∈ [Symbol("Y[alt]"),Symbol("Y[cep]"),Symbol("Y[ote]"),Symbol("Y[mmf]"),Symbol("Y[rec]"),Symbol("Y[min]"),Symbol("Y[pmt]"),Symbol("Y[sle]"),Symbol("Y[pet]"),Symbol("Y[oil]"),Symbol("Y[wst]"),Symbol("Y[agr]")], Compare))

print(filter(row -> row.var ∈ [Symbol("Y[alt]"),Symbol("Y[cep]"),Symbol("Y[ote]"),Symbol("Y[mmf]"),Symbol("Y[rec]"),Symbol("Y[min]"),Symbol("Y[pmt]"),Symbol("Y[sle]"),Symbol("Y[pet]"),Symbol("Y[oil]"),Symbol("Y[wst]"),Symbol("Y[agr]")], Compare))

set_silent(MultiNat)
unset_silent(MultiNat)

collect(1:.2:10) #turn iterator into vector
