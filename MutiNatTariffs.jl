if !@isdefined(MultiNat)
    include("./MultiGHGNat.jl")
end

set_value!(CH₄_tax,0.)
set_value!(CO₂_tax,0.)
 
set_value!(ta,ta_0DAA[Jp]); # make sure all WiNDC counterfactual taxes back to benchmark
set_value!(tm,tm_0DAA[Jp]); # make sure all WiNDC counterfactual taxes back to benchmark
set_value!(tx,tx_0DAA[Jp]); # make sure export tax set back to 0 
# for i in Ip 
#     set_value!(tx[i],0)
# end
solve!(MultiNat, cumulative_iteration_limit=10000)
sort(generate_report(MultiNat),:value) # just checking that we're back to benchmark

#### Add experiment for emissions impact of tarriff
for i in Ip
    if value(tm[i])<0.15
    set_value!(tm[i],0.15)
    end
end

solve!(MultiNat, cumulative_iteration_limit=10000)
imptar_results = generate_report(MultiNat)
imptar_results.var = Symbol.(imptar_results.var)
print(sort(imptar_results,:var))

CO2_tariff_15 = value(CO2TotEm)
CH4_tariff_15 = value(CH4TotEm)
GHG_tariff_15 = value(TotEm)

set_value!(tm,tm_0DAA[Jp]); # re-set import tariff to benchmark for 'export tariff' counterfactual

for i in Ip
    set_value!(tx[i],0.15)
end

#  ROWtaxrevenue = 302 #-sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:ROW])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])
solve!(MultiNat, cumulative_iteration_limit=10000)

exptar_results = generate_report(MultiNat)
exptar_results.var = Symbol.(exptar_results.var)
print(sort(exptar_results,:var))

CO2_Xtariff_15 = value(CO2TotEm)
CH4_Xtariff_15 = value(CH4TotEm)
GHG_Xtariff_15 = value(TotEm)

for i in Ip
    set_value!(tx[i],0.15)
end
for i in Ip
    if value(tm[i])<0.15
    set_value!(tm[i],0.15)
    end
end

solve!(MultiNat, cumulative_iteration_limit=10000)

CO2_ImXtariff_15 = value(CO2TotEm)
CH4_ImXtariff_15 = value(CH4TotEm)
GHG_ImXtariff_15 = value(TotEm)

### Reporting the emissions results from each counterfactual
EmissionsTariffs = DataFrame(
["CO2" "Gt" TotCO2bnchmk ;;
CO2_tariff_15;;
CO2_Xtariff_15;;
CO2_ImXtariff_15;
"CH4" "GtCO2eq" TotCH4bnchmk ;;
CH4_tariff_15 ;;
CH4_Xtariff_15;;
CH4_ImXtariff_15;
"GHGs" "GtCO2eq" TotGHGbnchmk ;;
 GHG_tariff_15;;
 GHG_Xtariff_15;;
 GHG_ImXtariff_15;
 "CO2" "%Δ" round((TotCO2bnchmk/TotCO2bnchmk-1)*100,digits=2) ;;
round((CO2_tariff_15/TotCO2bnchmk-1)*100,digits=2);;
round((CO2_Xtariff_15/TotCO2bnchmk-1)*100,digits=2);;
round((CO2_ImXtariff_15/TotCO2bnchmk-1)*100,digits=2);
"CH4" "%Δ" round((TotCH4bnchmk/TotCH4bnchmk-1) *100,digits=2);;
round((CH4_tariff_15/TotCH4bnchmk-1) *100,digits=2);;
round((CH4_Xtariff_15/TotCH4bnchmk-1)*100,digits=2);;
round((CH4_ImXtariff_15/TotCH4bnchmk-1)*100,digits=2);
"GHGs" "%Δ" round((TotGHGbnchmk/TotGHGbnchmk-1) *100,digits=2);;
 round((GHG_tariff_15/TotGHGbnchmk-1)*100,digits=2);;
 round((GHG_Xtariff_15/TotGHGbnchmk-1)*100,digits=2);;
 round((GHG_ImXtariff_15/TotGHGbnchmk-1)*100,digits=2);
],
["Emissions", "Unit", "Bnchmrk_Emissions", "Tariff", "Export tariffs", "Imp&Exp tariffs"])
EmissionsTariffs_Mt = vcat(hcat(EmissionsTariffs[1:3,1:2],EmissionsTariffs[1:3,3:end].*10^3),EmissionsTariffs[4:6,:]); EmissionsTariffs_Mt[:,2] = ["Mt" ,"MtCO2eq" ,"MtCO2eq", "%Δ", "%Δ", "%Δ"];
print(EmissionsTariffs_Mt)
# EmissionReductionResultsTariffs = DataFrame(
# ["CO2" "Gt" TotCO2bnchmk ;;
# TotCO2bnchmk - CO2_co2tax
# ; # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
# "CH4" "GtCO2eq" TotCH4bnchmk ;;
# TotCH4bnchmk - CH4_co2tax; # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
# "GHGs" "GtCO2eq" TotGHGbnchmk ;;
# TotGHGbnchmk - GHG_co2tax],
# ["Em_reductions", "Unit", "Bnchmrk_Emissions", "CO2tax_reduc", "CH4tax_reduc","Sum_of_each_tax", "both_taxes_combined" ,"Interactions"])
# EmissReductionsTariffs_Mt = hcat(EmissionReductionResults[:,1:2],EmissionReductionResults[:,3:end].*10^3); EmissionReductionResults_Mt[:,2] = ["Mt" ,"MtCO2eq" ,"MtCO2eq"];

