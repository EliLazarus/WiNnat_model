using WiNDC, DataFrames, MPSGE, NamedArrays, CSV, JuMP.Containers

##cd to top level of the data files
raw_data_directory = joinpath(@__DIR__,"./detailed_data/national")
# all_det_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :raw_detailed); # only for  2007, 2012, 2017
# Data is PROJECTED from summary => detailed bc detailed only 2007, 2012, 2017 
all_det_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :detailed);
all_summary_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :summary);
# For the summary level data, no issue using all years. Several years are failing
# calibration for the detailed data. 2022 works, as does 2017. Others do too, but
# those are the "important" ones for me at the moment.
# year = 2020
year = 2022
projected_det_national_data_yr = 
    NationalTable(
        get_table(all_det_national_data)|>
            x -> subset(x, :year => ByRow(==(year))) |>#,
## TEST: removing use and oth before callibration
 filter(x -> x[:commodities] !="S00401" && x[:commodities] !="S00402" && x[:commodities] !="S00300" && x[:commodities] !="S00900"), 
        all_det_national_data.sets
    )

raw_summary_national_data_yr = 
NationalTable(
    get_table(all_summary_national_data)|>
        x -> subset(x, :year => ByRow(==(year))) |>#,
        filter(x -> x[:commodities] !="S00401" && x[:commodities] !="S00402" && x[:commodities] !="S00300" && x[:commodities] !="S00900"), 
    all_summary_national_data.sets
)

# Calibrate the model. The JuMP model is also returned. 
# callibrated_det_national_data,M = calibrate(projected_det_national_data_yr)
# callibrated_summary_national_data,MS = calibrate(raw_summary_national_data_yr)

# checkCall = outerjoin(callibrated_summary_national_data.table, raw_summary_national_data_yr.table, on = [:commodities, :sectors, :year, :subtable], 
# renamecols = "" => "_raw")
# checkCall.diff = checkCall.value - checkCall.value_raw 
# sort!(checkCall, :diff)

# checkCall = outerjoin(callibrated_det_national_data.table, projected_det_national_data_yr.table, on = [:commodities, :sectors, :year, :subtable], 
# renamecols = "" => "_raw")
# checkCall.diff = checkCall.value - checkCall.value_raw 
# sort!(checkCall, :diff)

# checkCall = outerjoin(callibrated_det_national_data.table, projected_det_national_data_yr.table, on = [:commodities, :sectors, :year, :subtable], 
# renamecols = "" => "_raw")
# checkCall.diff = checkCall.value - checkCall.value_raw 
# sort!(checkCall, :diff)

## Aggregate and Compare
# Default detail to summary map
# summary_map = WiNDC.detailed_summary_map(raw_data_directory)
# aggregated_detailed_to_summary = aggregate(
#     callibrated_det_national_data,
#     :commodities => (summary_map, :detailed => :summary),
#     :sectors => (summary_map, :detailed => :summary),
#     :margin_demand => (summary_map, :detailed => :summary),
#     :margin_supply => (summary_map, :detailed => :summary),
#     :labor_demand => (summary_map, :detailed => :summary),
#     :exogenous_final_demand => (summary_map, :detailed => :summary),
#     :other_tax => (summary_map, :detailed => :summary),
#     :subsidies => (summary_map, :detailed => :summary),
#     :duty => (summary_map, :detailed => :summary),
#     :exports => (summary_map, :detailed => :summary),
#     :imports => (summary_map, :detailed => :summary),
#     :tax => (summary_map, :detailed => :summary),
#     :capital_demand => (summary_map, :detailed => :summary),
#     :personal_consumption => (summary_map, :detailed => :summary),
# )
# # Join the two DataFrames and look at the differences
# outerjoin(
#     get_table(aggregated_detailed_to_summary),
#     get_table(callibrated_summary_national_data),
#     on = [:commodities, :sectors, :year, :subtable],
#     renamecols = "_agg" => ""
# ) |>
# x -> coalesce.(x,0) |>
# x -> transform(x,
#     [:value_agg, :value] => ByRow((a,b) -> (a - b)) => :diff
# ) |>
# x -> sort(x, :diff)

# M1 = national_mpsge(aggregated_detailed_to_summary);
# solve!(M1, cumulative_iteration_limit = 0)

# M2 = national_mpsge(callibrated_summary_national_data);
# solve!(M2, cumulative_iteration_limit = 0)

# M3 = national_mpsge(callibrated_det_national_data);
# solve!(M3, cumulative_iteration_limit = 0)

# WiNDC plus = Multiple GHG update from national summary map (uti->[:uel,:ugs,:uwt], min->[:coa,:min])
summary_map = CSV.read(joinpath(@__DIR__,"summarywindcplus_detail.csv"), DataFrame)
# # Aggregate (WiNDC function) to WiNDC plus before oil/gas split 
# WplusAg = aggregate(
#     projected_det_national_data_yr,
# # callibrated_det_national_data, (for aggregate after calibration)
#     :commodities => (summary_map, :detailed => :summary),
#     :sectors => (summary_map, :detailed => :summary),
#     :margin_demand => (summary_map, :detailed => :summary),
#     :margin_supply => (summary_map, :detailed => :summary),
#     :labor_demand => (summary_map, :detailed => :summary),
#     :exogenous_final_demand => (summary_map, :detailed => :summary),
#     :other_tax => (summary_map, :detailed => :summary),
#     :subsidies => (summary_map, :detailed => :summary),
#     :duty => (summary_map, :detailed => :summary),
#     :exports => (summary_map, :detailed => :summary),
#     :imports => (summary_map, :detailed => :summary),
#     :tax => (summary_map, :detailed => :summary),
#     :capital_demand => (summary_map, :detailed => :summary),
#     :personal_consumption => (summary_map, :detailed => :summary),
# )

# # # Check benchmark solve.
# # M4 = national_mpsge(Wplus1);
# # solve!(M4, cumulative_iteration_limit = 0)
# # callibrated_det_national_data,M = calibrate(projected_det_national_data_yr)
# WplusAgC,M = calibrate(WplusAg)

###############################
## Double split process: 
# Run the first disag set, 
# run all the disag tables on the new NationalTable,
# generate new NAtionalTable, ready for Agg and then Calibrate
###############################

Gas_oil_splits= Dict(:prodgas=>0.3121, :prodoil=>0.6879, :exgas=>0.1707,:exoil=>0.8293,:imgas=>0.658,:imoil=>0.9342)
Elec_rnw_splits=Dict(:produel=>0.5768, :prodrnw=>0.4232, :exuel=>0.5768,:exrnw=>0.4232,:imuel=>0.1879,:imrnw=>0.8121)	

### Set up to Split
## Oil(&gas) to oil, gas; uel to uel, rnw; fen (fed elec) to uel, rnw; sle (state/loc elect) to uel, rnw
disagprod = DataFrame(
    old = ["211000", "211000","221100","221100","S00101","S00101","S00202","S00202"], # oil, oil, #uel,uel, #fen, fen, #sle, sle  
    new = ["211000", "211001","221100","221101","221100","221101","221100","221101"], # oil, gas, #uel, rnw, #uel, rnw, #uel, rnw
    shares = [Gas_oil_splits[:prodoil], Gas_oil_splits[:prodgas],
    Elec_rnw_splits[:produel], Elec_rnw_splits[:prodrnw],Elec_rnw_splits[:produel], Elec_rnw_splits[:prodrnw],Elec_rnw_splits[:produel], Elec_rnw_splits[:prodrnw]],
)
disagexp = DataFrame(
    old = ["211000", "211000","221100","221100","S00101","S00101","S00202","S00202"], # oil, oil, #uel,uel, #fen, fen, #sle, sle  
    new = ["211000", "211001","221100","221101","221100","221101","221100","221101"], # oil, gas, #uel, rnw, #uel, rnw, #uel, rnw
    shares = [Gas_oil_splits[:exoil], Gas_oil_splits[:exgas],
    Elec_rnw_splits[:exuel], Elec_rnw_splits[:exrnw],Elec_rnw_splits[:exuel], Elec_rnw_splits[:exrnw],Elec_rnw_splits[:exuel], Elec_rnw_splits[:exrnw]],
)
disagimp = DataFrame(
    old = ["211000", "211000","221100","221100","S00101","S00101","S00202","S00202"], # oil, oil, #uel,uel, #fen, fen, #sle, sle  
    new = ["211000", "211001","221100","221101","221100","221101","221100","221101"], # oil, gas, #uel, rnw, #uel, rnw, #uel, rnw
    shares = [Gas_oil_splits[:imoil], Gas_oil_splits[:imgas],
    Elec_rnw_splits[:imuel], Elec_rnw_splits[:imrnw],Elec_rnw_splits[:imuel], Elec_rnw_splits[:imrnw],Elec_rnw_splits[:imuel], Elec_rnw_splits[:imrnw]],
)


########################################
# projected_det_national_data_yr (split first), or WplusAgC (split last) 
# Split commodies and sectors by proportions in different disag DataFrames
# First the tables that split both ways, sector and commodity

# WplusIntprod1 = get_subtable(WplusAgC, ["intermediate_supply","intermediate_demand"]) |> # for agg calibrate first
WplusIntprod1 = get_subtable(projected_det_national_data_yr, ["intermediate_supply","intermediate_demand"]) |>
    x -> leftjoin(
            x,
            disagprod,
            on = :commodities => :old,
            renamecols = "" => "_com"
    ) |>
    x -> transform(x,
        [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
    ) |>
    x -> leftjoin(
            x,
            disagprod,
            on = :sectors => :old,
            renamecols = "" => "_sec"
    ) |>
    x -> transform(x,
        [:sectors, :new_sec] => ByRow((a,b) -> ismissing(b) ? a : b) => :sectors
    ) |>
    x -> select(x, Not(:new_com, :new_sec)) |>
    x -> coalesce.(x, 1) |>
    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x, 
        [:value, :shares_com, :shares_sec] => ((a,b,c) -> sum(a.*b.*c))=> :value
    )

# This is for value-added because it's only split by sector (l and k are the commods)
# WplusVAprod1 = get_subtable(WplusAgC, ["value_added"]) |>
WplusVAprod1 = get_subtable(projected_det_national_data_yr, ["value_added"]) |>
    x -> leftjoin(
            x,
            disagprod,
            on = :sectors => :old,
            renamecols = "" => "_sec"
    ) |>
    x -> transform(x,
        [:sectors, :new_sec] => ByRow((a,b) -> ismissing(b) ? a : b) => :sectors
    ) |>
    x -> select(x, Not(:new_sec)) |>
    x -> coalesce.(x, 1) |>
    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x, 
        [:value, :shares_sec] => ((a,b) -> sum(a.*b))=> :value
    )

# This is for the other tables/domains things that split only by commodity.    
# Wplusprod1 = get_subtable(WplusAgC, [
Wplusprod1 = get_subtable(projected_det_national_data_yr, [
    "exogenous_final_demand",
    "personal_consumption",
    "margin_demand",
    "margin_supply",
    "tax",
    "subsidies",
    ])|>
    x -> leftjoin(
            x,
            disagprod,
            on = :commodities => :old,
            renamecols = "" => "_com"
    ) |>
    x -> transform(x,
        [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
    ) |>
    x -> coalesce.(x, 1) |>
    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x, 
        [:value, :shares_com] => ((a,b) -> sum(a.*b))=> :value
    )
# Split exports with export specific ratio
# WplusEx1 = get_subtable(WplusAgC, [
WplusEx1 = get_subtable(projected_det_national_data_yr, [
        "exports"
        ])|>
        x -> leftjoin(
                x,
                disagexp,
                on = :commodities => :old,
                renamecols = "" => "_com"
        ) |>
        x -> transform(x,
            [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
        ) |>
        x -> coalesce.(x, 1) |>
        x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
        x -> combine(x, 
            [:value, :shares_com] => ((a,b) -> sum(a.*b))=> :value
        )
# Split imports with import specific ratio
# WplusIm1 = get_subtable(WplusAgC, [
WplusIm1 = get_subtable(projected_det_national_data_yr, ["imports", "duty"  ])|>
        x -> leftjoin(
                x,
                disagimp,
                on = :commodities => :old,
                renamecols = "" => "_com"
        ) |>
        x -> transform(x,
            [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
        ) |>
        x -> coalesce.(x, 1) |>
        x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
        x -> combine(x, 
            [:value, :shares_com] => ((a,b) -> sum(a.*b))=> :value
        )

## This is to build the new sets with the additional sectors/commodities, then vcat with tables to generate a NationalTable
new_sets1 = vcat(
    # get_set(WplusAgC) |>
    get_set(projected_det_national_data_yr) |>
        x -> subset(x, :set => ByRow(!in(["sectors", "commodities"]))),
        # get_set(WplusAgC, "sectors") |>  # For agg calibrate first
        get_set(projected_det_national_data_yr, "sectors") |>    
        x -> leftjoin(
            x,
            disagprod,
            on = :element => :old,
        ) |>
        x -> transform(x,
        [:element, :new] => ByRow((a,b) -> ismissing(b) ? a : b) => :element
        ) |>
        x -> select(x, Not(:new, :shares)),
        # get_set(WplusAgC, "commodities") |>    
        get_set(projected_det_national_data_yr, "commodities") |>    
        x -> leftjoin(
            x,
            disagprod,
            on = :element => :old,
        ) |>
        x -> transform(x,
        [:element, :new] => ByRow((a,b) -> ismissing(b) ? a : b) => :element
        ) |>
        x -> select(x, Not(:new, :shares))
)

"""
### Process to link fossil inputs with fossil electricity only (uel), in production ###
Step1: int_dem. commods pet, coa, pip, ugs, oil, gas => uel (0=>rnw)
Step2: sum difference for all (need that value to be adjusted in value added to balance)
Step3: int_sup. commod rnw and uel to sector rnw, => commod rnw, rnw and uel to sector uel => commod uel
Step4: int_sup. commod ugs to sector rnw + ugs to sector uel => commod ugs to sector uel (0=>rnw)
Step5: v_a rnw (surplus and compen) up by sum difference - old ugs commod=>sector rnw
Step6: v_a uel (surplus and compen) down by (sum difference - old ugs commod) =>sector uel"""

## Step1: int_dem. commods pet, coa, pip, ugs, oil, gas => uel (0=>rnw)
# Set up the codes and shares for the join
codes = DataFrame( sects = repeat(["221100","221101"],outer=[5]), #uel, rnw
     commods = ["211000","211000", #oil 1
            "211001", "211001", #gas 2
            "212100", "212100", #coa 3
            "221200", "221200", #ugs 4
            "324110", "324110", #pet 5 # "2122A0", "2122A0"]; #min 10 # For specific oil/gas splits for sectors # repeat(["211000","211001"], outer=[4])
],
     shares = [1,0,1,0,1,0,1,0,1,0,   #  0,1,  #  .93,.07, 0,1,0,1,0,1
    ], 
     group=[:oiluelrnw,:oiluelrnw,:gasuelrnw,:gasuelrnw,:coauelrnw,:coauelrnw,:ugsuelrnw,:ugsuelrnw,:petuelrnw,:petuelrnw
     #,:minuelrnw,:minuelrnw,
    #  :petoilgas,:petoilgas,:ugsoilgas,:ugsoilgas,:pipoilgas,:pipoilgas,:ueloilgas,:ueloilgas
     ])
# Join to set up for summing etc
setup = filter(x->x.subtable=="intermediate_demand",WplusIntprod1) |> 
    x -> leftjoin(x, codes, on = [:commodities => :commods , :sectors => :sects]) |>
    x -> transform(x, :group => ByRow(x -> coalesce(x, "unchanged")); renamecols=false)
        
## Step 2: sum of values for commods moving to uel (to distribute to the VA for rnw, and subtract from uel)
diff_forVA = only(combine(groupby(setup, [:shares])[1], :value=>sum)[!,:value_sum])
# Sums of each commodity for both sectors, set to pass to uel, and 0 for rnw
groupeddem2 =   coalesce.(setup,1)|>
    x -> groupby(x, :group) |>
    x -> combine(x, :value => y->sum(y))     
# if else to get all the right values, then back as orignal df structure        
WplusIntdem2 = leftjoin(setup, groupeddem2, on=[:group])|>
    x -> transform(x, [:value,:shares,:group, :value_function] => ByRow((a,b,c,d) -> c=="unchanged" ? a : b.*d) => :fvalue)  |>
    x -> transform(x, [:fvalue] => ByRow((a) -> a) => :value) |>
    x -> select!(x,[:commodities, :sectors, :year, :subtable, :value]) 

# Step3: int_sup. commod rnw and uel to sector rnw, => commod rnw, rnw and uel to sector eul => commod uel
# Step4: int_sup. commod ugs to sector rnw + ugs to sector uel => commod ugs to sector uel (0=>rnw)
# Update Intermediate Supply, as above, and so only oil commodity is produced by oil sector and only gas commodity is produced by gas sector
comstosecs = DataFrame(sects = [repeat(["221100","221101"],outer=[3]); repeat(["211000","211001"], inner=[2])], #; #uel, rnw x 3; oil, gas x 2
    commods = [repeat(["221100","221101"],inner=[2]);["221200","221200"]; repeat(["211000", "211001"], outer=[2])], #; #uel, rnw x 2, ugs; oil, gas x 2
    shares = [1,0,0,1,1,0,1,0,0,1], # Order super important. 1st 4 for 2 sec/com matches, then uel, rnw for ugs commod, then oil com for oil sec and gas com for gas sec
    group=[[:uel, :rnw,:uel, :rnw, :ugs, :ugs]; repeat([:oilgasoil, :oilgasgas], inner=[2])])
# Set up joined df with 'unchanged' for all other sector/commods
WplusIntsup2setup =  filter(x->x.subtable=="intermediate_supply",WplusIntprod1) |># just filtering from the first split intermediate tables
    x-> leftjoin(x, comstosecs, on=[:sectors=>:sects, :commodities=>:commods]) |>
    x -> transform(x, :group => ByRow(x -> coalesce(x, "unchanged")); renamecols=false)

## Step 2: sum of values for commods moving to uel (to distribute to the VA for rnw, and subtract from uel)
diff_forVA += 
  only(filter(x->(x.sectors=="221101" && x.commodities=="221200"), WplusIntsup2setup)[!,:value])

# Sums of each commodity for both sectors, set to pass to uel, and 0 for rnw
groupedsup2 =   coalesce.(WplusIntsup2setup,1)|>
    x -> groupby(x, :group) |>
    x -> combine(x, :value => y->sum(y)) 
# Multi step to link values with shares and groups, then re-assign and select back to value column    
WplusIntsup2  = leftjoin(WplusIntsup2setup, groupedsup2, on=[:group])|>
    x -> transform(x, [:value,:shares,:group, :value_function] => ByRow((a,b,c,d) -> c=="unchanged" ? a : b.*d) => :value)  |>
    x -> select!(x,[:commodities, :sectors, :year, :subtable, :value])  

### No longer needed? Previous to oil/gas etc reallocation, balance was off, and adjusting VA to counteract was needed.    
            # ## Step5: v_a rnw (surplus and compen) up by sum difference - old ugs commod=>sector rnw
            # ## Step6: v_a uel (surplus and compen) down by (sum difference - old ugs commod) =>sector uel
            # # Set up the codes and shares (for the join to update va to counterbalance the other inputs)
            # vaset = DataFrame(sects = repeat(["221100", "221101"], outer=[2]),
            # commods = repeat(["V00100","V00300"], inner=[2]), 
            # update = [-1,1,-1,-1],
            # group = repeat(["221100", "221101"], outer=[2]))
            # # Set up the joined df with 'unchanged' as the group for everything else
            # vasetup = WplusVAprod1|>
            # x-> leftjoin(x, vaset, on=[:sectors=>:sects, :commodities=>:commods]) |>
            # x -> transform(x, :group => ByRow(x -> coalesce(x, "unchanged")); renamecols=false)
            # # Sums of each sector for both commodities (labour and kapital), to use as the denominator for proportionally increasing/decreasing each commod to equal the diff
            # groupedVA2 =   coalesce.(vasetup,1)|>
            #     x -> groupby(x, :group) |>
            #     x -> combine(x, :value => y->sum(y)) 
            # # Multi step to link values with shares and groups, then re-assign and select back to value column    
            # WplusVAprod2 =  vasetup |> 
            #     x -> leftjoin(x, groupedVA2, on =[:sectors =>:group]) |>
            #     x-> transform(x,[:value,:update, :value_function, :group]=> ByRow((a,b,c,d) -> d=="unchanged" ? a : a + a*b/c * diff_forVA) =>:value)|>
            #     x -> select!(x,[:commodities, :sectors, :year, :subtable, :value])  

# Just initial look - delete next commit
# filter(x->(x.sectors in ["211000","211001","221200","221100","324110","324121","324122","324190"] && x.commodities in ["211000", "211001"]),WplusIntdem2) # oil, gas,  
# filter(x->(x.sectors in ["211000","211001"] && x.commodities in ["211000", "211001"]),WplusIntsup2)

## Re-allocate oil and gas demanded by relevant sectors
# Set up groups and shares
commshiftcodesintdem3 = DataFrame(
sects = repeat(["211000","211001","221200","221100","324110","324121","324122","324190"], inner=[2]), #oil, gas, ugs, uel, pet,pet,pet,pet
commods = repeat(["211000", "211001"], outer=[8]), # oil, gas
shares = [1,0,0,1,0,1,0,1,1.8777,.1223,1.8777,.1223,1.8777,.1223,1.8777,.1223], 	
group = repeat([:oilgasoil, :oilgasgas, :oilgasugs, :oilgasuel, :oilgaspet1, :oilgaspet2, :oilgaspet3, :oilgaspet4], inner=[2]) )
# Set up the joined df with 'unchanged' as the group for everything else
WplusIntdem3setup = leftjoin(WplusIntdem2, commshiftcodesintdem3, on = [:sectors=>:sects, :commodities=>:commods] ) |>
    x -> transform(x, :group => ByRow(x -> coalesce(x, "unchanged")); renamecols=false)
# Sums of each commodity for both sectors, as basis for shares
groupedIntdem3 =   coalesce.(WplusIntdem3setup,1)|>
    x -> groupby(x, :group) |>
    x -> combine(x, :value => y->sum(y)) 
# Multi step to link values with shares and groups, then re-assign and select back to value column    
WplusIntdem3 =  WplusIntdem3setup |> 
    x -> leftjoin(x, groupedIntdem3, on =[:group =>:group]) |>
x-> transform(x,[:value,:shares, :value_function, :group]=> ByRow((a,b,c,d) -> d=="unchanged" ? a : b * c) =>:value)|>
x -> select!(x,[:commodities, :sectors, :year, :subtable, :value])  

########################################
# TODONE!! fuction to split individual combinations, and then add a not in the rest of the splits
# oil and gas
## And other splits: "intermediate_demand": pet [.93,.07]
# "intermediate_demand": pip [.0,1 ]
# "intermediate_demand": ugs [.0,1 ]
# "intermediate_demand": uel [.0,1 ]

# "intermediate_supply": oil [.0,1]
# "intermediate_supply": gas [.0,1]
########################################

## All splits done, put it all back together
WplusSp1 = NationalTable(
    vcat(WplusIntsup2,WplusIntdem3,
    WplusVAprod1, Wplusprod1, WplusEx1, WplusIm1),
    new_sets1)

# WiNDC plus = Multiple GHG update from national summary map (uti->[:uel,:ugs,:uwt], min->[:coa,:min])
# summary_map = CSV.read("summarywindcplus_detail.csv", DataFrame)
# # Aggregate (WiNDC function) to WiNDC plus before oil/gas split 
WplusSpAg = aggregate(
    WplusSp1,
    :commodities => (summary_map, :detailed => :summary),
    :sectors => (summary_map, :detailed => :summary),
    :margin_demand => (summary_map, :detailed => :summary),
    :margin_supply => (summary_map, :detailed => :summary),
    :labor_demand => (summary_map, :detailed => :summary),
    :exogenous_final_demand => (summary_map, :detailed => :summary),
    :other_tax => (summary_map, :detailed => :summary),
    :subsidies => (summary_map, :detailed => :summary),
    :duty => (summary_map, :detailed => :summary),
    :exports => (summary_map, :detailed => :summary),
    :imports => (summary_map, :detailed => :summary),
    :tax => (summary_map, :detailed => :summary),
    :capital_demand => (summary_map, :detailed => :summary),
    :personal_consumption => (summary_map, :detailed => :summary),
)

WplusSpAgC,M = calibrate(WplusSpAg)

# M5 = national_mpsge(WplusSpAgC) 
# solve!(M5, cumulative_iteration_limit = 0)
# print(sort(generate_report(M5),:value))


# wp= [WplusCSpAg,WplusCAgSp,WplusSpCAg,WplusSpAgC,WplusAgCSp,WplusAgSpC]
# # a=3
# m = wp[a]
# M5 = national_mpsge(m) 
# solve!(M5, cumulative_iteration_limit = 0);a+=1
# print(sort(generate_report(M5),:margin))

#########################################
# Functions to transform from dfs to NamedArrays, each has own dimensions and elements
#########################################
function table_to_naCxS(
    X::WiNDCtable,
    subtable::String;
    columns_to_ignore = [:year])
    domain_columns = String.([a for a∈WiNDC.domain(X) if a∉columns_to_ignore])
    sets = get_set.(Ref(X), domain_columns) |>
            x -> map(y -> Symbol.(y[!,:element]), x)     
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns))
    for row in eachrow(get_subtable(X, subtable))
        idx = Tuple(row[domain_columns])
        out[Symbol.(idx)...] = row[:value]
    end
    return out
end

function table_to_naCxSecs(
    X::WiNDCtable,
    subtable::String,
    CorS::String;
    columns_to_ignore = [:year])
    domain_columns = String.([a for a∈domain(X) if a∉columns_to_ignore])
    sets = [Symbol.(filter(x->x ∉["use","oth"],get_set(X,CorS)[:,:element])),
    Symbol.(unique(get_subtable(X,subtable),:sectors)[:,:sectors])]
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns))
    for row = eachrow(get_subtable(X, subtable))
        idx = Tuple(row[domain_columns])
        out[Symbol.(idx)...] = row[:value]
    end
    return out
end

function table_to_nafn(
    X::WiNDCtable,
    subtable::DataFrame,
    CorS::String;
    columns_to_ignore = ["year"])
    domain_columns = String.([a for a in names(subtable) if a∉columns_to_ignore])
    sets = [Symbol.(filter(x->x ∉["use","oth"],get_set(X,CorS)[:,:element])), [:value]]
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns))
    for row in eachrow(subtable)
        idx = row[CorS]
        out[Symbol(idx),:value] = row[:value]
    end
    return out
end

function table_to_naSecsxC(
    X::WiNDCtable,
    subtable::String,
    CorS::String;
    columns_to_ignore = [:year])
    domain_columns = String.([a for a∈domain(X) if a∉columns_to_ignore])
    sets = [Symbol.(unique(get_subtable(X,subtable),:commodities)[:,:commodities]),Symbol.(filter(x->x ∉["use","oth"],get_set(X,CorS)[:,:element]))]
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns))
    for row = eachrow(get_subtable(X, subtable))
        idx = Tuple(row[domain_columns])
        out[Symbol.(idx)...] = row[:value]
    end
    return out
end
# a=1
# wp= [WplusCSpAg,WplusCAgSp,WplusSpCAg,WplusSpAgC,WplusAgCSp,WplusAgSpC] WplusSp2 WplusSpAgC
# table_to_naCxS(wp[a], "intermediate_demand")
# :($wp[a])
id_m0 = table_to_naCxS(WplusSpAgC, "intermediate_demand");
ys_m0 = table_to_naCxS(WplusSpAgC, "intermediate_supply");
fd_m0 = table_to_naCxSecs(WplusSpAgC,"exogenous_final_demand","commodities");
pce_m0 = table_to_naCxSecs(WplusSpAgC,"personal_consumption","commodities");
m_m0 = table_to_naCxSecs(WplusSpAgC,"imports","commodities");
x_m0 = table_to_naCxSecs(WplusSpAgC,"exports","commodities");
md_m0 = table_to_naCxSecs(WplusSpAgC,"margin_demand","commodities");
ms_m0 = table_to_naCxSecs(WplusSpAgC,"margin_supply","commodities");
y_m0    = table_to_nafn(WplusSpAgC,WiNDC.gross_output(WplusSpAgC),"commodities");
a_m0    = table_to_nafn(WplusSpAgC,WiNDC.armington_supply(WplusSpAgC),"commodities");
ty_m0   = table_to_nafn(WplusSpAgC,WiNDC.other_tax_rate(WplusSpAgC),"sectors"); # sectors only
ta_m0    = table_to_nafn(WplusSpAgC,WiNDC.absorption_tax_rate(WplusSpAgC),"commodities");
tm_m0   = table_to_nafn(WplusSpAgC,WiNDC.import_tariff_rate(WplusSpAgC),"commodities");
bopdef_m0  = NamedArray(WiNDC.balance_of_payments(WplusSpAgC)[:,:value],([Symbol(year)]),(:bop)) # not commodities or sectors
va_m0 = table_to_naSecsxC(WplusSpAgC, "value_added","sectors");

# WiNDC2022data2022 = Dict(
# WplusAgCSpdata2022 = Dict(
WplusSpAgCdata2022 = Dict(
  :a_0      => a_m0,
  :id_0     => id_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
  :ys_0     => ys_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
  :tm_0     => tm_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
  :va_0     => va_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
  :md_0     => md_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
  :fd_0     => fd_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
  :pce_0    => pce_m0,
  :m_0      => m_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
  :ty_0     => ty_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
  :ms_0     => ms_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
  :bopdef_0 => bopdef_m0, # 1-dimensional DenseAxisArray{Float64,1,...} 
  :x_0      => x_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
  :ta_0     => ta_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
  :y_0      => y_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
)

# ### MPSGE DataFrame format ###
# outerjoin(
#     WiNDC.armington_supply(WplusSpAgC; output = :as),
#     get_subtable(WplusSpAgC, "exports", output = :ex) |>
#         x -> select(x, Not(:sectors)),
#     get_subtable(WplusSpAgC, "margin_demand", output = :md) |>
#         x -> unstack(x, :sectors, :md), # :Trade, :Trans
#     WiNDC.gross_output(WplusSpAgC; output = :go),
#     get_subtable(WplusSpAgC, "imports", output = :im)|>
#         x -> select(x, Not(:sectors)),
#     WiNDC.absorption_tax_rate(WplusSpAgC, output = :atr),
#     WiNDC.import_tariff_rate(WplusSpAgC, output = :itr),
#     on = filter(y -> y!=:sectors, domain(WplusSpAgC))
# ) |>
# x -> coalesce.(x, 0) 


# Apparently there is an error here using the summary data. I have the value added
# categories hard coded... which is super smart of me! I should fix that.
# M = national_mpsge(national_data)
# Benchmark
# solve!(M, cumulative_iteration_limit = 0)

# There are three parameters you can experiment with
# Absorption_tax[commodities]
# Output_tax[sectors]
# Margin_tax[commodities]

#These default to their values in the dataset. 

# set_value!.(M[:Output_tax], .5*value.(M[:Output_tax]))

# solve!(M)

# generate_report(M)

########################
## Exploring the data ##
########################

# See the entire table. Most of the columns should be self-explanatory. The really
# important one is `subtable`, which is the human readable component of the table.
# get_table(callibrated_det_national_data)

# # You can extract a single subtable. Note the subtable column is dropped (should 
# # be fixed, I don't care for that).
# get_subtable(callibrated_det_national_data, "intermediate_demand")
# get_subtable(callibrated_det_national_data, "intermediate_demand")

# # Or multiple
# get_subtable(callibrated_det_national_data, ["intermediate_demand", "value_added"])

# # Turns out, that's in a branch... quick hack, but it drops subtables... 
# vcat(get_subtable.(Ref(callibrated_det_national_data), ["intermediate_demand", "value_added"])...)

# # You should notice there is no subtable for value added.
# get_table(callibrated_det_national_data) |>
#     x -> unique(x, :subtable)
# # This is because value added is a composite table, it consists of several smaller
# # subtables. To see which ones use get_set, the elements comprise value added.
# get_set(callibrated_det_national_data, "value_added")
# # If we look at labor_demand we will see two elements, labor demand and V00100. 
# # The V00100 is the NAICS code for labor compensation. 
# get_set(callibrated_det_national_data, "labor_demand")

# # example correspondence table
# # Currently has WiNDC_plus and WiNDC summary (theoretically to replicate summary but doesn't), but can be any abitraray mapping been BEA_detail, and any aggregation set
# Codes=CSV.read(joinpath(raw_data_directory,"./BEA_WiNDC_Detail-Summary_codes.csv"), DataFrame, header=1)

# #######################################################################################
# # X is the correspondence table between bea detail and whatever aggregation is set up #
# #######################################################################################
# # X = DataFrame([Codes[:,:BEA_detail], Codes[:,:WiNDC]], [:commodities, :name])
# X = DataFrame([Codes[:,:BEA_detail], Codes[:,:WiNDC_plus]], [:commodities, :name])
# aggcol = :WiNDC_plus # name of the columns just for a value comparison dataframe later
# # aggcol = :WiNDC # name of the columns just for a value comparison dataframe later

# ## To callibrate first=> get_table(callibrated_det_national_data): or to aggregate and then callibrate -> get_table(projected_det_national_data_yr) and callibrate below
# double_aggregate_subtables = ["intermediate_supply", "intermediate_demand","labor_demand","capital_demand","other_tax"]
# TMP = NationalTable( vcat(
#     get_table(callibrated_det_national_data) |>
#     x -> filter(:subtable => !(in(double_aggregate_subtables)), x) |> #ggregate the simple tables separately
#     x -> leftjoin(
#         x,
#         X,
#         on = :commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform(x,
#         [:commodities, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :commodities,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
#     x -> combine(x, :value => sum => :value),
# # Now aggregate the tables that also need to aggregate sectors, in 2 steps
#     get_table(callibrated_det_national_data) |>
#     x -> filter(:subtable => (in(double_aggregate_subtables)), x) |>
#     x -> leftjoin(
#         x,
#         X,
#         on = :commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform(x,
#         [:commodities, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :commodities,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
#     x -> combine(x, :value => sum => :value) |>
#     x -> leftjoin(
#     x,
#     X,
#     on = :sectors=>:commodities,
#     renamecols = "" => "_com") |>
#     x -> transform!(x,
#         [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
#     x -> combine(x, :value => sum => :value)
#     ),
#         callibrated_det_national_data.sets
# )


# subs_TMP = get_table(TMP) |>
#     x -> unique(x, :subtable)
# #print the numbers for raw and aggregated
# for i in subs_TMP[:,:subtable]
#     println(i,": ",length(unique(get_subtable(TMP,i)[!,:commodities])),": ",length(unique(get_subtable(projected_det_national_data_yr,i)[!,:commodities]))) # Compared to get_subtable(projected_det_national_data_yr,i)
# end

# # JLD2.save("WiNDCplus.jld2",MultiNatdata)
