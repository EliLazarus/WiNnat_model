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

Gas_oil_splits=Dict(:prodgas=>0.4003, :prodoil=>0.5997,
:imgas=>0.4308,:imoil=>0.5692,:exgas=>0.3961,:exoil=>0.6039)
## Split before Ag
# disag50_50 = DataFrame( #for testing
#     old = ["211000", "211000"], # oil, oil
#     new = ["211000", "211001"], # oil, gas
#     shares = [.5,.5],)
disagprod = DataFrame(
    old = ["211000", "211000"], # oil, oil
    new = ["211000", "211001"], # oil, gas
    shares = [Gas_oil_splits[:prodoil], Gas_oil_splits[:prodgas]],
)
disagexp = DataFrame(
    old = ["211000", "211000"], # oil, oil
    new = ["211000", "211001"], # oil, gas
    shares = [Gas_oil_splits[:exoil], Gas_oil_splits[:exgas]],
)
disagimp = DataFrame(
    old = ["211000", "211000"], # oil, oil
    new = ["211000", "211001"], # oil, gas
    shares = [Gas_oil_splits[:imoil], Gas_oil_splits[:imgas]],
)

## Secondary Split for uel to uel and rnw
disagprod2 = DataFrame(
    old = ["221100", "221100"], #uel,uel  
    new = ["221100", "221101"], #uel, rnw
    shares = [.2,.8], # just nominal so far
)
disagexp2 = DataFrame(
    old = ["221100", "221100"], #uel,uel
    new = ["221100", "221101"], #uel, rnw
    shares = [.4,.6], # just nominal so far
 
)
disagimp2 = DataFrame(
    old = ["221100", "221100"], #uel,uel
    new = ["221100", "221101"], #uel, rnw
    shares = [.6,.4], # just nominal so far
)

# ### For split after aggregations
# disagprod = DataFrame(
#     old = ["oil", "oil"],
#     new = ["oil", "gas"],
#     # shares = [Gas_oil_splits[:prodoil], Gas_oil_splits[:prodgas]],
# )
# disagexp = DataFrame(
#     old = ["oil", "oil"],
#     new = ["oil", "gas"],
# # shares = [Gas_oil_splits[:exoil], Gas_oil_splits[:exgas]],
# )
# disagimp = DataFrame(
#     old = ["oil", "oil"],
#     new = ["oil", "gas"],
#     shares = [.5,.5],
# # shares = [Gas_oil_splits[:imoil], Gas_oil_splits[:imgas]],
# )

# ## Secondary Split for after aggregation uel to uel and rnw
# disagprod2 = DataFrame(
#     old = ["uel", "uel"], #uel,uel  
#     new = ["uel", "rnw"], #uel, rnw
#     shares = [.5,.5],
# )
# disagexp2 = DataFrame(
#     old = ["uel", "uel"], #uel,uel
#     new = ["uel", "rnw"], #uel, rnw
#     shares = [.5,.5],
# )
# disagimp2 = DataFrame(
#     old = ["uel", "uel"], #uel,uel
#     new = ["uel", "rnw"], #uel, rnw
#     shares = [.5,.5],
# )

########################################
#!! TODO fuction to split individual combinations, and then add a not in the rest of the splits

## And other splits: "intermediate_demand": pet [.93,.07]
# "intermediate_demand": pip [.0,1 ]
# "intermediate_demand": ugs [.0,1 ]
# "intermediate_demand": uel [.0,1 ]
# "intermediate_demand" (all the rest prodoil and prodgas)

# intermediate_supply" : pet [.93,.07]
# "intermediate_supply" (all the rest prodoil and prodgas)
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

WplusSp1 = NationalTable(
    vcat(WplusIntprod1,WplusVAprod1, Wplusprod1, WplusEx1, WplusIm1),
    new_sets1
)

############################
### This is the secondary split, for uel into uel and rnw
############################
WplusIntprod2 = get_subtable(WplusSp1, ["intermediate_supply","intermediate_demand"]) |>
    x -> leftjoin( x,disagprod2,on = :commodities => :old,renamecols = "" => "_com"
    ) |>    x -> transform(x,        [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
    ) |>    x -> leftjoin(x,disagprod2,on = :sectors => :old,renamecols = "" => "_sec"
    ) |>    x -> transform(x,        [:sectors, :new_sec] => ByRow((a,b) -> ismissing(b) ? a : b) => :sectors
    ) |>    x -> select(x, Not(:new_com, :new_sec)) |>
    x -> coalesce.(x, 1) |>    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x,         [:value, :shares_com, :shares_sec] => ((a,b,c) -> sum(a.*b.*c))=> :value
    )
# This is for value-added because it's only split by sector (l and k are the commods)
WplusVAprod2 = get_subtable(WplusSp1, ["value_added"]) |>
    x -> leftjoin(x,disagprod2,on = :sectors => :old,renamecols = "" => "_sec"
    ) |>    x -> transform(x,        [:sectors, :new_sec] => ByRow((a,b) -> ismissing(b) ? a : b) => :sectors
    ) |>    x -> select(x, Not(:new_sec)) |>    x -> coalesce.(x, 1) |>    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x,         [:value, :shares_sec] => ((a,b) -> sum(a.*b))=> :value
    )
# This is for the other tables/domains things that split only by commodity.    
Wplusprod2 = get_subtable(WplusSp1, ["exogenous_final_demand","personal_consumption","margin_demand","margin_supply","tax","subsidies"    ])|>
    x -> leftjoin(x,disagprod2,on = :commodities => :old,renamecols = "" => "_com"
    ) |>    x -> transform(x,        [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
    ) |>    x -> coalesce.(x, 1) |>    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x,         [:value, :shares_com] => ((a,b) -> sum(a.*b))=> :value    )
# Split exports with export specific ratio
 WplusEx2 = get_subtable(WplusSp1, ["exports"])|>
        x -> leftjoin(    x,    disagexp2,    on = :commodities => :old,    renamecols = "" => "_com"        ) |>
        x -> transform(x,[:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
        ) |>        x -> coalesce.(x, 1) |>        x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
        x -> combine(x, [:value, :shares_com] => ((a,b) -> sum(a.*b))=> :value        )
# Split imports with import specific ratio
 WplusIm2 = get_subtable(WplusSp1, [ "imports","duty" ])|>
        x -> leftjoin(    x,    disagimp2,    on = :commodities => :old,    renamecols = "" => "_com"    ) |>
        x -> transform(x,   [:commodities, :new_com] => ByRow((a,b) -> ismissing(b) ? a : b) => :commodities
        ) |>        x -> coalesce.(x, 1) |>        x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
        x -> combine(x,             [:value, :shares_com] => ((a,b) -> sum(a.*b))=> :value        )
## This is to build the new sets with the additional sectors/commodities, then vcat with tables to generate a NationalTable
new_sets2 = vcat(
    get_set(WplusSp1) |>
        x -> subset(x, :set => ByRow(!in(["sectors", "commodities"]))),
    get_set(WplusSp1, "sectors") |>    
        x -> leftjoin(
            x,
            disagprod2,
            on = :element => :old,
        ) |>
        x -> transform(x,
        [:element, :new] => ByRow((a,b) -> ismissing(b) ? a : b) => :element
        ) |>
        x -> select(x, Not(:new, :shares)),
    get_set(WplusSp1, "commodities") |>    
        x -> leftjoin(
            x,
            disagprod2,
            on = :element => :old,
        ) |>
        x -> transform(x,
        [:element, :new] => ByRow((a,b) -> ismissing(b) ? a : b) => :element
        ) |>
        x -> select(x, Not(:new, :shares))
)

WplusSp2 = NationalTable(
    vcat(WplusIntprod2,WplusVAprod2, Wplusprod2, WplusEx2, WplusIm2),
    new_sets2
)

# WiNDC plus = Multiple GHG update from national summary map (uti->[:uel,:ugs,:uwt], min->[:coa,:min])
# summary_map = CSV.read("summarywindcplus_detail.csv", DataFrame)
# # Aggregate (WiNDC function) to WiNDC plus before oil/gas split 
WplusSpAg = aggregate(
    WplusSp2,
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
# [4.9571e+00,4.9571e+00,6.9223e-12,5.1176e-12,4.9571e+00,6.4337e-12]
#########################################
# Functions to transform from dfs to NamedArrays, each has own dimensions and elements
#########################################
function table_to_naCxS(
    X::WiNDCtable,
    subtable::String;
    columns_to_ignore = [:year]
)
    domain_columns = String.([a for a∈WiNDC.domain(X) if a∉columns_to_ignore])
    sets = get_set.(Ref(X), domain_columns) |>
            x -> map(y -> Symbol.(y[!,:element]), x)     
            
            
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns)
    )
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
    columns_to_ignore = [:year]
)
    domain_columns = String.([a for a∈domain(X) if a∉columns_to_ignore])
    sets = [Symbol.(filter(x->x ∉["use","oth"],get_set(X,CorS)[:,:element])),
    Symbol.(unique(get_subtable(X,subtable),:sectors)[:,:sectors])]
     
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns)
    )
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
    columns_to_ignore = ["year"]
)
    domain_columns = String.([a for a in names(subtable) if a∉columns_to_ignore])
    sets = [Symbol.(filter(x->x ∉["use","oth"],get_set(X,CorS)[:,:element])), [:value]]
     
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns)
    )
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
    columns_to_ignore = [:year]
)
    domain_columns = String.([a for a∈domain(X) if a∉columns_to_ignore])
    sets = [Symbol.(unique(get_subtable(X,subtable),:commodities)[:,:commodities]),Symbol.(filter(x->x ∉["use","oth"],get_set(X,CorS)[:,:element]))]
     
    out = NamedArray(
        zeros(length.(sets)...),
        Tuple(sets),
        Tuple(domain_columns)
    )
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

### MPSGE DataFrame format ###
outerjoin(
    WiNDC.armington_supply(WplusSpAgC; output = :as),
    get_subtable(WplusSpAgC, "exports", output = :ex) |>
        x -> select(x, Not(:sectors)),
    get_subtable(WplusSpAgC, "margin_demand", output = :md) |>
        x -> unstack(x, :sectors, :md), # :Trade, :Trans
    WiNDC.gross_output(WplusSpAgC; output = :go),
    get_subtable(WplusSpAgC, "imports", output = :im)|>
        x -> select(x, Not(:sectors)),
    WiNDC.absorption_tax_rate(WplusSpAgC, output = :atr),
    WiNDC.import_tariff_rate(WplusSpAgC, output = :itr),
    on = filter(y -> y!=:sectors, domain(WplusSpAgC))
) |>
x -> coalesce.(x, 0) 


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

#     ########################
#     ### Original hacky and messy oil => gas & oil ###
#     ########################
#     Gas_oil_splits=Dict(:imgas=>0.4308, :exgas=>0.3961, :prodgas=>0.4003, :imoil=>0.5692,:exoil=>0.6039,:prodoil=>0.5997)
#     ### This has to be done in the main dataframe, so it's pickup up with the tax functions etc.
#     # Add gas and give value of oil & gas x gas %
#     push!(TMP.table, (subtable="exports",commodities="gas", sectors=only(filter(x->(x.subtable=="exports") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="exports") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="exports").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:exgas]))
#     # update oil value to value of oil & gas x oil %
#     TMP.table[(TMP.table.subtable.=="exports").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="exports").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:exoil]
#     push!(TMP.table, (subtable="imports",commodities="gas", sectors=only(filter(x->(x.subtable=="imports") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="imports") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="imports").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:imgas]))
#     TMP.table[(TMP.table.subtable.=="imports").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="imports").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:imoil]
#     # -	“duty” w imports
#     push!(TMP.table, (subtable="duty",commodities="gas", sectors=only(filter(x->(x.subtable=="duty") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="duty") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="duty").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:imgas]))
#     TMP.table[(TMP.table.subtable.=="duty").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="duty").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:imoil]
    
#     # #### no oil&gas margin_supply, nothing to split 
#     # # push!(ms_tmp, (commodities=:gas, sectors=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trn),ms_tmp)[:,2]),year=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trn),ms_tmp)[:,3]),value=only(ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trn),:value])*Gas_oil_splits[:prodgas]))
#     # # push!(ms_tmp, (commodities=:gas, sectors=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trd),ms_tmp)[:,2]),year=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trd),ms_tmp)[:,3]),value=only(ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trd),:value])*Gas_oil_splits[:prodgas]))
#     # ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trn),:value]=ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trn),:value]*Gas_oil_splits[:prodoil]
#     # ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trd),:value]=ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trd),:value]*Gas_oil_splits[:prodoil]
#     push!(TMP.table, (subtable="margin_demand",commodities="gas", sectors=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRANS"),TMP.table)[:,2]),year=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRANS"),TMP.table)[:,3]),value=only(TMP.table[(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRANS"),:value])*Gas_oil_splits[:prodgas]))
#     push!(TMP.table, (subtable="margin_demand",commodities="gas", sectors=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRADE "),TMP.table)[:,2]),year=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRADE "),TMP.table)[:,3]),value=only(TMP.table[(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRADE "),:value])*Gas_oil_splits[:prodgas]))
#     TMP.table[(TMP.table.subtable.=="margin_demand").&(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRANS"),:value] =TMP.table[(TMP.table.subtable.=="margin_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRANS"),:value]*Gas_oil_splits[:prodoil]
#     TMP.table[(TMP.table.subtable.=="margin_demand").&(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRADE "),:value]=TMP.table[(TMP.table.subtable.=="margin_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRADE "),:value]*Gas_oil_splits[:prodoil]

#     # "personal_consumption" # 0 ie. not in TMP, nothing to split (this time)
#     # push!(TMP.table, (subtable="personal_consumption",commodities="gas", sectors=only(filter(x->(x.subtable=="personal_consumption") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="personal_consumption") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="personal_consumption").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:consgas]))
#     # TMP.table[(TMP.table.subtable.=="personal_consumption").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="personal_consumption").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:consoil]
    
#     # loop over final demand sectors in data, for each with data, add gas for each at oil x gas %, and multiply oil by oil %
#     for s in unique(TMP.table[TMP.table.subtable.=="exogenous_final_demand",:sectors])
#         println(s)
#         if size(TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s) .& (TMP.table.commodities.=="oil"),:],1)>0  # But actually just
#         push!(TMP.table, (subtable="exogenous_final_demand",commodities="gas", sectors=only(filter(x->(x.subtable=="exogenous_final_demand") & (x.sectors==s) & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="exogenous_final_demand") & (x.sectors==s) & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s) .& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:prodgas]))
#         TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s) .& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s).& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:prodoil]
#     end
#     end
#     push!(TMP.table, (subtable="labor_demand",sectors="gas", commodities=only(filter(x->(x.subtable=="labor_demand") & (x.sectors=="oil"),TMP.table)[:,1]),year=only(filter(x->(x.subtable=="labor_demand") & (x.sectors=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="labor_demand").& (TMP.table.sectors.=="oil"),:value])*Gas_oil_splits[:prodgas]))
#     TMP.table[(TMP.table.subtable.=="labor_demand").& (TMP.table.sectors.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="labor_demand").& (TMP.table.sectors.=="oil"),:value]*Gas_oil_splits[:prodoil]
#     push!(TMP.table, (subtable="capital_demand",sectors="gas", commodities=only(filter(x->(x.subtable=="capital_demand") & (x.sectors=="oil"),TMP.table)[:,1]),year=only(filter(x->(x.subtable=="capital_demand") & (x.sectors=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="capital_demand").& (TMP.table.sectors.=="oil"),:value])*Gas_oil_splits[:prodgas]))
#     TMP.table[(TMP.table.subtable.=="capital_demand").& (TMP.table.sectors.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="capital_demand").& (TMP.table.sectors.=="oil"),:value]*Gas_oil_splits[:prodoil]
#     # 👍   

#     # "intermediate_demand"
#     print(rightjoin(Sectors[:,1:2], TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="oil"),:], on=:index=>:commodities) )  # demand of commodities by oil sector: split by production
    
#     checksumsec = sum(filter(x->x.sectors=="oil",get_subtable(TMP,"intermediate_demand"))[:,:value])
#     checksumcom = sum(filter(x->x.commodities=="oil",get_subtable(TMP,"intermediate_demand"))[:,:value])
#     checksumall = sum(get_subtable(TMP,"intermediate_demand")[:,:value])
#     # Split sectors by demand by oil
#     for i in TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="oil"),:][!,:commodities];
#         push!(TMP.table, (commodities=i, sectors="gas", year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="oil") .& (TMP.table.commodities .==i),:][!,:value])*Gas_oil_splits[:prodgas]))
#         TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.sectors.=="oil") .& (TMP.table.commodities.==i),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.sectors.=="oil") .& (TMP.table.commodities.==i),:value]*Gas_oil_splits[:prodoil]
#     end
 
#     print(rightjoin(Sectors[:,1:2], TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="gas"),:], on=:index=>:commodities) )
    
#     # Check that oil is a commodity & how many, then gas should be 0
#     TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.commodities.=="oil"),:] # demand of gas/oil by sectors pet: .07/0.93, pip: 1/0.0, ugs: 1/0.0, uel:1/0.0, trn: 0.0/1, all other by production
#     TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.commodities.=="gas"),:] # demand of gas/oil by sectors pet: .07/0.93, pip: 1/0.0, ugs: 1/0.0, uel:1/0.0, trn: 0.0/1, all other by production
#     @assert(nrow(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.commodities.=="gas"),:])==0)
#     # Split commodities by demand for gas/oil
#     # Split the oil commodity to oil & gas for pet, pip, ugs, and uel by specific ratios
#     push!(TMP.table, (commodities = "gas", sectors="pet", year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="pet") .& (TMP.table.commodities .=="oil"),:][!,:value])*0.07)) #TODO add reference
#     TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="pet"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="pet"),:value]*0.93
#     push!(TMP.table, (commodities = "gas", sectors="pip", year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="pip") .& (TMP.table.commodities .=="oil"),:][!,:value])*1)) #TODO add reference
#     TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="pip"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="pip"),:value]*0.0
#     push!(TMP.table, (commodities = "gas", sectors="ugs", year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="ugs") .& (TMP.table.commodities .=="oil"),:][!,:value])*1)) #TODO add reference
#     TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="ugs"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="ugs"),:value]*0.0
#     push!(TMP.table, (commodities = "gas", sectors="uel", year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="uel") .& (TMP.table.commodities .=="oil"),:][!,:value])*1)) #TODO add reference
#     TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="uel"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="uel"),:value]*0.0
#     # Zero demand for gas/oil by trn transport
#     # push!(TMP.table, (commodities = "gas", sectors="trn", year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.sectors.=="trn") .& (TMP.table.commodities .=="oil"),:][!,:value])*0.0)) #TODO add reference
#     # TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="trn"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="trn"),:value]*1.0
# # Then split the oil commodity to oil & gas for all the other sectors    
#     for i in [i for i in TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.commodities.=="oil"),:][!,:sectors] if i ∉ ["pet", "pip", "ugs", "uel"]]
#         push!(TMP.table, (commodities="gas", sectors=i, year=parse(Int,string(year)), subtable="intermediate_demand", value=only(TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.commodities .=="oil") .& (TMP.table.sectors.==i),:][!,:value])*Gas_oil_splits[:prodgas]))
#         TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.sectors.==i) .& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.==i),:value]*Gas_oil_splits[:prodoil]
#     end

#     TMP.table[(TMP.table.subtable.=="intermediate_demand") .& (TMP.table.commodities.=="gas"),:] 
# # Make sure that totals for the subtable are still correct.
#     @assert(isapprox(checksumsec, sum(filter(x->x.sectors=="oil",get_subtable(TMP,"intermediate_demand"))[:,:value])+sum(filter(x->x.sectors=="gas",get_subtable(TMP,"intermediate_demand"))[:,:value]), atol=10^-12),"uh oh, post-split sum of sectors in int_dem doesn't match")
#     @assert(isapprox(checksumcom , sum(filter(x->x.commodities=="oil",get_subtable(TMP,"intermediate_demand"))[:,:value])+sum(filter(x->x.commodities=="gas",get_subtable(TMP,"intermediate_demand"))[:,:value]), atol=10^-12),"uh oh, post-split sum of commodities in int_dem doesn't match")
#     @assert(isapprox(checksumall , sum(get_subtable(TMP,"intermediate_demand")[:,:value]), atol=10^-11))

#     # "intermediate_supply"
#     checksumsec = sum(filter(x->x.sectors=="oil",get_subtable(TMP,"intermediate_supply"))[:,:value])
#     checksumcom = sum(filter(x->x.commodities=="oil",get_subtable(TMP,"intermediate_supply"))[:,:value])
# # Split oil sector into oil & gas sectors by supply of commodities
# # First just pet with a specific ratio
#     push!(TMP.table, (commodities = "pet", sectors="gas", year=parse(Int,string(year)), subtable="intermediate_supply", value=only(TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.sectors.=="oil") .& (TMP.table.commodities .=="pet"),:][!,:value])*0.07)) #TODO add reference
#     TMP.table[(TMP.table.subtable.=="intermediate_supply").& (TMP.table.commodities.=="pet") .& (TMP.table.sectors.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="intermediate_supply").& (TMP.table.commodities.=="pet") .& (TMP.table.sectors.=="oil"),:value]*0.93
# # Then the rest with ratio of production
#     for i in [i for i in TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.sectors.=="oil"),:][!,:commodities] if i ∉ ["pet"]]
#         push!(TMP.table, (commodities=i, sectors="gas", year=parse(Int,string(year)), subtable="intermediate_supply", value=only(TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.sectors.=="oil") .& (TMP.table.commodities .==i),:][!,:value])*Gas_oil_splits[:prodgas]))
#         TMP.table[(TMP.table.subtable.=="intermediate_supply").& (TMP.table.sectors.=="oil") .& (TMP.table.commodities.==i),:value]=TMP.table[(TMP.table.subtable.=="intermediate_supply").& (TMP.table.sectors.=="oil") .& (TMP.table.commodities.==i),:value]*Gas_oil_splits[:prodoil]
#     end
# # Split oil/gas commodity from all sectors into 2
# for i in TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.commodities.=="oil"),:][!,:sectors]
#     push!(TMP.table, (sectors=i, commodities="gas", year=parse(Int,string(year)), subtable="intermediate_supply", value=only(TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.commodities.=="oil") .& (TMP.table.sectors .==i),:][!,:value])*Gas_oil_splits[:prodgas]))
#     TMP.table[(TMP.table.subtable.=="intermediate_supply").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.==i),:value]=TMP.table[(TMP.table.subtable.=="intermediate_supply").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.==i),:value]*Gas_oil_splits[:prodoil]
# end

#    rightjoin(Sectors[:,1:2], TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.sectors.=="oil" .|| TMP.table.sectors.=="gas"),:], on=:index=>:commodities)  # supply of commodities from oil sector, what oil produces: pet: .07/0.93, all other by production
#    TMP.table[(TMP.table.subtable.=="intermediate_supply") .& (TMP.table.commodities.=="oil" .|| TMP.table.commodities.=="gas"),:] # supply of gas/oil, from each sector. all by production

#    @assert(isapprox(checksumsec , sum(filter(x->x.sectors=="oil",get_subtable(TMP,"intermediate_supply"))[:,:value])+sum(filter(x->x.sectors=="gas",get_subtable(TMP,"intermediate_supply"))[:,:value]), atol=10^-12),"uh oh, post-split sum of sectors in int_supp doesn't match")
#    @assert(isapprox(checksumcom , sum(filter(x->x.commodities=="oil",get_subtable(TMP,"intermediate_supply"))[:,:value])+sum(filter(x->x.commodities=="gas",get_subtable(TMP,"intermediate_supply"))[:,:value]), atol=10^12),"uh oh, post-split sum of commodities in int_supp doesn't match")

# #    -	“other_tax", w int_supply
#    push!(TMP.table, (subtable="other_tax",sectors="gas", commodities=only(filter(x->(x.subtable=="other_tax") & (x.sectors=="oil"),TMP.table)[:,1]),year=only(filter(x->(x.subtable=="other_tax") & (x.sectors=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="other_tax").& (TMP.table.sectors.=="oil"),:value])*Gas_oil_splits[:prodgas]))
#    TMP.table[(TMP.table.subtable.=="other_tax").& (TMP.table.sectors.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="other_tax").& (TMP.table.sectors.=="oil"),:value]*Gas_oil_splits[:prodoil]
    
# #    -	“tax” and “subsidies”, w arm_supply 
#     push!(TMP.table, (subtable="tax",commodities="gas", sectors=only(filter(x->(x.subtable=="tax") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="tax") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="tax").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:prodgas]))
#     TMP.table[(TMP.table.subtable.=="tax").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="tax").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:prodoil]
# # ###################################
# # # Joins to compare to BEA summary #
# # ###################################
# # # TODO fix the joins for exp/imp, bc the bea summary join doesn't work with the new uti split is mulitple joining on uti, and min/coa on min
# # ## Make DataFrames to compare values
# # x_0[yr,:]  #	    "Exports of goods and services",
# ex_tmp = get_subtable(TMP,"exports"); ex_tmp.commodities=Symbol.(ex_tmp.commodities)
# # ex = DataFrame([x_0[yr,:].axes[1] x_0[yr,:].data], [:Win, :value]); ex.Win = Symbol.(ex.Win)
# # ex.subtable .= "exports"# add subtable column for aggregation with all subtables to match and join .table versions
# # ex.variable .= "F04000"
# # ex_bea = leftjoin(filter(x->x.year==year,get_subtable(all_summary_national_data,"exports"))[:,[1,4]], unique(Codes[:,[:WiNDC_summary,aggcol]]), on=:commodities=>:WiNDC_summary); 
# # ex_bea[!,aggcol] = Symbol.(ex_bea[!,aggcol])

# # comp_x_1 = outerjoin(ex_tmp,ex, on=:commodities=>:Win, renamecols = ""=>"_JLD2") 
# # comp_x = outerjoin(comp_x_1,ex_bea, on=:commodities=>aggcol, renamecols = ""=>"_callBEA")
# # print(comp_x)# "BEA Summary doesn't split for 11CA for (agr &) fof, or 525	to fin hou ore"

# # m_0[yr,:]  #	    "Imports",
# im_tmp = get_subtable(TMP,"imports"); im_tmp.commodities=Symbol.(im_tmp.commodities)
# # im = DataFrame([m_0[yr,:].axes[1] m_0[yr,:].data], [:Win, :value]); #im.Win=string.(im.Win)
# # im.subtable .= "imports" # add subtable column for aggregation with all subtables to match and join .table versions
# # im.variable .= "MCIF"
# # im_bea = leftjoin(get_subtable(all_summary_national_data, "imports")[:,[1,4]], unique(Codes[:,[:WiNDC_summary,aggcol]]), on=:commodities=>:WiNDC_summary); im_bea.WiNDC=Symbol.(im_bea.WiNDC)
# # comp_m_1 = outerjoin(im_tmp,im, on=:commodities=>:Win, renamecols = ""=>"_JLD2") 
# # comp_m = outerjoin(comp_m_1,im_bea, on=:commodities=>aggcol, renamecols = ""=>"_callBEA")
# # print(comp_m)

# # ms_0[yr,:,:]  #	"Margin supply",
# ms_tmp = get_subtable(TMP,"margin_supply"); ms_tmp.commodities=Symbol.(ms_tmp.commodities)
# if "TRADE " in unique(get_subtable(TMP,"margin_supply")[!,:sectors])
#     ms_tmp[:,:sectors] = [item == "TRADE " ? "trd" : "trn" for item in ms_tmp[:,:sectors]]; ms_tmp.sectors=Symbol.(ms_tmp.sectors)
#     elseif "Trade" in unique(get_subtable(TMP,"margin_supply")[!,:sectors])
#     ms_tmp[:,:sectors] = [item == "Trade" ? "trd" : "trn" for item in ms_tmp[:,:sectors]]; ms_tmp.sectors=Symbol.(ms_tmp.sectors)
#     else exit
#     end
# # ms = stack(DataFrame([ms_0.axes[2] ms_0[yr,:,:].data], [:Win ; ms_0.axes[3]]),2:3); ms.variable=Symbol.(ms.variable); ms.Win=Symbol.(ms.Win)
# # ms.subtable .= "margin_supply" # add subtable column for aggregation with all subtables to match and join .table versions
# # comp_ms = outerjoin(ms_tmp,ms, on=[:commodities=>:Win, :sectors=>:variable], renamecols=""=>"_JLD2")
# # print(comp_ms)

# # md_0[yr,:,:]  #	"Margin demand",
# md_tmp = get_subtable(TMP,"margin_demand"); md_tmp.commodities=Symbol.(md_tmp.commodities)
# if "TRADE " in unique(get_subtable(TMP,"margin_demand")[!,:sectors])
# md_tmp[:,:sectors] = [item == "TRADE " ? "trd" : "trn" for item in md_tmp[:,:sectors]]; md_tmp.sectors=Symbol.(md_tmp.sectors)
# elseif "Trade" in unique(get_subtable(TMP,"margin_demand")[!,:sectors])
# md_tmp[:,:sectors] = [item == "Trade" ? "trd" : "trn" for item in md_tmp[:,:sectors]]; md_tmp.sectors=Symbol.(md_tmp.sectors)
# else exit
# end
# # md = stack(DataFrame([md_0.axes[3] transpose(md_0[yr,:,:].data)], [:Win ; md_0.axes[2]]),2:3); md.variable=Symbol.(md.variable); md.Win=Symbol.(md.Win)
# # md.subtable .= "margin_demand" # add subtable column for aggregation with all subtables to match and join .table versions

# # comp_md = outerjoin(md_tmp,md, on=[:commodities=>:Win,:sectors=>:variable], renamecols=""=>"_JLD2")
# # print(comp_md)

# # fd_0[yr,:,:pce]  #	"Final demand", # pce
# pce_tmp = get_subtable(TMP,"personal_consumption"); pce_tmp.commodities=Symbol.(pce_tmp.commodities)
# # pce = DataFrame([[i for x in  fd_0[yr,:,:pce].axes for i in x] fd_0[yr,:,:pce].data], [:Win, :value])
# # pce.subtable .= "personal_consumption" # add subtable column for aggregation with all subtables to match and join .table versions
# # pce.variable .= "F01000"
# # comp_fdpce = outerjoin(pce_tmp,pce, on=:commodities=>:Win, renamecols = ""=>"_JLD2")
# # print(comp_fdpce)

# # # TODO fix join for this too maybe, not showing the new disaggs: it's a 3-way
# # # Exog fd
# # fd = stack(DataFrame([fd_0.axes[2] fd_0[yr,:,[x for x in FD if x!=:pce]].data], [:Win ; [x for x in FD if x !=:pce]]),2:18); fd.variable=Symbol.(fd.variable); fd.Win=Symbol.(fd.Win)
# # fd.subtable .= "exogenous_final_demand" # add subtable column for aggregation with all subtables to match and join .table versions
# fd_temp = get_subtable(TMP,"exogenous_final_demand"); #fd_temp.sectors = Symbol.(fd_temp.sectors); 
# # fdbea = get_subtable(all_summary_national_data, "exogenous_final_demand")
# # fd_bea = leftjoin(fdbea, Codes, on=[:sectors=>:BEA_summary])#; fd_bea.WiNDC=Symbol.(fd_bea.WiNDC)
# # # comp_x_1 = outerjoin(ex_tmp,ex, on=:commodities=>:Win) 

# fd_tmp = leftjoin(fd_temp, Codes, on =:sectors=>:BEA_detail)[:,[:commodities,:sectors,:value,aggcol]];  fd_tmp[!,aggcol] = Symbol.(fd_tmp[!,aggcol]); fd_tmp.commodities = Symbol.(fd_tmp.commodities)
# #### TODONE THSI WAS FOR SUMMARY TEST, CHANGE BACK TO DETAIL FOR REAL AGGREGATION
# # fd_tmp = leftjoin(fd_temp, Codes, on =:sectors=>:BEA_summary)[:,[:commodities,:sectors,:value,aggcol]];  fd_tmp[!,aggcol] = Symbol.(fd_tmp[!,aggcol]); fd_tmp.commodities = Symbol.(fd_tmp.commodities)
# ####
# # comp_fd = outerjoin(fd_tmp, fd, on=[aggcol=>:variable, :commodities=>:Win], renamecols = "" => "_JLD2")[:1:125,:]
# # print(comp_fd)

# # va_0[yr,:,:]  #	"Value added",
# va_tmp = append!(get_subtable(TMP,"labor_demand"),get_subtable(TMP,"capital_demand")); va_tmp.commodities=Symbol.(va_tmp.commodities); va_tmp.sectors=Symbol.(va_tmp.sectors)
# # va = stack(DataFrame([va_0.axes[3] transpose(va_0[yr,[va for va in VA],:].data)], [:Win; [va for va in VA]]),2:3); va.Win = Symbol.(va.Win); va.variable = Symbol.(va.variable)
# # va.subtable .=  [var == :compen ? "labor_demand" : "capital_demand" for var in va.variable]# add subtable column for aggregation with all subtables to match and join .table versions
# # comp_va = outerjoin(va_tmp,va, on=[:commodities=>:variable, :sectors=>:Win], renamecols=""=>"_JLD2")
# # print(comp_va)

# # ys_0[yr,:,:]   #	"Intermediate Supply",
# ys_tmp = get_subtable(TMP,"intermediate_supply"); #ys_tmp.commodities=Symbol.(ys_tmp.commodities);# ys_tmp.sectors=Symbol.(ys_tmp.sectors)
# # ys = stack(DataFrame([ys_0.axes[3] ys_0[yr,:,:].data], [:variable_ ; [x for x in I]]),2:72); #ys.Win=string.(ys.Win)#ys.variable=Symbol.(ys.variable); 
# # ys.subtable .= "intermediate_supply" # add subtable column for aggregation with all subtables to match and join .table versions

# # comp_ys = outerjoin(ys_tmp, ys, on=[:commodities=>:variable, :sectors=>:Win], renamecols=""=>"_JLD2")
# # print(comp_ys[1:900,:])


# # id_0[yr,:,:]  #	"Intermediate demand",
# id_tmp = get_subtable(TMP,"intermediate_demand"); id_tmp.sectors = Symbol.(id_tmp.sectors); id_tmp.commodities = Symbol.(id_tmp.commodities)
# # id = stack(DataFrame([id_0.axes[3] id_0[yr,:,:].data], [:Win ; [x for x in I]]),2:72); id.variable=Symbol.(id.variable); id.Win=Symbol.(id.Win)
# # id.subtable .= "intermediate_demand" # add subtable column for aggregation with all subtables to match and join .table versions
# # comp_id = outerjoin(id_tmp, id, on=[:commodities=>:Win, :sectors=>:variable], renamecols=""=>"_JLD2")
# # print(comp_id[1:1000,:])

# # rename!(ys,[:variable_, :variable] .=> [:variable,:Win])

# # JLD2WinNat = vcat(ex,im,pce,ms,md,fd,va,ys,id); JLD2WinNat.Win = string.(JLD2WinNat.Win); JLD2WinNat.variable = string.(JLD2WinNat.variable)
# # compare = outerjoin(TMP.table, JLD2WinNat, on = [:subtable, :commodities=>:Win, :sectors=>:variable], renamecols=""=>"_JLD2")
# # # filter(x -> x.commodities=="sle",outerjoin(filter(x-> x.subtable=="intermediate_supply",TMP.table), filter(x -> x.subtable=="intermediate_supply", JLD2WinNat), on = [:subtable, :commodities=>:variable, :sectors=>:Win], renamecols=""=>"_JLD2"))
# # filter(x->x.commodities=="ins",filter(x-> x.subtable=="intermediate_demand", compare))

# # compare.diff .= compare.value - compare.value_JLD2
# # print(sort!(compare, :diff, by=abs)[3000:4250,:])
# # # x->(!ismissing(abs(x)), x) , rev=true

# bop_tmp = WiNDC.balance_of_payments(TMP)
# # bop = DataFrame([bopdef_0.axes[1][bopdef_0.axes[1].==[yr]] bopdef_0[yr]],[:year,:value_0]); bop.year = parse.(Int,string.(bop.year))  #	"Balance of payments deficit",

# # comp_bop = leftjoin(bop_tmp,bop, on=:year)
# # print(comp_bop)

# # a_0[yr,:]  #	    "Armington supply",
# a_tmp = WiNDC.armington_supply(TMP); a_tmp.commodities=Symbol.(a_tmp.commodities)
# # a = DataFrame([a_0[yr,:].axes[1] a_0[yr,:].data], [:Win, :a]); a.Win = Symbol.(a.Win)
# # comp_a = outerjoin(a_tmp,a, on=:commodities=>:Win)
# # print(comp_a)
# # y_0[yr,:]   #	"Gross output",
# y_tmp = WiNDC.gross_output(TMP); y_tmp.commodities=Symbol.(y_tmp.commodities)
# # y = DataFrame([y_0[yr,:].axes[1] y_0[yr,:].data], [:Win, :y]); y.Win = Symbol.(y.Win)
# # comp_y = outerjoin(y_tmp,y, on=:commodities=>:Win)
# # print(comp_y)

# # tm_0[yr,:]  #	"Import tariff"; Initial, for price 
# tm_tmp = WiNDC.import_tariff_rate(TMP); tm_tmp.commodities=Symbol.(tm_tmp.commodities)
# # tm = DataFrame([tm_0[yr,:].axes[1] tm_0[yr,:].data], [:Win, :tm]); tm.Win = Symbol.(tm.Win)
# # comp_tm = outerjoin(tm_tmp,tm, on=:commodities=>:Win)
# # print(comp_tm)

# # ta_0[yr,:]  #	"Tax net subsidy rate on intermediate demand", benchmark data also for price level
# ta_tmp = WiNDC.absorption_tax_rate(TMP); ta_tmp.commodities=Symbol.(ta_tmp.commodities)
# # ta = DataFrame([ta_0[yr,:].axes[1] ta_0[yr,:].data], [:Win, :ta]); ta.Win = Symbol.(ta.Win)
# # comp_ta = outerjoin(ta_tmp,ta, on=:commodities=>:Win)
# # print(comp_ta)

# # ty_0[yr,:]  #	"Output tax rate"
# ty_tmp = WiNDC.other_tax_rate(TMP); ty_tmp.sectors=Symbol.(ty_tmp.sectors)
# # ty = DataFrame([ty_0[yr,:].axes[1] ty_0[yr,:].data], [:Win, :ty]); ty.Win = Symbol.(ty.Win)
# # comp_ty = outerjoin(ty_tmp,ty, on=:sectors=>:Win)
# # print(comp_ty)


# # #############################################################################
# # # Generate DenseAxisArrays for WiNDC national/MultiNat model, one year only #
# # # Add 0s and missing indices (messy, but kinda works)
# # #############################################################################
# I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
# Ip = [[x for x in I if x∉[:uti]]; [:uel,:ugs, :uwt, :coa]]#, :gas]]
# # J = deepcopy(I) # Index for WiNDC BEA Sectors

# Iplus1 = deepcopy(Ip) #if I updated in MultiGHGNat           # [[x for x in I if x∉[:uti]]; :uel ; :ugs ; :uwt ; :coa ; :gas] # :use, :oth,

# othersectors_ex_tmp = [x for x in Iplus1 if x∉ex_tmp.commodities]
# x_m0 = DenseAxisArray(vcat(ex_tmp.value,zeros(length(othersectors_ex_tmp))), vcat(ex_tmp.commodities,othersectors_ex_tmp))
# othersectors_im_tmp = [x for x in Iplus1 if x∉im_tmp.commodities]
# m_m0 = DenseAxisArray(vcat(im_tmp.value,zeros(length(othersectors_im_tmp))), vcat(im_tmp.commodities,othersectors_im_tmp))
# ms_m0df = unstack(ms_tmp,:commodities, :sectors,:value, fill=0)
# othersectors_ms_tmp = [x for x in Iplus1 if x∉ms_tmp.commodities]
# ms_m0 = DenseAxisArray(vcat(hcat(ms_m0df[:,:trn], ms_m0df[:,:trd]),zeros(length(othersectors_ms_tmp),2)), vcat(ms_m0df.commodities,othersectors_ms_tmp), [:trn, :trd])
# md_m0df = unstack(md_tmp,:commodities, :sectors,:value, fill=0)
# othersectors_md_tmp = [x for x in Iplus1 if x∉md_tmp.commodities]
# md_m0 = DenseAxisArray(transpose(vcat(hcat(md_m0df[:,:trn], md_m0df[:,:trd]),zeros(length(othersectors_md_tmp),2))), [:trn, :trd], vcat(md_m0df.commodities,othersectors_md_tmp))
# fd_m0exog = unstack(fd_tmp,:commodities, aggcol,:value, fill=0)
# fd_m0df = coalesce.(outerjoin(fd_m0exog,rename(pce_tmp[:,[:commodities,:value]],:value=>:pce), on=:commodities),0)
# othersectors_fddf_tmp = [x for x in Iplus1 if x∉fd_m0df.commodities]
# fd_m0 = DenseAxisArray(vcat(Matrix(fd_m0df[:,2:end]),zeros(length(othersectors_fddf_tmp),size(fd_m0df,2)-1)), vcat(fd_m0df.commodities,othersectors_fddf_tmp), Symbol.(names(fd_m0df)[2:end]))
# bop_m0 = DenseAxisArray(bop_tmp.value, [Symbol(year)])
# othersectors_a_tmp = [x for x in Iplus1 if x∉a_tmp.commodities]
# a_m0 = DenseAxisArray(vcat(a_tmp.value,zeros(length(othersectors_a_tmp))), vcat(a_tmp.commodities,othersectors_a_tmp))
# othersectors_y_tmp = [x for x in Iplus1 if x∉y_tmp.commodities]# Nice this has everything, but still works
# y_m0 = DenseAxisArray(vcat(y_tmp.value,zeros(length(othersectors_y_tmp))), vcat(y_tmp.commodities,othersectors_y_tmp))
# # tiny values generated for these sectors casue mismatch with 0s in output, so change back to 0
# y_m0[:fbt]= 0; y_m0[:gmt]= 0; y_m0[:mvt]= 0

# othersectors_tm_tmp = [x for x in Iplus1 if x∉tm_tmp.commodities]
# tm_m0 = DenseAxisArray(vcat(tm_tmp.value,zeros(length(othersectors_tm_tmp))), vcat(tm_tmp.commodities,othersectors_tm_tmp))
# othersectors_ta_tmp = [x for x in Iplus1 if x∉ta_tmp.commodities]
# ta_m0 = DenseAxisArray(vcat(ta_tmp.value,zeros(length(othersectors_ta_tmp))), vcat(ta_tmp.commodities,othersectors_ta_tmp))

# ys_tmp = get_subtable(TMP,"intermediate_supply"); #ys_tmp.commodities=Symbol.(ys_tmp.commodities);# ys_tmp.sectors=Symbol.(ys_tmp.sectors)


# othersectors_ys_tmp = [x for x in Iplus1 if x∉Symbol.(unique(ys_tmp.commodities))]# None, so no set up YET to add sectors
# for i in [x for x in Iplus1 if x∉ Symbol.(unique(ys_tmp.sectors))]
# push!(ys_tmp,[string(i) string(i) ys_tmp.year[1] 0])
# end
# ys_m0df = unstack(ys_tmp,:commodities, :sectors,:value, fill=0)

# # for i in (1:length(othersectors_ys_tmp2))
# #     ys_m0df[!,othersectors_ys_tmp2[i]] .= 0
# # end
# ys_m0 = DenseAxisArray(Matrix(ys_m0df[!,2:end]), Symbol.(ys_m0df.commodities), Symbol.(names(ys_m0df)[2:end]))
# id_m0df = unstack(id_tmp,:commodities, :sectors,:value, fill=0)
# othersectors_id_tmp = [x for x in Iplus1 if x∉id_tmp.commodities]
# othersectors_id_tmp2 = [x for x in Iplus1 if x∉ unique(id_tmp.sectors)]
# for i in (1:length(othersectors_id_tmp2))
#     id_m0df[!,othersectors_id_tmp2[i]] .= 0
# end
# id_m0 = DenseAxisArray(vcat(Matrix(id_m0df[:,2:end]),zeros(length(othersectors_id_tmp),size(id_m0df,2)-1))
# , vcat(id_m0df.commodities,othersectors_id_tmp), Symbol.(names(id_m0df)[2:end]))
# va_m0df = unstack(va_tmp,:commodities, :sectors,:value, fill=0)
# othersectors_va_tmp = [x for x in Iplus1 if x∉va_tmp.sectors]
# for i in (1:length(othersectors_va_tmp))
#     va_m0df[!,othersectors_va_tmp[i]] .= 0
# end
# va_m0 = DenseAxisArray(Matrix(va_m0df[!,2:end]), va_m0df[:,:commodities], Symbol.(names(va_m0df)[2:end]))
# othersectors_ty_tmp = [x for x in Iplus1 if x∉ty_tmp.sectors]
# ty_m0 = DenseAxisArray(vcat(ty_tmp.value,zeros(length(othersectors_ty_tmp))), vcat(ty_tmp.sectors,othersectors_ty_tmp))


# # CSV.write("id_plus.csv",Tables.table(transpose(id_m0.data)), header=axes(id_m0)[1])
# # id_m0[row,col] id_m0[secs,coms]

# macro varname(arg)
#     string(arg)
# end
# using XLSX
# filename="CompareSummary.xlsx"
# # filename="Combos.xlsx"

# # XLSX.openxlsx(filename, mode="w") do creation_of_the_excelfile
# # end

# # function exp_to_xlsx(df, dfname, filename="Combos.xlsx")
# function exp_to_xlsx(df, dfname, filename="CompareSummary.xlsx")
#     XLSX.openxlsx(filename,mode="rw") do xf
#     XLSX.addsheet!(xf,dfname)
#     XLSX.writetable!(xf[dfname], df)
# end
# end



# MultiNatdata = Dict(
#   :a_0      => a_m0,
#   :id_0     => id_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
#   :ys_0     => ys_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
#   :tm_0     => tm_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
#   :va_0     => va_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
#   :md_0     => md_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
#   :fd_0     => fd_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
#   :m_0      => m_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
#   :ty_0     => ty_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
#   :ms_0     => ms_m0, # 3-dimensional DenseAxisArray{Float64,3,...} 
#   :bopdef_0 => bop_m0, # 1-dimensional DenseAxisArray{Float64,1,...} 
#   :x_0      => x_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
#   :ta_0     => ta_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
#   :y_0      => y_m0, # 2-dimensional DenseAxisArray{Float64,2,...} 
# )

# # JLD2.save("WiNDCplus.jld2",MultiNatdata)
# # CSV.write("WiNDCplusDF.csv",TMP.table)

# # comp_x.commodities=string.(comp_x.commodities)
# # exp_to_xlsx(comp_x, @varname(comp_x))
# # comp_m.commodities=string.(comp_m.commodities)
# # exp_to_xlsx(comp_m, @varname(comp_m))
# # exp_to_xlsx(comp_bop, @varname(comp_bop))
# # comp_a.commodities=string.(comp_a.commodities)
# # exp_to_xlsx(comp_a, @varname(comp_a))
# # comp_y.commodities=string.(comp_y.commodities)
# # exp_to_xlsx(comp_y, @varname(comp_y))
# # comp_ms.commodities=string.(comp_ms.commodities); comp_ms.sectors=string.(comp_ms.sectors)
# # exp_to_xlsx(comp_ms, @varname(comp_ms))
# # comp_md.commodities=string.(comp_md.commodities); comp_md.sectors=string.(comp_md.sectors)
# # exp_to_xlsx(comp_md, @varname(comp_md))
# # comp_fdpce.commodities=string.(comp_fdpce.commodities)
# # exp_to_xlsx(comp_fdpce, @varname(comp_fdpce))
# # comp_fd.commodities=string.(comp_fd.commodities); comp_fd[!,aggcol]=string.(comp_fd[!,aggcol])
# # exp_to_xlsx(comp_fd, @varname(comp_fd))
# # comp_tm.commodities=string.(comp_tm.commodities)
# # exp_to_xlsx(comp_tm, @varname(comp_tm))
# # comp_ta.commodities=string.(comp_ta.commodities)
# # exp_to_xlsx(comp_ta, @varname(comp_ta))
# # comp_ys.commodities=string.(comp_ys.commodities); comp_ys.sectors=string.(comp_ys.sectors)
# # exp_to_xlsx(comp_ys, @varname(comp_ys))
# # comp_id.commodities=string.(comp_id.commodities); comp_id.sectors=string.(comp_id.sectors)
# # exp_to_xlsx(comp_id, @varname(comp_id))
# # comp_va.commodities=string.(comp_va.commodities); comp_va.sectors=string.(comp_va.sectors)
# # exp_to_xlsx(comp_va, @varname(comp_va))
# # comp_ty.sectors = string.(comp_ty.sectors) 
# # exp_to_xlsx(comp_ty, @varname(comp_ty))

# # Doesn't work because of BEAcodes -> WiNDC_plus codes
# outerjoin(
#     get_table(TMP),
#     get_table(callibrated_summary_national_data),
#     on = [:commodities, :sectors, :year, :subtable],
#     renamecols = "_agg" => ""
# ) |>
# x -> coalesce.(x,0) |>
# x -> transform(x,
#     [:value_agg, :value] => ByRow((a,b) -> (a - b)) => :diff
# ) |>
# x -> sort(x, :diff)



# function table_to_na(
#     X::WiNDCtable,
#     subtable::String;
#     columns_to_ignore = [:year]
# )
#     domain_columns = String.([a for a∈domain(X) if a∉columns_to_ignore])
#     sets = get_set.(Ref(X), domain_columns) |>
#             x -> map(y -> y[!,:element], x)       
#     out = NamedArray(
#         zeros(length.(sets)...),
#         sets
#     )
#     for row in eachrow(get_subtable(X, subtable))
#         idx = Tuple(row[domain_columns])
#         out[idx...] = row[:value]
#     end
#     return out
# end

# TestTMPArray = table_to_na(callibrated_summary_national_data, "intermediate_demand")
# # Troubleshoot function


# domain_columns = String.([a for a∈domain(callibrated_summary_national_data) if a∉[:year]])
#      sets = 
#     get_set.(Ref(callibrated_summary_national_data), domain_columns) |>
#             x -> map(y -> y[!,:element], x)       
#     out = NamedArray(
#         zeros(length.(sets)...),
#         sets
#     )
#     for row in eachrow(get_subtable(callibrated_summary_national_data, "imports"))
#         idx = Tuple(row[domain_columns])
#         out[idx...] = row[:value]
#     end