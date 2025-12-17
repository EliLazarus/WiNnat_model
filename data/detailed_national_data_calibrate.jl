#############################
### Generate detailed sector/commodity for each year from mapping to annual summary level
### This is a set of functions from a deprecated branch of WiNDC. Mapping to current versions of the module may eventually be possible, but complicated
### as many functions have changed significantly in development since this version (though the main process is ultimately paralllel)
#############################
using XLSX, Ipopt

"""
    WiNDCtable

Base abstract type for all structures. Subtypes should have field names

- `table` - A DataFrame with columns [`domain`](@ref), `:subtable`, and `:value`.
- `sets` - A DataFrame with columns `:set`, `:element`, and `:description`.

And implement the function `domain(data::T) where T<:WiNDCtable` which should 
    return a vector of symbols representing the domain of the table.
"""
abstract type WiNDCtable end;

abstract type AbstractNationalTable <: WiNDCtable end

abstract type AbstractRegionalTable <: WiNDCtable end


"""
    domain(data::T) where T<:WiNDCtable

Return the domain of the table as a vector of symbols. Must be implemented for
any subtype of a WiNDCtable.
"""
domain(data::WiNDCtable) = throw(ArgumentError("domain not implemented for WiNDCtable"))

"""
    get_table(data::T) where T<:WiNDCtable

Return the main table of the WiNDCtable object as a DataFrame

## Required Arguments

1. `data` - A WiNDCtable-like object.

## Output

Returns a DataFrame with columns `domain(data)`, `subtable`, and `value`.
"""
function get_table(data::T) where T<:WiNDCtable
    return data.table
end

"""
    get_set(data::T) where T<:WiNDCtable  

    get_set(data::T, set_name::String) where T<:WiNDCtable

    get_set(data::T, set_name::Vector{String}) where T<:WiNDCtable

Return the elements of the given sets. If no set is given, return all sets.

## Required Arguments
1. `data` - A WiNDCtable-like object. 
2. `set_name` - A string or vector of strings representing the set names to be extracted.

## Returns

Returns a DataFrame with three columns, `:element`, `:set` and `:description`
"""
function get_set(data::T, set_name::String) where T<:WiNDCtable
    data.sets |>
        x -> subset(x, 
            :set => ByRow(==(set_name))
        )    
end

function get_set(data::T, set_name::Vector{String}) where T<:WiNDCtable
    data.sets |>
        x -> subset(x, 
            :set => ByRow(x -> in(x, set_name))
        )    
end

function get_set(data::T) where T<:WiNDCtable
    return data.sets
end

"""
    get_subtable(data::T, subtable::String, column::Vector{Symbol}; negative::Bool = false, keep_all_columns = false) where T<:WiNDCtable

    get_subtable(data::T, subtable::String; column::Symbol = :value, output::Symbol = :value, negative = false) where T<:WiNDCtable

    get_subtable(data::T, subtable::Vector{String}) where T<:WiNDCtable

Return the subtables requested as a DataFrame

## Required Arguments
1. `data` - A WiNDCtable-like object.
2. `subtable` - A string or vector of strings representing the subtable names to be extracted.

## Optional Arguments
- `column` - A symbol representing the column to be extracted. Default is `:value`.
- `output` - A symbol representing the output column name. Default is `:value`.
- `negative` - A boolean representing whether the values should be negated. Default is `false`.

## Returns

Returns a DataFrame with the requested subtables and columns.

"""
function get_subtable(
    data::WiNDCtable,
    subtable::String,
    column::Vector{Symbol};
    negative::Bool = false,
    keep_all_columns = false
)  

    columns = domain(data)
    append!(columns, column)

    elements = get_set(data, subtable) |>
        x -> select(x, :element)

    @assert(size(elements, 1) > 0, "Error: No elements found in subtable $subtable")

    return get_table(data) |>
        x -> innerjoin(
                x,
                elements,
            on = [:subtable => :element]
        ) |>
        x -> ifelse(keep_all_columns, x, select(x, columns))

end

function get_subtable(
        data::WiNDCtable, 
        subtable::String;
        column::Symbol = :value,
        output::Symbol = :value,
        negative = false
    )

    return get_subtable(data, subtable, [column]) |>
        x -> rename(x, column => output) |>
        x -> transform(x, output => ByRow(y -> negative ? -y : identity(y)) => output)
end

function get_subtable(
        data::WiNDCtable,
        subtable::Vector{String}
)
    
    return reduce(
        (x,y) -> append!(x, get_subtable(data, y, [:value]; keep_all_columns = true)),
        subtable,
        init = DataFrame()
    )

end

"""
    create_national_sets(
        use::XLSX.Worksheet, 
        supply::XLSX.Worksheet,
        set_regions)

This function creates the sets for the detailed national data.

set regions for detailed table

    Dict(
        "commodities" => ("use", ["A7:B408"], false, :commodities),
        "labor_demand" => ("use", ["A410:B410"], false, :commodities),
        "other_tax" => ("use", ["A411:B411"], false, :commodities),
        "capital_demand" => ("use", ["A412:B412"], false, :commodities),
        "sectors" => ("use", ["C5:ON6"], true, :sectors),
        "personal_consumption" => ("use", ["OP5:OP6"], true, :sectors),
        "household_supply" => ("use", ["OP5:OP6"], true, :sectors),
        "exports" => ("use", ["OV5:OV6"], true, :sectors),
        "exogenous_final_demand" => ("use", ["OQ5:OU6","OW5:PH6"], true, :sectors),
        "imports" => ("supply", ["OP5:OP6"], true, :sectors),
        "margin_demand" => ("supply", ["OS5:OT6"], true, :sectors),
        "margin_supply" => ("supply", ["OS5:OT6"], true, :sectors),
        "duty" => ("supply", ["OV5:OV6"], true, :sectors),
        "tax" => ("supply", ["OW5:OW6"], true, :sectors),
        "subsidies" => ("supply", ["OX5:OX6"], true, :sectors)
    )
"""
function create_national_sets(
        use::XLSX.Worksheet, 
        supply::XLSX.Worksheet,
        set_regions;
        table_type = :detailed)

    aggregate_sets = DataFrame(
        [
        ("labor_demand", "Labor Demand", "value_added", :composite),
        ("capital_demand", "Capital Demand", "value_added", :composite),
        ("other_tax", "Other taxes on production", "value_added", :composite),
        ("exogenous_final_demand", "Exogenous portion of final demand", "final_demand", :composite),
        ("exports", "Exports of goods and services", "final_demand", :composite),
        ("personal_consumption", "Personal consumption expenditures", "final_demand", :composite),
        ("intermediate_demand", "", "intermediate_demand", :composite),
        ("labor_demand", "", "labor_demand", :composite),
        ("capital_demand", "", "capital_demand", :composite),
        ("other_tax", "", "other_tax", :composite),
        ("personal_consumption", "Personal consumption expenditures", "personal_consumption", :composite),
        ("exogenous_final_demand","", "exogenous_final_demand", :composite),
        ("exports", "", "exports", :composite),
        ("intermediate_supply", "", "intermediate_supply", :composite),
        ("imports", "", "imports", :composite),
        ("household_supply", "Personal consumption expenditures, values <0", "household_supply", :composite),
        ("margin_demand", "", "margin_demand", :composite),
        ("margin_supply", "", "margin_supply", :composite),
        ("duty", "", "duty", :composite),
        ("tax", "", "tax", :composite),
        ("subsidies", "", "subsidies", :composite),
        ],
    [:element, :description, :set, :domain]
    )

    S = [aggregate_sets]
    for (key, (sheet, range, flip, dom)) in set_regions
        if sheet == "use"
            table = use
        else
            table= supply
        end
        if flip
            region = hcat([table[x] for x in range]...)
            region = permutedims(region, (2,1))
            if table_type == :detailed
                region[:,1], region[:,2] = region[:,2], region[:,1]
            end
        else
            region = vcat([table[x] for x in range]...)
        end
        df = DataFrame(region, [:element, :description]) |>
            x -> transform!(x, 
                :element => ByRow(x -> key) => :set,
                :element => ByRow(x -> dom) => :domain
            )
        push!(S, df)
    end

    return vcat(
            S...
        ) |>
        x -> transform(x,
            [:element, :set] .=> ByRow(x -> string(x)) .=> [:element, :set]
        )
end

function get_column_from_set(
        data::WiNDCtable,
        set_names::Vector{Symbol},
)
    D = []
    for a in set_names
        S = get_set(data) |>
            x -> subset(x, 
                :set => ByRow(==(String(a))),
                :domain => ByRow(!=(:composite))
            ) 
        elements = S[:, :element]
        new_d = unique(S[:, :domain])

        #return new_d
        @assert length(new_d) <= 1 "More than one domain found for set $(a)"
        @assert length(new_d) !=0 "No domain found for set $(a)" 
        push!(D, (new_d[1], elements))
    end
    return D

end

"""
    make_subtable(sets, rows, columns, table, subtable)

A helper function for extracting subtables.
"""
function make_subtable(sets, rows, columns, table, subtable)
    return crossjoin(
        sets |>
            x -> subset(x,:set => ByRow(==(rows))) |>
            x -> select(x, :element) |>
            x -> rename(x, :element => :commodities),
        sets |>
            x -> subset(x,:set => ByRow(==(columns))) |>
            x -> select(x, :element) |>
            x -> rename(x, :element => :sectors)
    ) |>
    x -> transform(x, :commodities => ByRow(x -> (table,subtable)) => [:table,:subtable]) 
end

"""
    create_national_subtables(sets)

This function creates the subtables for the detailed national data.
"""
function create_national_subtables(sets)
    return vcat(
        make_subtable(sets, "commodities", "sectors", "use", "intermediate_demand"),
        make_subtable(sets, "labor_demand", "sectors", "use", "labor_demand"),
        make_subtable(sets, "other_tax", "sectors", "use", "other_tax"),
        make_subtable(sets, "capital_demand", "sectors", "use", "capital_demand"),
        make_subtable(sets, "commodities", "personal_consumption", "use", "personal_consumption"),
        make_subtable(sets, "commodities", "household_supply", "use", "household_supply"),
        make_subtable(sets, "commodities", "exogenous_final_demand", "use", "exogenous_final_demand"),
        make_subtable(sets, "commodities", "exports", "use", "exports"),
        make_subtable(sets, "commodities", "sectors", "supply", "intermediate_supply"),
        make_subtable(sets, "commodities", "imports", "supply", "imports"),
        make_subtable(sets, "commodities", "cif", "supply", "cif"),
        make_subtable(sets, "commodities", "margin_demand", "supply", "margin_demand"),
        make_subtable(sets, "commodities", "margin_supply", "supply", "margin_supply"),
        make_subtable(sets, "commodities", "duty", "supply", "duty"),
        make_subtable(sets, "commodities", "tax", "supply", "tax"),
        make_subtable(sets, "commodities", "subsidies", "supply", "subsidies"),
    )
end

function load_national_year(
            X::XLSX.XLSXFile,
            year,
            range,
            table_name::String;
            scale = 1_000,
            replace_missing = false,
            data_start_row = 2)

    U = X[year][range]

    U[1,1] = :commodities
    U[1,2] = :drop

    if replace_missing
        U[U.=="..."] .= missing
    end

    return DataFrame(U[data_start_row:end,1:end], string.(U[1,:])) |>
                    x -> select(x, Not(:drop)) |>
                    x -> stack(x, Not("commodities"), variable_name = :sectors) |>
                    x -> coalesce.(x, 0) |>
                    x -> subset(x,
                        :value => ByRow(!=(0))
                    ) |>
                    x -> transform(x,
                        :value => (y -> y/scale) => :value,
                        [:commodities,:sectors] .=> ByRow(string) .=> [:commodities,:sectors],
                        :commodities => ByRow(y -> parse(Int, year)) => :year,
                        :commodities => ByRow(y -> table_name) => :table
                    )
end

function load_national_tables(
            use, 
            supply, 
            year::String;
            table_type = :detailed,
            use_range = "A6:PI417",
            supply_range = "A6:OZ409"
        )

    insurance_codes = table_type == :detailed ? ["524113","5241XX","524200"] : ["524"]    
    replace_missing = table_type != :detailed
    data_start_row = table_type == :detailed ? 2 : 3


    detailed_use = load_national_year(
        use,
        year,
        use_range,
        "use";
        replace_missing = replace_missing,
        data_start_row = data_start_row
    )

    trans_col = table_type == :detailed ? :TRANS : :Trans

    detailed_supply = load_national_year(
            supply,
            year,
            supply_range,
            "supply";
            replace_missing = replace_missing,
            data_start_row = data_start_row
        ) |>
        x -> unstack(x, :sectors, :value) |>
        x -> coalesce.(x, 0) |>
        x -> transform(x, 
            # adjust transport margins for transport sectors according to CIF/FOP 
            # adjustments. Insurance imports are specified as net of adjustments.
        [:commodities, trans_col, :MADJ] => ByRow((c,t,f) -> c∈insurance_codes ? t : t+f) => trans_col,
        [:commodities, :MCIF, :MADJ] => ByRow((c,i,f) -> c∈insurance_codes ? i+f : i) => :MCIF,
        ) |>
        x -> select(x, Not(:MADJ)) |>
            x -> stack(x, Not(:commodities, :year,:table), variable_name = :sectors, value_name = :value) |>
        x -> dropmissing(x) |>
        x -> subset(x, :value => ByRow(x -> x!=0)) 


    return vcat(
        detailed_use,
        detailed_supply
    ) |>
    x -> transform(x,
        [:commodities, :sectors] .=> ByRow(x -> string(x)) .=> [:commodities, :sectors]
    ) 
end

function load_detailed_national_tables(data_path::String)
    use = XLSX.readxlsx(joinpath(data_path, "use_detailed.xlsx"))
    supply = XLSX.readxlsx(joinpath(data_path, "supply_detailed.xlsx"))
    
    detailed_sets = Dict(
        "commodities" => ("use", ["A7:B408"], false, :commodities),
        "labor_demand" => ("use", ["A410:B410"], false, :commodities),
        "other_tax" => ("use", ["A411:B411"], false, :commodities),
        "capital_demand" => ("use", ["A412:B412"], false, :commodities),
        "sectors" => ("use", ["C5:ON6"], true, :sectors),
        "personal_consumption" => ("use", ["OP5:OP6"], true, :sectors),
        "household_supply" => ("use", ["OP5:OP6"], true, :sectors),
        "exports" => ("use", ["OV5:OV6"], true, :sectors),
        "exogenous_final_demand" => ("use", ["OQ5:OU6","OW5:PH6"], true, :sectors),
        "imports" => ("supply", ["OP5:OP6"], true, :sectors),
        "margin_demand" => ("supply", ["OS5:OT6"], true, :sectors),
        "margin_supply" => ("supply", ["OS5:OT6"], true, :sectors),
        "duty" => ("supply", ["OV5:OV6"], true, :sectors),
        "tax" => ("supply", ["OW5:OW6"], true, :sectors),
        "subsidies" => ("supply", ["OX5:OX6"], true, :sectors)
    )

    sets = create_national_sets(use["2017"], supply["2017"], detailed_sets)
    #detailed_subtables = WiNDC.create_national_subtables(sets)
    
    tables = []
    for year in [f for f in XLSX.sheetnames(use) if f!="NAICS Codes"]
        push!(
            tables, 
            load_national_tables(
                use, 
                supply, 
                year::String
            )
        )
    end

    return (vcat(tables...), sets)

end

"""
    function aggregate(
        data::T,
        aggregations...
    ) where T<:WiNDCtable


## Required Arguments

- `data::T`: The `WiNDCtable` to aggregate.
- `aggregations`: Takes the form `set_name => (X, original => new))`. 
    - `set_name` is a symbol, the name of the set to aggregate.
    - `X` a dataframe with the columns `original` and `new`.
    - `original` is the name of the column in with the elements to be aggregated.
    - `new` is the name of the column with the aggregated names.

"""
function aggregate(
    data::T,
    aggregations...
) where T<:WiNDCtable

    df = get_table(data)
    sets = get_set(data)

    for (set, (X, (original, new))) in aggregations
        (column, elements) = get_column_from_set(data, [set])[1]
        # (column, elements) = WiNDC.get_column_from_set(data, [set])[1]

        aggr = X |>
            y -> subset(y,
                original => ByRow(e -> in(e, elements))
            ) |>
            x -> select(x, [original, new])

        df = df |>
            x -> leftjoin(
                x,
                aggr,
                on = [column => original]
            ) |>
            x -> transform(x,
                [column, new] => ByRow((c, n) -> ismissing(n) ? c : n) => column
            ) |>
            x -> select(x, Not(new)) 

        sets = get_set(data, string(set)) |>
            x -> leftjoin(
                x,
                aggr,
                on = [:element => original]
            ) |>
            x -> transform(x,
                [:element, new] => ByRow((e, n) -> ismissing(n) ? e : n) => :element
            ) |>
            x -> select(x, Not(new))  |>
            x -> unique(x, :element) |>
            x -> vcat(
                sets |>
                    y -> subset(y, 
                        :set => ByRow(!=(string(set)))
                    ),
                    x
            )
        
    end
            
    df = df |>
        x -> groupby(x, [:subtable; domain(data)]) |>
        x -> combine(x, :value => sum => :value)



    return T(df, sets)

end

function load_summary_national_tables(data_path::String)
    summary_use = XLSX.readxlsx(joinpath(data_path, "use_summary.xlsx"))
    summary_supply = XLSX.readxlsx(joinpath(data_path, "supply_summary.xlsx"))
    

    summary_set_regions = Dict(
        "commodities" => ("use", ["A8:B80"], false, :commodities),
        "labor_demand" => ("use", ["A82:B82"], false, :commodities),
        "other_tax" => ("use", ["A83:B83"], false, :commodities),
        "capital_demand" => ("use", ["A85:B85"], false, :commodities),
        "sectors" => ("use", ["C6:BU7"], true, :sectors),
        "personal_consumption" => ("use", ["BW6:BW7"], true, :sectors),
        "household_supply" => ("use", ["BW6:BW7"], true, :sectors),
        "exports" => ("use", ["CC6:CC7"], true, :sectors),
        "exogenous_final_demand" => ("use", ["BX6:CB7","CD6:CO7"], true, :sectors),
        "imports" => ("supply", ["BW6:BW7"], true, :sectors),
        "margin_demand" => ("supply", ["BZ6:CA7"], true, :sectors),
        "margin_supply" => ("supply", ["BZ6:CA7"], true, :sectors),
        "duty" => ("supply", ["CC6:CC7"], true, :sectors),
        "tax" => ("supply", ["CD6:CD7"], true, :sectors),
        "subsidies" => ("supply", ["CE6:CE7"], true, :sectors)
    )

    summary_sets = create_national_sets(summary_use["2017"], summary_supply["2017"], summary_set_regions; table_type=:summary)

    summary_tables = []
    for year in XLSX.sheetnames(summary_use)
        push!(
            summary_tables, 
            load_national_tables(
                summary_use, 
                summary_supply, 
                year;
                table_type = :summary,
                use_range = "A6:CP90",
                supply_range = "A6:CG81"
            )
        )
    end

    return (vcat(summary_tables...), summary_sets)
end


function apply_national_subtables(data::DataFrame, subtables::DataFrame)

    return data |>
        x -> innerjoin(
            x,
            subtables,
            on = [:commodities, :sectors, :table]
        ) |>
        x -> select(x, :commodities, :sectors, :year, :subtable, :value) |>
        x -> unstack(x, :subtable, :value) |>
        x -> coalesce.(x,0) |>
        x -> transform(x, 
            [:intermediate_demand, :intermediate_supply] => ByRow(
                (d,s) -> (max(0, d - min(0, s)), max(0, s - min(0, d)))) => [:intermediate_demand, :intermediate_supply], #negative flows are reversed
            :subsidies => ByRow(y -> -y) => :subsidies,
            :margin_demand => ByRow(y ->  max(0,y)) => :margin_demand,
            :margin_supply => ByRow(y -> -min(0,y)) => :margin_supply,
            :personal_consumption => ByRow(y -> max(0,y)) => :personal_consumption,
            :household_supply => ByRow(y -> -min(0,y)) => :household_supply,
        ) |>
        x -> stack(x, Not(:commodities, :sectors, :year), variable_name = :subtable, value_name = :value) |>
        x -> subset(x, :value => ByRow(!=(0)))

end


function national_tables(data_path::String; aggregation = :detailed)
    
    @assert(aggregation∈[:summary,:detailed,:raw_detailed], "Error: aggregation must be either :summary or :detailed")

    if aggregation == :summary
        X, sets = load_summary_national_tables(data_path)
        subtables = create_national_subtables(sets)

        X = apply_national_subtables(X, subtables)

        return NationalTable(X, sets)
    end

    if aggregation == :raw_detailed
        X,sets = load_detailed_national_tables(data_path)
        subtables = create_national_subtables(sets)
        X = apply_national_subtables(X, subtables)
        return NationalTable(X, sets)
    end

    if aggregation == :detailed
        
        summary,_ = load_summary_national_tables(data_path)
        detailed,sets = load_detailed_national_tables(data_path)
        subtables = create_national_subtables(sets)
        
        summary_map = detailed_summary_map(data_path)

        X = national_disaggragate_summary_to_detailed(detailed, summary, summary_map)
        X = apply_national_subtables(X, subtables)
        return NationalTable(X,sets)
    end

end

############################
#### Disagregation Code ####
############################
"""
    down_fill(X)

This function fills in the missing values in a column with the last non-missing value.
"""
function down_fill(X)
    output = Vector{Any}()
    current_value = X[1]
    for row in X
        if !ismissing(row)
            current_value = row
        end
        push!(output, current_value)
    end
    return output
end

"""
    detailed_summary_map(detailed_path)

This function reads the detailed table and returns a DataFrame that maps the detailed
sectors to the summary sectors. The first sheet of the detailed table is a map between
the detailed sectors and the summary sectors. In addition this maps value added, final
demand and supply extras to the summary sectors.
"""
function detailed_summary_map(detailed_path::String)

    detailed_xlsx = XLSX.readdata(joinpath(detailed_path,"use_detailed.xlsx"), "NAICS Codes", "B5:E1022")

    df = detailed_xlsx |>
        x -> DataFrame(x[4:end,:], [x[1,1:end-1]..., "description"]) |>
        x -> select(x, Not("U.Summary")) |>
        x -> transform(x,
            :Summary => down_fill => :Summary,
            :Detail => down_fill => :Detail
        ) |>
        x -> dropmissing(x) |>
        x -> rename(x, :Summary => :summary, :Detail => :detailed) |>
        x -> transform(x,
            [:summary,:detailed] .=> ByRow(string) .=> [:summary,:detailed]
        )

        df = vcat(df, 
            DataFrame([
                (summary = "T00OTOP", detailed = "T00OTOP", description = "Other taxes on production"),
                (summary = "V003", detailed = "V00300", description = "Gross operating surplus"),
                (summary = "V001", detailed = "V00100", description = "Compensation of employees"),
                (summary = "F010", detailed = "F01000", description = "Personal consumption expenditures"),
                (summary = "F02E", detailed = "F02E00", description = "Nonresidential private fixed investment in equipment"),
                (summary = "F02N", detailed = "F02N00", description = "Nonresidential private fixed investment in intellectual property products"),
                (summary = "F02R", detailed = "F02R00", description = "Residential private fixed investment"),
                (summary = "F02S", detailed = "F02S00", description = "Nonresidential private fixed investment in structures"),
                (summary = "F030", detailed = "F03000", description = "Change in private inventories"),
                (summary = "F040", detailed = "F04000", description = "Exports of goods and services"),
                (summary = "F06C", detailed = "F06C00", description = "National defense: Consumption expenditures"),
                (summary = "F06E", detailed = "F06E00", description = "Federal national defense: Gross investment in equipment"),
                (summary = "F06N", detailed = "F06N00", description = "Federal national defense: Gross investment in intellectual property products"),
                (summary = "F06S", detailed = "F06S00", description = "Federal national defense: Gross investment in structures"),
                (summary = "F07C", detailed = "F07C00", description = "Nondefense: Consumption expenditures"),
                (summary = "F07E", detailed = "F07E00", description = "Federal nondefense: Gross investment in equipment"),
                (summary = "F07N", detailed = "F07N00", description = "Federal nondefense: Gross investment in intellectual property products"),
                (summary = "F07S", detailed = "F07S00", description = "Federal nondefense: Gross investment in structures"),
                (summary = "F10C", detailed = "F10C00", description = "State and local government consumption expenditures"),
                (summary = "F10E", detailed = "F10E00", description = "State and local: Gross investment in equipment"),
                (summary = "F10N", detailed = "F10N00", description = "State and local: Gross investment in intellectual property products"),
                (summary = "F10S", detailed = "F10S00", description = "State and local: Gross investment in structures"),
                (summary = "MCIF", detailed = "MCIF", description = "Imports"),
                (summary = "MADJ", detailed = "MADJ", description = "CIF/FOB Adjustments on Imports"),
                (summary = "Trade", detailed = "TRADE ", description = "Trade margins"),
                (summary = "Trans", detailed = "TRANS", description = "Transport margins"),
                (summary = "MDTY", detailed = "MDTY", description = "Import duties"),
                (summary = "TOP", detailed = "TOP", description = "Tax on products"),
                (summary = "SUB", detailed = "SUB", description = "Subsidies on products"),
                ])
        )
    return df
end


"""
    weight_function(year_detail, year_summary, minimum_detail, maximum_detail)

Create the weight function for the interpolation of the detailed table to the summary table
based solely on the year.
"""
function weight_function(year_detail::Int, year_summary::Int, minimum_detail::Int, maximum_detail::Int)
    if year_detail == minimum_detail && year_summary < year_detail
        return 1
    elseif year_detail == maximum_detail && year_summary > year_detail
        return 1
    elseif abs(year_detail.-year_summary) < 5
        return 1 .-abs.(year_detail.-year_summary)/5
    else
        return 0
    
    end
end


function national_disaggragate_summary_to_detailed(detailed::DataFrame, summary::DataFrame, summary_map::DataFrame)

    min_detail_year = minimum(detailed[!, :year])
    max_detail_year = maximum(detailed[!, :year])

#1111A0 111CA
#1111B0 111CA

    detailed_value_share = detailed |>
        x -> innerjoin(x, summary_map, on = :commodities => :detailed, renamecols = "" => "_commodities") |>
        x -> innerjoin(x, summary_map, on = :sectors => :detailed, renamecols = "" => "_sectors") |>
        x -> groupby(x, [:summary_sectors,:summary_commodities,:year, :table]) |>
        x -> combine(x,
            :value => (y -> y./sum(y)) => :value_share,
            [:commodities, :sectors] .=> identity .=> [:commodities, :sectors]
        ) |>
        x -> select(x, :commodities, :summary_commodities, :sectors, :summary_sectors, :year, :table, :value_share)

    df = innerjoin(
            detailed_value_share,
            summary,
            on = [:summary_commodities => :commodities, :summary_sectors => :sectors, :table],
            renamecols = "" => "_summary"
        ) |>
        x -> transform(x,
            [:year, :year_summary, :value_share, :value_summary] => 
                ByRow((dy,y,share,summary) -> weight_function(dy,y,min_detail_year,max_detail_year)*share*summary) => 
                :value
        ) |>
        x -> groupby(x, [:commodities, :sectors, :year_summary, :table]) |>
        x -> combine(x, :value => sum => :value) |>
        x -> rename(x, :year_summary => :year) |>
        x -> subset(x, :value => ByRow(!=(0)))    


    return df

end

#############################
### Begin Calibrate section
#############################
function WiNDC_zero_profit(data::AbstractRegionalTable; column = :value, output = :zero_profit)

    ag_columns = filter(y -> y!=:commodities, domain(data))
    return vcat(
            get_subtable(data, "intermediate_demand", column = column, output = output),
            get_subtable(data, "value_added", column = column, output = output), 
            get_subtable(data, "intermediate_supply", column = column, output = output, negative = true) 
        ) |> 
        x -> groupby(x, ag_columns) |>
        x -> combine(x, output => (y -> sum(y;init=0)) => output)    
end

function WiNDC_zero_profit(data::AbstractNationalTable; column = :value, output = :zero_profit)

    ag_columns = filter(y -> y!=:commodities, domain(data))
    return vcat(
            get_subtable(data, "intermediate_demand", column = column, output = output),
            get_subtable(data, "value_added", column = column, output = output), 
            get_subtable(data, "intermediate_supply", column = column, output = output, negative = true) 
        ) |> 
        x -> groupby(x, ag_columns) |>
        x -> combine(x, output => (y -> sum(y;init=0)) => output)    
end


function WiNDC_market_clearance(data::AbstractRegionalTable; column = :value, output = :market_clearance) 
    ag_columns = filter(y -> y!=:sectors, domain(data))

        return vcat(
            get_subtable(data, "intermediate_demand", column = column, output = output) ,
            get_subtable(data, "final_demand", column = column, output = output),
            get_subtable(data, "household_supply", column = column, output = output, negative = true),

            get_subtable(data, "intermediate_supply", column = column, output = output, negative = true),
            get_subtable(data, "imports", column = column, output = output, negative = true),
            get_subtable(data, "margin_demand", column = column, output = output, negative = true),
            get_subtable(data, "margin_supply", column = column, output = output), # Made negative earlier
            get_subtable(data, "duty", column = column, output = output, negative = true),
            get_subtable(data, "tax", column = column, output = output, negative = true),
            get_subtable(data, "subsidies", column = column, output = output) # Made negative earlier
        ) |>
        x -> groupby(x, ag_columns) |>
        x -> combine(x, output => (y -> sum(y;init=0)) => output)
end

function WiNDC_market_clearance(data::AbstractNationalTable; column = :value, output = :market_clearance) 
    ag_columns = filter(y -> y!=:sectors, domain(data))

        return vcat(
            get_subtable(data, "intermediate_demand", column = column, output = output) ,
            get_subtable(data, "final_demand", column = column, output = output),
            get_subtable(data, "household_supply", column = column, output = output, negative = true),

            get_subtable(data, "intermediate_supply", column = column, output = output, negative = true),
            get_subtable(data, "imports", column = column, output = output, negative = true),
            get_subtable(data, "margin_demand", column = column, output = output, negative = true),
            get_subtable(data, "margin_supply", column = column, output = output), # Made negative earlier
            get_subtable(data, "duty", column = column, output = output, negative = true),
            get_subtable(data, "tax", column = column, output = output, negative = true),
            get_subtable(data, "subsidies", column = column, output = output) # Made negative earlier
        ) |>
        x -> groupby(x, ag_columns) |>
        x -> combine(x, output => (y -> sum(y;init=0)) => output)
end

WiNDC_margin_balance(data::AbstractRegionalTable; column = :value, output = :margin_balance) =
    vcat(
        # WiNDC.get_subtable(data, "margin_supply", column = column, output = output, negative = true),
        # WiNDC.get_subtable(data, "margin_demand", column = column, output = output)
        get_subtable(data, "margin_supply", column = column, output = output, negative = true),
        get_subtable(data, "margin_demand", column = column, output = output)
    ) |>
    x -> groupby(x, filter(y -> y!=:commodities, domain(data))) |>
    x -> combine(x, output => (y -> sum(y;init=0)) => output)

WiNDC_margin_balance(data::AbstractNationalTable; column = :value, output = :margin_balance) =
    vcat(
        # WiNDC.get_subtable(data, "margin_supply", column = column, output = output, negative = true),
        # WiNDC.get_subtable(data, "margin_demand", column = column, output = output)
        get_subtable(data, "margin_supply", column = column, output = output, negative = true),
        get_subtable(data, "margin_demand", column = column, output = output)
    ) |>
    x -> groupby(x, filter(y -> y!=:commodities, domain(data))) |>
    x -> combine(x, output => (y -> sum(y;init=0)) => output)

using JuMP
"""
    WiNDCtable

Base abstract type for all structures. Subtypes should have field names

- `table` - A DataFrame with columns [`domain`](@ref), `:subtable`, and `:value`.
- `sets` - A DataFrame with columns `:set`, `:element`, and `:description`.

And implement the function `domain(data::T) where T<:WiNDCtable` which should 
    return a vector of symbols representing the domain of the table.
"""
abstract type WiNDCtable end;

"""
    domain(data::T) where T<:WiNDCtable

Return the domain of the table as a vector of symbols. Must be implemented for
any subtype of a WiNDCtable.
"""
domain(data::WiNDCtable) = throw(ArgumentError("domain not implemented for WiNDCtable"))

abstract type AbstractNationalTable <: WiNDCtable end

domain(data::AbstractNationalTable) = [:commodities, :sectors, :year]


struct NationalTable <: AbstractNationalTable
    table::DataFrame
    sets::DataFrame
end

"""
    calibrate(data::AbstractNationalTable)

This is currently geared toward calibrating the national dataset.
I'll be working to make this be a general calibration function.

Returns a new AbstractNationalTable with the calibrated values and the model.

There are three primary balancing operations:

1. Zero Profit - Column sums are equal
2. Market Clearance - Row sums are equal
3. Margin Balance - The margins balance

The three tax rates are fixed. The tax rates are:

1. Output Tax Rate
2. Absorption Tax Rate
3. Import Tariff Rate

The following are fixed:

1. Labor Compensation
2. Imports
3. Exports
4. Household Supply

Any zero values will remain zero. 

"""
function calibrate(data::T; silent = false) where T<:AbstractNationalTable

    
    M = Model(Ipopt.Optimizer)

    if silent
        set_silent(M)
    end

    @variable(M, 
        x[1:size(get_table(data),1)]
    )


    # Attach variables to dataframe
    get_table(data) |>
    x -> transform!(x,
        :value => (y -> M[:x]) => :variable
    )

    lob = .01
    upb = 100

    # set bounds and start values
    for row in eachrow(get_table(data))
        set_start_value(row[:variable], row[:value])
        lower_bound = row[:value]>0 ? row[:value]*lob : row[:value]*upb
        upper_bound = row[:value]>0 ? row[:value]*upb : row[:value]*lob
        set_lower_bound(row[:variable], lower_bound)
        set_upper_bound(row[:variable], upper_bound)
        if row[:value] == 0
            fix(row[:variable], 0; force=true)
        end
    end


    # Fix certain parameters -- exogenous portions of final demand,
    # value added, imports, exports and household supply
    vcat(
        get_subtable(data, "imports", [:value, :variable]),
        get_subtable(data, "exports", [:value, :variable]),
        get_subtable(data, "labor_demand", [:value, :variable]),
        get_subtable(data, "household_supply", [:value, :variable]),
    ) |>
    x -> transform(x,
        [:value, :variable] => ByRow((val, var) -> fix(var, val; force=true))
    ) 

    @objective(
        M, 
        Min, 
        get_table(data) |> 
            x -> transform(x,
                [:value, :variable] => ByRow((val, var) -> 
                    abs(val) * (var/val - 1)^2) => :objective
            ) |>
            x -> combine(x, :objective => sum => :objective) |>
            x -> x[1,:objective]
    )

    WiNDC_zero_profit(data; column = :variable) |> 
    x -> @constraint(M, 
        WiNDC_zero_profit[i=1:size(x,1)],
        x[i,:zero_profit] == 0
    )

    WiNDC_market_clearance(data; column = :variable) |>
    x -> @constraint(M,
        mkt[i=1:size(x,1)],
        x[i,:market_clearance] == 0
    )

    WiNDC_margin_balance(data; column = :variable) |>
    x -> @constraint(M,
        WiNDC_margin_balance[i=1:size(x,1)],
        x[i,:margin_balance] == 0
    )
    
    
    # Bound gross output
    outerjoin(
        gross_output(data; column = :variable, output = :expr),
        gross_output(data; column = :value),
        on = [:commodities, :year]
    ) |> 
    x -> transform(x,
            :value => ByRow(v -> v>0 ? floor(lob*v)-5 : upb*v) => :lower, 
            :value => ByRow(v -> v>0 ? upb*v : ceil(lob*v)+5) => :upper, 
    )|>
    x -> @constraint(M,
        gross_output[i=1:size(x,1)],
        x[i,:lower] <= x[i,:expr] <= x[i,:upper]
    )
        
    
    # Bound armington supply
    outerjoin(
        armington_supply(data; column = :variable, output = :expr),
        armington_supply(data; column = :value),
        on = [:commodities, :year]
    ) |>
    x -> @constraint(M,
        armington_supply[i=1:size(x,1)],
        max(0,lob * x[i,:value]) <= x[i,:expr] <= abs(upb * x[i,:value])
    )
    
    
    # Fix tax rates
    outerjoin(
        get_subtable(data, "intermediate_supply", column = :variable, output = :is) |>
            x -> groupby(x, filter(y -> y!=:commodities, domain(data))) |>
            x -> combine(x, :is => sum => :is),
        get_subtable(data, "other_tax", column = :variable, output = :ot) |>
            x -> select(x, Not(:commodities)),
        other_tax_rate(data, column = :value, output = :otr),
        on = filter(y -> y!=:commodities, domain(data))
    ) |>
    x -> dropmissing(x) |>
    x -> @constraint(M, 
        Output_Tax_Rate[i=1:size(x,1)],
        x[i,:ot] == x[i,:is] * x[i,:otr]
    )
    
    outerjoin(
        absorption_tax(data, column = :variable, output = :at),
        armington_supply(data, column = :variable, output = :as),
        absorption_tax_rate(data, output = :atr),
        on = filter(y -> y!=:sectors, domain(data))
    ) |>
    x -> dropmissing(x) |>
    x -> @constraint(M,
        Absorption_Tax_Rate[i=1:size(x,1)],
        x[i,:at] == x[i,:as] * x[i,:atr]
    )
    
    outerjoin(
        get_subtable(data, "duty", column = :variable, output = :it) |>
            x -> select(x, Not(:sectors)),
        get_subtable(data, "imports", column = :variable, output = :imports) |>
            x -> select(x, Not(:sectors)),
        import_tariff_rate(data, output = :itr),
        on = filter(y -> y!=:sectors, domain(data))
    ) |>
    x -> dropmissing(x) |>
    x -> @constraint(M,
        Import_Tariff_Rate[i=1:size(x,1)],
        x[i,:it] == x[i,:imports] * x[i,:itr]
    )
    
    optimize!(M)

    @assert is_solved_and_feasible(M) "Error: The model was not solved to optimality."

    df = get_table(data) |>
        x -> transform(x,
            :variable => ByRow(value) => :value
        ) |>
        x -> select(x, Not(:variable))


    get_table(data) |>
        x -> select!(x, Not(:variable))

    return (T(df, data.sets), M)

end

######################
## Aggregate Tables ##
######################

gross_output(data::AbstractNationalTable; column::Symbol = :value, output::Symbol = :value) =
    vcat(
        get_subtable(data, "intermediate_supply", column = column, output = output),
        get_subtable(data, "household_supply", column = column, output = output),
        get_subtable(data, "margin_supply", column = column, output = output, negative = true)
    ) |>
    x -> groupby(x, filter(y -> y!=:sectors, domain(data))) |>
    x -> combine(x, output => (y -> sum(y;init=0)) => output)


armington_supply(data::AbstractNationalTable; column = :value, output = :value) = 
    vcat(
        get_subtable(data, "intermediate_demand", column = column, output = output),
        get_subtable(data, "exogenous_final_demand", column = column, output = output),
        get_subtable(data, "personal_consumption", column = column, output = output),
    ) |>
    x -> groupby(x, filter(y -> y!=:sectors, domain(data))) |>
    x -> combine(x, output => (y -> sum(y;init=0)) => output)


other_tax_rate(data::AbstractNationalTable; column = :value, output = :value) = 
    outerjoin(
        get_subtable(data, "intermediate_supply", column = column, output = :is) |>
            x -> groupby(x, filter(y -> y!=:commodities, domain(data))) |>
            x -> combine(x, :is => sum => :is),
        get_subtable(data, "other_tax", column = column, output = :ot) |>
            x -> select(x, Not(:commodities)),
        on = filter(y -> y!=:commodities, domain(data))
    ) |>
    x -> coalesce.(x,0) |>
    x -> transform(x,
        [:is, :ot] => ByRow((i,o) -> i == 0 ? 0 : o/i) => output
    ) |>
    x -> select(x, Not(:is,:ot))


absorption_tax(data::AbstractNationalTable; column = :value, output = :value) = 
    outerjoin(
        get_subtable(data, "tax", column = column, output = :tax) |>
            x -> select(x, Not(:sectors)),
        get_subtable(data, "subsidies", column = column, output = :subsidies) |>
            x -> select(x, Not(:sectors)),
        on = filter(y -> y!=:sectors, domain(data))
    ) |>
    x -> coalesce.(x, 0) |>
    x -> transform(x,
        [:tax, :subsidies] => ByRow((t,s) -> t-s) => output
    ) |>
    x -> select(x, Not(:tax, :subsidies))

absorption_tax_rate(data::AbstractNationalTable; column = :value, output = :value) =     
    outerjoin(
        absorption_tax(data; column = column, output = :total_tax),
        armington_supply(data; column = column, output = :arm_sup),
        on = filter(y -> y!=:sectors, domain(data))
    )|>
    x -> coalesce.(x, 0) |>
    x -> transform(x,
        [:arm_sup, :total_tax] => ByRow((v,t) -> v == 0 ? 0 : t/v) => output
    ) |>
    x -> select(x, Not(:total_tax, :arm_sup)) |>
    x -> subset(x, output => ByRow(!=(0)))


import_tariff_rate(data::AbstractNationalTable; column = :value, output = :value) = 
    outerjoin(
        get_subtable(data, "duty", column = column, output = :duty) |>
            x -> select(x, Not(:sectors)),
        get_subtable(data, "imports", column = column, output = :imports) |>
            x -> select(x, Not(:sectors)),
        on = filter(y -> y!=:sectors, domain(data))
    ) |>
    x -> coalesce.(x, 0) |>
    x -> transform(x,
        [:duty, :imports] => ByRow((d,i) -> i==0 ? 0 : d/i) => output
    ) |>
    x -> select(x, Not(:duty, :imports)) |>
    x -> subset(x, output => ByRow(!=(0)))


balance_of_payments(data::AbstractNationalTable; column = :value, output = :value) = 
    outerjoin(
        get_subtable(data, "imports", output = :im) |>
            x -> select(x, Not(:sectors)),
        get_subtable(data, "exports", output = :ex) |>
            x -> select(x, Not(:sectors)),
        armington_supply(data, output = :as),
        on = filter(y -> y!=:sectors, domain(data))
    ) |>
    x -> coalesce.(x,0) |>
    x -> transform(x,
        #[:im, :ex, :as] => ByRow((im, ex, a) -> a!= 0 ? im - ex : 0) => output
        [:im, :ex, :as] => ByRow((im, ex, a) -> im - ex) => output
    ) |>
    x -> groupby(x, :year) |>
    x -> combine(x, output => sum => output)