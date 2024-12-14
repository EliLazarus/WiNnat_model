using WiNDC
using JuMP

using GamsStructure
using DataFrames

data_dir = "C:\\Users\\Eli\\.julia\\dev\\WiNnat_model\\data\\windc_2021_julia"
data_dir = "..\\..\\..\\OneDrive - Stanford\\Box_2024-11-15\\CGE\\windc_2021_julia"


GU = WiNDC.load_state_data(data_dir)
GU = WiNDC.load_national_data(data_dir)
GU = WiNDC.load_detailed_national_data(data_dir)
WiNDC.create_national_sets()
test = WiNDC.NationalTable
year = Symbol(2021)

# m = state_disaggregation_model_mcp_year(GU,year)
m = national_disaggregation_model_mcp_year(GU,year)
WiNDC.create_national_subtables(GU)
import WiNDC.create_national_subtables
#function from "model.jl" in 'national'.
national_mpsge(data::NationalTable)

# Fix an income level to normalize prices in the MCP model 
fix(m[:RA][:CA],GU[:c0_][[year],[:CA]],force=true)

set_attribute(m, "cumulative_iteration_limit", 10_000)

optimize!(m)