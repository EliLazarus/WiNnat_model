using WiNDC, DataFrames, MPSGE, GamsStructure, CSV

##cd to top level of the data files
raw_data_directory = "./detailed_data/national"

# aggregation in [:summary, :detailed]. I'd like to have this be a wider range
# of options.
all_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :detailed);
all_summary_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :summary);
# For the summary level data, no issue using all years. Several years are failing
# calibration for the detailed data. 2022 works, as does 2017. Others do too, but
# those are the "important" ones for me at the moment.

raw_national_data = 
    NationalTable(
        get_table(all_national_data)|>
            x -> subset(x, :year => ByRow(==(2020))), 
        all_national_data.sets
    )

raw_national_summary_data = 
NationalTable(
    get_table(all_summary_national_data)|>
        x -> subset(x, :year => ByRow(==(2020))), 
    all_summary_national_data.sets
)

# Calibrate the model. The JuMP model is also returned. 
callibrated_national_data,M = calibrate(raw_national_data)
callibrated_national_data_summary,MS = calibrate(raw_national_summary_data)


# Apparently there is an error here using the summary data. I have the value added
# categories hard coded... which is super smart of me! I should fix that.
# M = national_mpsge(national_data)


# Benchmark
# solve!(M, cumulative_iteration_limit = 0)


# There are three parameters you can experiment with
#
# Absorption_tax[commodities]
# Output_tax[sectors]
# Margin_tax[commodities]
#
#These default to their values in the dataset. 

# set_value!.(M[:Output_tax], .5*value.(M[:Output_tax]))

# solve!(M)

# generate_report(M)

########################
## Exploring the data ##
########################

# See the entire table. Most of the columns should be self-explanatory. The really
# important one is `subtable`, which is the human readable component of the table.
get_table(callibrated_national_data)

# You can extract a single subtable. Note the subtable column is dropped (should 
# be fixed, I don't care for that).
get_subtable(callibrated_national_data, "intermediate_demand")
get_subtable(callibrated_national_data, "intermediate_demand")

# Or multiple
get_subtable(callibrated_national_data, ["intermediate_demand", "value_added"])

# Turns out, that's in a branch... quick hack, but it drops subtables... 
vcat(get_subtable.(Ref(callibrated_national_data), ["intermediate_demand", "value_added"])...)

# You should notice there is no subtable for value added.
get_table(callibrated_national_data) |>
    x -> unique(x, :subtable)

# This is because value added is a composite table, it consists of several smaller
# subtables. To see which ones use get_set, the elements comprise value added.
get_set(callibrated_national_data, "value_added")

# If we look at labor_demand we will see two elements, labor demand and V00100. 
# The V00100 is the NAICS code for labor compensation. 
get_set(callibrated_national_data, "labor_demand")

# example correspondence table
# X = DataFrame([["111A0","1111B0"], ["A","A"]], [:commodities, :name])
Codes=CSV.read(joinpath(raw_data_directory,"./BEA_WiNDC_Detail-Summary_codes.csv"), DataFrame, header=1)
# Codes=[Symbol.(Codes[:,1:4]) Codes[:,5]] NOT GOOD -> raw_national_data has Strings

X = DataFrame([Codes[:,:BEA_detail], Codes[:,:WiNDC]], [:commodities, :name])

## callibrate first=> get_table(callibrated_national_data)
test_subtables = ["intermediate_supply", "intermediate_demand","labor_demand","capital_demand","other_tax"]
TMP = NationalTable( vcat(
    get_table(callibrated_national_data) |>
    x -> filter(:subtable => !(in(test_subtables)), x) |>
    x -> leftjoin(
        x,
        X,
        on = :commodities,
        renamecols = "" => "_com"
    ) |>
    x -> transform(x,
        [:commodities, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :commodities,
    ) |>
    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x, :value => sum => :value),


    get_table(callibrated_national_data) |>
    x -> filter(:subtable => (in(test_subtables)), x) |>
    x -> leftjoin(
        x,
        X,
        on = :commodities,
        renamecols = "" => "_com"
    ) |>
    x -> transform(x,
        [:commodities, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :commodities,
    ) |>
    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x, :value => sum => :value) |>
    x -> leftjoin(
    x,
    X,
    on = :sectors=>:commodities,
    renamecols = "" => "_com") |>
    x -> transform!(x,
        [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
    ) |>
    x -> groupby(x, [:commodities, :sectors, :year, :subtable]) |>
    x -> combine(x, :value => sum => :value)
    ),
        callibrated_national_data.sets
)

# TMP.table[findall(in(test_subtables), TMP.table.subtable), :] =
# print(filter(:subtable => !(in(test_subtables)), TMP.table))
# # TMP.table = vcat(filter(:subtable => !(in(test_subtables)), TMP.table),
#  filter(:subtable => in(test_subtables), TMP.table) |>
#         x -> leftjoin(
#         x,
#         X,
#         on = :sectors=>:commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform!(x,
#         [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :subtable, :year]) |>
#     x -> combine(x, :value => sum => :value))
# # get_table(callibrated_national_data)[findall(in(test_subtables), TMP.table.subtable), :]
# get_subtable(TMP, "intermediate_supply") |>
# df[findall(in(["stringA", "stringB", "stringC"]), df.col), :]
# TMP = NationalTable(
#     TMP.table[findall(in(test_subtables), TMP.table.subtable), :] |> 
#         x -> leftjoin(
#         x,
#         X,
#         on = :sectors=>:commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform!(x,
#         [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :subtable, :year]) |>
#     x -> combine(x, :value => sum => :value),
# TMP.sets)
    
# testys.commodities = Symbol.(testys.commodities); testys.sectors = Symbol.(testys.sectors)   

# TMP,M = calibrate(TMPraw)



# [i for i in sort!(unique(TMP.table[:,:sectors])) ]
# [i for i in unique(TMP.table[:,:commodities])]






get_subtable(TMP, "labor_demand")
get_subtable(TMP, "intermediate_demand")
get_subtable(TMP, "exports")
get_subtable(callibrated_national_data, "intermediate_demand")
get_subtable(callibrated_national_data_summary, "intermediate_demand")

subs_TMP = get_table(TMP) |>
    x -> unique(x, :subtable)

# function bea_io_notations() in core_data_defines.jl from WiNDC.jl itelf has BEA to WiNDC.

for i in subs_TMP[:,:subtable]
    println(i,": ",length(unique(get_subtable(TMP,i)[!,:commodities])),": ",length(unique(get_subtable(raw_national_data,i)[!,:commodities]))) # Compared to get_subtable(raw_national_data,i)
end
TMP
subs_TMP

x_0[yr,:]  #	    "Exports of goods and services",
ex_tmp = get_subtable(TMP,"exports"); ex_tmp.commodities=Symbol.(ex_tmp.commodities)
ex = DataFrame([x_0[yr,:].axes[1] x_0[yr,:].data], [:Win, :exp]); ex.Win = Symbol.(ex.Win)
ex_bea = leftjoin(get_subtable(national_data_summary, "exports")[:,[1,4]], unique(Codes[:,[1,4]]), on=:commodities=>:WiNDC_summary); ex_bea.WiNDC=Symbol.(ex_bea.WiNDC)
comp_x_1 = outerjoin(ex_tmp,ex, on=:commodities=>:Win) 
comp_x = outerjoin(comp_x_1,ex_bea, on=:commodities=>:WiNDC, makeunique=true)
print(comp_x)# "Summary doesn't split for 11CA for (agr &) fof, or 525	to fin hou ore"

m_0[yr,:]  #	    "Imports",
im_tmp = get_subtable(TMP,"imports")#; im_tmp.commodities=Symbol.(im_tmp.commodities)
im = DataFrame([m_0[yr,:].axes[1] m_0[yr,:].data], [:Win, :im]); im.Win=string.(im.Win)
im_bea = leftjoin(get_subtable(national_data_summary, "imports")[:,[1,4]], unique(Codes[:,[1,4]]), on=:commodities=>:WiNDC_summary);# ex_bea.WiNDC=Symbol.(ex_bea.WiNDC)
comp_m_1 = outerjoin(im_tmp,im, on=:commodities=>:Win) 
comp_m = outerjoin(comp_m_1,im_bea, on=:commodities=>:WiNDC, makeunique=true)
print(comp_m)

ms_0[yr,:,:]  #	"Margin supply",
ms_tmp = get_subtable(TMP,"margin_supply"); ms_tmp.commodities=Symbol.(ms_tmp.commodities)
ms_tmp[:,:sectors] = [item == "TRADE " ? "trd" : "trn" for item in ms_tmp[:,:sectors]]; ms_tmp.sectors=Symbol.(ms_tmp.sectors)
ms = stack(DataFrame([ms_0.axes[2] ms_0[yr,:,:].data], [:Win ; ms_0.axes[3]]),2:3); ms.variable=Symbol.(ms.variable); ms.Win=Symbol.(ms.Win)
comp_ms = outerjoin(ms_tmp,ms, on=[:commodities=>:Win, :sectors=>:variable], renamecols=""=>"_ms_0")
print(comp_ms)
md_0[yr,:,:]  #	"Margin demand",
md_tmp = get_subtable(TMP,"margin_demand"); md_tmp.commodities=Symbol.(md_tmp.commodities)
md_tmp[:,:sectors] = [item == "TRADE " ? "trd" : "trn" for item in md_tmp[:,:sectors]]; md_tmp.sectors=Symbol.(md_tmp.sectors)
md = stack(DataFrame([md_0.axes[3] transpose(md_0[yr,:,:].data)], [:Win ; md_0.axes[2]]),2:3); md.variable=Symbol.(md.variable); md.Win=Symbol.(md.Win)
comp_md = outerjoin(md_tmp,md, on=[:commodities=>:Win,:sectors=>:variable], renamecols=""=>"_md_0")
print(comp_md)

fd_0[yr,:,:pce]  #	"Final demand", # pce
pce_tmp = get_subtable(TMP,"personal_consumption"); pce_tmp.commodities=Symbol.(pce_tmp.commodities)
pce = DataFrame([[i for x in  fd_0[yr,:,:pce].axes for i in x] fd_0[yr,:,:pce].data], [:Win, :pce])
comp_fdpce = outerjoin(pce_tmp,pce, on=:commodities=>:Win)
print(comp_fdpce)
# Exog fd
fd = stack(DataFrame([fd_0.axes[2] fd_0[yr,:,[x for x in FD if x!=:pce]].data], [:Win ; [x for x in FD if x !=:pce]]),2:18); fd.variable=Symbol.(fd.variable); fd.Win=Symbol.(fd.Win)
fd_temp = get_subtable(TMP,"exogenous_final_demand"); #fd_temp.sectors = Symbol.(fd_temp.sectors); 
fdbea = get_subtable(national_data_summary, "exogenous_final_demand")
fd_bea = leftjoin(fdbea, Codes, on=[:sectors=>:BEA_summary])#; fd_bea.WiNDC=Symbol.(fd_bea.WiNDC)
comp_x_1 = outerjoin(ex_tmp,ex, on=:commodities=>:Win) 
fd_tmp = leftjoin(fd_temp, Codes, on =:sectors=>:BEA_detail)[:,[1,2,4,7]];  fd_tmp.WiNDC = Symbol.(fd_tmp.WiNDC); fd_tmp.commodities = Symbol.(fd_tmp.commodities)
comp_fd = outerjoin(fd_tmp, fd, on=[:WiNDC=>:variable, :commodities=>:Win], renamecols = "" => "_fd_0")[:1:125,:]
print(comp_fd)

bop_tmp = WiNDC.balance_of_payments(TMP)
bop = DataFrame([bopdef_0.axes[1][bopdef_0.axes[1].==[yr]] bopdef_0[yr]],[:year,:value_0]); bop.year = parse.(Int,string.(bop.year))  #	"Balance of payments deficit",
comp_bop = leftjoin(bop_tmp,bop, on=:year)
print(comp_bop)

a_0[yr,:]  #	    "Armington supply",
a_tmp = WiNDC.armington_supply(TMP); a_tmp.commodities=Symbol.(a_tmp.commodities)
a = DataFrame([a_0[yr,:].axes[1] a_0[yr,:].data], [:Win, :a]); a.Win = Symbol.(a.Win)
comp_a = outerjoin(a_tmp,a, on=:commodities=>:Win)
print(comp_a)
y_0[yr,:]   #	"Gross output",
y_tmp = WiNDC.gross_output(TMP); y_tmp.commodities=Symbol.(y_tmp.commodities)
y = DataFrame([y_0[yr,:].axes[1] y_0[yr,:].data], [:Win, :y]); y.Win = Symbol.(y.Win)
comp_y = outerjoin(y_tmp,y, on=:commodities=>:Win)
print(comp_y)

tm_0[yr,:]  #	"Import tariff"; Initial, for price 
tm_tmp = WiNDC.import_tariff_rate(TMP); tm_tmp.commodities=Symbol.(tm_tmp.commodities)
tm = DataFrame([tm_0[yr,:].axes[1] tm_0[yr,:].data], [:Win, :tm]); tm.Win = Symbol.(tm.Win)
comp_tm = outerjoin(tm_tmp,tm, on=:commodities=>:Win)
print(comp_tm)
ta_0[yr,:]  #	"Tax net subsidy rate on intermediate demand", benchmark data also for price level
ta_tmp = WiNDC.absorption_tax_rate(TMP); ta_tmp.commodities=Symbol.(ta_tmp.commodities)
ta = DataFrame([ta_0[yr,:].axes[1] ta_0[yr,:].data], [:Win, :ta]); ta.Win = Symbol.(ta.Win)
comp_ta = outerjoin(ta_tmp,ta, on=:commodities=>:Win)
print(comp_ta)
# testys = get_subtable(TMP, "intermediate_supply") |>
#     x -> leftjoin(
#         x,
#         X,
#         on = :sectors=>:commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform(x,
#         [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :year]) |>
#     x -> combine(x, :value => sum => :value)
# testys.commodities = Symbol.(testys.commodities); testys.sectors = Symbol.(testys.sectors)   

# "Complicated"
ys_0[yr,:,:]   #	"Intermediate Supply",
ys_tmp = get_subtable(TMP,"intermediate_supply"); ys_tmp.commodities=Symbol.(ys_tmp.commodities); ys_tmp.sectors=Symbol.(ys_tmp.sectors)
# ys_tmp_win = leftjoin(ys_tmp,X, on= :sectors=>:commodities, renamecols=""=>"_com"); ys_tmp_win.commodities=Symbol.(ys_tmp_win.commodities); ys_tmp_win.name_com=Symbol.(ys_tmp_win.name_com)
ys = stack(DataFrame([ys_0.axes[3] ys_0[yr,:,:].data], [:Win ; [x for x in I]]),2:72); ys.variable=Symbol.(ys.variable); ys.Win=Symbol.(ys.Win)
comp_ys = rightjoin(ys_tmp, ys, on=[:commodities=>:variable, :sectors=>:Win], renamecols=""=>"_win")
print(comp_ys[1:900,:])
# ys = DataFrame([[i for x in  ys_0[yr,:].axes for i in x] ys_0[yr,:].data], [:Win, :ys]); ys.Win = Symbol.(ys.Win)
# print(outerjoin(ys_tmp,y, on=:commodities=>:Win))

# testid = get_subtable(TMP, "intermediate_demand") |>
#     x -> leftjoin(
#         x,
#         X,
#         on = :sectors=>:commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform(x,
#         [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :year]) |>
#     x -> combine(x, :value => sum => :value)
# testid.commodities = Symbol.(testid.commodities); testid.sectors = Symbol.(testid.sectors)   

id_0[yr,:,:]  #	"Intermediate demand",
id_tmp = get_subtable(TMP,"intermediate_demand"); id_tmp.sectors = Symbol.(id_tmp.sectors); id_tmp.commodities = Symbol.(id_tmp.commodities)
# leftjoin(id_tmp,Symbol.(unique(Codes[:,[3,4]])),on=:sectors=>:BEA_detail)
id = stack(DataFrame([id_0.axes[3] id_0[yr,:,:].data], [:Win ; [x for x in I]]),2:72); id.variable=Symbol.(id.variable); id.Win=Symbol.(id.Win)
comp_id = leftjoin(id_tmp, id, on=[:commodities=>:Win, :sectors=>:variable], renamecols=""=>"_win")
print(comp_id)

# testva = append!(get_subtable(TMP,"labor_demand"),get_subtable(TMP,"capital_demand")) |>
#     x -> leftjoin(
#         x,
#         X,
#         on = :sectors=>:commodities,
#         renamecols = "" => "_com"
#     ) |>
#     x -> transform(x,
#         [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
#     ) |>
#     x -> groupby(x, [:commodities, :sectors, :year]) |>
#     x -> combine(x, :value => sum => :value)
# testva.commodities = Symbol.(testva.commodities); testva.sectors = Symbol.(testva.sectors)   

va_0[yr,:,:]  #	"Value added",
va_tmp = append!(get_subtable(TMP,"labor_demand"),get_subtable(TMP,"capital_demand")); va_tmp.commodities=Symbol.(va_tmp.commodities); va_tmp.sectors=Symbol.(va_tmp.sectors)
# va_tmp = leftjoin(va_tmp, fd_codes, on=:commodities =>:detailed)
# va_tmp = leftjoin(va_tmp, Codes, on=:sectors =>:BEA_detail, renamecols=""=>"_win")
# va_tmp = va_tmp[:,[1 ; 2 ;3;4;5;10]]
va = stack(DataFrame([va_0.axes[3] transpose(va_0[yr,[va for va in VA],:].data)], [:Win; [va for va in VA]]),2:3); va.Win = Symbol.(va.Win); va.variable = Symbol.(va.variable)
comp_va = outerjoin(testva,va, on=[:commodities=>:variable, :sectors=>:Win], renamecols=""=>"_win")
print(comp_va)

# other_tax_rate(data::AbstractNationalTable; column = :value, output = :value) = 
# testothtax = get_subtable(TMP, "other_tax") |>
# x -> leftjoin(
#     x,
#     X,
#     on = :sectors=>:commodities,
#     renamecols = "" => "_com"
# ) |>
# x -> transform(x,
#     [:sectors, :name_com] => ByRow((x,y) -> ismissing(y) ? x : y) => :sectors,
# ) |>
# x -> groupby(x, [:commodities, :sectors, :year]) |>
# x -> combine(x, :value => sum => :value)
# testothtax.commodities = Symbol.(testothtax.commodities); testothtax.sectors = Symbol.(testothtax.sectors)   

# # testystest = copy(testys)    
# testother_tax_rate = outerjoin(
#     rename!(testys, :value=>:is) |>
#             x -> groupby(x, filter(y -> y!=:commodities, [:commodities, :sectors, :year])) |>
#             x -> combine(x, :is => sum => :is),
#         rename!(testothtax,:value => :ot) |>
#             x -> select(x, Not(:commodities)) |>
#             x -> groupby(x, [:sectors, :year]) |>
#             x -> combine(x, :ot => sum => :ot),
#         on = filter(y -> y!=:commodities, [:commodities, :sectors, :year])
#     ) |>
#     x -> coalesce.(x,0) |>
#     x -> transform(x,
#         [:is, :ot] => ByRow((i,o) -> i == 0 ? 0 : o/i) => :value
#     ) |>
#     x -> select(x, Not(:is,:ot))
# testother_tax_rate.sectors = Symbol.(testother_tax_rate.sectors)   
## Existing Taxes
ty_0[yr,:]  #	"Output tax rate"
ty_tmp = WiNDC.other_tax_rate(TMP); ty_tmp.sectors=Symbol.(ty_tmp.sectors)
ty = DataFrame([ty_0[yr,:].axes[1] ty_0[yr,:].data], [:Win, :ty]); ty.Win = Symbol.(ty.Win)
comp_ty = outerjoin(testother_tax_rate,ty, on=:sectors=>:Win)
print(comp_ty)

macro varname(arg)
    string(arg)
end
using XLSX
filename="CompareSummary.xlsx"
XLSX.openxlsx(filename, mode="w") do creation_of_the_excelfile
end
function exp_to_xlsx(df, dfname, filename="CompareSummary.xlsx")
    # XLSX.openxlsx(filename, mode="w") do xf
    XLSX.openxlsx(filename,mode="rw") do xf
    XLSX.addsheet!(xf,dfname)
    XLSX.writetable!(xf[dfname], df)#, anchor_cell=XLSX.CellRef("A1"))
    # xf[dfname]=df
    # sheet = xf[dfname]
    # for r in 1:size(df,1), c in 1:size(df,2)
        #  sheet[XLSX.CellRef(r , c )] = df[r,c]
    # end
end
end

comp_x.commodities=string.(comp_x.commodities)
exp_to_xlsx(comp_x, @varname(comp_x))
comp_m.commodities=string.(comp_m.commodities)
exp_to_xlsx(comp_m, @varname(comp_m))
exp_to_xlsx(comp_bop, @varname(comp_bop))
comp_a.commodities=string.(comp_a.commodities)
exp_to_xlsx(comp_a, @varname(comp_a))
comp_y.commodities=string.(comp_y.commodities)
exp_to_xlsx(comp_y, @varname(comp_y))
comp_ms.commodities=string.(comp_ms.commodities); comp_ms.sectors=string.(comp_ms.sectors)
exp_to_xlsx(comp_ms, @varname(comp_ms))
comp_md.commodities=string.(comp_md.commodities); comp_md.sectors=string.(comp_md.sectors)
exp_to_xlsx(comp_md, @varname(comp_md))
comp_fdpce.commodities=string.(comp_fdpce.commodities)
exp_to_xlsx(comp_fdpce, @varname(comp_fdpce))
comp_fd.commodities=string.(comp_fd.commodities); comp_fd.WiNDC=string.(comp_fd.WiNDC)
exp_to_xlsx(comp_fd, @varname(comp_fd))
comp_tm.commodities=string.(comp_tm.commodities)
exp_to_xlsx(comp_tm, @varname(comp_tm))
comp_ta.commodities=string.(comp_ta.commodities)
exp_to_xlsx(comp_ta, @varname(comp_ta))
comp_ys.commodities=string.(comp_ys.commodities); comp_ys.sectors=string.(comp_ys.sectors)
exp_to_xlsx(comp_ys, @varname(comp_ys))
comp_id.commodities=string.(comp_id.commodities); comp_id.sectors=string.(comp_id.sectors)
exp_to_xlsx(comp_id, @varname(comp_id))
comp_va.commodities=string.(comp_va.commodities); comp_va.sectors=string.(comp_va.sectors)
exp_to_xlsx(comp_va, @varname(comp_va))
comp_ty.sectors = string.(comp_ty.sectors) 
exp_to_xlsx(comp_ty, @varname(comp_ty))
