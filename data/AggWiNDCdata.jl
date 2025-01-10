using WiNDC, DataFrames, MPSGE, GamsStructure, CSV
##cd to top level of the data files
raw_data_directory = "./detailed_data/national"
# all_det_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :raw_detailed); # only for  2007, 2012, 2017
# Data is PROJECTED from summary => detailed bc detailed only 2007, 2012, 2017 
all_det_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :detailed);
all_summary_national_data = WiNDC.national_tables(raw_data_directory; aggregation = :summary);
# For the summary level data, no issue using all years. Several years are failing
# calibration for the detailed data. 2022 works, as does 2017. Others do too, but
# those are the "important" ones for me at the moment.

projected_det_national_data_yr = 
    NationalTable(
        get_table(all_det_national_data)|>
            x -> subset(x, :year => ByRow(==(2020))), # |>#,
## TEST: removing use and oth before callibration
##  # filter(x -> x[:commodities] !="S00401" && x[:commodities] !="S00402" && x[:commodities] !="S00300" && x[:commodities] !="S00900"), 
        all_det_national_data.sets
    )


raw_summary_national_data_yr = 
NationalTable(
    get_table(all_summary_national_data)|>
        x -> subset(x, :year => ByRow(==(2020))), 
    all_summary_national_data.sets
)

# projected_det_national_data_yr.table = filter(x->!(x.commodities in(["S00401", "S00402","S00300","S00900"])),projected_det_national_data_yr_full.table)

# Calibrate the model. The JuMP model is also returned. 
callibrated_det_national_data,M = calibrate(projected_det_national_data_yr)
# TMP,MS = calibrate(raw_summary_national_data_yr)
callibrated_summary_national_data,MS = calibrate(raw_summary_national_data_yr)

# check match from callibrated summary to raw summary
checkCall = outerjoin(callibrated_summary_national_data.table, raw_summary_national_data_yr.table, on = [:commodities, :sectors, :year, :subtable], 
renamecols = "" => "_raw")
checkCall.diff = checkCall.value - checkCall.value_raw 
sort!(checkCall, :diff)

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

#These default to their values in the dataset. 

# set_value!.(M[:Output_tax], .5*value.(M[:Output_tax]))

# solve!(M)

# generate_report(M)

########################
## Exploring the data ##
########################

# See the entire table. Most of the columns should be self-explanatory. The really
# important one is `subtable`, which is the human readable component of the table.
get_table(callibrated_det_national_data)

# You can extract a single subtable. Note the subtable column is dropped (should 
# be fixed, I don't care for that).
get_subtable(callibrated_det_national_data, "intermediate_demand")
get_subtable(callibrated_det_national_data, "intermediate_demand")

# Or multiple
get_subtable(callibrated_det_national_data, ["intermediate_demand", "value_added"])

# Turns out, that's in a branch... quick hack, but it drops subtables... 
vcat(get_subtable.(Ref(callibrated_det_national_data), ["intermediate_demand", "value_added"])...)

# You should notice there is no subtable for value added.
get_table(callibrated_det_national_data) |>
    x -> unique(x, :subtable)
# This is because value added is a composite table, it consists of several smaller
# subtables. To see which ones use get_set, the elements comprise value added.
get_set(callibrated_det_national_data, "value_added")
# If we look at labor_demand we will see two elements, labor demand and V00100. 
# The V00100 is the NAICS code for labor compensation. 
get_set(callibrated_det_national_data, "labor_demand")

# example correspondence table
# Currently has WiNDC_plus and WiNDC summary (theoretically to replicate summary but doesn't), but can be any abitraray mapping been BEA_detail, and any aggregation set
Codes=CSV.read(joinpath(raw_data_directory,"./BEA_WiNDC_Detail-Summary_codes.csv"), DataFrame, header=1)
X = DataFrame([unique(Codes[:,:WiNDC_summary]), unique(Codes[:,:WiNDC])], [:commodities, :name])
double_aggregate_subtables = ["intermediate_supply", "intermediate_demand","labor_demand","capital_demand","other_tax"]
TMP = NationalTable( vcat(
    get_table(raw_summary_national_data_yr) |>
    x -> filter(:subtable => !(in(double_aggregate_subtables)), x) |> #ggregate the simple tables separately
    x -> leftjoin(
        x,
        X,
        on = :commodities,
        # renamecols = "" => "_com"
    ) |>
    x ->hcat(x, x.sec_codes.=missing) |>
    x ->x[:,1:7] |>
    x -> rename!(x, :commodities => :codes, :name =>:commodities)
,
# # Now aggregate the tables that also need to aggregate sectors, in 2 steps
    get_table(raw_summary_national_data_yr) |>
    x -> filter(:subtable => (in(double_aggregate_subtables)), x) |>
    x -> leftjoin(
        x,
        X,
        on = :commodities,
        # renamecols = "" => "_com"
    ) |> 
    # |>
    x -> leftjoin(
    x,
    X,
    on = :sectors=>:commodities,
    # )
    # ,
    renamecols = "" => "_sec") |>
    x -> rename!(x, :commodities => :codes, :name =>:commodities, :sectors=>:sec_codes, :name_sec=>:sectors)
    )
    ,
        callibrated_det_national_data.sets
)

# TMP = NationalTable( 
#         get_table(callibrated_summary_national_data) |>
#         x -> leftjoin(
#             x,
#             X,
#             on = :commodities,
            
#             # renamecols .= ["commodities","name"] .=> ["codes", "commodities"]),
#         ) 
#         |>
#         x -> rename!(x, :commodities => :codes, :name =>:commodities)
#         |>
#         x -> leftjoin(
#             x,
#             X,
#             on = :sectors => :commodities
#             ,)
#             |>
#             x -> rename!(x, :sectors => :seccodes, :name =>:sectors)
#         # )
#         ,
#         callibrated_summary_national_data.sets
# )

#######################################################################################
# X is the correspondence table between bea detail and whatever aggregation is set up #
#######################################################################################
# X = DataFrame([Codes[:,:BEA_detail], Codes[:,:WiNDC]], [:commodities, :name])
# X = DataFrame([Codes[:,:BEA_detail], Codes[:,:WiNDC_plus]], [:commodities, :name])
aggcol = :WiNDC_plus # name of the columns just for a value comparison dataframe later
aggcol = :WiNDC # name of the columns just for a value comparison dataframe later

## To callibrate first=> get_table(callibrated_det_national_data): or to aggregate and then callibrate -> get_table(projected_det_national_data_yr) and callibrate below
double_aggregate_subtables = ["intermediate_supply", "intermediate_demand","labor_demand","capital_demand","other_tax"]
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

#Take a look to check and look at the difference
get_subtable(TMP, "labor_demand")
get_subtable(TMP, "intermediate_demand")
get_subtable(TMP, "exports")
get_subtable(callibrated_det_national_data, "intermediate_demand")
get_subtable(callibrated_summary_national_data, "intermediate_demand")

subs_TMP = get_table(TMP) |>
    x -> unique(x, :subtable)
#print the numbers for raw and aggregated
for i in subs_TMP[:,:subtable]
    println(i,": ",length(unique(get_subtable(TMP,i)[!,:commodities])),": ",length(unique(get_subtable(projected_det_national_data_yr,i)[!,:commodities]))) # Compared to get_subtable(projected_det_national_data_yr,i)
end

    # #oil => gas & oil
    # Gas_oil_splits=Dict(:imgas=>0.4308, :exgas=>0.3961, :prodgas=>0.4003, :imoil=>0.5692,:exoil=>0.6039,:prodoil=>0.5997)
    # ### This HAS to be done in the main dataframe, so it's pickup up with the tax functions etc.
    # push!(TMP.table, (subtable="exports",commodities="gas", sectors=only(filter(x->(x.subtable=="exports") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="exports") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="exports").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:exgas]))
    # TMP.table[(TMP.table.subtable.=="exports").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="exports").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:exoil]
    # push!(TMP.table, (subtable="imports",commodities="gas", sectors=only(filter(x->(x.subtable=="imports") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="imports") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="imports").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:imgas]))
    # TMP.table[(TMP.table.subtable.=="imports").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="imports").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:imoil]
    # # #### no oil&gas import, nothing to split 
    # # # push!(ms_tmp, (commodities=:gas, sectors=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trn),ms_tmp)[:,2]),year=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trn),ms_tmp)[:,3]),value=only(ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trn),:value])*Gas_oil_splits[:prodgas]))
    # # # push!(ms_tmp, (commodities=:gas, sectors=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trd),ms_tmp)[:,2]),year=only(filter(x->(x.commodities ==:oil) .& (x.sectors ==:trd),ms_tmp)[:,3]),value=only(ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trd),:value])*Gas_oil_splits[:prodgas]))
    # # ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trn),:value]=ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trn),:value]*Gas_oil_splits[:prodoil]
    # # ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trd),:value]=ms_tmp[(ms_tmp.commodities.==:oil) .& (ms_tmp.sectors.==:trd),:value]*Gas_oil_splits[:prodoil]
    # push!(TMP.table, (subtable="margin_demand",commodities="gas", sectors=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRANS"),TMP.table)[:,2]),year=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRANS"),TMP.table)[:,3]),value=only(TMP.table[(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRANS"),:value])*Gas_oil_splits[:prodgas]))
    # push!(TMP.table, (subtable="margin_demand",commodities="gas", sectors=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRADE "),TMP.table)[:,2]),year=only(filter(x->(x.subtable =="margin_demand") .&(x.commodities =="oil") .& (x.sectors =="TRADE "),TMP.table)[:,3]),value=only(TMP.table[(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRADE "),:value])*Gas_oil_splits[:prodgas]))
    # TMP.table[(TMP.table.subtable.=="margin_demand").&(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRANS"),:value] =TMP.table[(TMP.table.subtable.=="margin_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRANS"),:value]*Gas_oil_splits[:prodoil]
    # TMP.table[(TMP.table.subtable.=="margin_demand").&(TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRADE "),:value]=TMP.table[(TMP.table.subtable.=="margin_demand").& (TMP.table.commodities.=="oil") .& (TMP.table.sectors.=="TRADE "),:value]*Gas_oil_splits[:prodoil]

    # # "personal_consumption" # 0 ie. not in TMP, nothing to split (this time)
    # # push!(TMP.table, (subtable="personal_consumption",commodities="gas", sectors=only(filter(x->(x.subtable=="personal_consumption") & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="personal_consumption") & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="personal_consumption").& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:consgas]))
    # # TMP.table[(TMP.table.subtable.=="personal_consumption").& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="personal_consumption").& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:consoil]
    # for s in unique(TMP.table[TMP.table.subtable.=="exogenous_final_demand",:sectors])
    #     println(s)
    #     if size(TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s) .& (TMP.table.commodities.=="oil"),:],1)>0  # But actually just
    #     push!(TMP.table, (subtable="exogenous_final_demand",commodities="gas", sectors=only(filter(x->(x.subtable=="exogenous_final_demand") & (x.sectors==s) & (x.commodities=="oil"),TMP.table)[:,2]),year=only(filter(x->(x.subtable=="exogenous_final_demand") & (x.sectors==s) & (x.commodities=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s) .& (TMP.table.commodities.=="oil"),:value])*Gas_oil_splits[:prodgas]))
    #     TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s) .& (TMP.table.commodities.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="exogenous_final_demand").& (TMP.table.sectors.==s).& (TMP.table.commodities.=="oil"),:value]*Gas_oil_splits[:prodoil]
    # end
    # end
    # push!(TMP.table, (subtable="labor_demand",sectors="gas", commodities=only(filter(x->(x.subtable=="labor_demand") & (x.sectors=="oil"),TMP.table)[:,1]),year=only(filter(x->(x.subtable=="labor_demand") & (x.sectors=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="labor_demand").& (TMP.table.sectors.=="oil"),:value])*Gas_oil_splits[:prodgas]))
    # TMP.table[(TMP.table.subtable.=="labor_demand").& (TMP.table.sectors.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="labor_demand").& (TMP.table.sectors.=="oil"),:value]*Gas_oil_splits[:prodoil]
    # push!(TMP.table, (subtable="capital_demand",sectors="gas", commodities=only(filter(x->(x.subtable=="capital_demand") & (x.sectors=="oil"),TMP.table)[:,1]),year=only(filter(x->(x.subtable=="capital_demand") & (x.sectors=="oil"),TMP.table)[:,3]), value=only(TMP.table[(TMP.table.subtable.=="capital_demand").& (TMP.table.sectors.=="oil"),:value])*Gas_oil_splits[:prodgas]))
    # TMP.table[(TMP.table.subtable.=="capital_demand").& (TMP.table.sectors.=="oil"),:value]=TMP.table[(TMP.table.subtable.=="capital_demand").& (TMP.table.sectors.=="oil"),:value]*Gas_oil_splits[:prodoil]

    # "intermediate_supply"

    # "intermediate_demand"

###################################
# Joins to compare to BEA summary #
###################################
# TODO fix the joins for exp/imp, bc the bea summary join doesn't work with the new uti split is mulitple joining on uti, and min/coa on min
## Make DataFrames to compare values
x_0[yr,:]  #	    "Exports of goods and services",
ex_tmp = get_subtable(TMP,"exports"); ex_tmp.commodities=Symbol.(ex_tmp.commodities)

# ex_tmp = leftjoin(ex_tmp1, unique(Codes[:,[:WiNDC_summary,aggcol]]), on=:commodities=>:WiNDC_summary);  ex_tmp.WiNDC=Symbol.(ex_tmp.WiNDC)
ex = DataFrame([x_0[yr,:].axes[1] x_0[yr,:].data], [:Win, :exp]); ex.Win = Symbol.(ex.Win)

ex_bea = leftjoin(filter(x->x.year==2020,get_subtable(all_summary_national_data,"exports"))[:,[1,4]], unique(Codes[:,[:WiNDC_summary,aggcol]]), on=:commodities=>:WiNDC_summary); 
ex_bea[!,aggcol] = Symbol.(ex_bea[!,aggcol])

comp_x_1 = outerjoin(ex_tmp,ex, on=:commodities=>:Win) 
comp_x = outerjoin(comp_x_1,ex_bea, on=:commodities=>aggcol, renamecols=""=>"_beasum")#, makeunique=true)
print(comp_x)# "BEA Summary doesn't split for 11CA for (agr &) fof, or 525	to fin hou ore"

m_0[yr,:]  #	    "Imports",
# im_tmp = filter(x->x.subtable=="imports",TMP.table); im_tmp.commodities=Symbol.(im_tmp.commodities)
im_tmp = get_subtable(TMP,"imports"); #im_tmp.commodities=Symbol.(im_tmp.commodities)
im = DataFrame([m_0[yr,:].axes[1] m_0[yr,:].data], [:Win, :im]); im.Win=string.(im.Win)
im_bea = leftjoin(filter(x->x.year==2020,get_subtable(all_summary_national_data, "imports"))[:,[1,4]], unique(Codes[:,[:WiNDC_summary,aggcol]]), on=:commodities=>:WiNDC_summary);# ex_bea.WiNDC=Symbol.(ex_bea.WiNDC)
# comp_m_1 = outerjoin(im_tmp,im, on=:name_com=>:Win) 
comp_m_1 = outerjoin(im_tmp,im, on=:commodities=>:Win) 
comp_m = outerjoin(comp_m_1,im_bea, on=:commodities=>aggcol, makeunique=true)
print(comp_m)

ms_0[yr,:,:]  #	"Margin supply",
ms_tmp = get_subtable(TMP,"margin_supply"); ms_tmp.commodities=Symbol.(ms_tmp.commodities)
ms_tmp[:,:sectors] = [item == "Trans" ? "trn" : "trd" for item in ms_tmp[:,:sectors]]; ms_tmp.sectors=Symbol.(ms_tmp.sectors)
# ms_tmp[:,:sectors] = [item == "TRADE " ? "trd" : "trn" for item in ms_tmp[:,:sectors]]; ms_tmp.sectors=Symbol.(ms_tmp.sectors)
ms = stack(DataFrame([ms_0.axes[2] ms_0[yr,:,:].data], [:Win ; ms_0.axes[3]]),2:3); ms.variable=Symbol.(ms.variable); ms.Win=Symbol.(ms.Win)
comp_ms = outerjoin(ms_tmp,ms, on=[:commodities=>:Win, :sectors=>:variable], renamecols=""=>"_ms_0")
print(comp_ms)

md_0[yr,:,:]  #	"Margin demand",
md_tmp = get_subtable(TMP,"margin_demand"); md_tmp.commodities=Symbol.(md_tmp.commodities)
md_tmp[:,:sectors] = [item == "Trans" ? "trn" : "trd" for item in md_tmp[:,:sectors]]; md_tmp.sectors=Symbol.(md_tmp.sectors)
# md_tmp[:,:sectors] = [item == "TRADE " ? "trd" : "trn" for item in md_tmp[:,:sectors]]; md_tmp.sectors=Symbol.(md_tmp.sectors)
md = stack(DataFrame([md_0.axes[3] transpose(md_0[yr,:,:].data)], [:Win ; md_0.axes[2]]),2:3); md.variable=Symbol.(md.variable); md.Win=Symbol.(md.Win)
comp_md = outerjoin(md_tmp,md, on=[:commodities=>:Win,:sectors=>:variable], renamecols=""=>"_md_0")
print(comp_md)

fd_0[yr,:,:pce]  #	"Final demand", # pce
pce_tmp = get_subtable(TMP,"personal_consumption"); pce_tmp.commodities=Symbol.(pce_tmp.commodities)
pce = DataFrame([[i for x in  fd_0[yr,:,:pce].axes for i in x] fd_0[yr,:,:pce].data], [:Win, :pce])
comp_fdpce = outerjoin(pce_tmp,pce, on=:commodities=>:Win)
print(comp_fdpce)

# TODO fix join for this too maybe, not showing the new disaggs: it's a 3-way
# Exog fd
fd = stack(DataFrame([fd_0.axes[2] fd_0[yr,:,[x for x in FD if x!=:pce]].data], [:Win ; [x for x in FD if x !=:pce]]),2:18); fd.variable=Symbol.(fd.variable); fd.Win=Symbol.(fd.Win)
fd_temp = get_subtable(TMP,"exogenous_final_demand"); #fd_temp.sectors = Symbol.(fd_temp.sectors); 
fdbea = filter(x->x.year==2020,get_subtable(all_summary_national_data, "exogenous_final_demand"))
# fd_bea = leftjoin(fdbea, Codes, on=[:sectors=>:BEA_summary])#; fd_bea.WiNDC=Symbol.(fd_bea.WiNDC)
# comp_x_1 = outerjoin(ex_tmp,ex, on=:commodities=>:Win) 
# fd_tmp = leftjoin(fd_temp, Codes, on =:sectors=>:BEA_summary)[:,[:commodities,:sectors,:value,aggcol]];  fd_tmp[!,aggcol] = Symbol.(fd_tmp[!,aggcol]); fd_tmp.commodities = Symbol.(fd_tmp.commodities)
# fd_tmp = leftjoin(fd_temp, Codes, on =:sectors=>:BEA_detail)[:,[:commodities,:sectors,:value,aggcol]];  fd_tmp[!,aggcol] = Symbol.(fd_tmp[!,aggcol]); fd_tmp.commodities = Symbol.(fd_tmp.commodities)
comp_fd = outerjoin(fd_tmp, fd, on=[aggcol=>:variable, :commodities=>:Win], renamecols = "" => "_fd_0")[:1:125,:]
print(comp_fd)

va_0[yr,:,:]  #	"Value added",
va_tmp = append!(get_subtable(TMP,"labor_demand"),get_subtable(TMP,"capital_demand")); va_tmp.commodities=Symbol.(va_tmp.commodities); va_tmp.sectors=Symbol.(va_tmp.sectors)
va = stack(DataFrame([va_0.axes[3] transpose(va_0[yr,[va for va in VA],:].data)], [:Win; [va for va in VA]]),2:3); va.Win = Symbol.(va.Win); va.variable = Symbol.(va.variable)
comp_va = outerjoin(va_tmp,va, on=[:commodities=>:variable, :sectors=>:Win], renamecols=""=>"_win")
print(comp_va)

ys_0[yr,:,:]   #	"Intermediate Supply",
ys_tmp = get_subtable(TMP,"intermediate_supply"); ys_tmp.commodities=Symbol.(ys_tmp.commodities); ys_tmp.sectors=Symbol.(ys_tmp.sectors)
ys = stack(DataFrame([ys_0.axes[3] ys_0[yr,:,:].data], [:Win ; [x for x in I]]),2:72); ys.variable=Symbol.(ys.variable); ys.Win=Symbol.(ys.Win)
comp_ys = outerjoin(ys_tmp, ys, on=[:commodities=>:variable, :sectors=>:Win], renamecols=""=>"_win")
print(comp_ys[1:900,:])


id_0[yr,:,:]  #	"Intermediate demand",
id_tmp = get_subtable(TMP,"intermediate_demand"); id_tmp.sectors = Symbol.(id_tmp.sectors); id_tmp.commodities = Symbol.(id_tmp.commodities)
id = stack(DataFrame([id_0.axes[3] id_0[yr,:,:].data], [:Win ; [x for x in I]]),2:72); id.variable=Symbol.(id.variable); id.Win=Symbol.(id.Win)
comp_id = outerjoin(id_tmp, id, on=[:commodities=>:Win, :sectors=>:variable], renamecols=""=>"_win")
print(comp_id[1:1000,:])

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

ty_0[yr,:]  #	"Output tax rate"
ty_tmp = WiNDC.other_tax_rate(TMP); ty_tmp.sectors=Symbol.(ty_tmp.sectors)
ty = DataFrame([ty_0[yr,:].axes[1] ty_0[yr,:].data], [:Win, :ty]); ty.Win = Symbol.(ty.Win)
comp_ty = outerjoin(ty_tmp,ty, on=:sectors=>:Win)
print(comp_ty)


# # #############################################################################
# # # Generate DenseAxisArrays for WiNDC PLUS national/MultiNat model, one year only #
# # # Add 0s and missing indices (messy, but kinda works)
# # #############################################################################
# # Iplus1 = copy(I) To generate DAAs for WiNDC aggregation
# Iplus1 = [[x for x in I if x∉[:uti]]; [:uel,:ugs, :uwt, :coa, :use, :oth, :gas]]

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
# fd_m0exog = unstack(fd_tmp,:commodities, :WiNDC_plus,:value, fill=0)
# fd_m0df = coalesce.(outerjoin(fd_m0exog,rename(pce_tmp[:,[:commodities,:value]],:value=>:pce), on=:commodities),0)
# othersectors_fddf_tmp = [x for x in Iplus1 if x∉fd_m0df.commodities]
# fd_m0 = DenseAxisArray(vcat(Matrix(fd_m0df[:,2:end]),zeros(length(othersectors_fddf_tmp),size(fd_m0df,2)-1)), vcat(fd_m0df.commodities,othersectors_fddf_tmp), Symbol.(names(fd_m0df)[2:end]))
# bop_m0 = DenseAxisArray(bop_tmp.value, [:2020])
# othersectors_a_tmp = [x for x in Iplus1 if x∉a_tmp.commodities]
# a_m0 = DenseAxisArray(vcat(a_tmp.value,zeros(length(othersectors_a_tmp))), vcat(a_tmp.commodities,othersectors_a_tmp))
# othersectors_y_tmp = [x for x in Iplus1 if x∉y_tmp.commodities]# Nice this has everything, but still works
# y_m0 = DenseAxisArray(vcat(y_tmp.value,zeros(length(othersectors_y_tmp))), vcat(y_tmp.commodities,othersectors_y_tmp))
# othersectors_tm_tmp = [x for x in Iplus1 if x∉tm_tmp.commodities]
# tm_m0 = DenseAxisArray(vcat(tm_tmp.value,zeros(length(othersectors_tm_tmp))), vcat(tm_tmp.commodities,othersectors_tm_tmp))
# othersectors_ta_tmp = [x for x in Iplus1 if x∉ta_tmp.commodities]
# ta_m0 = DenseAxisArray(vcat(ta_tmp.value,zeros(length(othersectors_ta_tmp))), vcat(ta_tmp.commodities,othersectors_ta_tmp))
# ys_m0df = unstack(ys_tmp,:commodities, :sectors,:value, fill=0)
# othersectors_ys_tmp = [x for x in Iplus1 if x∉ys_tmp.commodities]# None, so no set up YET to add sectors
# othersectors_ys_tmp2 = [x for x in Iplus1 if x∉ unique(ys_tmp.sectors)]
# for i in (1:length(othersectors_ys_tmp2))
#     ys_m0df[!,othersectors_ys_tmp2[i]] .= 0
# end
# ys_m0 = DenseAxisArray(Matrix(ys_m0df[!,2:end]), ys_m0df.commodities, Symbol.(names(ys_m0df)[2:end]))
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



# macro varname(arg)
#     string(arg)
# end
# using XLSX
# filename="CompareSummary.xlsx"
# XLSX.openxlsx(filename, mode="w") do creation_of_the_excelfile
# end
# function exp_to_xlsx(df, dfname, filename="CompareSummary.xlsx")
#     XLSX.openxlsx(filename,mode="rw") do xf
#     XLSX.addsheet!(xf,dfname)
#     XLSX.writetable!(xf[dfname], df)
# end
# end







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
comp_fd.commodities=string.(comp_fd.commodities); comp_fd[!,aggcol]=string.(comp_fd[!,aggcol])
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
