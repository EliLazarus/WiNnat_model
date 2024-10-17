# Adapted from WiNDC National Model
using MPSGE
using DataFrames, JLD2
using JuMP
using MPSGE.JuMP.Containers
using CSV, Plots
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S
Sectors = CSV.read("Sectors.csv", DataFrame);

I = [i for i∈S[:i] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
J = [i for i∈S[:j] if i∉[:use,:oth]] # Index for WiNDC BEA Sectors
CH4sectors = [:agr,:min,:pip,:oil,:wst] #:uti? # subset index for relevant CH4 mitigation sectors (VA slack in benchmark)
VA = [va for va∈S[:va] if va!=:othtax] # Index Value Added (compen = returns to labour/wage, 'surplus' = returns to Kapital)
FD = S[:fd] # Index for final demand categories (pce, and investment types)
TS = S[:ts] #index for taxes/subsidies
YR = S[:yr] # Index for years for potential multi year runs
M = S[:m] # Index for margins (transport and trade)
a_0 = P[:a_0] #	    "Armington supply",
id_0 = P[:id_0] #	"Intermediate demand",
ys_0 = P[:ys_0]#	"Sectoral supply",
va_0 = P[:va_0] #	"Value added",
md_0 = P[:md_0] #	"Margin demand",
fd_0 = P[:fd_0] #	"Final demand",
m_0 = P[:m_0] #	    "Imports",
ms_0 = P[:ms_0] #	"Margin supply",
bopdef_0 = P[:bopdef_0] #	"Balance of payments deficit",
x_0 = P[:x_0] #	    "Exports of goods and services",
y_0 = P[:y_0]  #	"Gross output",
## Existing Taxes
ty_0 = P[:ty_0] #	"Output tax rate"
tm_0 = P[:tm_0] #	"Import tariff"; Initial, for price 
ta_0 = P[:ta_0] #	"Tax net subsidy rate on intermediate demand", benchmark data also for price level


# ##### - *** - Split :min into :min and :mnr as first pass VERY rough separation, splitting everything 50/50 until I have time to go thorugh all data - *** - ####
#             push!(Sectors, ["mnr" "mining_mineral" "mineral, non-coal mining, split as 1/2 of mining" Sectors[71,4]/2 Sectors[71,5]/2 "$(round((Sectors[71,4]*0.5)/((Sectors[71,4]*0.5)+Sectors[71,5]*0.5),digits=2))%"])
#             #	"Intermediate demand",
#             #  Split min into min and mnr along 2nd dimenstion (using permutdims), with values as min x 0.5 
#             id_0mnr1 = DenseAxisArray(permutedims([id_0[:,:ppd,:];;;id_0[:,:res,:];;;id_0[:,:com,:];;;id_0[:,:amb,:];;;id_0[:,:fbp,:];;;id_0[:,:rec,:];;;id_0[:,:con,:];;;id_0[:,:agr,:];;;id_0[:,:eec,:];;;id_0[:,:fnd,:];;;id_0[:,:pub,:];;;id_0[:,:hou,:];;;id_0[:,:fbt,:];;;id_0[:,:ins,:];;;id_0[:,:tex,:];;;id_0[:,:leg,:];;;id_0[:,:fen,:];;;id_0[:,:uti,:];;;id_0[:,:nmp,:];;;id_0[:,:brd,:];;;id_0[:,:bnk,:];;;id_0[:,:ore,:];;;id_0[:,:edu,:];;;id_0[:,:ote,:];;;id_0[:,:man,:];;;id_0[:,:mch,:];;;id_0[:,:dat,:];;;id_0[:,:amd,:];;;id_0[:,:oil,:];;;id_0[:,:hos,:];;;id_0[:,:rnt,:];;;id_0[:,:pla,:];;;id_0[:,:fof,:];;;id_0[:,:fin,:];;;id_0[:,:tsv,:];;;id_0[:,:nrs,:];;;id_0[:,:sec,:];;;id_0[:,:art,:];;;id_0[:,:mov,:];;;id_0[:,:fpd,:];;;id_0[:,:slg,:];;;id_0[:,:pri,:];;;id_0[:,:grd,:];;;id_0[:,:pip,:];;;id_0[:,:sle,:];;;id_0[:,:osv,:];;;id_0[:,:trn,:];;;id_0[:,:smn,:];;;id_0[:,:fmt,:];;;id_0[:,:pet,:];;;id_0[:,:mvt,:];;;id_0[:,:cep,:];;;id_0[:,:wst,:];;;id_0[:,:mot,:];;;id_0[:,:adm,:];;;id_0[:,:soc,:];;;id_0[:,:alt,:];;;id_0[:,:pmt,:];;;id_0[:,:trk,:];;;id_0[:,:fdd,:];;;id_0[:,:gmt,:];;;id_0[:,:wtt,:];;;id_0[:,:wpd,:];;;id_0[:,:wht,:];;;id_0[:,:wrh,:];;;id_0[:,:ott,:];;;id_0[:,:che,:];;;id_0[:,:air,:];;;id_0[:,:mmf,:];;;id_0[:,:otr,:];;;id_0[:,:min,:]*0.5;;;id_0[:,:min,:]*0.5],[1,3,2]), axes(id_0)[1], [I ; :mnr], axes(id_0)[3])
#             # Then also Split min into min and mnr along 2nd dimenstion (using permutdims), with values as min x 0.5 
#             id_0 = DenseAxisArray([id_0mnr1[:,:,:ppd];;;id_0mnr1[:,:,:res];;;id_0mnr1[:,:,:com];;;id_0mnr1[:,:,:amb];;;id_0mnr1[:,:,:fbp];;;id_0mnr1[:,:,:rec];;;id_0mnr1[:,:,:con];;;id_0mnr1[:,:,:agr];;;id_0mnr1[:,:,:eec];;;id_0mnr1[:,:,:fnd];;;id_0mnr1[:,:,:pub];;;id_0mnr1[:,:,:hou];;;id_0mnr1[:,:,:fbt];;;id_0mnr1[:,:,:ins];;;id_0mnr1[:,:,:tex];;;id_0mnr1[:,:,:leg];;;id_0mnr1[:,:,:fen];;;id_0mnr1[:,:,:uti];;;id_0mnr1[:,:,:nmp];;;id_0mnr1[:,:,:brd];;;id_0mnr1[:,:,:bnk];;;id_0mnr1[:,:,:ore];;;id_0mnr1[:,:,:edu];;;id_0mnr1[:,:,:ote];;;id_0mnr1[:,:,:man];;;id_0mnr1[:,:,:mch];;;id_0mnr1[:,:,:dat];;;id_0mnr1[:,:,:amd];;;id_0mnr1[:,:,:oil];;;id_0mnr1[:,:,:hos];;;id_0mnr1[:,:,:rnt];;;id_0mnr1[:,:,:pla];;;id_0mnr1[:,:,:fof];;;id_0mnr1[:,:,:fin];;;id_0mnr1[:,:,:tsv];;;id_0mnr1[:,:,:nrs];;;id_0mnr1[:,:,:sec];;;id_0mnr1[:,:,:art];;;id_0mnr1[:,:,:mov];;;id_0mnr1[:,:,:fpd];;;id_0mnr1[:,:,:slg];;;id_0mnr1[:,:,:pri];;;id_0mnr1[:,:,:grd];;;id_0mnr1[:,:,:pip];;;id_0mnr1[:,:,:sle];;;id_0mnr1[:,:,:osv];;;id_0mnr1[:,:,:trn];;;id_0mnr1[:,:,:smn];;;id_0mnr1[:,:,:fmt];;;id_0mnr1[:,:,:pet];;;id_0mnr1[:,:,:mvt];;;id_0mnr1[:,:,:cep];;;id_0mnr1[:,:,:wst];;;id_0mnr1[:,:,:mot];;;id_0mnr1[:,:,:adm];;;id_0mnr1[:,:,:soc];;;id_0mnr1[:,:,:alt];;;id_0mnr1[:,:,:pmt];;;id_0mnr1[:,:,:trk];;;id_0mnr1[:,:,:fdd];;;id_0mnr1[:,:,:gmt];;;id_0mnr1[:,:,:wtt];;;id_0mnr1[:,:,:wpd];;;id_0mnr1[:,:,:wht];;;id_0mnr1[:,:,:wrh];;;id_0mnr1[:,:,:ott];;;id_0mnr1[:,:,:che];;;id_0mnr1[:,:,:air];;;id_0mnr1[:,:,:mmf];;;id_0mnr1[:,:,:otr];;;id_0mnr1[:,:,:min]*0.5;;;id_0mnr1[:,:,:min]*0.5], axes(id_0mnr1)[1], [I ; :mnr],[I ; :mnr])

#             #	"Sectoral supply",
#             #  Split min into min and mnr along 2nd dimenstion (using permutdims), with values as min x 0.5 
#             ys_0mnr1 = DenseAxisArray(permutedims([ys_0[:,:ppd,:];;;ys_0[:,:res,:];;;ys_0[:,:com,:];;;ys_0[:,:amb,:];;;ys_0[:,:fbp,:];;;ys_0[:,:rec,:];;;ys_0[:,:con,:];;;ys_0[:,:agr,:];;;ys_0[:,:eec,:];;;ys_0[:,:fnd,:];;;ys_0[:,:pub,:];;;ys_0[:,:hou,:];;;ys_0[:,:fbt,:];;;ys_0[:,:ins,:];;;ys_0[:,:tex,:];;;ys_0[:,:leg,:];;;ys_0[:,:fen,:];;;ys_0[:,:uti,:];;;ys_0[:,:nmp,:];;;ys_0[:,:brd,:];;;ys_0[:,:bnk,:];;;ys_0[:,:ore,:];;;ys_0[:,:edu,:];;;ys_0[:,:ote,:];;;ys_0[:,:man,:];;;ys_0[:,:mch,:];;;ys_0[:,:dat,:];;;ys_0[:,:amd,:];;;ys_0[:,:oil,:];;;ys_0[:,:hos,:];;;ys_0[:,:rnt,:];;;ys_0[:,:pla,:];;;ys_0[:,:fof,:];;;ys_0[:,:fin,:];;;ys_0[:,:tsv,:];;;ys_0[:,:nrs,:];;;ys_0[:,:sec,:];;;ys_0[:,:art,:];;;ys_0[:,:mov,:];;;ys_0[:,:fpd,:];;;ys_0[:,:slg,:];;;ys_0[:,:pri,:];;;ys_0[:,:grd,:];;;ys_0[:,:pip,:];;;ys_0[:,:sle,:];;;ys_0[:,:osv,:];;;ys_0[:,:trn,:];;;ys_0[:,:smn,:];;;ys_0[:,:fmt,:];;;ys_0[:,:pet,:];;;ys_0[:,:mvt,:];;;ys_0[:,:cep,:];;;ys_0[:,:wst,:];;;ys_0[:,:mot,:];;;ys_0[:,:adm,:];;;ys_0[:,:soc,:];;;ys_0[:,:alt,:];;;ys_0[:,:pmt,:];;;ys_0[:,:trk,:];;;ys_0[:,:fdd,:];;;ys_0[:,:gmt,:];;;ys_0[:,:wtt,:];;;ys_0[:,:wpd,:];;;ys_0[:,:wht,:];;;ys_0[:,:wrh,:];;;ys_0[:,:ott,:];;;ys_0[:,:che,:];;;ys_0[:,:air,:];;;ys_0[:,:mmf,:];;;ys_0[:,:otr,:];;;ys_0[:,:min,:]*0.5;;;ys_0[:,:min,:]*0.5],[1,3,2]), axes(ys_0)[1], [I ; :mnr], axes(ys_0)[3])
#             # Then also Split min into min and mnr along 2nd dimenstion (using permutdims), with values as min x 0.5 
#             ys_0 = DenseAxisArray([ys_0mnr1[:,:,:ppd];;;ys_0mnr1[:,:,:res];;;ys_0mnr1[:,:,:com];;;ys_0mnr1[:,:,:amb];;;ys_0mnr1[:,:,:fbp];;;ys_0mnr1[:,:,:rec];;;ys_0mnr1[:,:,:con];;;ys_0mnr1[:,:,:agr];;;ys_0mnr1[:,:,:eec];;;ys_0mnr1[:,:,:fnd];;;ys_0mnr1[:,:,:pub];;;ys_0mnr1[:,:,:hou];;;ys_0mnr1[:,:,:fbt];;;ys_0mnr1[:,:,:ins];;;ys_0mnr1[:,:,:tex];;;ys_0mnr1[:,:,:leg];;;ys_0mnr1[:,:,:fen];;;ys_0mnr1[:,:,:uti];;;ys_0mnr1[:,:,:nmp];;;ys_0mnr1[:,:,:brd];;;ys_0mnr1[:,:,:bnk];;;ys_0mnr1[:,:,:ore];;;ys_0mnr1[:,:,:edu];;;ys_0mnr1[:,:,:ote];;;ys_0mnr1[:,:,:man];;;ys_0mnr1[:,:,:mch];;;ys_0mnr1[:,:,:dat];;;ys_0mnr1[:,:,:amd];;;ys_0mnr1[:,:,:oil];;;ys_0mnr1[:,:,:hos];;;ys_0mnr1[:,:,:rnt];;;ys_0mnr1[:,:,:pla];;;ys_0mnr1[:,:,:fof];;;ys_0mnr1[:,:,:fin];;;ys_0mnr1[:,:,:tsv];;;ys_0mnr1[:,:,:nrs];;;ys_0mnr1[:,:,:sec];;;ys_0mnr1[:,:,:art];;;ys_0mnr1[:,:,:mov];;;ys_0mnr1[:,:,:fpd];;;ys_0mnr1[:,:,:slg];;;ys_0mnr1[:,:,:pri];;;ys_0mnr1[:,:,:grd];;;ys_0mnr1[:,:,:pip];;;ys_0mnr1[:,:,:sle];;;ys_0mnr1[:,:,:osv];;;ys_0mnr1[:,:,:trn];;;ys_0mnr1[:,:,:smn];;;ys_0mnr1[:,:,:fmt];;;ys_0mnr1[:,:,:pet];;;ys_0mnr1[:,:,:mvt];;;ys_0mnr1[:,:,:cep];;;ys_0mnr1[:,:,:wst];;;ys_0mnr1[:,:,:mot];;;ys_0mnr1[:,:,:adm];;;ys_0mnr1[:,:,:soc];;;ys_0mnr1[:,:,:alt];;;ys_0mnr1[:,:,:pmt];;;ys_0mnr1[:,:,:trk];;;ys_0mnr1[:,:,:fdd];;;ys_0mnr1[:,:,:gmt];;;ys_0mnr1[:,:,:wtt];;;ys_0mnr1[:,:,:wpd];;;ys_0mnr1[:,:,:wht];;;ys_0mnr1[:,:,:wrh];;;ys_0mnr1[:,:,:ott];;;ys_0mnr1[:,:,:che];;;ys_0mnr1[:,:,:air];;;ys_0mnr1[:,:,:mmf];;;ys_0mnr1[:,:,:otr];;;ys_0mnr1[:,:,:min]*0.5;;;ys_0mnr1[:,:,:min]*0.5], axes(ys_0mnr1)[1], [I ; :mnr],[I ; :mnr])

#             #	"Final demand",
#             # 3D, but not i x j, so append along 2nd dimension (using permudims) copy of :min array for :mnr and multiply both :min and :mnr by 0.5
#             fd_0 = DenseAxisArray(permutedims([fd_0[:,:ppd,:];;;fd_0[:,:res,:];;;fd_0[:,:com,:];;;fd_0[:,:amb,:];;;fd_0[:,:fbp,:];;;fd_0[:,:rec,:];;;fd_0[:,:con,:];;;fd_0[:,:agr,:];;;fd_0[:,:eec,:];;;fd_0[:,:fnd,:];;;fd_0[:,:pub,:];;;fd_0[:,:hou,:];;;fd_0[:,:fbt,:];;;fd_0[:,:ins,:];;;fd_0[:,:tex,:];;;fd_0[:,:leg,:];;;fd_0[:,:fen,:];;;fd_0[:,:uti,:];;;fd_0[:,:nmp,:];;;fd_0[:,:brd,:];;;fd_0[:,:bnk,:];;;fd_0[:,:ore,:];;;fd_0[:,:edu,:];;;fd_0[:,:ote,:];;;fd_0[:,:man,:];;;fd_0[:,:mch,:];;;fd_0[:,:dat,:];;;fd_0[:,:amd,:];;;fd_0[:,:oil,:];;;fd_0[:,:hos,:];;;fd_0[:,:rnt,:];;;fd_0[:,:pla,:];;;fd_0[:,:fof,:];;;fd_0[:,:fin,:];;;fd_0[:,:tsv,:];;;fd_0[:,:nrs,:];;;fd_0[:,:sec,:];;;fd_0[:,:art,:];;;fd_0[:,:mov,:];;;fd_0[:,:fpd,:];;;fd_0[:,:slg,:];;;fd_0[:,:pri,:];;;fd_0[:,:grd,:];;;fd_0[:,:pip,:];;;fd_0[:,:sle,:];;;fd_0[:,:osv,:];;;fd_0[:,:trn,:];;;fd_0[:,:smn,:];;;fd_0[:,:fmt,:];;;fd_0[:,:pet,:];;;fd_0[:,:mvt,:];;;fd_0[:,:cep,:];;;fd_0[:,:wst,:];;;fd_0[:,:mot,:];;;fd_0[:,:adm,:];;;fd_0[:,:soc,:];;;fd_0[:,:alt,:];;;fd_0[:,:pmt,:];;;fd_0[:,:trk,:];;;fd_0[:,:fdd,:];;;fd_0[:,:gmt,:];;;fd_0[:,:wtt,:];;;fd_0[:,:wpd,:];;;fd_0[:,:wht,:];;;fd_0[:,:wrh,:];;;fd_0[:,:ott,:];;;fd_0[:,:che,:];;;fd_0[:,:air,:];;;fd_0[:,:mmf,:];;;fd_0[:,:otr,:];;;fd_0[:,:min,:]*0.5;;;fd_0[:,:min,:]*0.5],[1,3,2]), axes(fd_0)[1], [I ; :mnr], axes(fd_0)[3])
#             #	"Value added",
#             # 3D, but not i x j, so append along 3rd dimensions (no permudims needed) copy of :min array for :mnr and multiply both :min and :mnr by 0.5
#             va_0 = DenseAxisArray([va_0[:,:,:ppd];;;va_0[:,:,:res];;;va_0[:,:,:com];;;va_0[:,:,:amb];;;va_0[:,:,:fbp];;;va_0[:,:,:rec];;;va_0[:,:,:con];;;va_0[:,:,:agr];;;va_0[:,:,:eec];;;va_0[:,:,:fnd];;;va_0[:,:,:pub];;;va_0[:,:,:hou];;;va_0[:,:,:fbt];;;va_0[:,:,:ins];;;va_0[:,:,:tex];;;va_0[:,:,:leg];;;va_0[:,:,:fen];;;va_0[:,:,:uti];;;va_0[:,:,:nmp];;;va_0[:,:,:brd];;;va_0[:,:,:bnk];;;va_0[:,:,:ore];;;va_0[:,:,:edu];;;va_0[:,:,:ote];;;va_0[:,:,:man];;;va_0[:,:,:mch];;;va_0[:,:,:dat];;;va_0[:,:,:amd];;;va_0[:,:,:oil];;;va_0[:,:,:hos];;;va_0[:,:,:rnt];;;va_0[:,:,:pla];;;va_0[:,:,:fof];;;va_0[:,:,:fin];;;va_0[:,:,:tsv];;;va_0[:,:,:nrs];;;va_0[:,:,:sec];;;va_0[:,:,:art];;;va_0[:,:,:mov];;;va_0[:,:,:fpd];;;va_0[:,:,:slg];;;va_0[:,:,:pri];;;va_0[:,:,:grd];;;va_0[:,:,:pip];;;va_0[:,:,:sle];;;va_0[:,:,:osv];;;va_0[:,:,:trn];;;va_0[:,:,:smn];;;va_0[:,:,:fmt];;;va_0[:,:,:pet];;;va_0[:,:,:mvt];;;va_0[:,:,:cep];;;va_0[:,:,:wst];;;va_0[:,:,:mot];;;va_0[:,:,:adm];;;va_0[:,:,:soc];;;va_0[:,:,:alt];;;va_0[:,:,:pmt];;;va_0[:,:,:trk];;;va_0[:,:,:fdd];;;va_0[:,:,:gmt];;;va_0[:,:,:wtt];;;va_0[:,:,:wpd];;;va_0[:,:,:wht];;;va_0[:,:,:wrh];;;va_0[:,:,:ott];;;va_0[:,:,:che];;;va_0[:,:,:air];;;va_0[:,:,:mmf];;;va_0[:,:,:otr];;;va_0[:,:,:min]*0.5;;;va_0[:,:,:min]*0.5], axes(va_0)[1], axes(va_0)[2], [I ; :mnr])
#             #	"Margin demand",
#             # 3D, but not i x j, so append along 3rd dimensions (no permudims needed) copy of :min array for :mnr and multiply both :min and :mnr by 0.5
#             md_0 = DenseAxisArray([md_0[:,:,:ppd];;;md_0[:,:,:res];;;md_0[:,:,:com];;;md_0[:,:,:amb];;;md_0[:,:,:fbp];;;md_0[:,:,:rec];;;md_0[:,:,:con];;;md_0[:,:,:agr];;;md_0[:,:,:eec];;;md_0[:,:,:fnd];;;md_0[:,:,:pub];;;md_0[:,:,:hou];;;md_0[:,:,:fbt];;;md_0[:,:,:ins];;;md_0[:,:,:tex];;;md_0[:,:,:leg];;;md_0[:,:,:fen];;;md_0[:,:,:uti];;;md_0[:,:,:nmp];;;md_0[:,:,:brd];;;md_0[:,:,:bnk];;;md_0[:,:,:ore];;;md_0[:,:,:edu];;;md_0[:,:,:ote];;;md_0[:,:,:man];;;md_0[:,:,:mch];;;md_0[:,:,:dat];;;md_0[:,:,:amd];;;md_0[:,:,:oil];;;md_0[:,:,:hos];;;md_0[:,:,:rnt];;;md_0[:,:,:pla];;;md_0[:,:,:fof];;;md_0[:,:,:fin];;;md_0[:,:,:tsv];;;md_0[:,:,:nrs];;;md_0[:,:,:sec];;;md_0[:,:,:art];;;md_0[:,:,:mov];;;md_0[:,:,:fpd];;;md_0[:,:,:slg];;;md_0[:,:,:pri];;;md_0[:,:,:grd];;;md_0[:,:,:pip];;;md_0[:,:,:sle];;;md_0[:,:,:osv];;;md_0[:,:,:trn];;;md_0[:,:,:smn];;;md_0[:,:,:fmt];;;md_0[:,:,:pet];;;md_0[:,:,:mvt];;;md_0[:,:,:cep];;;md_0[:,:,:wst];;;md_0[:,:,:mot];;;md_0[:,:,:adm];;;md_0[:,:,:soc];;;md_0[:,:,:alt];;;md_0[:,:,:pmt];;;md_0[:,:,:trk];;;md_0[:,:,:fdd];;;md_0[:,:,:gmt];;;md_0[:,:,:wtt];;;md_0[:,:,:wpd];;;md_0[:,:,:wht];;;md_0[:,:,:wrh];;;md_0[:,:,:ott];;;md_0[:,:,:che];;;md_0[:,:,:air];;;md_0[:,:,:mmf];;;md_0[:,:,:otr];;;md_0[:,:,:min]*0.5;;;md_0[:,:,:min]*0.5], axes(md_0)[1], axes(md_0)[2], [I ; :mnr])
#             #	"Margin supply",
#             # 3D, but not i x j, so append along 2nd dimension (using permudims) copy of :min array for :mnr and multiply both :min and :mnr by 0.5
#             ms_0 = DenseAxisArray(permutedims([ms_0[:,:ppd,:];;;ms_0[:,:res,:];;;ms_0[:,:com,:];;;ms_0[:,:amb,:];;;ms_0[:,:fbp,:];;;ms_0[:,:rec,:];;;ms_0[:,:con,:];;;ms_0[:,:agr,:];;;ms_0[:,:eec,:];;;ms_0[:,:fnd,:];;;ms_0[:,:pub,:];;;ms_0[:,:hou,:];;;ms_0[:,:fbt,:];;;ms_0[:,:ins,:];;;ms_0[:,:tex,:];;;ms_0[:,:leg,:];;;ms_0[:,:fen,:];;;ms_0[:,:uti,:];;;ms_0[:,:nmp,:];;;ms_0[:,:brd,:];;;ms_0[:,:bnk,:];;;ms_0[:,:ore,:];;;ms_0[:,:edu,:];;;ms_0[:,:ote,:];;;ms_0[:,:man,:];;;ms_0[:,:mch,:];;;ms_0[:,:dat,:];;;ms_0[:,:amd,:];;;ms_0[:,:oil,:];;;ms_0[:,:hos,:];;;ms_0[:,:rnt,:];;;ms_0[:,:pla,:];;;ms_0[:,:fof,:];;;ms_0[:,:fin,:];;;ms_0[:,:tsv,:];;;ms_0[:,:nrs,:];;;ms_0[:,:sec,:];;;ms_0[:,:art,:];;;ms_0[:,:mov,:];;;ms_0[:,:fpd,:];;;ms_0[:,:slg,:];;;ms_0[:,:pri,:];;;ms_0[:,:grd,:];;;ms_0[:,:pip,:];;;ms_0[:,:sle,:];;;ms_0[:,:osv,:];;;ms_0[:,:trn,:];;;ms_0[:,:smn,:];;;ms_0[:,:fmt,:];;;ms_0[:,:pet,:];;;ms_0[:,:mvt,:];;;ms_0[:,:cep,:];;;ms_0[:,:wst,:];;;ms_0[:,:mot,:];;;ms_0[:,:adm,:];;;ms_0[:,:soc,:];;;ms_0[:,:alt,:];;;ms_0[:,:pmt,:];;;ms_0[:,:trk,:];;;ms_0[:,:fdd,:];;;ms_0[:,:gmt,:];;;ms_0[:,:wtt,:];;;ms_0[:,:wpd,:];;;ms_0[:,:wht,:];;;ms_0[:,:wrh,:];;;ms_0[:,:ott,:];;;ms_0[:,:che,:];;;ms_0[:,:air,:];;;ms_0[:,:mmf,:];;;ms_0[:,:otr,:];;;ms_0[:,:min,:]*0.5;;;ms_0[:,:min,:]*0.5],[1,3,2]), axes(ms_0)[1], [I ; :mnr], axes(ms_0)[3])

#             # 2D arrays so append a vector for mnr with 0.5 x values of :min
#             m_0 = DenseAxisArray([m_0.data deepcopy(m_0[:,:min]/2)], axes(m_0)[1], [axes(m_0)[2]; :mnr])  #	    "Imports",
#             # 0.5 times original value for :min
#             m_0[:,:min]=deepcopy(m_0[:,:min]/2)
#             # append a vector for mnr with 0.5 x values of :min
#             a_0 = DenseAxisArray([a_0.data deepcopy(a_0[:,:min]/2)], axes(a_0)[1], [axes(a_0)[2]; :mnr])
#             # 0.5 times original value for :min
#             a_0[:,:min]=deepcopy(a_0[:,:min]/2) #	    "Armington supply",
#             # append a vector for mnr with 0.5 x values of :min
#             x_0  = DenseAxisArray([x_0.data deepcopy(x_0[:,:min].data)/2], axes(x_0)[1], [axes(x_0)[2]; :mnr])
#             # 0.5 times original value for :min
#             x_0[:,:min]=deepcopy(x_0[:,:min]/2) #	    "Exports of goods and services",
#             # append a vector for mnr with 0.5 x values of :min
#             y_0 = DenseAxisArray([y_0.data deepcopy(y_0[:,:min]/2)], axes(y_0)[1], [axes(y_0)[2]; :mnr]) #	"Output tax rate"  Copied bc rate
#             # 0.5 times original value for :min
#             y_0[:,:min]=deepcopy(y_0[:,:min]/2) #	"Gross output",

#             ## Existing Taxes
#             # append a vector for mnr with *copy* of values of :min (because taxes are a ratio)
#             ty_0 = DenseAxisArray([ty_0.data deepcopy(ty_0[:,:min])], axes(ty_0)[1], [axes(ty_0)[2]; :mnr]) #	"Output tax rate"  Copied bc rate 
#             tm_0 = DenseAxisArray([tm_0.data deepcopy(tm_0[:,:min])], axes(tm_0)[1], [axes(tm_0)[2]; :mnr]) #	"Output tax rate"  Copied bc rate 
#             ta_0 = DenseAxisArray([ta_0.data deepcopy(ta_0[:,:min])], axes(ta_0)[1], [axes(ta_0)[2]; :mnr]) #	"Output tax rate"  Copied bc rate
#             I = [I ; :mnr]
#             J = [J ; :mnr]
# ###### - End Split min section - #######

yr = Symbol(2020)

## Base Marginal Abatement Cost EPA data (2020)
MAC_CH4_data=CSV.read("./data/EPA_CH4_MAC_2020_data.csv", DataFrame, header=2, limit=14)
MAC_CH4_totemiss=CSV.read("./data/EPA_CH4_MAC_2020_data.csv", DataFrame, header=2, skipto=17)
# % split for pipelines and oil from MAC for 'GAS' by proportion of combined economic output
pip_of_GAS = sum(ys_0[yr,:pip,:])/((sum(ys_0[yr,:pip,:])+sum(ys_0[yr,:oil,:])))
oil_of_GAS = sum(ys_0[yr,:oil,:])/((sum(ys_0[yr,:pip,:])+sum(ys_0[yr,:oil,:])))
# Aggregate/disaggregate for WiNDC sectors 
MAC_CH4_WiNDC=DataFrame([MAC_CH4_data[:,:cost_per_t], MAC_CH4_data[:,:agr_livestock]+MAC_CH4_data[:,:agr_rice],
 MAC_CH4_data[:,:min], 
 MAC_CH4_data[:,:GAS]*pip_of_GAS,
 MAC_CH4_data[:,:GAS]*oil_of_GAS,
 MAC_CH4_data[:,:wst_land]+MAC_CH4_data[:,:wst_water]],
 [:cost_per_t; CH4sectors])
# Aggregate/disaggregate Total CH4 Emissions for WiNDC sectors, and include total VA for convenience 
MAC_CH4_WiNDC_tot=DataFrame([MAC_CH4_totemiss[:,:cost_per_t], MAC_CH4_totemiss[:,:agr_livestock]+MAC_CH4_totemiss[:,:agr_rice],
 MAC_CH4_totemiss[:,:min], 
 MAC_CH4_totemiss[:,:GAS]*pip_of_GAS,
 MAC_CH4_totemiss[:,:GAS]*oil_of_GAS,
 MAC_CH4_totemiss[:,:wst_land]+MAC_CH4_totemiss[:,:wst_water]],
 [:cost_per_t; CH4sectors])
 push!(MAC_CH4_WiNDC_tot, ["Total va_cost"; [sum(va_0[yr,:,Symbol(sector)]) for sector in names(MAC_CH4_WiNDC[:,2:end])]])

# Initialise df with 0s to fill with calculations
CH4_cost_per_tier = copy(MAC_CH4_WiNDC[:,2:end]); CH4_cost_per_tier[:,:].=0
# Mulitply costs by tons of abatement potential for the first row
    for c in CH4sectors
        CH4_cost_per_tier[1,c] = MAC_CH4_WiNDC[1,:cost_per_t]*MAC_CH4_WiNDC[1,c]
    end
# Calculate Costs per tier of marginal abatement: subtract all emissions with abatement potential below tier, and mulitply by cost at that tier
    for c in CH4sectors
        for i in 2:length(CH4_cost_per_tier[:,1])
            CH4_cost_per_tier[i,c] = MAC_CH4_WiNDC[i,:cost_per_t]*(MAC_CH4_WiNDC[i,c]-MAC_CH4_WiNDC[i-1,c])
        end
    end
### Calculate the weighted cost for cumulative abatement at each level ###
    CH4_cumul_costAll = copy(MAC_CH4_WiNDC[:,2:end]); CH4_cumul_costAll.=0 # Intitialize df
    CH4_cumul_costAll[1,:]=copy(CH4_cost_per_tier[1,:])  # 1st Row copied as just tier costs
    # Rest of dataframe sums cost of each tier with all tiers below
    for c in CH4sectors
        for i in 2:length(CH4_cumul_costAll[:,1])
            CH4_cumul_costAll[i,c] = CH4_cost_per_tier[i,c]+CH4_cumul_costAll[i-1,c]
        end
    end
# Value-Added block names for DenseAxisArray indexing
VAMset = [:VAM5,:VAM10,:VAM15,:VAM20,:VAM30,:VAM40,:VAM50,:VAM100,:VAM500,:VAM1000]
# dfs to DenseAxisArrays and filter to cumulative Marginal costs which are positive for at least one sector i.e. 5th row, $5/t and up
    CH4_cumul_cost = DenseAxisArray(Matrix(CH4_cumul_costAll[5:end,:]),VAMset, CH4sectors)
    CH4_EmMitigated = DenseAxisArray(Matrix(MAC_CH4_WiNDC[5:end,2:end]), VAMset, CH4sectors)
## CH4 Emissions intensity of Standard value-added activities: total emissions/total value-added cost 
    VASInt =  [MAC_CH4_WiNDC_tot[1,c]*10^-3/MAC_CH4_WiNDC_tot[2,c] for c in CH4sectors]
## CH4 emission intensity (emissions/total va cost) for each sector at each abatement tier, subtracts cumulatively abated emissions and adds additional cost of abatement activities
    VAM_CH4EmInt = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Initializing DenseAxisArray
    [VAM_CH4EmInt[v,c] = (MAC_CH4_WiNDC_tot[1,c]*10^-3-CH4_EmMitigated[v,c]*10^-3)/(sum(va_0[yr,:,c])+CH4_cumul_cost[v,c]*10^-3) for v in VAMset for c in CH4sectors]
## Relative cost of VA: (standard va cost + mitigation cost)/(standard cost) - used to multiply BOTH the va[:surplus] and va[:compen] equally in the blocks
    VAM_costover = DenseAxisArray(zeros(length(VAMset),length(CH4sectors)),VAMset,CH4sectors) # Initialise DenseAxisArray
    [VAM_costover[cost,c] = (sum(va_0[yr,:,c])+CH4_cumul_cost[cost,c]*10^-3)/sum(va_0[yr,:,c]) for cost in VAMset for c in CH4sectors]

#####--------Single Mitigation step data set up-----------#####
#### - TODO Also NEEDS UPDATE FOR min split when that's done - ####
            ## Base Data for reference
                # EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDs : EPAnonCO2-report-data-annex-9-30-19_0\NonCO2 MACs-Summary Tables.xlsx
                CH4emissdatadf = DataFrame(Wsector = [:agr,:agr,:min,:oil,:wst,:wst],
                ## EPA Total Emissions per sector, MMt CO2eq 
                EPAemiss =[260.483532*10^-3,13.70952225*10^-3,59.31302643*10^-3,224.8979059*10^-3,111.5049515*10^-3,20.36144996*10^-3],
                ## EPA maximum % of abatement per sector at <$1000/t
                MaxpercMit = [.304,.280,.645,.475,.050,.350],
                PAsector = ["AGRICULTURE, LIVE","AGRICULTURE, RICE","ENERGY, COL",
                "ENERGY, GAS","WASTE, LAN","WASTE, WWR"])
            ## Calculate CH4 Intensity factors from EPA data
            # EPA Non-CO2 Marginal Abatment Curve data, 2019, dataframe because non-unique row IDs
            # Emissions MMt, sum and weighted average with current aggregation (before potential disaggregation of BEA/WiNDC data at some point)
            CH4emiss = DenseAxisArray([260.483532*10^-3+13.70952225*10^-3 59.31302643*10^-3 224.8979059*10^-3*pip_of_GAS 224.8979059*10^-3*oil_of_GAS 111.5049515*10^-3+20.36144996*10^-3
            ## Weighted average mitigated potential per sector
            (260.483532*0.3039385+13.70952225*0.2803694)/(260.483532+13.70952225)  0.6452054 0.4749547 0.4749547 (111.5049515*0.049615808+20.36144996*0.122234272)/(111.5049515+20.36144996)
            ## Mitigation cost per sector in Million $US 2019, for the maximum % mititation potential
            ## EPA Non-CO2 MAC Sum of each $s/ton mit x tons mitigated at that wedge of abatement cost potential - calculated in Excel
            ## cost x sum(MMT for that sector) + cost x sum(MMT additional at that cost for that sector) + etc.
            ## :oil and :pip split EPA "Gas" by proportion of combined output from ys_0
                7090.604857 269.0188218 6862.194574*pip_of_GAS 6862.194574*oil_of_GAS 1795.700322],
            # DenseAxisArray Indices, 2D
            [:EPAemiss :MaxpercMit :MitCostTot],
            [:agr,:min,:pip,:oil,:wst])

            CH4calc = DenseAxisArray([
            ## CH4 Emissions (MMt), 2019 / $US Billion (2017) Value Added inputs (kapital and labor, i.e. productive actiity)
                [CH4emiss[:EPAemiss,c]/sum(va_0[yr,:,c]) for c in CH4sectors];;
            ## subtract (maximum) mitigated CH4, so CH4 of remaining emissions after maximum abatement (at <$1000/t)
                [CH4emiss[:EPAemiss,c]*(1-CH4emiss[:MaxpercMit,c])/(sum(va_0[yr,:,c])+CH4emiss[:MitCostTot,c]*10^-3) for c in CH4sectors];;
            ## Total Cost of Mitigation is the standard VA inputs + the additional cost of the mitigation in US$Bill
                [CH4emiss[:MitCostTot,c]*10^-3+sum(va_0[yr,:,c]) for c in CH4sectors]],
            CH4sectors, # dimension 1 indexed by sector
            [:CH4Intens :CH4MitIntens :TotCostwMit ]) # dimension 2 indexed by values
            
            ## Relative cost of VA including max mitigation
            MitCostoverVA = DenseAxisArray([CH4calc[c,:TotCostwMit]/sum(va_0[yr,:,c]) for c in CH4sectors], CH4sectors)

            vam_0 = deepcopy(va_0) #copy for slack mitigating activties
            # # benchmark value added input levels for with additional % of costs for mitigation
            for c in CH4sectors
                vam_0[yr,:,c] = va_0[yr,:,c].data .* MitCostoverVA[c]
            end
            ## default of 0 for non-emitting and less significant emitting sectors
            ch4VASInt = DenseAxisArray(zeros(length(J)),J)
            ## Set vector of standard CH4 intensities: CH4 (in CO2eq)/value-added factor inputs
            for c in CH4sectors
                ch4VASInt[c] = CH4calc[c,:CH4Intens]
            end
            ## Still default of 0 for non-emitting and less significant emitting sectors
            ch4VAMInt = DenseAxisArray(zeros(length(J)),J)
            ## Set vector of CH4 *mitigated* intensities: Mitigated CH4 (in CO2eq)/value-added factor inputs
            for c in CH4sectors
                ch4VAMInt[c] = CH4calc[c,:CH4MitIntens]
            end
#####--------End single Mitigation step data set up-----------#####

CO2Int = DenseAxisArray(zeros(length(J)),J)
# 2024 GHG Inv Table 3-5: CO2 Emissions from Fossil Fuel Combustion by Fuel Type and Sector (MMT CO2 Eq.)
TotalCO2EmMMt_coal = 835.6 # 2020 895.9 # EPA Inventory CO2 Stationary Combustion - Coal sum(Electricity, Industrial, Commercial, & Residential=0)
# Option with no transport emissions in the model: assumption that direct
# TotalCO2EmMMt_gas_oil = 2104.80
Natural_gasCO2 = 1615.7 #2020 (incl 58.7 transport)
PetroleumCO2 = 1890.0 # 2020 (incl 1,514.2 transport)
# Option with all Transport CO2 emissions attributed to oil inputs: assumption that forms of oil fuel all transport that has direct CO2 emissions, and so taxing CO2 is total emissions from all oil as an input 
TotalCO2EmMMt_gas_oil = Natural_gasCO2 + PetroleumCO2 # EPA inventory all CO2 transport + CO2 Stationary Combustion - both Oil & Natural Gas sum(Electricity, Industrial, Commercial, Residential & oil U.S. Territories)Sta
# (MMtCO2eq/$Bill x 10-^3 -> Gt/$B=t/$) EPA Inventory 2022 sum of CO2 MMt for coal, and for gas & oil per Billion of total benchmark intermediate input from sector
CO2Int[:min] =  TotalCO2EmMMt_coal*10^-3/sum(id_0[yr,:min,:]) 
CO2Int[:oil] =   TotalCO2EmMMt_gas_oil*10^-3/sum(id_0[yr,:oil,:])  

TotCO2bnchmk =  TotalCO2EmMMt_coal*10^-3 + TotalCO2EmMMt_gas_oil*10^-3
# TotCH4bnchmk = sum(CH4emiss[:EPAemiss,:]) # from single step data version, same data, same value
TotCH4bnchmk = sum(MAC_CH4_WiNDC_tot[1,2:end])*10^-3
TotGHGbnchmk =  TotCO2bnchmk + TotCH4bnchmk
## End data preparations

## Set tax rates
CO2_taxrate = 190 # SC CO2 EPA 2023 SCGHG report, 2020 year, 2020US$, central 2% discount rate
CH4_taxrate = 190 # using SC CO2 because CH4 data is in MtCO2eq

MultiNat = MPSGEModel()

@parameters(MultiNat, begin
    ta[J], ta_0[yr,J]
    ty[J], ty_0[yr,J]
    tm[J], tm_0[yr,J]
    CH4_tax, 0.
    CO2_tax, 0.
end)

@sectors(MultiNat,begin
    Y[J],      (description = "Sectoral Production",)
    A[I],      (description = "Armington Supply",)
    VAS[J],    (description = "Value Added, standard")
    # VAM[J],    (description = "Value Added, with additional max mitigating activity") # For use to compared to the previous, single step set-up
    MS[M],     (description = "Margin Supply",)
    VAM5[CH4sectors],   (description = "Value Added, with mitigating activity up to \$5/t")
    VAM10[CH4sectors],  (description = "Value Added, with mitigating activity up to \$10/t")
    VAM15[CH4sectors],  (description = "Value Added, with mitigating activity up to \$15/t")
    VAM20[CH4sectors],  (description = "Value Added, with mitigating activity up to \$20/t")
    VAM30[CH4sectors],  (description = "Value Added, with mitigating activity up to \$30/t")
    VAM40[CH4sectors],  (description = "Value Added, with mitigating activity up to \$40/t")
    VAM50[CH4sectors],  (description = "Value Added, with mitigating activity up to \$50/t")
    VAM100[CH4sectors], (description = "Value Added, with mitigating activity up to \$100/t")
    VAM500[CH4sectors], (description = "Value Added, with mitigating activity up to \$500/t")
    VAM1000[CH4sectors],(description = "Value Added, with max mitigating activity, up to \$1000/t")
end)

@commodities(MultiNat,begin
    PA[I],   (description = "Armington Price")
    PY[J],   (description = "Supply",)
    PVA[VA], (description = "Value-added Input to VA blocks",)
    PVAM[J], (description = "Value-added output - Input to Y",)
    PM[M],   (description = "Margin Price",)
    PFX,     (description = "Foreign Exachange",)
end)

# Variables to track and report levels of CO2 emissions
@auxiliary(MultiNat, CO2em, index = [[:min, :oil]])
@auxiliary(MultiNat, CO2TotEm, description = "Total CO2 emissions from fossil fuels")
# Variables to track and report levels of CH4 emissions
@auxiliary(MultiNat, CH4em, index = [[:agr,:min,:oil,:pip,:wst]])
@auxiliary(MultiNat, CH4TotEm, description = "Total CH4 emissions")
@auxiliary(MultiNat, TotEm, description = "Total both emissions")


@consumer(MultiNat, RA, description = "Representative Agent")

# Domestic production for all sectors
for j∈J
    @production(MultiNat, Y[j], [t=0, s = 0], begin
        [@output(PY[i],ys_0[yr,j,i], t, taxes = [Tax(RA,ty[j])]) for i∈I]... 
        [@input(PA[i], id_0[yr,i,j], s, taxes = [Tax(RA,CO2_tax * CO2Int[i])]) for i∈I]...
         @input(PVAM[j], sum(va_0[yr,VA,j]), s)
    end)
end

# Total value added cost as a function labor (compen) and kapital (surplus), standard (no mitigation)
for j∈J
    @production(MultiNat, VAS[j], [t=0, s = 0, va => s = 1], begin
        [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
        [@input(PVA[va], va_0[yr,va,j], va, taxes = [Tax(RA,CH4_tax* ch4VASInt[j])]) for va∈VA]...
    end)
end

## Loop over all the Marginal Abatement tiers as Value-Added production blocks
VAMcommodSet = [VAM5,VAM10,VAM15,VAM20,VAM30,VAM40,VAM50,VAM100,VAM500,VAM1000]
for vam in VAMcommodSet
    for j∈CH4sectors
        if VAM_costover[vam.name,j]>1 # Some sectors are still cumulatively -negative costs at $5/t, so filtering those out.
            @production(MultiNat, vam[j], [t=0, s = 0, va => s = 1], begin
                [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
                [@input(PVA[va], va_0[yr,va,j]*VAM_costover[vam.name,j], va, taxes = [Tax(RA, CH4_tax*VAM_CH4EmInt[vam.name,j])]) for va∈VA]...
            end)
        end
    end
end

## Alternate structure, all EPA mitigation potential up to $1,000/t
# Slack mitigating VA activities for main CH4 producing sectors: total weighted Marginal Abatment version
# for j∈CH4sectors
#     @production(MultiNat, VAM[j], [t=0, s = 0, va => s = 1], begin
#         [@output(PVAM[j],sum(va_0[yr,:,j]), t)]... 
#         [@input(PVA[va], vam_0[yr,va,j], va, taxes = [Tax(RA, CH4_tax* ch4VAMInt[j])]) for va∈VA]...
#     end)
# end

for m∈M
    @production(MultiNat, MS[m], [t = 0, s = 0], begin
        [@output(PM[m], sum(ms_0[yr,i,m] for i∈I), t)]...
        [@input(PY[i], ms_0[yr,i,m], s) for i∈I]...
    end)
end

for i∈I
    @production(MultiNat, A[i], [t = 2, s = 0, dm => s = 2], begin
        [@output(PA[i], a_0[yr,i], t, taxes=[Tax(RA,ta[i])],reference_price=1-ta_0[yr,i])]...
        [@output(PFX, x_0[yr,i], t)]...
        [@input(PM[m], md_0[yr,m,i], s) for m∈M]...
        @input(PY[i], y_0[yr,i], dm)
        @input(PFX, m_0[yr,i], dm, taxes = [Tax(RA,tm[i])],reference_price=1+tm_0[yr,i])
    end)
end

@demand(MultiNat, RA, begin
    [@final_demand(PA[i], fd_0[yr,i,:pce]) for i∈I]...
    @endowment(PFX, bopdef_0[yr])
    [@endowment(PA[i], -sum(fd_0[yr,i,xfd] for xfd∈FD if xfd!=:pce)) for i∈I]...
    [@endowment(PVA[va], sum(va_0[yr,va,j] for j∈J)) for va∈VA]...
end, elasticity = 1)

## CO2 emissions for fossil fuel sectors are the activity levels times the (base) total emissions intensity 
@aux_constraint(MultiNat, CO2em[:min],  CO2em[:min] - Y[:min]*TotalCO2EmMMt_coal*10^-3)
@aux_constraint(MultiNat, CO2em[:oil],  CO2em[:oil] - Y[:oil]*TotalCO2EmMMt_gas_oil*10^-3)
## Total CO2 emissions are the sum of emissions from the 2 fossil fuel sectors (constraint expressed as equantion = 0)
@aux_constraint(MultiNat, CO2TotEm, CO2TotEm - (CO2em[:min] + CO2em[:oil]))
## CH4 emissions for each CH4 emitting sector are the sum of (either): VA Standard activity levels x standard CH4 emissions intensity (benchmark = 1 x base emissions)
## +/or VA Mitigating activities at each tier, activity level x base emissions x mitigated emissions factor (mitigated intensity/baseline intensity)
for c in CH4sectors
    # VAM, the old one, for testing.
    # @aux_constraint(MultiNat, CH4em[c],  CH4em[c] - (VAS[c]*CH4emiss[:EPAemiss,c]+VAM[c]*CH4emiss[:EPAemiss,c]*ch4VAMInt[c]/ch4VASInt[c]))
    @aux_constraint(MultiNat, CH4em[c],  CH4em[c] - (VAS[c]*CH4emiss[:EPAemiss,c] +#  )) # +
    ifelse(VAM_costover[:VAM5,c]>1, VAM5[c]*CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM5,c]/ch4VASInt[c] , 0) + # The $5/t tier includes sectors with negative costs, so have to filter those out here (AS WELL as in the production block)
    VAM10[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM10,c]/ch4VASInt[c]+
    VAM15[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM15,c]/ch4VASInt[c]+
    VAM20[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM20,c]/ch4VASInt[c]+
    VAM30[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM30,c]/ch4VASInt[c]+
    VAM40[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM40,c]/ch4VASInt[c]+
    VAM50[c]  *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM50,c]/ch4VASInt[c]+
    VAM100[c] *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM100,c]/ch4VASInt[c]+
    VAM500[c] *CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM500,c]/ch4VASInt[c]+
    VAM1000[c]*CH4emiss[:EPAemiss,c]*VAM_CH4EmInt[:VAM1000,c]/ch4VASInt[c]
    ))
end
 
## Total CH4 Emissions are the sum of emissions from CH4 emitting sectors
@aux_constraint(MultiNat, CH4TotEm, CH4TotEm - (CH4em[:agr] + CH4em[:min] + CH4em[:oil] + CH4em[:pip] + CH4em[:wst] ))
## Total GHG (CO2 & CH4) emissions in Mt CO2eq
@aux_constraint(MultiNat, TotEm, TotEm - (CH4TotEm + CO2TotEm))
set_silent(MultiNat)

# Benchmark 
# fix(RA, sum(fd_0[yr,I,:pce])
## Note: Benchmark doesn't solve at 0 interation because of margins of slack activity. Does balance with interactions or slack vars and production commented out.
solve!(MultiNat)
#; cumulative_iteration_limit = 0)

fullvrbnch = generate_report(MultiNat);
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
# print(sort(fullvrbnch, :bmkmarg, by= abs))#, rev=true))
# fullvrbnch[!,:var] = Symbol.(fullvrbnch[:,:var])

# Initialize a Dataframe to save final demand results
FDemand = DataFrame(index=Vector{Symbol}(undef, length(I)),
desc=Vector{Symbol}(undef, length(I)), 
bnch=Vector{Float64}(undef, length(I)), 
WiNcntfact=Vector{Float64}(undef, length(I)), 
ch4tax=Vector{Float64}(undef, length(I)),
cO2tax=Vector{Float64}(undef, length(I)),
both=Vector{Float64}(undef, length(I)),
ch4Qdelta=Vector{Float64}(undef, length(I)),
CO2Qdelta=Vector{Float64}(undef, length(I)),
bothQdelta=Vector{Float64}(undef, length(I)),
bncQpc=Vector{Float64}(undef, length(I)),
ch4Qpc=Vector{Float64}(undef, length(I)),
CO2Qpc=Vector{Float64}(undef, length(I)),
bothQpc=Vector{Float64}(undef, length(I))
)
for (n,i) in enumerate(I)
    FDemand[n,:index]= i
    FDemand[n,:desc] = Symbol(Sectors[Sectors.index.==string(i),2][1])
    FDemand[n,:bnch] = value(demand(RA,PA[i]))
end
# unfix(RA)
## Check against WiNDC standard counterfactual
set_value!(ta,0)
set_value!(tm,0)

solve!(MultiNat)

fullvrWiNcntfact = generate_report(MultiNat)
rename!(fullvrWiNcntfact, :value => :WiNcntfact, :margin => :WiNcntfactmarg)

for (n,i) in enumerate(I)
    FDemand[n,:WiNcntfact] = value(demand(RA,PA[i]))
end

## DataFrame to hold the industry indices with descriptions
Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

## Re-set to benchmark taxes
set_value!(ta,ta_0[yr,J])
set_value!(tm,tm_0[yr,J])

# tax are at $/t of CH4(CO2eq)
## "EPA SC CH4 is $1600/t. Possibly set as CO2eq * CO2 conversion rate back to per ton of CH4? But 1600/190 is 8.42, not 29.8 conversion...  
## Or alternatively, re-work data to calculate CH4 in actual tons"
set_value!(CH4_tax, CH4_taxrate)
# Set CO2 tax to 0 for running separately.
set_value!(CO2_tax,0.)
solve!(MultiNat)

for (n,i) in enumerate(I)
    FDemand[n,:ch4tax] = value(demand(RA,PA[i]))
end

Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Sorted, to report highest and lowest 4 output activity levels for this policy simulation
Rs[:,2][1:4]
Rs[:,2][68:71]

fullvrch4 = generate_report(MultiNat)
rename!(fullvrch4, :value => :ch4tax, :margin => :ch4marg)

# print(sort(df, :margin, by= abs, rev=true))
# print(df)

# # Counterfactual Fossil fuel extraction is ALSO taxed at emissions intensitiy of input x tax in $/ton
set_value!(CO2_tax, CO2_taxrate)

solve!(MultiNat, cumulative_iteration_limit=10000) #;

for (n,i) in enumerate(I)
    FDemand[n,:both] = value(demand(RA,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:5]
# Rs[:,2][68:71]

fullvrboth = generate_report(MultiNat)
rename!(fullvrboth, :value => :bothtaxes, :margin => :bothmarg)

## Then, set CH4 taxes back to 0 to generate CO2 tax ONLY
    set_value!(CH4_tax, 0.0)

solve!(MultiNat, cumulative_iteration_limit=10000) #;

for (n,i) in enumerate(I)
    FDemand[n,:cO2tax] = value(demand(RA,PA[i]))
end

# Rs = DataFrame([Y value.(Y) last.(first.(string.(Y),6),3)][sortperm([Y value.(Y)][:,2], rev= true), :], [:var, :val, :index])
# Rs = innerjoin(Sectors[:,[1,2]], Rs[:,[2,3]], on = :index)
# Rs[:,2][1:4]
# Rs[:,2][68:71]
fullvrCO2 = generate_report(MultiNat)
rename!(fullvrCO2, :value => :CO2tax, :margin => :CO2marg)

#Generate Dataframe with all results (including names expressions)
FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrCO2, fullvrboth, fullvrWiNcntfact, on = [:var], makeunique=true);
Compare = FullResults[1:end,[1,2,4,6,8,10]];
## Sum the difference of each tax applied individually
Compare.sum = Compare.ch4tax .- 1 + Compare.CO2tax .-1

## sum difference between the sum of the individual taxes to the combined taxes
Compare.diff = Vector{Union{Missing, Float64}}(undef, length(Compare[:,1]))
for n in 1:length(Compare[:,1])
    if Compare.ch4tax[n] == 0 || Compare.CO2tax[n] == 0
        Compare.diff[n] = missing
    elseif Compare.bothtaxes[n] >= 1 &&  Compare.sum[n] >= 0
        Compare.diff[n] = (Compare.bothtaxes[n] - 1) - Compare.sum[n]
    elseif Compare.bothtaxes[n] >= 1 &&  Compare.sum[n] <= 0
        Compare.diff[n] = 1 - Compare.bothtaxes[n] + Compare.sum[n]
    elseif Compare.bothtaxes[n] <= 1 &&  Compare.sum[n] <= 0
        Compare.diff[n] = 1 - Compare.bothtaxes[n] + Compare.sum[n]
    elseif Compare.bothtaxes[n] <= 1 &&  Compare.sum[n] >= 0
        Compare.diff[n] = 1 - Compare.bothtaxes[n] - Compare.sum[n]
    end
end

Compare[!,:var] = Symbol.(Compare[:,:var]);
# print(sort!(Compare, :var));
# println(sort!(Compare, :diff, by = abs, rev=true))#[1:25,:])
# CSV.write("C:\\Users\\Eli\\Box\\CGE\\MPSGE-JL\\First Mulit GHG taxes Paper\\MultiResults$(Dates.format(now(),"yyyy-mm-d_HhM")).csv", Compare, missingstring="missing", bom=true)
compCO2em = filter(:var => ==(:CO2TotEm), Compare);
compCH4em = filter(:var => ==(:CH4TotEm), Compare);
compTotem = filter(:var => ==(:TotEm), Compare);
EmissionReductionResults = DataFrame(
["CO2" "Gt" TotCO2bnchmk TotCO2bnchmk - only(compCO2em[:,:CO2tax]) TotCO2bnchmk - only(compCO2em[:,:ch4tax]) TotCO2bnchmk - only(compCO2em[:,:ch4tax]) + TotCO2bnchmk - only(compCO2em[:,:CO2tax]) TotCO2bnchmk - only(compCO2em[:,:bothtaxes])  TotCO2bnchmk - only(compCO2em[:,:ch4tax]) + TotCO2bnchmk - only(compCO2em[:,:CO2tax]) - (TotCO2bnchmk - only(compCO2em[:,:bothtaxes])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"CH4" "GtCO2eq" TotCH4bnchmk TotCH4bnchmk - only(compCH4em[:,:CO2tax]) TotCH4bnchmk - only(compCH4em[:,:ch4tax]) TotCH4bnchmk - only(compCH4em[:,:ch4tax]) + TotCH4bnchmk - only(compCH4em[:,:CO2tax]) TotCH4bnchmk - only(compCH4em[:,:bothtaxes]) TotCH4bnchmk - only(compCH4em[:,:ch4tax]) + TotCH4bnchmk - only(compCH4em[:,:CO2tax]) - (TotCH4bnchmk - only(compCH4em[:,:bothtaxes])); # Interactions = the amount not reduced with taxes combined, compared to expections from individual taxes
"GHGs" "GtCO2eq" TotGHGbnchmk TotGHGbnchmk - only(compTotem[:,:CO2tax]) TotGHGbnchmk - only(compTotem[:,:ch4tax]) TotGHGbnchmk - only(compTotem[:,:ch4tax]) + TotGHGbnchmk - only(compTotem[:,:CO2tax]) TotGHGbnchmk - only(compTotem[:,:bothtaxes]) TotGHGbnchmk - only(compTotem[:,:ch4tax]) + TotGHGbnchmk - only(compTotem[:,:CO2tax]) - (TotGHGbnchmk - only(compTotem[:,:bothtaxes]))], ["Emissions", "Unit", "Bnchmrk_Emissions", "CO2tax_reduc", "CH4tax_reduc","Sum_of_both_taxes", "taxes_combined" ,"Interactions"])

## Generate subset DataFrame with just the Value-Added activity for the emitting sectors, show those results
#filter(row -> row.var ∈ [Symbol("VAM[agr]"),Symbol("VAM[min]"),Symbol("VAM[pip]"),Symbol("VAS[oil]"),Symbol("VAM[oil]"),Symbol("VAS[min]"),Symbol("VAS[pip]"),Symbol("VAS[agr]"),Symbol("VAS[wst]"),Symbol("VAM[wst]"),], Compare)

function plottaxemisscurve(tax1, tax2, start, interval, finish, RAval, isfixed, cnst=1)
    """# runs a loop increasing each tax by \$1/t and then plotting Total GHG (CO2 & CH4) **incorporated** emissions 
    # Arguments are: which tax to change, other tax to either change simultaneously OR keep at 0, st=initial \$ tax value, fin= final \$ tax value,
    # and final (optional) argument can be set to 0 to remove other tax, by default """
    margemiss = DataFrame(tax=Float64[], Emissions=Float64[], CH4Emissions=Float64[],CO2Emissions=Float64[])
    Testvars = DataFrame(taxrt=Float64[], 
    # Yagr=Float64[],Ymin=Float64[],Ypip=Float64[],Yoil=Float64[],Apip=Float64[],Aoil=Float64[],Ywst=Float64[],
    # CompDYPApip=Float64[],CompDApipPApip=Float64[],DemRAPApip=Float64[],
    # PAagr=Float64[],PAmin=Float64[],PApip=Float64[],PAoil=Float64[],PAwst=Float64[],PAuti=Float64[],compDApipPAoil=Float64[],compdDAoilPAoil=Float64[],
    # VASagr=Float64[],VAMagr=Float64[],VASmin=Float64[],VAMmin=Float64[],VASpip=Float64[],VAMpip=Float64[],VASoil=Float64[],VAMoil=Float64[],VASwst=Float64[],VAMwst=Float64[],
    VASagr=Float64[],VAM10agr=Float64[],VAM100agr=Float64[],VAM500agr=Float64[],VAMkagr=Float64[],VASmin=Float64[],VAM10min=Float64[],VAM100min=Float64[],VAM500min=Float64[],VAMkmin=Float64[],VASpip=Float64[],VAM10pip=Float64[],VAM100pip=Float64[],VAM500pip=Float64[],VAMkpip=Float64[],VASoil=Float64[],VAM10oil=Float64[],VAM100oil=Float64[],VAM500oil=Float64[],VAMkoil=Float64[],VASwst=Float64[],VAM10wst=Float64[],VAM100wst=Float64[],VAM500wst=Float64[],VAMkwst=Float64[],
    TotEm=Float64[],CH4TotEm=Float64[],CO2TotEm=Float64[]
    # CH4emoil=Float64[],CH4empip=Float64[],CO2emin=Float64[],CO2emoil=Float64[]
    )
    ResultsTroubleshoot = DataFrame(var=[], value=Float64[], margin=Float64[], x1=Float64[]) 
    for i in start:interval:finish
        set_value!(tax1, i)
        set_value!(tax2, cnst*i)
        solve!(MultiNat, output="no");
        Results = generate_report(MultiNat)
        Results[!,:var] = Symbol.(Results[:,:var]);
        ResultsTroubleshoot =vcat(ResultsTroubleshoot, [Results fill(i,length(Results[:,1]))])
        push!(margemiss, [i only(filter(:var => ==(:TotEm), Results)[:, :value]) only(filter(:var => ==(:CH4TotEm), Results)[:, :value]) only(filter(:var => ==(:CO2TotEm), Results)[:, :value])])
        push!(Testvars, [i,                
        # value(Y[:agr]),value(Y[:min]),value(Y[:pip]),value(Y[:oil]),value(A[:pip]),value(A[:oil]),value(Y[:wst]),
        # value(compensated_demand(MultiNat[:Y][:pip],MultiNat[:PA][:pip])),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:pip])),value(demand(MultiNat[:RA],MultiNat[:PA][:pip])),
        # value(PA[:agr]),value(PA[:min]),value(PA[:pip]),value(PA[:oil]),value(PA[:wst]),value(PA[:uti]),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:oil])),value(compensated_demand(MultiNat[:A][:oil],MultiNat[:PA][:oil])),
        # value(VAS[:agr]),value(VAM[:agr]),value(VAS[:min]),value(VAM[:min]),value(VAS[:pip]),value(VAM[:pip]),value(VAS[:oil]),value(VAM[:oil]),value(VAS[:wst]),value(VAM[:wst]),
        value(VAS[:agr]),value(VAM10[:agr]),value(VAM100[:agr]),value(VAM500[:agr]),value(VAM1000[:agr]),value(VAS[:min]),value(VAM10[:min]),value(VAM100[:min]),value(VAM500[:min]),value(VAM1000[:min]),value(VAS[:pip]),value(VAM10[:pip]),value(VAM100[:pip]),value(VAM500[:pip]),value(VAM1000[:pip]),value(VAS[:oil]),value(VAM10[:oil]),value(VAM100[:oil]),value(VAM500[:oil]),value(VAM1000[:oil]),value(VAS[:wst]),value(VAM10[:wst]),value(VAM100[:wst]),value(VAM500[:wst]),value(VAM1000[:wst]),
        value(TotEm),value(CH4TotEm),value(CO2TotEm)
        # value(CH4em[:oil]),value(CH4em[:pip]),value(CO2em[:min]),value(CO2em[:oil])
        ] 
            )
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    plt = plot(margemiss[!,:tax], margemiss[!,:Emissions].*10^3, title= "RA:\$$RAval fxd:$isfixed", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,1600))
    plt = plot!([190], seriestype=:vline, label="SCCO2", ylim=(0,TotGHGbnchmk*10^3))
    plch4 = plot(margemiss[!,:tax], margemiss[!,:CH4Emissions].*10^3, label=false, title="CH4 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,1600))
    plch4 = plot!([190], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
    plco2 = plot(margemiss[!,:tax], margemiss[!,:CO2Emissions].*10^3, label=false, title="CO2 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,1600))
    plco2 = plot!([190], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
    plt3 = plot(margemiss[!,:tax], margemiss[!,:Emissions].*10^3, title= "Emissions with $tax1 $tax2in", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="\$/t", xlims=(0,1600), color=:black, linewidth=3)
    plt3 = plot!(margemiss[!,:tax], margemiss[!,:CH4Emissions].*10^3, label="CH4 Emissions", color=:green)
    plt3 = plot!(margemiss[!,:tax], margemiss[!,:CO2Emissions].*10^3, label="CO2 Emissions", color=:blue)
    plt3 = plot!([190], seriestype=:vline, label="SCCO2", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
    # plcdyp = plot(margemiss[!,:tax],Testvars[!,:CompDYPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(Y:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDYPApip]),maximum(Testvars[!,:CompDYPApip])), xlabel="$tax1 $tax2in \$/t")
    # plcdap = plot(margemiss[!,:tax],Testvars[!,:CompDApipPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(A:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDApipPApip]),maximum(Testvars[!,:CompDApipPApip])), xlabel="$tax1 $tax2in \$/t")
    # plfdrap = plot(margemiss[!,:tax],Testvars[!,:DemRAPApip], title= "RA:\$$RAval fxd:$isfixed", label="final_dem(RA,PA:pip)", ylim=(minimum(Testvars[!,:DemRAPApip]),maximum(Testvars[!,:DemRAPApip])), xlabel="$tax1 $tax2in \$/t")
    # pla = plot(margemiss[!,:tax],Testvars[!,:Yagr], label=false, title="Y:agr", ylim=(minimum(Testvars[!,:Yagr]),maximum(Testvars[!,:Yagr])), xlabel="$tax1 $tax2in \$/t")
    # plm = plot(margemiss[!,:tax],Testvars[!,:Ymin], title= "RA:\$$RAval fxd:$isfixed", label="Y:min", ylim=(minimum(Testvars[!,:Ymin]),maximum(Testvars[!,:Ymin])), xlabel="$tax1 $tax2in \$/t")
    # plp = plot(margemiss[!,:tax],Testvars[!,:Ypip], title= "RA:\$$RAval fxd:$isfixed", label="Y:pip", ylim=(minimum(Testvars[!,:Ypip]),maximum(Testvars[!,:Ypip])), xlabel="$tax1 $tax2in \$/t")
    # plo = plot(margemiss[!,:tax],Testvars[!,:Yoil], title= "RA:\$$RAval fxd:$isfixed", label="Y:oil", ylim=(minimum(Testvars[!,:Yoil]),maximum(Testvars[!,:Yoil])), xlabel="$tax1 $tax2in \$/t")
    # plw = plot(margemiss[!,:tax],Testvars[!,:Ywst], title= "RA:\$$RAval fxd:$isfixed", label="Y:wst", ylim=(minimum(Testvars[!,:Ywst]),maximum(Testvars[!,:Ywst])), xlabel="$tax1 $tax2in \$/t")
    # plpa = plot(margemiss[!,:tax],Testvars[!,:PAagr], title= "RA:\$$RAval fxd:$isfixed", label="PA:agr", ylim=(minimum(Testvars[!,:PAagr]),maximum(Testvars[!,:PAagr])), xlabel="$tax1 $tax2in \$/t")
    # plpm = plot(margemiss[!,:tax],Testvars[!,:PAmin], title= "RA:\$$RAval fxd:$isfixed", label="PA:min", ylim=(minimum(Testvars[!,:PAmin]),maximum(Testvars[!,:PAmin])), xlabel="$tax1 $tax2in \$/t")
    # plpp = plot(margemiss[!,:tax],Testvars[!,:PApip], title= "RA:\$$RAval fxd:$isfixed", label="PA:pip", ylim=(minimum(Testvars[!,:PApip]),maximum(Testvars[!,:PApip])), xlabel="$tax1 $tax2in \$/t")
    # plpo = plot(margemiss[!,:tax],Testvars[!,:PAoil], title= "RA:\$$RAval fxd:$isfixed", label="PA:oil", ylim=(minimum(Testvars[!,:PAoil]),maximum(Testvars[!,:PAoil])), xlabel="$tax1 $tax2in \$/t")
    # plpw = plot(margemiss[!,:tax],Testvars[!,:PAwst], title= "RA:\$$RAval fxd:$isfixed", label="PA:wst", ylim=(minimum(Testvars[!,:PAwst]),maximum(Testvars[!,:PAwst])), xlabel="$tax1 $tax2in \$/t")
    # plAp = plot(margemiss[!,:tax],Testvars[!,:Apip], title= "RA:\$$RAval fxd:$isfixed", label="A:pip", ylim=(minimum(Testvars[!,:Apip]),maximum(Testvars[!,:Apip])), xlabel="$tax1 $tax2in \$/t")
    # plAo = plot(margemiss[!,:tax],Testvars[!,:Aoil], title= "RA:\$$RAval fxd:$isfixed", label="A:oil", ylim=(minimum(Testvars[!,:Aoil]),maximum(Testvars[!,:Aoil])), xlabel="$tax1 $tax2in \$/t")     # Or label=false, title="price of oil commodity"
    return margemiss, plt, Testvars, ResultsTroubleshoot, plch4, plco2, plt3 #,pla, plm, plp, plo, plw, plpa, plpm, plpp, plpo, plpw, 
    # plcdyp, plcdap, plfdrap, plAp, plAo
end

## Add Quant diff and % of consumption columns for Final Demand report dataframe
FDemand[:,:ch4Qdelta]=FDemand[:,:ch4tax].-FDemand[:,:bnch]
FDemand[:,:CO2Qdelta]=FDemand[:,:cO2tax].-FDemand[:,:bnch]
FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
FDemand[:,:ch4Qpc]=FDemand[:,:ch4tax]./sum(FDemand[:,:ch4tax])*100
FDemand[:,:CO2Qpc]=FDemand[:,:cO2tax]./sum(FDemand[:,:cO2tax])*100
FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

set_silent(MultiNat)
fix(RA,16426.2) # RA value at $190/t
# checkch4CO2 = plottaxemisscurve(CH4_tax, CO2_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# checkch4CO2[7]
# png(checkch4CO2[7], "./Results/Bothtax-Allemiss")
# checkch4CO2[2]
# print("checkch4CO2",checkch4CO2[3])

# fix(RA, sum(fd_0[yr,I,:pce]))
set_upper_bound(MultiNat[:A][:pip], 10)
# # MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pip]))
# # fix(RA, 14008.668551652801)
# # unfix(RA)
# fix(RA,16030.7) # RA value at $190/t
# # checkCO2 = plottaxemisscurve(CO2_tax, CH4_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2 = plottaxemisscurve(CO2_tax, CH4_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkCO2[2] # Total emissions
# # # println("PApip up lim =", upper_bound(MultiNat[:A][:pip]))

# checkCO2[5] # CH4 Emissions, was Y agr
# checkCO2[6] # CO2 emissions, was Y min
# checkCO2[7]
# png(checkCO2[7], "./Results/CO2tax-Allemiss")

# checkCO2[7] # Y pip
# checkCO2[8] # Y oil
# checkCO2[9] # Y wst
# checkCO2[10] # PA agr
# checkCO2[11] # PA min
# checkCO2[12] # PA pip
# checkCO2[13] # PA oil
# checkCO2[14] # PA wst
# checkCO2[15] # ch4 emissions
# checkCO2[16] # co2 emissions
# checkCO2[17] # comp demand Y pip
# checkCO2[18] # comp demand A pip
# checkCO2[19] # Final demand pip
# checkCO2[20] # A pip
# checkCO2[21] # A oil

# # print("checkCO2",checkCO2[3])
# # # checkCO2[4]
# fix(RA,16426.2) # RA value at $190/t
# checkch4CO2 = plottaxemisscurve(CH4_tax, CO2_tax, 0, 1, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# # checkch4CO2[2]

# fix(RA,fix(RA,16030.7) # RA value at $190/t
# checkch4 = plottaxemisscurve(CH4_tax,CO2_tax, 0, 10, 1600, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# checkch4[2]
# checkch4[5]
# checkch4[6]
# checkch4[7]
# png(checkch4[7], "./Results/CH4tax-Allemiss")

# checkch4[8]
# checkch4[9]
# checkch4[10]
# checkch4[11]
# checkch4[12]
# checkch4[13]
# checkch4[14]
# checkch4[15] #ch3 emissions
# checkch4[16] # Co2 emissions
# print("checkch4",checkch4[3])

# png(checkch4[2], "./Results/CH4to500-Totemiss")
# png(checkch4[15], "./Results/CH4to500-ch4emiss")
# png(checkch4[16], "./Results/CH4to500-co2emiss")
# png(checkch4[8], "./Results/CH4to500-oilactivity")
# png(checkch4[13], "./Results/CH4to500-oilprice")

# print("checkch4",checkch4[3])
# print(check[1])

# CH4sectors = [:agr,:min,:pip,:oil,:wst] #:uti? # subset index for relevant CH4 mitigation sectors (VA slack in benchmark)

# png(checkCO2[2], "./Results/CO2to1600RAUnfxd")
println(length(Sectors[:,1]))
EmissUnits_mt = DataFrame();
EmissUnits_mt.Unit=["Mt"; "MtCO2eq"; "MtCO2eq"];
EmissionReductionResults_Mt =[EmissionReductionResults[:,1:1] EmissUnits_mt EmissionReductionResults[:,3:end].*10^3] 
