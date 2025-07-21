"""
Function to calculate CES ustility. Currently set specifically to "FDem" sector and "PA" indexed commodity
Also very slow bc dictionaries.
"""
function ces_node(
    bnch_quants::Dict{Symbol,Float64}, # hold the name of the inputs and benchmark data quantities (can be leaves or nests, with the total benchmark - data for leaves, or from nest returns
    quantities::Dict{Symbol,Float64},#; # hold the name of the inputs and utilities (can be leaves or nests, with the total benchmark - data for leaves, or from nest returns
)   
    nests = []
    for i in keys(quantities)
        push!(nests, ch_pr[i]) # Dict that has children as keys and parents as values
    end
    if !all( ==(nests[1]), nests)
        error("Error: not the same nest")
    end
    nest = nests[1] # this is the parent, the one they're going into (just one, to use to get elasticity)
    σ = elasticity(FDem,nest)
    # ρ = 1 - 1 / σ   # calculate once, so \rho is a constant
    total_benchmark = sum(values(bnch_quants)) # values from the benchmark data, used for shares
    shares = Dict{Symbol, Float64}(q => bnch_quants[q] / total_benchmark for q in keys(quantities)) 
    utility =  value(FDem) * sum(shares[i] * quantities[i]^((σ-1)/σ) for i in keys(quantities) if quantities[i]>0)^(σ/(σ-1)) #value(FDem) is the activity level of the utility 'sector' quantities are the comp_dem [leaves], or the utility of input composite
    # utility =  value(FDem) * sum(shares[i]^σ * quantities[i]^(1-σ) for i in keys(quantities) if quantities[i]>0)^(1 / (σ-1)) #value(FDem) is the activity level of the utility 'sector' quantities are the comp_dem [leaves], or the utility of input composite

    return (
        Dict{Symbol, Float64}(nest => total_benchmark),
        Dict{Symbol, Float64}(nest => utility),)
end

function CESutility_multinat()

### For leaves
leaf_ugs = Dict(:ugs => pce_0[:ugs,:pce])
leaf_coa = Dict(:coa => pce_0[:coa,:pce])
# Leaves [:ugs,:coa] into :homefuels
homefuel_i,homefuel_ute = ces_node(Dict{Symbol, Float64}(:ugs => bnchdata[:ugs,:pce],:coa => bnchdata[:coa,:pce]), 
    Dict{Symbol, Float64}(:coa=>value(compensated_demand(FDem,PA[:coa])),:ugs=>value(compensated_demand(FDem,PA[:ugs]))))
# homefuel_i,homefuel_ute = ces_node(merge(leaf_coa,leaf_ugs), merge(leaf_coa,leaf_ugs))
# bnch_quants = merge(leaf_coa,leaf_ugs)
# quantities =   Dict{Symbol, Float64}([x=>value(compensated_demand(FDem,PA[x])) for x in [:coa,:ugs]])

# Leaves [:rnw, :uel] into :elect
elect_i,elect_ute = ces_node(Dict{Symbol, Float64}(:rnw_elect => bnchdata[:rnw,:pce] * (1-veh_chrg_pc),:uel_elect => bnchdata[:uel,:pce] * (1-veh_chrg_pc)),
    Dict{Symbol, Float64}(:rnw_elect=>value(compensated_demand(FDem,PA[:rnw],:elect)),:uel_elect=>value(compensated_demand(FDem,PA[:uel],:elect))))
# elect_i,elect_ute = ces_node(merge(leaf_rnw,leaf_uel), merge(leaf_rnw,leaf_uel))
# bnch_quants = merge(leaf_rnw,leaf_uel), Dict{Symbol, Float64}(Electq_names[x]=>value(compensated_demand(FDem,PA[x],:elect)) for x in [:rnw, :uel])
# quantities = 

# Nests [:homefuels, :elect] into [:home_nrg_expd]
home_nrg_expd_i,home_nrg_expd_ute = ces_node(merge(homefuel_i,elect_i),merge(homefuel_ute,elect_ute))
# bnch_quants = merge(homefuel_i,elect_i),merge(homefuel_ute,elect_ute))
# quantities = 

# # Leaves [:hou,:ore] into :own_occ_exp
own_occ_exp_i, own_occ_exp_ute= ces_node(Dict{Symbol, Float64}(:hou => bnchdata[:hou,:pce],:ore => bnchdata[:ore,:pce]),  
    Dict{Symbol, Float64}([x=>value(compensated_demand(FDem,PA[x])) for x in [:hou,:ore]]))
# own_occ_exp_i, own_occ_exp_ute= ces_node(merge(leaf_hou,leaf_ore), merge(leaf_hou,leaf_ore))
# bnch_quants = merge(leaf_hou,leaf_ore)
# quantities = Dict{Symbol, Float64}([x=>value(compensated_demand(FDem,PA[x])) for x in [:hou,:ore]])


# Nests [:own_occ_exp, :home_nrg_expd] into [:housing_exp]
housing_exp_i, housing_exp_ute = ces_node(merge(home_nrg_expd_i,own_occ_exp_i),merge(home_nrg_expd_ute,own_occ_exp_ute))
# bnch_quants = merge(home_nrg_expd_i,own_occ_exp_i)
# quantities = merge(home_nrg_expd_ute,own_occ_exp_ute)

# leaves_goods = Dict{Symbol, Float64}(x => bnchdata[x,:pce] for x in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa :oil ]),Ip))
goods_i, goods_ute = ces_node(Dict{Symbol, Float64}(x => bnchdata[x,:pce] for x in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa :oil ]),Ip)),
    Dict{Symbol, Float64}(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa :oil]),Ip)))
# goods_i, goods_ute = ces_node(leaves_goods,leaves_goods)
# bnch_quants = leaves_goods
# quantities = Dict{Symbol, Float64}(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa ]),Ip))

# Nest :goods , :housing_exp into :non_transp, 
non_transp_i, non_transp_ute = ces_node(merge(housing_exp_i,goods_i),merge(housing_exp_ute,goods_ute))
# bnch_quants = merge(housing_exp_i,goods_i)
# quantities = merge(housing_exp_ute,goods_ute)

# # # Leaves[[:rnw, :uel]] into [:veh_elect]
# leaf_rnw_veh = Dict{Symbol, Float64}(:rnw_veh_elect => bnchdata[:rnw,:pce] * veh_chrg_pc)
# leaf_uel_veh = Dict{Symbol, Float64}(:uel_veh_elect => bnchdata[:uel,:pce] * veh_chrg_pc)
veh_elect_i,veh_elect_ute = ces_node(Dict{Symbol, Float64}(:rnw_veh_elect => bnchdata[:rnw,:pce] * veh_chrg_pc, :uel_veh_elect => bnchdata[:uel,:pce] * veh_chrg_pc), 
    Dict{Symbol, Float64}(:rnw_veh_elect=>value(compensated_demand(FDem,PA[:rnw],:veh_elect)),:uel_veh_elect=>value(compensated_demand(FDem,PA[:uel],:veh_elect))))
# bnch_quants = merge(leaf_rnw_veh,leaf_uel_veh)
# quantities = Dict{Symbol, Float64}(Veh_electq_names[x]=>value(compensated_demand(FDem,PA[x],:veh_elect)) for x in [:rnw, :uel])

# # # Leaves[:pet, :oil] into :petr
# leaf_pet = Dict{Symbol, Float64}(:pet => bnchdata[:pet,:pce])
# leaf_oil = Dict{Symbol, Float64}(:oil => bnchdata[:oil,:pce])
petr_i, petr_ute= ces_node(Dict{Symbol, Float64}(:pet => bnchdata[:pet,:pce], :oil => bnchdata[:oil,:pce]),  
    Dict{Symbol, Float64}([x=>value(compensated_demand(FDem,PA[x])) for x in [:pet,:oil]]))
# petr_i, petr_ute= ces_node(merge(leaf_pet,leaf_oil), merge(leaf_pet,leaf_oil))
# bnch_quants = merge(leaf_pet,leaf_oil)
# quantities = Dict{Symbol, Float64}([x=>value(compensated_demand(FDem,PA[x])) for x in [:pet,:oil]])

# Nest :petr, :veh_elect into :fuels, 
fuels_i, fuels_ute = ces_node(merge(veh_elect_i,petr_i),merge(veh_elect_ute,petr_ute))

# Leaves [:mot :ote :mvt] =>veh_exp
veh_exp_i, veh_exp_ute = ces_node(Dict{Symbol, Float64}(x => bnchdata[x,:pce] for x in [:mot, :ote, :mvt]),
    Dict{Symbol, Float64}(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in [:mot :ote :mvt]))
# veh_exp_i, veh_exp_ute = ces_node(leaves_vehexp,leaves_vehexp)
# bnch_quants = leaves_vehexp
# quantities = Dict{Symbol, Float64}(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in [:mot :ote :mvt])

# Nest fuels, veh_exp=>pers_transp=0.25,                     
pers_transp_i, pers_transp_ute = ces_node(merge(fuels_i, veh_exp_i),merge(fuels_ute, veh_exp_ute))

# Leaves [:air :trn :wtt :trk :grd :otr] => non_pers_transp
non_pers_transp_i, non_pers_transp_ute = ces_node(Dict{Symbol, Float64}(x => bnchdata[x,:pce] for x in [:air,:trn,:wtt,:trk,:grd,:otr]),
    Dict{Symbol, Float64}(x=>value(compensated_demand(FDem,PA[x])) for x in [:air,:trn,:wtt,:trk,:grd,:otr]))
# bnch_quants = leaves_np_trans
# quantities = Dict{Symbol, Float64}(x=>value(compensated_demand(FDem,PA[x])) for x in [:air,:trn,:wtt,:trk,:grd,:otr])

# Nest non_pers_transp=>transp=0.25,                  pers_transp=>transp=0.4, 
transp_i, transp_ute = ces_node(merge(non_pers_transp_i, pers_transp_i),merge(non_pers_transp_ute, pers_transp_ute))

# Nest transp, non_transp=>s,
s_i, s_ute = ces_node(merge(transp_i, non_transp_i),merge(transp_ute, non_transp_ute))

    return only(values(s_ute))
end

# # Edited from elastiticity vector
# ### Keys are the ends, and values are the nodes above
ch_pr = Dict{Symbol, Symbol}(
:transp => :s,  :non_transp => :s,
                            
:non_pers_transp => :transp,       :pers_transp => :transp, 
:veh_exp => :pers_transp,           :fuels => :pers_transp,
:air => :non_pers_transp,          :trn => :non_pers_transp, :wtt => :non_pers_transp,
:trk => :non_pers_transp,          :grd => :non_pers_transp, :otr => :non_pers_transp, 
:mot => :veh_exp,                  :ote => :veh_exp, :mvt => :veh_exp,
:veh_elect => :fuels,               :petr => :fuels, 
:uel_veh_elect => :veh_elect,             :rnw_veh_elect => :veh_elect, #These have different names to distinguish both times that electricity appears
:pet => :petr, :oil => :petr,

:goods => :non_transp,                        :housing_exp => :non_transp, 
:pla  => :goods, :soc  => :goods, :leg  => :goods, :rnt  => :goods, :fof  => :goods, :mov  => :goods, :adm  => :goods, :slg  => :goods,
 :fdd  => :goods, :man  => :goods, :fen  => :goods, :pri  => :goods, :tex  => :goods, :sle  => :goods, :pip  => :goods, :pub  => :goods,
  :brd  => :goods, :bnk  => :goods, :cep  => :goods, :rec  => :goods, :uwt  => :goods, :eec  => :goods, :nmp  => :goods, :osv  => :goods,
   :con  => :goods, :ppd  => :goods, :amd  => :goods, :edu  => :goods, :wst  => :goods, :smn  => :goods, :fpd  => :goods, :pmt  => :goods,
    :sec  => :goods, :wrh  => :goods, :hos  => :goods, :gas  => :goods, :min  => :goods, :fbp  => :goods, :fbt  => :goods, :tsv  => :goods,
     :res  => :goods, :dat  => :goods, :fnd  => :goods, :ins  => :goods, :wht  => :goods, :nrs  => :goods, :art  => :goods, :fmt  => :goods,
      :che  => :goods,  :alt  => :goods, :agr  => :goods, :gmt  => :goods, :fin  => :goods, :mmf  => :goods, :amb  => :goods,
       :mch  => :goods, :wpd  => :goods, :ott  => :goods, :com  => :goods, #:oil  => :goods, oil which is 0, is combined with pet to make a nest bc leaf+nest doesn't work
:own_occ_exp => :housing_exp,              :home_nrg_expd => :housing_exp,
:hou => :own_occ_exp,                      :ore => :own_occ_exp,
:elect => :home_nrg_expd,  :homefuels => :home_nrg_expd,
:ugs => :homefuels,                       :coa => :homefuels,
:uel_elect => :elect,             :rnw_elect => :elect,  #These have different names to distinguish both times that electricity appears
)

# ### Alternate version for testing a single nest level of the same structure, to see if it matches the old demand elasticity
function ces_nodeAlt(
    bnch_quants::Dict{Symbol,Float64}, # hold the name of the inputs and benchmark data quantities (can be leaves or nests, with the total benchmark - data for leaves, or from nest returns
    quantities::Dict{Symbol,Float64},#; # hold the name of the inputs and utilities (can be leaves or nests, with the total benchmark - data for leaves, or from nest returns
)   
    
    σ = 0.999999 #elasticity(FDem,0.999)
    # ρ = 1 - 1 / σ   # calculate once, so \rho is a constant
    total_benchmark = sum(values(bnch_quants))
    weights = Dict(q => bnch_quants[q] / total_benchmark for q in keys(quantities))
    utility = value(FDem) *sum(weights[i] *  quantities[i]^((σ-1)/σ) for i in keys(quantities) if quantities[i]>0)^(σ / (σ-1))
    # utility = sum(weights[i] * quantities[i]^ρ for i in keys(quantities) if quantities[i]>0)^(1 / ρ)

    return (
        Dict(nest => total_benchmark),
        Dict(nest => utility),)
end

bnch_quants = leaves
quantities = Dict(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in Ip)

function CESutility_multinatAlt()

leaves = Dict(x => pce_0[x,:pce] for x in Ip)
goods_i, goods_ute = ces_nodeAlt(leaves,Dict(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in Ip))
    return only(values(goods_ute))
end