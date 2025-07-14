# # CES Aggregator Function
# function ces_aggregate(quantities::Vector{Float64}, weights::Vector{Float64}, rho::Float64)
#     # if abs(rho - 0.0) < 1e-8
#         # Cobb-Douglas limit
#         # return prod(q^w for (q, w) in zip(quantities, weights))
#     # else
#         return (sum(w * q^rho for (q, w) in zip(quantities, weights)))^(1/rho)
#     # end
# end

# # Want to return quantities, weights, parent_elasticity so they can commect as arguments at the level above. 
# ## BUT, then have to be joined, like
# ces_agg([q2 ces_agg(s,q,benchmark_q_combo,parent_rho)[1]],[w2 ces_agg(s,q,benchmark_q_combo,parent_rho)[2]],[bnch2 ces_agg(s,q,benchmark_q_combo,parent_rho)[3]],parent_rho)

# ##args #1 ces_util (name quantity) #2 share based on either index, OR combination of indexes (can I do that recursively? )

# # Shares, instead of a vector, another Dict --> terrible idea

# Dict(:pet=>pce_0[:pet,:pce], :veh_elect=>sum(pce_0[[:rnw, :uel],:pce]))

# # For shares: for each I
# for (q,w) in
# zip(x, params)
# println(q,"  ",w)
# end

# # Nested CES Utility Function
# function nested_ces_utility(x::Dict{String, Float64}, params) # x is the elasticities, params has the name of shares (keys), and shares (values)
#     # Unpack parameters
#     rho_s = params["rho_s"]
#     rho_a = params["rho_a"]
#     rho_b = params["rho_b"]
#     weights_s = params["weights_s"]
#     weights_a = params["weights_a"]
#     weights_b = params["weights_b"]

#     # Lower level nests: a and b
#     C_a = ces_aggregate([x["x1"], x["x2"]], weights_a, rho_a)#for the leaves
#     C_b = ces_aggregate([x["x3"], x["x4"]], weights_b, rho_b)

#     # Top level nest: s
#     U = ces_aggregate([C_a, C_b], weights_s, rho_s)  # top level, shares between C_a and Cb, elasticity of 
#     return U
# end

# ### Arguments => also returns 
# sector, node_name->nest->Symbol, bnchmark_sum, quantity, 

# CES(FDem, [:homefuels :elect],[] )

# p = production(FDem)
# n = MPSGE.netputs(p)
# MPSGE.parent(n[PA[:uwt]][1])
# elasticity(MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:pet]][1])[1])
# # get own next name
# MPSGE.name(MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:pet]][1])[1])
# ## Extract parent name
# MPSGE.name((MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:pet]][1])[1]).parent)
# elasticity((MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:pet]][1])[1]).parent)

# value(compensated_demand(FDem,PA[i]))

# value(compensated_demand(sector,PA[i]))
#     # quantities = 
#     # quantities for leaves
#     Dict(x=>value(compensated_demand(sector,PA[x])) for x in [:ugs,:coa])
# CESutil(FDem, Dict(:ugs=>pce_0[:ugs,:pce],:coa=>pce_0[:coa,:pce]),:homefuels, Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:ugs,:coa]))[1]
# CESutil(FDem, Dict(:ugs=>pce_0[:ugs,:pce],:coa=>pce_0[:coa,:pce]),:homefuels, Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:ugs,:coa]))[2]
#     # quantities for leaf and CES branch
# # Dict(:pet => value(compensated_demand(sector,PA[:pet], CESutil(FDem, [:rnw, :uel],:veh_elect,))))
# CESutil(FDem, Dict(:pet=>pce_0[:ugs,:pce], :fuels=>sum(values(CESutil(FDem, Dict(:rnw=>pce_0[:rnw,:pce],:uel=>pce_0[:uel,:pce]),:veh_elect,Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:rnw, :uel]))[2]))),:fuels, Dict(:pet=>value(compensated_demand(FDem,PA[:pet])),:veh_elect=>CESutil(FDem,Dict(:rnw=>pce_0[:rnw,:pce],:uel=>pce_0[:uel,:pce]) ,:veh_elect,Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:rnw, :uel]))[1]))

# I = [:pet; CESutil(FDem, [:rnw, :uel],:veh_elect,Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:rnw, :uel]))[2]]
# MPSGE.name(MPSGE.parent(MPSGE.netputs(production(FDem))[PA[i]][1])[1])



### This now theoretically works for leaves
# function CESutil(sector::ScalarSector,I::Dict, nest::Symbol,quantities::Dict)
#     # shares = Dict(x=>pce_0[x,:pce] for x in I)
#     # OR it's the value of the CESutil(s,II, n)
#     return sum([(I[i]/sum(values(I)))*(value(sector)*quantities[i])^((elasticity(sector,nest)-1)/elasticity(sector,nest)) for i in keys(I) if quantities[i]>0])^(elasticity(sector,nest)/(elasticity(sector,nest)-1)), I
# end
# using MPSGE
# indcommods = [a, c]

"""
At the actual leaf levels of each electricity production, it doesn't matter because the weights are the same.
It matters at every level above because oterhwise the electricity would be over-weighted (by a tiny bit for house, a massive amount for veh_elect)
AND at the top level that would lead to complete double counting of the electricity

if uel in indexes && rnw in indexes && electritity = elect
 subtract the amount at the leaf level of pce_0 that is not for elect

if uel in indexes && rnw in indexes && electritity = veh_elect

if uel in indexes && rnw in indexes && electritity = all
Ummm, here the electricty is finally broaght together. So we need uel and rnw to be incorporated twice and separately in each of the (transp and non_transp) branches.

!!!!Maybe this suggests a different way of getting the shares. Maybe instead of saving the indices, I save the total benchmark values!! Then I just need to adjust the leaves, and the rest will be inherited.!!!

else error

"""


"""
Issue with new plan is just that I can't use the index values (quanitites) to look
"""
function ces_node(
    # benchmark_data::Dict{Symbol, Dict{Symbol, Float64}},
    # indcommods::Vector{Dict{Symbol, Vector{Symbol}}},
    bnch_quants::Dict{Symbol,Float64}, # hold the name of the inputs and benchmark data quantities (can be leaves or nests, with the total benchmark - data for leaves, or from nest returns
    # sector::MPSGE.ScalarSector,
    # commodity::MPSGE.IndexedCommodity,
    quantities::Dict{Symbol,Float64},#; # hold the name of the inputs and utilities (can be leaves or nests, with the total benchmark - data for leaves, or from nest returns
    # electricity=nothing
    # nests::Vector{Symbol}
)   
    # ar = collect(Iterators.flatten(values(indcommods[1])))
    # br= collect(Iterators.flatten(values(indcommods[2])))
    # indexes = vcat(ar,br)
    # indexes = vcat(ar,br)
    nests = []
    elases = []
    for i in keys(quantities)
        # push!(nests, MPSGE.name(MPSGE.parent(MPSGE.netputs(production(FDem))[PA[i]][1])[1])) # Oh, this only works for leaves
        push!(nests, ch_pr[i][1])
        push!(elases, ch_pr[i][2])
        # push!(nests, nest_name][1])
    end
    if !all( ==(nests[1]), nests)
        error("Error: not the same nest")
    end
    if !all( ==(elases[1]), elases)
        error("Error: These elasticity values do not match")
    end
    # if all( ==(:electr), nests)
    #     nests .= :elect
    # end
    nest = nests[1] # this is the parent, the one they're going into 
    σ = elases[1]
    
    # bnch = Dict(i=>pce_0[i,:pce] for i in indexes)
    # total_benchmark = sum(values(bnch))
    total_benchmark = sum(values(bnch_quants))
    # prepare benchmark values for use
    # df to connect the values to the nest names, and then groupby
#     bnch_nestqs = DataFrame(nest=Symbol[], bnval=[])
#     for k in keys(bnch)
#         if k∈keys(quantities) # check if it's a leaf, then leave the name
#             push!(bnch_nestqs,[k, bnch[k]])
#         else # otherwise get the parent name for grouping
#             push!(bnch_nestqs,[ch_pr[k][1], bnch[k]])
#         end
#     end
#     bnch_nestqs = combine(groupby(bnch_nestqs, :nest),:bnval => sum =>:bnval)
#     transform!(bnch_nestqs, :nest => ByRow(a -> a ==:electr ? electricity : a)=> :nest)
#     regroup = DataFrame(nest=Symbol[], bnval=[])
#     for n in eachrow(bnch_nestqs)
#         if n[:nest] ∉keys(quantities)
#             println("nest level is ",n[:nest],", nest above, not in ",keys(quantities)," is",ch_pr[n[:nest]][1])
#             # push!(regroup,[ch_pr[n[:nest]][1], n[:bnval]])
#             # bnch_nestqs = filter(r -> r.nest !=n[:nest], bnch_nestqs)
#         end
#     end
#     if !isempty(regroup)
#         regroup = combine(groupby(regroup, :nest),:bnval => sum =>:bnval)
#     end
#     bnch_nestqs2 = vcat(regroup,bnch_nestqs)
# transform!(bnch_nestqs2, :nest => ByRow(a -> a ==:electr ? electricity : a)=> :nest)
#     regroup2 = DataFrame(nest=Symbol[], bnval=[])
#     for n in eachrow(bnch_nestqs2)
#         if n[:nest] ∉keys(quantities)
#             println(n[:nest],ch_pr[n[:nest]][1])
#             push!(regroup2,[ch_pr[n[:nest]][1], n[:bnval]])
#             bnch_nestqs2 = filter(r -> r.nest !=n[:nest], bnch_nestqs2)
#         end
#     end
#     if !isempty(regroup2)
#         regroup2 = combine(groupby(regroup2, :nest),:bnval => sum =>:bnval)
#     end
#     # Turn into NArray so the slicing is easier
#     nest_qs = NamedArray(vcat(regroup2,bnch_nestqs2)[:,:bnval],vcat(regroup2,bnch_nestqs2)[:,:nest])
    # nest_qs = NamedArray(bnch_nestqs[:,:bnval_sum],bnch_nestqs[:,:nest])

        weights = Dict(q => bnch_quants[q] / total_benchmark for q in keys(quantities)) # TODO How to get the weights categorised by the quantities, not just the single level up.
        # weights = Dict(q => nest_qs[q] / total_benchmark for q in keys(quantities)) # TODO How to get the weights categorised by the quantities, not just the single level up.
    # weights = Dict(q => nest_qs[q] / total_benchmark for q in names(nest_qs))

    # --- Get σ and compute ρ ---
    # σ = elasticity((MPSGE.parent(MPSGE.netputs(production(FDem))[PA[indexes[end]]][1])[1]).parent)  #TODO, is this right? Good enough just to take the last?
    # σ = ch_pr[nest][2]
    ρ = 1 - 1 / σ
    # --- Get compensated demand from model results ---
    # quantities = Dict(i => comp_demand(sector, commodity, nest, i) for i in indexes)

    # --- CES aggregation ---
    utility = sum(weights[i] * quantities[i]^ρ for i in keys(quantities) if quantities[i]>0)^(1 / ρ)
    # --- Return in standardized form ---
    return (
        # Dict(nest => [i for i in indexes]),
        Dict(nest => total_benchmark),
        # FDem,
        # commodity,
        Dict(nest => utility),
        # electricity
    )
end
## Multiple dispatch, this one just for the two electricity leaf nodes
# # function ces_node(
#     indexes::Vector{Symbol},
#     # sector::MPSGE.ScalarSector,
#     # commodity::MPSGE.IndexedCommodity,
#     quantities::Dict{Symbol,Float64},
#     electricity::Symbol)

#     nests = []
# # uel and rnw appear twice, so hack is to add a keyword for them
#     for i in keys(quantities)
#         if i ∉ [:uel, :rnw] 
#             error("This method is only for electricity at the leaves")
#         end
#         push!(nests, electricity)
#     end
#     if !all( ==(nests[1]), nests)
#         error("Error: not the same nest")
#     end
#     nest = nests[1]
#     bnch = Dict()
#     for i in indexes
#         if nest == :elect
#             bnch[i] = pce_0[i,:pce]*(1-veh_chrg_pc)
#         elseif nest== :veh_elect
#             bnch[i] = pce_0[i,:pce]*veh_chrg_pc
#         end
#     end
#     total_benchmark = sum(values(bnch))
#         # prepare benchmark values for use
#         # df to connect the values to the nest names, and then groupby
#     bnch_nestqs = DataFrame(nest=[], bnval=[])
#         for k in keys(bnch)
#             if k∈keys(quantities) # check if it's a leaf, then leave the name
#                 push!(bnch_nestqs,[k, bnch[k]])
#             else # otherwise get the parent name for grouping
#                 push!(bnch_nestqs,[ch_pr[k][1], bnch[k]])
#             end
#         end
#     bnch_nestqs = combine(groupby(bnch_nestqs, :nest),:bnval => sum)
#     # Turn into NArray so the slicing is easier
#     nest_qs = NamedArray(bnch_nestqs[:,:bnval_sum],bnch_nestqs[:,:nest])
#     weights = Dict(q => nest_qs[q] / total_benchmark for q in keys(quantities))

#     # --- Get σ and compute ρ ---
#     # σ = elasticity((MPSGE.parent(MPSGE.netputs(production(FDem))[PA[indexes[1]]][1])[1]).parent)
#     σ = ch_pr[nest][2]
#     ρ = 1 - 1 / σ

#     if nest == :elect
#         utility = sum(weights[i] * (quantities[i]*(1-veh_chrg_pc))^ρ for i in keys(quantities))^(1 / ρ)
#     elseif nest== :veh_elect
#         utility = sum(weights[i] * (quantities[i]*veh_chrg_pc)^ρ for i in keys(quantities))^(1 / ρ)
#     end
#     # --- CES aggregation ---
#     # --- Return in standardized form ---
#     return (
#         Dict(nest => [i for i in indexes]),
#         # FDem,
#         # commodity,
#         Dict(nest => utility),
#         electricity
#     )
# end

# indexes = [collect(Iterators.flatten(values(a)));collect(Iterators.flatten(values(c)))]
# quantities = merge(b,d)
# indexes = [collect(Iterators.flatten(values(a)));collect(Iterators.flatten(values(c)))]
# quantities = merge(b,d)
# electricity=z

# indexes = [collect(Iterators.flatten(values(a))); collect(Iterators.flatten(values(e)))]
# quantities = merge(d,h)
# indexes = [:coa, :ugs]
# indexes = [:rnw, :uel]
# quantities = Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:rnw, :uel])
# quantities = Dict([x=>value(compensated_demand(FDem,PA[x])) for x in [:coa,:ugs]])
# include("MultiGHGNat.jl")

# TODO 4) Email Larry - include job app, workshops, MEEW... maybe note re slow, re work feeding into implementing new MPSGE.jl functionality

function CESutility_multinat()
### For leaves
leaf_ugs = Dict(:ugs => pce_0[:ugs,:pce])
leaf_coa = Dict(:coa => pce_0[:coa,:pce])
# Leaves [:ugs,:coa] into :homefuels
homefuel_i,homefuel_ute = ces_node(merge(leaf_coa,leaf_ugs),  Dict([x=>value(compensated_demand(FDem,PA[x])) for x in [:coa,:ugs]]))#,electricity=:elect)
homefuel_i 
homefuel_ute


# Leaves [:rnw, :uel] into :elect
leaf_rnw = Dict(:rnw_elect => pce_0[:rnw,:pce] * (1-veh_chrg_pc))
leaf_uel = Dict(:uel_elect => pce_0[:uel,:pce] * (1-veh_chrg_pc))
Electq_names = Dict(:rnw=>:rnw_elect,:uel=>:uel_elect)
elect_i,elect_ute = ces_node(merge(leaf_rnw,leaf_uel), Dict(Electq_names[x]=>value(compensated_demand(FDem,PA[x],:elect)) for x in [:rnw, :uel]))
elect_i
elect_ute

# Nests [:homefuels, :elect] into [:home_nrg_expd]
home_nrg_expd_i,home_nrg_expd_ute = ces_node(merge(homefuel_i,elect_i),merge(homefuel_ute,elect_ute))

# Leaves [:hou,:ore] into :own_occ_exp
leaf_hou = Dict(:hou => pce_0[:hou,:pce])
leaf_ore = Dict(:ore => pce_0[:ore,:pce])
own_occ_exp_i, own_occ_exp_ute= ces_node(merge(leaf_hou,leaf_ore),  Dict([x=>value(compensated_demand(FDem,PA[x])) for x in [:hou,:ore]]))
    
# Nests [:own_occ_exp, :home_nrg_expd] into [:housing_exp]
housing_exp_i, housing_exp_ute = ces_node(merge(home_nrg_expd_i,own_occ_exp_i),merge(home_nrg_expd_ute,own_occ_exp_ute))

leaves_goods = Dict(x => pce_0[x,:pce] for x in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa ]),Ip))
goods_i, goods_ute = ces_node(leaves_goods,Dict(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in filter(∉([:air :trn :wtt :trk :grd :otr :mot :ote :mvt :uel :rnw :pet :hou :ore :ugs :coa ]),Ip)))

# Nest :goods , :housing_exp into :non_transp, 
non_transp_i, non_transp_ute = ces_node(merge(housing_exp_i,goods_i),merge(housing_exp_ute,goods_ute))


# Leaves [:rnw, :uel] into :elect
leaf_rnw = Dict(:rnw_elect => pce_0[:rnw,:pce] * (1-veh_chrg_pc))
leaf_uel = Dict(:uel_elect => pce_0[:uel,:pce] * (1-veh_chrg_pc))
Electq_names = Dict(:rnw=>:rnw_elect,:uel=>:uel_elect)
elect_i,elect_ute = ces_node(merge(leaf_rnw,leaf_uel), Dict(Electq_names[x]=>value(compensated_demand(FDem,PA[x],:elect)) for x in [:rnw, :uel]))
elect_i
elect_ute

# Leaves[[:rnw, :uel]] into [:veh_elect]
leaf_rnw_veh = Dict(:rnw_veh_elect => pce_0[:rnw,:pce] * veh_chrg_pc)
leaf_uel_veh = Dict(:uel_veh_elect => pce_0[:uel,:pce] * veh_chrg_pc)
Veh_electq_names = Dict(:rnw=>:rnw_veh_elect,:uel=>:uel_veh_elect)
veh_elect_i,veh_elect_ute = ces_node(merge(leaf_rnw_veh,leaf_uel_veh), Dict(Veh_electq_names[x]=>value(compensated_demand(FDem,PA[x],:veh_elect)) for x in [:rnw, :uel]))

# Nest/Leaf :pet, :veh_elect into :fuels, 
leaf_pet = Dict(:pet => pce_0[:pet,:pce])
fuels_i, fuels_ute = ces_node(merge(veh_elect_i,leaf_pet),merge(veh_elect_ute,Dict(:pet=>value(compensated_demand(FDem,PA[:pet])))))

# Leaves [:mot :ote :mvt] =>veh_exp
leaves_vehexp = Dict(x => pce_0[x,:pce] for x in [:mot, :ote, :mvt])
veh_exp_i, veh_exp_ute = ces_node(leaves_vehexp,Dict(x=>Float64(value(compensated_demand(FDem,PA[x]))) for x in [:mot :ote :mvt]))

# Nest fuels, veh_exp=>pers_transp=0.25,                     
pers_transp_i, pers_transp_ute = ces_node(merge(fuels_i, veh_exp_i),merge(fuels_ute, veh_exp_ute))

# Leaves [:air :trn :wtt :trk :grd :otr] => non_pers_transp
leaves_np_trans = Dict(x => pce_0[x,:pce] for x in [:air,:trn,:wtt,:trk,:grd,:otr])
non_pers_transp_i, non_pers_transp_ute = ces_node(leaves_np_trans,Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:air,:trn,:wtt,:trk,:grd,:otr]))

# Nest non_pers_transp=>transp=0.25,                  pers_transp=>transp=0.4, 
transp_i, transp_ute = ces_node(merge(non_pers_transp_i, pers_transp_i),merge(non_pers_transp_ute, pers_transp_ute))

# Nest transp, non_transp=>s,
s_i, s_ute = ces_node(merge(transp_i, non_transp_i),merge(transp_ute, non_transp_ute))

    return only(values(s_ute))
end

# quantities = merge(d,Dict(:elect=>value(compensated_demand(FDem,PA[:uel]))))
# quantities = Dict([:nat_gas, :fuel_oil] .=>[value(compensated_demand(FDem,PA[x])) for x in [:ugs,:coa]])
# elasticity(MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:ugs]][1])[1])
# # get own next name
# MPSGE.name(MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:ugs]][1])[1])
# ## Extract parent name
# MPSGE.name((MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:ugs]][1])[1]).parent)
# elasticity((MPSGE.parent(MPSGE.netputs(production(FDem))[PA[:ugs]][1])[1]).parent)
# MPSGE.name(MPSGE.parent(MPSGE.netputs(production(FDem))[:veh_fuel][1])[1])

# CESutil(FDem, Dict(:ugs=>pce_0[:ugs,:pce],:coa=>pce_0[:coa,:pce]),:homefuels, Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:ugs,:coa]))[1]
# CESutil(FDem, Dict(:pet=>pce_0[:ugs,:pce], :fuels=>sum(values(CESutil(FDem, Dict(:rnw=>pce_0[:rnw,:pce],:uel=>pce_0[:uel,:pce]),:veh_elect,Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:rnw, :uel]))[2]))),:fuels, Dict(:pet=>value(compensated_demand(FDem,PA[:pet])),:veh_elect=>CESutil(FDem,Dict(:rnw=>pce_0[:rnw,:pce],:uel=>pce_0[:uel,:pce]) ,:veh_elect,Dict(x=>value(compensated_demand(FDem,PA[x])) for x in [:rnw, :uel]))[1]))

# # Edited from elastiticity vector
# ### Keys are the ends, and values are the nodes and the elasticity from the node
ch_pr = Dict(
:transp=>[:s,0.5],  :non_transp=>[:s,0.5],
                            
:non_pers_transp=>[:transp,0.35],       :pers_transp=>[:transp,0.35], 
:veh_exp=>[:pers_transp,0.4],           :fuels=>[:pers_transp,0.4],
:air=>[:non_pers_transp,0.25],          :trn=>[:non_pers_transp,0.25], :wtt=>[:non_pers_transp,0.25],
:trk=>[:non_pers_transp,0.25],          :grd=>[:non_pers_transp,0.25], :otr=>[:non_pers_transp,0.25], 
:mot=>[:veh_exp,0.25],                  :ote=>[:veh_exp,0.25], :mvt=>[:veh_exp,0.25],
:veh_elect=>[:fuels,0.25],               :pet=>[:fuels,0.25], 
:uel_veh_elect=>[:veh_elect,Elas[:uel,:SAGE_en]],             :rnw_veh_elect=>[:veh_elect, Elas[:uel,:SAGE_en]], #This is for aggregating to get the share wieghts

:goods=>[:non_transp,0.5],                        :housing_exp=>[:non_transp,0.5], 
:pla => [:goods, 0.87], :soc => [:goods, 0.87], :leg => [:goods, 0.87], :rnt => [:goods, 0.87], :fof => [:goods, 0.87], :mov => [:goods, 0.87], :adm => [:goods, 0.87], :slg => [:goods, 0.87],
 :fdd => [:goods, 0.87], :man => [:goods, 0.87], :fen => [:goods, 0.87], :pri => [:goods, 0.87], :tex => [:goods, 0.87], :sle => [:goods, 0.87], :pip => [:goods, 0.87], :pub => [:goods, 0.87],
  :brd => [:goods, 0.87], :bnk => [:goods, 0.87], :cep => [:goods, 0.87], :rec => [:goods, 0.87], :uwt => [:goods, 0.87], :eec => [:goods, 0.87], :nmp => [:goods, 0.87], :osv => [:goods, 0.87],
   :con => [:goods, 0.87], :ppd => [:goods, 0.87], :amd => [:goods, 0.87], :edu => [:goods, 0.87], :wst => [:goods, 0.87], :smn => [:goods, 0.87], :fpd => [:goods, 0.87], :pmt => [:goods, 0.87],
    :sec => [:goods, 0.87], :wrh => [:goods, 0.87], :hos => [:goods, 0.87], :gas => [:goods, 0.87], :min => [:goods, 0.87], :fbp => [:goods, 0.87], :fbt => [:goods, 0.87], :tsv => [:goods, 0.87],
     :res => [:goods, 0.87], :dat => [:goods, 0.87], :fnd => [:goods, 0.87], :ins => [:goods, 0.87], :wht => [:goods, 0.87], :nrs => [:goods, 0.87], :art => [:goods, 0.87], :fmt => [:goods, 0.87],
      :che => [:goods, 0.87], :oil => [:goods, 0.87], :alt => [:goods, 0.87], :agr => [:goods, 0.87], :gmt => [:goods, 0.87], :fin => [:goods, 0.87], :mmf => [:goods, 0.87], :amb => [:goods, 0.87],
       :mch => [:goods, 0.87], :wpd => [:goods, 0.87], :ott => [:goods, 0.87], :com => [:goods, 0.87],
:own_occ_exp=>[:housing_exp,0.3],              :home_nrg_expd=>[:housing_exp,0.3],
:hou=>[:own_occ_exp,0.25],                      :ore=>[:own_occ_exp,0.25],
:elect=>[:home_nrg_expd,0.75],  :homefuels=>[:home_nrg_expd,0.75],
# homefuelhome_fossil=>[:homefuels,0.999],
:ugs=>[:homefuels,0.5],                       :coa=>[:homefuels,0.5],
# :uel=>[:homefuels,],                         :rnw=>[:veh_elect,],
:uel_elect=>[:elect,Elas[:uel,:SAGE_en]],             :rnw_elect=>[:elect, Elas[:uel,:SAGE_en]], #This is for aggregating to get the share wieghts
)

## Benchmark
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
set_value!(CO₂_tax, 0.0)
set_value!(ta,ta_0DAA[Jp])
set_value!(tm,tm_0DAA[Jp])
solve!(MultiNat)# Temp measure to address residual price changes
bnch_util = CESutility_multinat()
println(" Consumption utility is : ",CESutility_multinat())
println("% change is" ,(CESutility_multinat()-bnch_util)/bnch_util*100)

## Counterfactual
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
set_value!(CO₂_tax, 0.0)
set_value!(ta,0)
set_value!(tm,0)
solve!(MultiNat)# Temp measure to address residual price changes
println(" Consumption utility is : ",CESutility_multinat())
println("% change is" ,(CESutility_multinat()-bnch_util)/bnch_util*100)

## Set other taxes back to benchmark
set_value!(ta,ta_0DAA[Jp])
set_value!(tm,tm_0DAA[Jp])

## CH4 tax
set_value!(CO₂_tax, 0)
set_value!(CH₄_tax, CH4_taxrate)
solve!(MultiNat, cumulative_iteration_limit=10000) #;
println(" Consumption utility is : ",CESutility_multinat())
println("% change is" ,(CESutility_multinat()-bnch_util)/bnch_util*100)

# CO2 tax
set_value!(CH₄_tax, 0.0) ## Set CH4 taxes back to 0 to generate CO2 tax ONLY
set_value!(CO₂_tax, CO2_taxrate)
solve!(MultiNat, cumulative_iteration_limit=10000) #;
println(" Consumption utility is : ",CESutility_multinat())
println("% change is" ,(CESutility_multinat()-bnch_util)/bnch_util*100)

# Both
set_value!(CO₂_tax, CO2_taxrate)
set_value!(CH₄_tax, CH4_taxrate)
solve!(MultiNat, cumulative_iteration_limit=10000) #;
println(" Consumption utility is : ",CESutility_multinat())
println("% change is" ,(CESutility_multinat()-bnch_util)/bnch_util*100)
