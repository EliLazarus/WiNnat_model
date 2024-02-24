# Replication of the WiNDC national MGE model
using MPSGE, JLD2, CSV

using JuMP,PATHSolver
# For sensitivity tests, timing, and plotting
using DataFrames, Plots, Tables, Dates, Distributions

""" 
This function generates a GAMS-like report detailing both the values of the variables,
but also the evaluation of the complementary constraint. Very useful for debugging as the
complementary constraint should always have value 0.
"""
# for i in 1:length(all_constraints(WiNnat._jump_model; include_variable_in_set_constraints = false))
# 	 try println(i,"::",extract_variable_ref(constraint_object(all_constraints(WiNnat._jump_model;include_variable_in_set_constraints = false)[i]).func[2]),": ", value(constraint_object(all_constraints(WiNnat._jump_model;include_variable_in_set_constraints 
# 	= false)[i]).func[1]))
#  	catch 
# 	println(i, "::",extract_variable_ref(constraint_object(all_constraints(WiNnat._jump_model; include_variable_in_set_constraints= false)[i]).func[2]))
# 	end
# 	   end

extract_variable_ref(v::NonlinearExpr) = v.args[1]
extract_variable_ref(v::AffExpr) = collect(keys(v.terms))[1]
extract_variable_ref(v::QuadExpr) = extract_variable_ref(v.aff)

function generate_report(m::JuMP.Model; decimals::Int = 15, mdecimals::Int = 4)
	#mcp_data = Complementarity.get_MCP_data(m)
	#vars = all_variables(m)
	#sols = Dict(zip(vars,value.(vars)));

	mapping = Dict()
	for ci in all_constraints(m; include_variable_in_set_constraints = false)
		c = constraint_object(ci)
		mapping[extract_variable_ref(c.func[2])] = c.func[1]
	end

	out = "var_name,value,margin\n"
	for elm in all_variables(m)

		# val = round(value(elm),digits = decimals)
		val = JuMP.is_parameter(elm) ? round(JuMP.parameter_value(elm), digits = decimals) : round(JuMP.value(elm), digits=decimals)

		margin = "."
		try
			margin = round(value(mapping[elm]),digits = mdecimals)
		catch
			margin = "."
		end
		

		out = out*"\"$elm\",$val,$margin\n"
	end

	return(out)
end


cd(dirname(Base.source_path()))
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in date from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in date from saved Notebook output Dict, named S
# Alternate, Julia WiNDC generated data
PJ= load(joinpath(@__DIR__,"./data/nationaldata_julia/JDAAData.jld2"))["data"] # load in date from saved Notebook output Dict, named P
SJ= load(joinpath(@__DIR__,"./data/nationaldata_julia/JIndices.jld2"))["data"] # load in date from saved Notebook output Dict, named S

y_ = filter!(x -> x != :oth && x!= :use, S[:i][:]) # These 2 sectors 'use' & 'oth' are in the indices list, but have no data (and therefore cause problems)
a_ = filter!(x -> x != :fbt && x != :mvt && x != :gmt, copy(y_))
# a_ = copy(y_)

# Indexes (set from the data files, via the notebook)
sectorsi = S[:i]#[:] # "BEA Goods and sectors categories", is "i" in GAMS
sectorsj = copy(sectorsi) # "BEA Goods and sectors categories", is "j" in GAMS, for iterating over double index
xfd = filter!(x -> x != :pce, S[:fd]) # "BEA Final demand categories",
ts = S[:ts] # "BEA Taxes and subsidies categories",
valueadded = filter!(s -> s != :othtax, S[:va]) # "BEA Value added categories excluding othtax", va in GAMS
margin  = S[:m] # "Margins (trade or transport)"; m in GAMS

# for datayear in 1997:2017
datayear = 2017
yr = Symbol(datayear)
# PARAMETERS

# Data For a single year, knock out one dimension
y_0 = P[:y_0][yr,:][y_] #	"Gross output",
ys_0 = P[:ys_0][yr,y_,y_] #	"Sectoral supply",
ty_0 = P[:ty_0][yr,:][y_] #	"Output tax rate"
fs_0 = P[:fs_0][yr,:][y_] #	"Household supply", # All zeros
id_0 = P[:id_0][yr,y_,y_] #	"Intermediate demand",
fd_0 = P[:fd_0][yr,y_,:] #	"Final demand",
va_0 = P[:va_0][yr,:,y_] #	"Value added",
# ts_0 = P[:ts_0][yr,:,:] #	"Taxes and subsidies", Not in this model
m_0 = P[:m_0][yr,:][y_] #	"Imports",
x_0 = P[:x_0][yr,:][y_] #	"Exports of goods and services",
# mrg_0 = P[:mrg_0][yr,:] #	"Trade margins", Not in this model
# trn_0 = P[:trn_0][yr,:] #	"Transportation costs",  Not in this model
# duty_0 = P[:ty_0][yr,:] #	"Import duties", Not in this model
# sbd_0 = P[:sbd_0][yr,:] #	"Subsidies on products", Not in this model
# tax_0 = P[:tax_0][yr,:] #	"Taxes on products", Not in this model
ms_0 = P[:ms_0][yr,y_,:] #	"Margin supply",
md_0 = P[:md_0][yr,:,y_] #	"Margin demand",
# s_0 = P[:s_0][yr,:] #	"Aggregate supply", Not in this model
a_0 = P[:a_0][yr,:][a_]  #	"Armington supply",
bopdef_0 = P[:bopdef_0][yr] #	"Balance of payments deficit",
ta_0 = P[:ta_0][yr,:][y_] #	"Tax net subsidy rate on intermediate demand", Initial, for price
tm_0 = P[:tm_0][yr,:][y_] #	"Import tariff"; Initial, for price 

# ty_0 = add!(WiNnat, Parameter(:ty, indices = (sectorsj,), value=P[:ty_0][year,:].data)) #	"Output tax rate", Not in this model

# Option to set model build and solve as function for time tests
# function timeWiNnat(n::Int64)
WiNnat = MPSGE.Model()

	# parameters
	ta = add!(WiNnat, MPSGE.Parameter(:ta, indices = (a_,), value=P[:ta_0][yr,a_].data)) #	"Tax net subsidy rate on intermediate demand",
	tm = add!(WiNnat, MPSGE.Parameter(:tm, indices = (a_,), value=P[:tm_0][yr,a_].data)) #	"Import tariff";

	# Elasticity parameters
	t_elas_y =  add!(WiNnat, MPSGE.Parameter(:t_elas_y,  value=0.))
	elas_y =    add!(WiNnat, MPSGE.Parameter(:elas_y,    value=0.))
	elas_va =   add!(WiNnat, MPSGE.Parameter(:elas_va,   value=1.))
	t_elas_m =  add!(WiNnat, MPSGE.Parameter(:t_elas_m,  value=0.))
	elas_m =    add!(WiNnat, MPSGE.Parameter(:elas_m,    value=0.))
	t_elas_a =  add!(WiNnat, MPSGE.Parameter(:t_elas_a,  value=2.))
	elas_a =    add!(WiNnat, MPSGE.Parameter(:elas_a,    value=0.))
	elas_dm =   add!(WiNnat, MPSGE.Parameter(:elas_dm,   value=2.))
	d_elas_ra = add!(WiNnat, MPSGE.Parameter(:d_elas_ra, value=1.))

	# sectors:
	Y = add!(WiNnat, Sector(:Y, indices=(sectorsj,)))
	A = add!(WiNnat, Sector(:A, indices=(a_,)))

	MS = add!(WiNnat, Sector(:MS, indices=(margin,)))

	# commodities:
	PA  = add!(WiNnat, Commodity(:PA, indices=(a_, ))) #	Armington price
	PY  = add!(WiNnat, Commodity(:PY, indices=(sectorsi,))) #	Supply
	PVA = add!(WiNnat, Commodity(:PVA, indices=(valueadded,))) #		Value-added
	PM  = add!(WiNnat, Commodity(:PM, indices=(margin,))) #		Margin
	PFX = add!(WiNnat, Commodity(:PFX))	#	Foreign exchnage

	# consumers:
	RA = add!(WiNnat, Consumer(:RA, benchmark = sum(fd_0[:,:pce]) ))

	# production functions
	for j in y_
	# Commenting in/out for options to use Floats or Parameters for elasticities (parameters for sensitivity tests)	
		@production(WiNnat, Y[j], 0., 0.,
		# @production(WiNnat, Y[j], :($(t_elas_y)*1), 0.,
		# @production(WiNnat, Y[j], :($(t_elas_y)*1), :($(elas_y)*1),
		# @production(WiNnat, Y[j], 0., :($(elas_y)*1),
		[	
			Output(PY[i], ys_0[j,i], taxes=[Tax(ty_0[j], RA)]) for i in sectorsi if ys_0[j,i]>0
		], 
		[
			[Input(PA[i], id_0[i,j]) for i in a_ if id_0[i,j]>0];  # filtered to A
			[Input(Nest(
					Symbol("VA$j"),
					1.,
					# :($(elas_va)*1),
					sum(va_0[:,j]),
							[Input(PVA[va], va_0[va,j]) for va in valueadded if va_0[va,j]>0.] 
						),
						sum(va_0[:,j] )
				  )
			]
		]
	)
	end

	for m in margin
		add!(WiNnat, Production(MS[m], 0., 0., 
		# add!(WiNnat, Production(MS[m], 0, :($(elas_m)*1), 
		# add!(WiNnat, Production(MS[m], :($(t_elas_m)*1), :($(elas_m)*1), 
			[Output(PM[m], sum(ms_0[:,m]) ) ],
			[Input(PY[i], ms_0[i,m]) for i in sectorsi if ms_0[i,m]>0])) 
	end

	for i in a_  
		# @production(WiNnat, A[i], 2, :($(elas_a)*1),
		# @production(WiNnat, A[i], :($(t_elas_a)*1), :($(elas_a)*1),
		@production(WiNnat, A[i], 2., 0.,
			[
				[
				Output(PA[i], a_0[i], taxes=[Tax(:($(ta[i])*1), RA)], price=(1-ta_0[i]) )
				];
				[
					Output(PFX, x_0[i])
				]
			]
			,
				[
					[	
						Input(Nest(Symbol("dm$i"),
						2.,
						# :($(elas_dm)*1),
						(y_0[i]+m_0[i]+m_0[i]*get_value(tm[tm[i].subindex])),
						if m_0[i]>0 && y_0[i]>0
							[
								Input(PY[i], y_0[i] ),
								Input(PFX, m_0[i], taxes=[Tax(:($(tm[i])*1), RA)],  price=(1+tm_0[i]*1)  )
							]
						elseif y_0[i]>0
							[
								Input(PY[i], y_0[i] )
							]
						end
								),
						(y_0[i]+m_0[i]+m_0[i]*get_value(tm[tm[i].subindex])))
					];
					[Input(PM[m], md_0[m,i]) for m in margin if md_0[m,i]>0]
				]
				)
	end

	add!(WiNnat, DemandFunction(RA, 1.,
	# add!(WiNnat, DemandFunction(RA, :($(d_elas_ra)*1),
		[Demand(PA[i], fd_0[i,:pce]) for i in a_],
		[
			[Endowment(PY[i], fs_0[i]) for i in a_];
			[Endowment(PA[i], -sum(fd_0[i,xfd])) for i in a_];  
			[Endowment(PVA[va], sum(va_0[va,sectorsi])) for va in valueadded];
			Endowment(PFX, bopdef_0)
		]
		))

# End of optional function version
	# return WiNnat
# end

# @time solve!(WiNnat, cumulative_iteration_limit=0)
# MPSGE.build(WiNnat)


# Counterfactual solve
# set_value((A[(:gmt)]), 1.0)
# set_value((A[(:mvt)]), 1.0)
# set_value((A[(:fbt)]), 1.0)
# set_fixed!(A[(:gmt)], true)
# set_fixed!(A[(:mvt)], true)
# set_fixed!(A[(:fbt)], true)
# set_value((PA[(:gmt)]), 1.0)
# set_value((PA[(:mvt)]), 1.0)
# set_value((PA[(:fbt)]), 1.0)
# set_fixed!(PA[(:gmt)], true)
# set_fixed!(PA[(:mvt)], true)
# set_fixed!(PA[(:fbt)], true)

set_value(RA, 13138.7573)
set_fixed!(RA, true)

for i in a_
	set_value(ta[i], P[:ta_0][yr,a_][i])
	set_value(tm[i], P[:tm_0][yr,a_][i])
end
solve!(WiNnat, cumulative_iteration_limit=0)
println("$year ","bnchmk")

# Report = CSV.File(IOBuffer(generate_report(WiNnat._jump_model, mdecimals = 6)))
# CSV.write("FullCounterTest2023-11-22.csv", Report, missingstring="missing", bom=true)

# Counterfactual
for i in a_
	set_value(ta[i], 0.)
	set_value(tm[i], 0.)
end
# set_value(RA,  11846.425) #So far, this updated default normalization value needs to be set, value from GAMS output. 12453.8963
# set_value(RA,  12453.8963154469) #So far, this updated default normalization value needs to be set, value from GAMS output. 12453.8963
set_fixed!(RA, false)
# set_value(Y[:ppd],1.0187954); set_fixed!(Y[:ppd], false)
# set_value(Y[:res],1.03916451); set_fixed!(Y[:res], false)
# set_value(Y[:com],0.99921351); set_fixed!(Y[:com], false)
# set_value(Y[:amb],0.96924161); set_fixed!(Y[:amb], false)
# set_value(Y[:fbp],1.04401988); set_fixed!(Y[:fbp], false)
# set_value(Y[:rec],1.02557666); set_fixed!(Y[:rec], false)
# set_value(Y[:con],0.99872785); set_fixed!(Y[:con], false)
# set_value(Y[:agr],1.02650938); set_fixed!(Y[:agr], false)
# set_value(Y[:eec],0.99342307); set_fixed!(Y[:eec], false)
# set_value(Y[:fnd],1.); set_fixed!(Y[:fnd], false)
# set_value(Y[:pub],0.99505822); set_fixed!(Y[:pub], false)
# set_value(Y[:hou],0.94669728); set_fixed!(Y[:hou], false)
# set_value(Y[:fbt],1.02333194); set_fixed!(Y[:fbt], false)
# set_value(Y[:ins],0.99526435); set_fixed!(Y[:ins], false)
# set_value(Y[:tex],0.98775501); set_fixed!(Y[:tex], false)
# set_value(Y[:leg],1.00528428); set_fixed!(Y[:leg], false)
# set_value(Y[:fen],1.00420712); set_fixed!(Y[:fen], false)
# set_value(Y[:uti],1.0281457); set_fixed!(Y[:uti], false)
# set_value(Y[:nmp],0.99766877); set_fixed!(Y[:nmp], false)
# set_value(Y[:brd],1.02314763); set_fixed!(Y[:brd], false)
# set_value(Y[:bnk],0.98197656); set_fixed!(Y[:bnk], false)
# set_value(Y[:ore],1.00433844); set_fixed!(Y[:ore], false)
# set_value(Y[:edu],0.96132564); set_fixed!(Y[:edu], false)
# set_value(Y[:ote],1.00279349); set_fixed!(Y[:ote], false)
# set_value(Y[:man],1.01575578); set_fixed!(Y[:man], false)
# set_value(Y[:mch],1.00583442); set_fixed!(Y[:mch], false)
# set_value(Y[:dat],0.9970277); set_fixed!(Y[:dat], false)
# set_value(Y[:amd],1.05731717); set_fixed!(Y[:amd], false)
# set_value(Y[:oil],1.076019); set_fixed!(Y[:oil], false)
# set_value(Y[:hos],0.96974167); set_fixed!(Y[:hos], false)
# set_value(Y[:rnt],1.0200621); set_fixed!(Y[:rnt], false)
# set_value(Y[:pla],1.00803412); set_fixed!(Y[:pla], false)
# set_value(Y[:fof],1.01200162); set_fixed!(Y[:fof], false)
# set_value(Y[:fin],0.97302438); set_fixed!(Y[:fin], false)
# set_value(Y[:tsv],1.00152715); set_fixed!(Y[:tsv], false)
# set_value(Y[:nrs],0.98639151); set_fixed!(Y[:nrs], false)
# set_value(Y[:sec],0.98179106); set_fixed!(Y[:sec], false)
# set_value(Y[:art],1.00647162); set_fixed!(Y[:art], false)
# set_value(Y[:mov],1.00701579); set_fixed!(Y[:mov], false)
# set_value(Y[:fpd],1.01933465); set_fixed!(Y[:fpd], false)
# set_value(Y[:slg],1.); set_fixed!(Y[:slg], false)
# set_value(Y[:pri],1.00376859); set_fixed!(Y[:pri], false)
# set_value(Y[:grd],0.99154134); set_fixed!(Y[:grd], false)
# set_value(Y[:pip],1.02794336); set_fixed!(Y[:pip], false)
# set_value(Y[:sle],0.99645864); set_fixed!(Y[:sle], false)
# set_value(Y[:osv],0.99268578); set_fixed!(Y[:osv], false)
# set_value(Y[:trn],1.02187139); set_fixed!(Y[:trn], false)
# set_value(Y[:smn],0.97047217); set_fixed!(Y[:smn], false)
# set_value(Y[:fmt],1.00183519); set_fixed!(Y[:fmt], false)
# set_value(Y[:pet],1.08463382); set_fixed!(Y[:pet], false)
# set_value(Y[:mvt],1.02309901); set_fixed!(Y[:mvt], false)
# set_value(Y[:cep],0.98526698); set_fixed!(Y[:cep], false)
# set_value(Y[:wst],1.00291955); set_fixed!(Y[:wst], false)
# set_value(Y[:mot],1.02404048); set_fixed!(Y[:mot], false)
# set_value(Y[:adm],1.00240441); set_fixed!(Y[:adm], false)
# set_value(Y[:soc],0.97751794); set_fixed!(Y[:soc], false)
# set_value(Y[:alt],0.84929079); set_fixed!(Y[:alt], false)
# set_value(Y[:pmt],1.01858673); set_fixed!(Y[:pmt], false)
# set_value(Y[:trk],1.02647471); set_fixed!(Y[:trk], false)
# set_value(Y[:fdd],1.); set_fixed!(Y[:fdd], false)
# set_value(Y[:gmt],1.02311042); set_fixed!(Y[:gmt], false)
# set_value(Y[:wtt],1.01549527); set_fixed!(Y[:wtt], false)
# set_value(Y[:wpd],1.00651815); set_fixed!(Y[:wpd], false)
# set_value(Y[:wht],1.02303098); set_fixed!(Y[:wht], false)
# set_value(Y[:wrh],1.01943032); set_fixed!(Y[:wrh], false)
# set_value(Y[:ott],1.02349455); set_fixed!(Y[:ott], false)
# set_value(Y[:che],1.00544525); set_fixed!(Y[:che], false)
# set_value(Y[:air],1.08534791); set_fixed!(Y[:air], false)
# set_value(Y[:mmf],0.99696959); set_fixed!(Y[:mmf], false)
# set_value(Y[:otr],1.02245485); set_fixed!(Y[:otr], false)
# set_value(Y[:min],1.0168042); set_fixed!(Y[:min], false)
# set_value(A[:ppd],1.01598686); set_fixed!(A[:ppd], false)
# set_value(A[:res],1.03671031); set_fixed!(A[:res], false)
# set_value(A[:com],1.00094989); set_fixed!(A[:com], false)
# set_value(A[:amb],0.97034602); set_fixed!(A[:amb], false)
# set_value(A[:fbp],1.04186461); set_fixed!(A[:fbp], false)
# set_value(A[:rec],1.02538056); set_fixed!(A[:rec], false)
# set_value(A[:con],0.99850451); set_fixed!(A[:con], false)
# set_value(A[:agr],1.02359577); set_fixed!(A[:agr], false)
# set_value(A[:eec],1.0072111); set_fixed!(A[:eec], false)
# set_value(A[:fnd],1.); set_fixed!(A[:fnd], false)
# set_value(A[:pub],0.99539065); set_fixed!(A[:pub], false)
# set_value(A[:hou],0.94734024); set_fixed!(A[:hou], false)
# set_value(A[:ins],0.99501231); set_fixed!(A[:ins], false)
# set_value(A[:tex],1.02455568); set_fixed!(A[:tex], false)
# set_value(A[:leg],1.0053916); set_fixed!(A[:leg], false)
# set_value(A[:fen],1.0041957); set_fixed!(A[:fen], false)
# set_value(A[:uti],1.02048042); set_fixed!(A[:uti], false)
# set_value(A[:nmp],1.00572162); set_fixed!(A[:nmp], false)
# set_value(A[:brd],1.02313972); set_fixed!(A[:brd], false)
# set_value(A[:bnk],0.98224769); set_fixed!(A[:bnk], false)
# set_value(A[:ore],1.00425055); set_fixed!(A[:ore], false)
# set_value(A[:edu],0.97190425); set_fixed!(A[:edu], false)
# set_value(A[:ote],1.00325542); set_fixed!(A[:ote], false)
# set_value(A[:man],1.01575578); set_fixed!(A[:man], false)
# set_value(A[:mch],1.00647437); set_fixed!(A[:mch], false)
# set_value(A[:dat],0.99709514); set_fixed!(A[:dat], false)
# set_value(A[:amd],1.04227943); set_fixed!(A[:amd], false)
# set_value(A[:oil],1.07357625); set_fixed!(A[:oil], false)
# set_value(A[:hos],0.97584341); set_fixed!(A[:hos], false)
# set_value(A[:rnt],1.0129865); set_fixed!(A[:rnt], false)
# set_value(A[:pla],1.01315757); set_fixed!(A[:pla], false)
# set_value(A[:fof],1.01507934); set_fixed!(A[:fof], false)
# set_value(A[:fin],0.97723212); set_fixed!(A[:fin], false)
# set_value(A[:tsv],1.00197671); set_fixed!(A[:tsv], false)
# set_value(A[:nrs],0.98630195); set_fixed!(A[:nrs], false)
# set_value(A[:sec],0.98184695); set_fixed!(A[:sec], false)
# set_value(A[:art],1.00644844); set_fixed!(A[:art], false)
# set_value(A[:mov],1.00730735); set_fixed!(A[:mov], false)
# set_value(A[:fpd],1.01083025); set_fixed!(A[:fpd], false)
# set_value(A[:slg],1.); set_fixed!(A[:slg], false)
# set_value(A[:pri],1.00264223); set_fixed!(A[:pri], false)
# set_value(A[:grd],0.99257781); set_fixed!(A[:grd], false)
# set_value(A[:pip],1.06222612); set_fixed!(A[:pip], false)
# set_value(A[:sle],0.99699646); set_fixed!(A[:sle], false)
# set_value(A[:osv],0.99904312); set_fixed!(A[:osv], false)
# set_value(A[:trn],0.98696555); set_fixed!(A[:trn], false)
# set_value(A[:smn],1.00731963); set_fixed!(A[:smn], false)
# set_value(A[:fmt],1.00894214); set_fixed!(A[:fmt], false)
# set_value(A[:pet],1.07877535); set_fixed!(A[:pet], false)
# set_value(A[:cep],0.9986603); set_fixed!(A[:cep], false)
# set_value(A[:wst],1.00298117); set_fixed!(A[:wst], false)
# set_value(A[:mot],1.01595927); set_fixed!(A[:mot], false)
# set_value(A[:adm],1.00241194); set_fixed!(A[:adm], false)
# set_value(A[:soc],0.97811812); set_fixed!(A[:soc], false)
# set_value(A[:alt],1.07760725); set_fixed!(A[:alt], false)
# set_value(A[:pmt],1.01403892); set_fixed!(A[:pmt], false)
# set_value(A[:trk],1.0181068); set_fixed!(A[:trk], false)
# set_value(A[:fdd],1.); set_fixed!(A[:fdd], false)
# set_value(A[:wtt],1.00693768); set_fixed!(A[:wtt], false)
# set_value(A[:wpd],1.0063465); set_fixed!(A[:wpd], false)
# set_value(A[:wht],1.01980127); set_fixed!(A[:wht], false)
# set_value(A[:wrh],1.01969958); set_fixed!(A[:wrh], false)
# set_value(A[:ott],0.98433973); set_fixed!(A[:ott], false)
# set_value(A[:che],1.00690582); set_fixed!(A[:che], false)
# set_value(A[:air],1.07797921); set_fixed!(A[:air], false)
# set_value(A[:mmf],1.00711561); set_fixed!(A[:mmf], false)
# set_value(A[:otr],1.0210803); set_fixed!(A[:otr], false)
# set_value(A[:min],1.01630156); set_fixed!(A[:min], false)
# set_value(MS[:trn],1.02748421); set_fixed!(MS[:trn], false)
# set_value(MS[:trd],1.02278522); set_fixed!(MS[:trd], false)
# set_value(PA[:ppd],0.94563437); set_fixed!(PA[:ppd], false)
# set_value(PA[:res],0.90543868); set_fixed!(PA[:res], false)
# set_value(PA[:com],0.97436569); set_fixed!(PA[:com], false)
# set_value(PA[:amb],0.97709585); set_fixed!(PA[:amb], false)
# set_value(PA[:fbp],0.90558227); set_fixed!(PA[:fbp], false)
# set_value(PA[:rec],0.92394951); set_fixed!(PA[:rec], false)
# set_value(PA[:con],0.96193526); set_fixed!(PA[:con], false)
# set_value(PA[:agr],0.96906201); set_fixed!(PA[:agr], false)
# set_value(PA[:eec],0.93784143); set_fixed!(PA[:eec], false)
# set_value(PA[:fnd],0.98007818); set_fixed!(PA[:fnd], false)
# set_value(PA[:pub],0.95583419); set_fixed!(PA[:pub], false)
# set_value(PA[:hou],1.00056422); set_fixed!(PA[:hou], false)
# set_value(PA[:ins],0.95538235); set_fixed!(PA[:ins], false)
# set_value(PA[:tex],0.91272373); set_fixed!(PA[:tex], false)
# set_value(PA[:leg],0.93514899); set_fixed!(PA[:leg], false)
# set_value(PA[:fen],0.97748696); set_fixed!(PA[:fen], false)
# set_value(PA[:uti],0.91926868); set_fixed!(PA[:uti], false)
# set_value(PA[:nmp],0.93936815); set_fixed!(PA[:nmp], false)
# set_value(PA[:brd],0.9128476); set_fixed!(PA[:brd], false)
# set_value(PA[:bnk],0.97900839); set_fixed!(PA[:bnk], false)
# set_value(PA[:ore],0.9651338); set_fixed!(PA[:ore], false)
# set_value(PA[:edu],0.97861036); set_fixed!(PA[:edu], false)
# set_value(PA[:ote],0.96144191); set_fixed!(PA[:ote], false)
# set_value(PA[:man],0.97556216); set_fixed!(PA[:man], false)
# set_value(PA[:mch],0.94763907); set_fixed!(PA[:mch], false)
# set_value(PA[:dat],0.96795255); set_fixed!(PA[:dat], false)
# set_value(PA[:amd],0.89282254); set_fixed!(PA[:amd], false)
# set_value(PA[:oil],0.94498367); set_fixed!(PA[:oil], false)
# set_value(PA[:hos],0.97130593); set_fixed!(PA[:hos], false)
# set_value(PA[:rnt],0.93206828); set_fixed!(PA[:rnt], false)
# set_value(PA[:pla],0.93675232); set_fixed!(PA[:pla], false)
# set_value(PA[:fof],0.96616522); set_fixed!(PA[:fof], false)
# set_value(PA[:fin],0.97051374); set_fixed!(PA[:fin], false)
# set_value(PA[:tsv],0.97425832); set_fixed!(PA[:tsv], false)
# set_value(PA[:nrs],0.96122427); set_fixed!(PA[:nrs], false)
# set_value(PA[:sec],0.97361367); set_fixed!(PA[:sec], false)
# set_value(PA[:art],0.94578077); set_fixed!(PA[:art], false)
# set_value(PA[:mov],0.95053652); set_fixed!(PA[:mov], false)
# set_value(PA[:fpd],0.92686303); set_fixed!(PA[:fpd], false)
# set_value(PA[:slg],0.97120258); set_fixed!(PA[:slg], false)
# set_value(PA[:pri],0.95332269); set_fixed!(PA[:pri], false)
# set_value(PA[:grd],0.96118423); set_fixed!(PA[:grd], false)
# set_value(PA[:pip],0.77986139); set_fixed!(PA[:pip], false)
# set_value(PA[:sle],0.95551988); set_fixed!(PA[:sle], false)
# set_value(PA[:osv],0.95176253); set_fixed!(PA[:osv], false)
# set_value(PA[:trn],1.14242988); set_fixed!(PA[:trn], false)
# set_value(PA[:smn],0.96675554); set_fixed!(PA[:smn], false)
# set_value(PA[:fmt],0.94700924); set_fixed!(PA[:fmt], false)
# set_value(PA[:pet],0.82092499); set_fixed!(PA[:pet], false)
# set_value(PA[:cep],0.95765729); set_fixed!(PA[:cep], false)
# set_value(PA[:wst],0.95400451); set_fixed!(PA[:wst], false)
# set_value(PA[:mot],0.92995882); set_fixed!(PA[:mot], false)
# set_value(PA[:adm],0.96772238); set_fixed!(PA[:adm], false)
# set_value(PA[:soc],0.96915909); set_fixed!(PA[:soc], false)
# set_value(PA[:alt],0.87606611); set_fixed!(PA[:alt], false)
# set_value(PA[:pmt],0.95970015); set_fixed!(PA[:pmt], false)
# set_value(PA[:trk],0.94204935); set_fixed!(PA[:trk], false)
# set_value(PA[:fdd],0.97539708); set_fixed!(PA[:fdd], false)
# set_value(PA[:wtt],0.94749965); set_fixed!(PA[:wtt], false)
# set_value(PA[:wpd],0.94499167); set_fixed!(PA[:wpd], false)
# set_value(PA[:wht],0.97501489); set_fixed!(PA[:wht], false)
# set_value(PA[:wrh],0.96857117); set_fixed!(PA[:wrh], false)
# set_value(PA[:ott],0.97525797); set_fixed!(PA[:ott], false)
# set_value(PA[:che],0.94464432); set_fixed!(PA[:che], false)
# set_value(PA[:air],0.85437139); set_fixed!(PA[:air], false)
# set_value(PA[:mmf],0.93363776); set_fixed!(PA[:mmf], false)
# set_value(PA[:otr],0.96096079); set_fixed!(PA[:otr], false)
# set_value(PA[:min],0.9334001); set_fixed!(PA[:min], false)
# set_value(PY[:ppd],0.95786829); set_fixed!(PY[:ppd], false)
# set_value(PY[:res],0.96851404); set_fixed!(PY[:res], false)
# set_value(PY[:com],0.98199622); set_fixed!(PY[:com], false)
# set_value(PY[:amb],0.97757205); set_fixed!(PY[:amb], false)
# set_value(PY[:fbp],0.95537463); set_fixed!(PY[:fbp], false)
# set_value(PY[:rec],0.9694208); set_fixed!(PY[:rec], false)
# set_value(PY[:con],0.96306156); set_fixed!(PY[:con], false)
# set_value(PY[:agr],0.95842558); set_fixed!(PY[:agr], false)
# set_value(PY[:eec],0.96761881); set_fixed!(PY[:eec], false)
# set_value(PY[:fnd],0.98007818); set_fixed!(PY[:fnd], false)
# set_value(PY[:pub],0.97987017); set_fixed!(PY[:pub], false)
# set_value(PY[:hou],0.98212598); set_fixed!(PY[:hou], false)
# set_value(PY[:fbt],0.97687732); set_fixed!(PY[:fbt], false)
# set_value(PY[:ins],0.97110252); set_fixed!(PY[:ins], false)
# set_value(PY[:tex],0.956419); set_fixed!(PY[:tex], false)
# set_value(PY[:leg],0.97904606); set_fixed!(PY[:leg], false)
# set_value(PY[:fen],0.97746561); set_fixed!(PY[:fen], false)
# set_value(PY[:uti],0.96287344); set_fixed!(PY[:uti], false)
# set_value(PY[:nmp],0.96390791); set_fixed!(PY[:nmp], false)
# set_value(PY[:brd],0.96273768); set_fixed!(PY[:brd], false)
# set_value(PY[:bnk],0.9780216); set_fixed!(PY[:bnk], false)
# set_value(PY[:ore],0.96458078); set_fixed!(PY[:ore], false)
# set_value(PY[:edu],0.97927884); set_fixed!(PY[:edu], false)
# set_value(PY[:ote],0.97034296); set_fixed!(PY[:ote], false)
# set_value(PY[:man],0.97807942); set_fixed!(PY[:man], false)
# set_value(PY[:mch],0.96595707); set_fixed!(PY[:mch], false)
# set_value(PY[:dat],0.97130917); set_fixed!(PY[:dat], false)
# set_value(PY[:amd],0.97021248); set_fixed!(PY[:amd], false)
# set_value(PY[:oil],0.97018593); set_fixed!(PY[:oil], false)
# set_value(PY[:hos],0.97289506); set_fixed!(PY[:hos], false)
# set_value(PY[:rnt],0.97590621); set_fixed!(PY[:rnt], false)
# set_value(PY[:pla],0.95899025); set_fixed!(PY[:pla], false)
# set_value(PY[:fof],0.97627359); set_fixed!(PY[:fof], false)
# set_value(PY[:fin],0.97071421); set_fixed!(PY[:fin], false)
# set_value(PY[:tsv],0.97667404); set_fixed!(PY[:tsv], false)
# set_value(PY[:nrs],0.97276553); set_fixed!(PY[:nrs], false)
# set_value(PY[:sec],0.97581762); set_fixed!(PY[:sec], false)
# set_value(PY[:art],0.97695408); set_fixed!(PY[:art], false)
# set_value(PY[:mov],0.97459103); set_fixed!(PY[:mov], false)
# set_value(PY[:fpd],0.96162601); set_fixed!(PY[:fpd], false)
# set_value(PY[:slg],0.97120258); set_fixed!(PY[:slg], false)
# set_value(PY[:pri],0.96262109); set_fixed!(PY[:pri], false)
# set_value(PY[:grd],0.97169525); set_fixed!(PY[:grd], false)
# set_value(PY[:pip],0.98129409); set_fixed!(PY[:pip], false)
# set_value(PY[:sle],0.95551988); set_fixed!(PY[:sle], false)
# set_value(PY[:osv],0.97432865); set_fixed!(PY[:osv], false)
# set_value(PY[:trn],0.96153926); set_fixed!(PY[:trn], false)
# set_value(PY[:smn],0.96686952); set_fixed!(PY[:smn], false)
# set_value(PY[:fmt],0.96954022); set_fixed!(PY[:fmt], false)
# set_value(PY[:pet],0.95076004); set_fixed!(PY[:pet], false)
# set_value(PY[:mvt],0.97812318); set_fixed!(PY[:mvt], false)
# set_value(PY[:cep],0.98341694); set_fixed!(PY[:cep], false)
# set_value(PY[:wst],0.96838549); set_fixed!(PY[:wst], false)
# set_value(PY[:mot],0.9503919); set_fixed!(PY[:mot], false)
# set_value(PY[:adm],0.97437075); set_fixed!(PY[:adm], false)
# set_value(PY[:soc],0.97291565); set_fixed!(PY[:soc], false)
# set_value(PY[:alt],0.95896045); set_fixed!(PY[:alt], false)
# set_value(PY[:pmt],0.96219127); set_fixed!(PY[:pmt], false)
# set_value(PY[:trk],0.95299057); set_fixed!(PY[:trk], false)
# set_value(PY[:fdd],0.97539708); set_fixed!(PY[:fdd], false)
# set_value(PY[:gmt],0.97774744); set_fixed!(PY[:gmt], false)
# set_value(PY[:wtt],0.95750015); set_fixed!(PY[:wtt], false)
# set_value(PY[:wpd],0.96244699); set_fixed!(PY[:wpd], false)
# set_value(PY[:wht],0.97498827); set_fixed!(PY[:wht], false)
# set_value(PY[:wrh],0.9697362); set_fixed!(PY[:wrh], false)
# set_value(PY[:ott],0.97525797); set_fixed!(PY[:ott], false)
# set_value(PY[:che],0.96176982); set_fixed!(PY[:che], false)
# set_value(PY[:air],0.95238392); set_fixed!(PY[:air], false)
# set_value(PY[:mmf],0.97032749); set_fixed!(PY[:mmf], false)
# set_value(PY[:otr],0.96283448); set_fixed!(PY[:otr], false)
# set_value(PY[:min],0.96259614); set_fixed!(PY[:min], false)
# set_value(PVA[:compen],0.99159958); set_fixed!(PVA[:compen], false)
# set_value(PVA[:surplus],0.98420975); set_fixed!(PVA[:surplus], false)
# set_value(PM[:trn],0.95763115); set_fixed!(PM[:trn], false)
# set_value(PM[:trd],0.97554125); set_fixed!(PM[:trd], false)
# set_value(PFX,0.97385956); set_fixed!(PFX, false)
# set_value(RA,12453.8963154469); set_fixed!(RA, false)
solve!(WiNnat, cumulative_iteration_limit=10000)
println("$year ","counter")
# end
#, convergence_tolerance=1e-3);
# @time timetest()
# timearray=[]
# for in in 1:n
# tmel = @elapsed timetest()
# push!(timearray, tmel)
# end
# println(timearray, sum(timearray/n))
# @profview solve!(WiNnat)
# @time MPSGE.build(WiNnat);
# @elapsed solve!(WiNnat, cumulative_iteration_limit=0)

##  Write the full algebraic model to a file for viewing
	# open("WiNnat_Algebraic2.txt", "w") do file
	# 	show(file, algebraic_version(WiNnat))
	# end

	# open("Report.txt", "w") do file
	# 	write(file, generate_report(WiNnat._jump_model))
	# end

	# Report = CSV.File(IOBuffer(generate_report(WiNnat._jump_model, decimals=8, mdecimals=6)));
	# CSV.write("./Results/FullReportCounterTazpoint1PJfd_0"*Dates.format(now(),"yyyy-mm-dd,H-M")*".csv", (Report), missingstring="missing", bom=true)
	# println("done ",value(WiNnat._jump_model[:RA]))
# ## For testing with variable numbers of sectors	
# 	# timeWiNnat()
# 	# [@elapsed timeWiNnat(t) for t in [2 2 2 8 16]]
# 	# []@time timeWiNnat(t) for t in [2 2 2 8 16]]
# 	# [@time timeWiNnat(t) for t in [2 2 8 73]]
# 	# @profview timeWiNnat(111)

# 	###  SENSITIVITY TESTS
# ## First, save counterfactual baseline. Run this once, save and then comment out so it doesn't overwrite with values from random sampling
# # basecountervalues = [:t_elas_y           :elas_y           :elas_va           :t_elas_m           :elas_m           :t_elas_a           :elas_a           :elas_dm           :d_elas_ra              transpose(all_variables(WiNnat._jump_model)[[152:295;435:578],:])
# 			# get_value(t_elas_y) get_value(elas_y) get_value(elas_va) get_value(t_elas_m) get_value(elas_m) get_value(t_elas_a) get_value(elas_a) get_value(elas_dm) get_value(d_elas_ra) transpose(JuMP.value.(all_variables(WiNnat._jump_model)[[152:295;435:578],:]))]
# # CSV.write("CounterVals.csv", Tables.table(basecountervalues), missingstring="missing", bom=true)

# ## Before 1st run, set and write headers and baseline values into CSV for sensitivity test data output
# # basecountervalues = DataFrame(CSV.File(joinpath(@__DIR__,"CounterVals.csv"), header = 2, skipto=3))
# # CSV.write("ExpectedVal.csv", Tables.table(basecountervalues), missingstring="missing", bom=true)

# # Visualise Distributions by taking a random sample and plotting
# # histogram(rand(Gamma(1,2),10^6), label = "Gamma(1,2)")

# ## For this sensetivity test, we alter all elasticities simultaneously for each Monte Carlo run- except demand (any alternate value for demand within the pdf pushed all variable values back to the benchmark).
# function ExpectedVal(n)
# 	starttime = now()
# 	out = Array{Float64}(undef, n, 297) # 9 elasticity values plus the 288 variables I'll report
# 	for i in 1:n
# 		t_elas_y_val = rand(Gamma(1,2))
# 		elas_y_val = rand(Gamma(1,2)) 
# 		elas_va_val = rand(truncated(Normal(1,2),lower=0)) #1
# 		t_elas_m_val = rand(Gamma(1,2)) 
# 		elas_m_val = rand(Gamma(1,2)) 
# 		t_elas_a_val = rand(truncated(Normal(2,2),lower=0)) #2
# 		elas_a_val = rand(Gamma(1,2)) 
# 		elas_dm_val = rand(truncated(Normal(2,2),lower=0)) #2
# 		# d_elas_ra_val = rand(truncated(Normal(1,2),lower=0)) #1
# 		set_value(d_elas_ra, 1.)# Any value other than 1 generates benchmark results, re-set to 1 to make sure.
	
# 		set_value(t_elas_y, t_elas_y_val)
# 		set_value(elas_y, elas_y_val)
# 		set_value(elas_va, elas_va_val)
# 		set_value(t_elas_m, t_elas_m_val)
# 		set_value(elas_m, elas_m_val)
# 		set_value(t_elas_a, t_elas_a_val)
# 		set_value(elas_a, elas_a_val)
# 		set_value(elas_dm, elas_dm_val)
# 		# set_value(d_elas_ra, d_elas_ra_val)

# 		solve!(WiNnat, cumulative_iteration_limit=10000);#, convergence_tolerance=1e-0);
# 		# save sampled elasticity values, and indices for Y through to MS, and PA through PFX in results, variables we're most interested in
# 	out[i,:] = [get_value(t_elas_y) get_value(elas_y) get_value(elas_va) get_value(t_elas_m) get_value(elas_m) get_value(t_elas_a) get_value(elas_a) get_value(elas_dm) get_value(d_elas_ra) transpose(JuMP.value.(all_variables(WiNnat._jump_model)[[152:295;435:578],:]))]
# 	println(i,": So far ",now()-starttime)
# 	end
# 	println("It took",now()-starttime)
# 	return out
# end

# ## Run Monte Carlo n times, and append to csv (can run many small batches)
# # outall = ExpectedVal(10)
# ## Write the results from the sensitivity test to csv
# # CSV.write("./SensitivityTests/ExpectedVal.csv", Tables.table(outall), missingstring="missing", bom=true, append=true)

# ## get all results back for plotting
# # ExVal = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTests/ExpectedVal.csv"), header = 2, skipto=3))

# ## df data frame; fv first variable index from SensElas; lv last variable index, for loop, 21 variables at a time for plotting in 3x7
# function plotgroups(df, fv, lv)
# 	pltvarsAll = df[:,fv:lv]
# 	saveplots =[]
# 	# Loop over 21 variables for plotting
# 	for i in 1:length(pltvarsAll[1,:])
# 	# 	# Start with original value from the model, and set the axis ranges from the values (including the originals)
# 		plt =	scatter([pltvarsAll[1,i]],[pltvarsAll[1,i]], size=(600,1400),color = "green",
# 		 label = false, markersize = 9,
# 		ylimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(maximum(maximum(eachcol(pltvarsAll)))/8)),
# 		xlimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(minimum(minimum(eachcol(pltvarsAll)))/8)), 
# 		 title=replace(names(pltvarsAll)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>""))
# 		 plt =	plot!(pltvarsAll[2:end,i],pltvarsAll[2:end,i], seriestype=:scatter,
# 				label=false, xaxis = ("",font(7)), yaxis=("",font(7)),
# 				titlefontsize=10, markersize=3, markerstrokewidth=0, color="red", alpha=0.4)
# 		# display(plt)
# 		# sleep(.5)
# 		push!(saveplots, plt)
# 	end
# 	pltAll = plot((saveplots...), layout=(7,3), size=(900,1000),  left_margin=1*Plots.mm,right_margin=1*Plots.mm,
# 	 bottom_margin=-1*Plots.mm, top_margin=1*Plots.mm, plot_title="Counterfactual Base Value,\nand the impact on each variable from varying all elasticities",
# 	  plot_titlefontsize=9, plot_titlefonthalign=:left, fontfamily="Palatino Roman")
# 	  display(pltAll)
# 	savefig(replace("EvVal"*names(ExVal)[fv]*"-"*names(ExVal)[lv]*".png",":" => ""))
# end
# # plotgroups(ExVal, 10,30) # fv = 10; lv=30; # Y
# # plotgroups(ExVal, 31,51) # Y
# # plotgroups(ExVal, 52,72) # Y
# # plotgroups(ExVal, 73,80) # Y
# # plotgroups(ExVal, 81,101) # A
# # plotgroups(ExVal, 102,122) # A
# # plotgroups(ExVal, 123,143) # A
# # plotgroups(ExVal, 144,151) # A
# # plotgroups(ExVal, 152,153) # MS
# # plotgroups(ExVal, 154,174) # PA
# # plotgroups(ExVal, 175,195) # PA
# # plotgroups(ExVal, 196,216) # PA
# # plotgroups(ExVal, 217,221) # PA
# # plotgroups(ExVal, 222,242) # PY
# # plotgroups(ExVal, 243,263) # PY
# # plotgroups(ExVal, 264,284) # PY
# # plotgroups(ExVal, 285,292) # PY
# # plotgroups(ExVal, 293,297) # PVA -> PFX


# ## Selection of variables to show as a sample in one collection
# # plotgroups(ExVal,10,21,24,66,81,92,95,137,154,218,212,222,289,283,293,294,295,296,152,153,297)
# ## I didn't figure out to make it work using the above function, so just set it separately
# pltvarsAll = ExVal[:,[10,21,24,66,81,92,95,137,154,218,212,222,289,283,293,294,295,296,152,153,297]]
# 	saveplots =[]
# 	for i in 1:length(pltvarsAll[1,:])
# 		plt =	scatter([pltvarsAll[1,i]],[pltvarsAll[1,i]], size=(600,1400),color = "green",
# 		 label = false, markersize = 9,
# 		ylimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(maximum(maximum(eachcol(pltvarsAll)))/8)),
# 		xlimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(minimum(minimum(eachcol(pltvarsAll)))/8)), 
# 		 title=replace(names(pltvarsAll)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>""))
# 		 plt =	plot!(pltvarsAll[2:end,i],pltvarsAll[2:end,i], seriestype=:scatter,
# 				label=false, xaxis = ("",font(7)), yaxis=("",font(7)),
# 				titlefontsize=10, markersize=3, markerstrokewidth=0, color="red", alpha=0.4)
# 		# display(plt)
# 				# sleep(.5)
# 		push!(saveplots, plt)
# 	end
# 	pltAll = plot((saveplots...), layout=(7,3), size=(900,1000),  left_margin=1*Plots.mm,right_margin=1*Plots.mm,
# 	 bottom_margin=-1*Plots.mm, top_margin=1*Plots.mm, plot_title="Counterfactual Base Value,\nand the impact on each variable from varying all elasticities",
# 	  plot_titlefontsize=9, plot_titlefonthalign=:left, fontfamily="Palatino Roman")
# 	  display(pltAll)
# 	savefig("ExValSupSet.png")

# ## SENSETIVITY TEST FOR INDIVIDUAL ELASTICITY PARAMETERS
# ## First, function to set header and baseline counterfactual values rows.
# function SensitivityTestHead(param)
# 	out = Array{Any}(undef, 2, 289)
# 	set_value(t_elas_y, 0.)
# 	set_value(elas_y, 0.)
# 	set_value(elas_va, 1.)
# 	set_value(t_elas_m, 0.)
# 	set_value(elas_m, 0.)
# 	set_value(t_elas_a, 2.)
# 	set_value(elas_a, 0.)
# 	set_value(elas_dm, 2.)
# 	set_value(d_elas_ra, 1.)
# 	solve!(WiNnat, cumulative_iteration_limit=10000);
# 	paramname = param.model._parameters[param.index].name
# 	out[1,:] = [paramname        transpose(all_variables(WiNnat._jump_model)[[152:295;435:578],:])]
# 	out[2,:] = [get_value(param)  					 transpose(JuMP.value.(all_variables(WiNnat._jump_model)[[152:295;435:578],:]))]
# 	CSV.write("Sens_"*String(paramname)*".csv", Tables.table(out), missingstring="missing", bom=true)
# end

# ## Here, we set a function to select a single elasticity and vary that according to the distribution argument, n times in each Monte Carlo run, and append values to csv
# function SensitivityTest(param, n, dist)
# 	out = Array{Float64}(undef, n, 289)
# 	paramname = param.model._parameters[param.index].name
# 	set_value(t_elas_y, 0.)
# 	set_value(elas_y, 0.)
# 	set_value(elas_va, 1.)
# 	set_value(t_elas_m, 0.)
# 	set_value(elas_m, 0.)
# 	set_value(t_elas_a, 2.)
# 	set_value(elas_a, 0.)
# 	set_value(elas_dm, 2.)
# 	set_value(d_elas_ra, 1.)
# 	for i in 1:n
# 		pdr = rand(dist)
# 		set_value(param, pdr)
# 		println(i,"/",n,". ",param)
# 		solve!(WiNnat, cumulative_iteration_limit=10000);#, convergence_tolerance=1e-0);
# 	##indexes for Y through to MS, and PA through PFX in results, variables we're most interested in
# 		out[i,:] = [get_value(param)  transpose(JuMP.value.(all_variables(WiNnat._jump_model)[[152:295;435:578],:]))]
# 	end
# 	CSV.write("Sens_"*String(paramname)*".csv", Tables.table(out), missingstring="missing", bom=true, append=true)
# 	# return out
# end

# ## Run heading function once, and then can run ST in batches
# # Example: for the first elasticity parameter
# # SensitivityTestHead(t_elas_y)
# # SensitivityTest(t_elas_y, 5, Gamma(1,2))


# ### Get All Sensitivity results for plotting together
# # SenseoutputT_Y = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_t_elas_y.csv"), header = 2, skipto=3))
# # SenseoutputY = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_elas_y.csv"), header = 2, skipto=3))
# # SenseoutputVA = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_elas_va.csv"), header = 2, skipto=3))
# # SenseoutputT_M = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_t_elas_m.csv"), header = 2, skipto=3))
# # SenseoutputM = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_elas_m.csv"), header = 2, skipto=3))
# # SenseoutputT_A = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_t_elas_a.csv"), header = 2, skipto=3))
# # SenseoutputA = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_elas_a.csv"), header = 2, skipto=3))
# # SenseoutputDM = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_elas_dm.csv"), header = 2, skipto=3))
# # SenseoutputD_RA = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest/Sens_d_elas_ra.csv"), header = 2, skipto=3))

# # # #fv first variable index from SensElas; lv last variable index, 21 variables, for plotting loop and 3x7 plot collection
# function SensePlot(fv, lv)
# 	pltvarsT_Y = [SenseoutputT_Y[:,:t_elas_y] SenseoutputT_Y[:,fv:lv]]
# 	pltvarsY = [SenseoutputY[:,:elas_y] SenseoutputY[:,fv:lv]]
# 	pltvarsVA = [SenseoutputVA[:,:elas_va] SenseoutputVA[:,fv:lv]] #1
# 	pltvarsT_M = [SenseoutputT_M[:,:t_elas_m] SenseoutputT_M[:,fv:lv]]
# 	pltvarsM = [SenseoutputM[:,:elas_m] SenseoutputM[:,fv:lv]] 
# 	pltvarsT_A = [SenseoutputT_A[:,:t_elas_a] SenseoutputT_A[:,fv:lv]] #2
# 	pltvarsA = [SenseoutputA[:,:elas_a] SenseoutputA[:,fv:lv]] 
# 	pltvarsDM = [SenseoutputDM[:,:elas_dm] SenseoutputDM[:,fv:lv]] #2
# 	pltvarsD_RA = [SenseoutputD_RA[:,:d_elas_ra] SenseoutputD_RA[:,fv:lv]] #1
# 	groupplots =[]
# 	plt = scatter(1,markersize=12, markerstrokewidth=2, size=(300,400), xaxis=false, yaxis=false,
# 	label="   Baseline\n   Counterfactual\n   Value", color="green", legendfontsize=9,legend_title="Elasticity\n  Parameters",
# 	legend_title_font_size=11, fg_legend=:transparent,labelspacing=0.8, framestyle=:none,
# 	fontfamily="Computer Modern")
# 	plot!(1,1, label="   t_elas_y", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightgreen")
# 	plot!(1,1, label="   elas_y", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange")
# 	plot!(1,1, label="   elas_va", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55")
# 	plot!(1,1, label="   t_elas_m", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightblue")
# 	plot!(1,1, label="   elas_m", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue")
# 	plot!(1,1, label="   t_elas_a", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple")
# 	plot!(1,1, label="   elas_a", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red" )
# 	plot!(1,1, label="   elas_dm", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow")
# 	plot!(1,1, label="   d_elas_ra", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="cyan") # Any value other than 1 makes all values 1...
# 	push!(groupplots, plt)
# 	for i in 2:length(pltvarsA[1,:])
# 	# 	# Start with original value from the model
# 		plt =	scatter([pltvarsM[1,1]],[pltvarsM[1,i]], color = "darkgreen", label = false, markersize = 9, 
# 		ylimits=(min(pltvarsA[1,i],minimum(pltvarsA[2:end,i]),minimum(pltvarsM[2:end,i]),minimum(pltvarsY[2:end,i]),minimum(pltvarsDM[2:end,i]),minimum(pltvarsVA[2:end,i]),minimum(pltvarsT_A[2:end,i]),minimum(pltvarsT_M[2:end,i]),minimum(pltvarsT_Y[2:end,i]),minimum(pltvarsD_RA[2:end,i]))-.1,
# 		max(pltvarsA[1,i],maximum(pltvarsA[2:end,i]),maximum(pltvarsM[2:end,i]),maximum(pltvarsY[2:end,i]),maximum(pltvarsDM[2:end,i]),maximum(pltvarsVA[2:end,i]),maximum(pltvarsT_A[2:end,i]),maximum(pltvarsT_M[2:end,i]),maximum(pltvarsT_Y[2:end,i]),maximum(pltvarsD_RA[2:end,i]))+.1))
# 		plt =	plot!(pltvarsT_Y[2:end,i].-1,pltvarsT_Y[2:end,i], seriestype=:scatter, title=(replace(names(pltvarsT_Y)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>"")),
# 				label=false, xaxis = ("Change in each Elasticity of Substitution",font(7)), yaxis=("Diff in Sector activity",font(7)),# ylim=(minimum(Senseoutput[:,i]), maximum(Senseoutput[:,i])),
# 				titlefontsize=10, markersize=3, markerstrokewidth=0, color="lightgreen", alpha=0.3)
# 		plt =	plot!(pltvarsY[2:end,1],pltvarsY[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange", alpha=0.3)
# 		plt =	plot!(pltvarsVA[2:end,1].-1,pltvarsVA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55", alpha=0.3)
# 		plt =	plot!(pltvarsT_M[2:end,1],pltvarsT_M[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="light blue", alpha=0.3)
# 		plt =	plot!(pltvarsM[2:end,1],pltvarsM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue", alpha=0.3)
# 		plt =	plot!(pltvarsT_A[2:end,1].-2,pltvarsT_A[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple", alpha=0.3)
# 		plt =	plot!(pltvarsA[2:end,1],pltvarsA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red", alpha=0.3)
# 		plt =	plot!(pltvarsDM[2:end,1].-2,pltvarsDM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow", alpha=0.3)
# 		plt =	plot!(pltvarsD_RA[2:end,1].-1,pltvarsD_RA[2:end,i], label=false, seriestype=:scatter, markersize=2, markerstrokewidth=0, color="cyan", alpha=0.1) # Any value other than 1 makes all values 1...
# 			#   display(plt)		
# 	# sleep(.5)
# 		push!(groupplots, plt)
# 	end
# 	l = @layout([a{0.1w} [a d e; f g h ; i j k ; l m n; o p q; r s t; u v w]])
# 	plot((groupplots...), layout=l, size=(1000,1200), left_margin=2*Plots.mm, bottom_margin=-1*Plots.mm, top_margin=3*Plots.mm, show=true)
# 	savefig(replace("Sense_"*names(pltvarsT_Y)[2]*"-"*names(pltvarsT_Y)[end]*".png",":" => ""))
# end	

# # SensePlot(2, 22) # Y
# # SensePlot(23, 43) # Y
# # SensePlot(44, 64) # Y
# # SensePlot(65, 72) # Y
# # SensePlot(73, 93) # A
# # SensePlot(94, 114) # A
# # SensePlot(115, 135) # A
# # SensePlot(136, 143) # A
# # SensePlot(144, 145) # MS
# # SensePlot(146, 166) # PA
# # SensePlot(167, 187) # PA
# # SensePlot(188, 208) # PA
# # SensePlot(209, 213) # PA
# # SensePlot(214, 234) # PY
# # SensePlot(235, 255) # PY
# # SensePlot(256, 276) # PY
# # SensePlot(277, 284) # PY
# # SensePlot(285, 289) # PM -> PFX

# ## Function to plot non-contiguous set of 21 variables, to show as a sample
# function SensePlotInd(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)
# 	pltvarsT_Y = [SenseoutputT_Y[:,:t_elas_y] SenseoutputT_Y[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]]
# 	pltvarsY = [SenseoutputY[:,:elas_y] SenseoutputY[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]]
# 	pltvarsVA = [SenseoutputVA[:,:elas_va] SenseoutputVA[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #1
# 	pltvarsT_M = [SenseoutputT_M[:,:t_elas_m] SenseoutputT_M[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]]
# 	pltvarsM = [SenseoutputM[:,:elas_m] SenseoutputM[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] 
# 	pltvarsT_A = [SenseoutputT_A[:,:t_elas_a] SenseoutputT_A[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #2
# 	pltvarsA = [SenseoutputA[:,:elas_a] SenseoutputA[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] 
# 	pltvarsDM = [SenseoutputDM[:,:elas_dm] SenseoutputDM[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #2
# 	pltvarsD_RA = [SenseoutputD_RA[:,:d_elas_ra] SenseoutputD_RA[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #1
# 	# Generate a legend for the whole group of plots as a separate 'plot'
# 	groupplots =[]
# 	plt = scatter(1,markersize=12, markerstrokewidth=2, size=(300,400), xaxis=false, yaxis=false,
# 	label="   Baseline\n   Counterfactual\n   Value", color="green", legendfontsize=9,legend_title="Elasticity\n  Parameters",
# 	legend_title_font_size=11, fg_legend=:transparent,labelspacing=0.8, framestyle=:none,
# 	fontfamily="Computer Modern")
# 	plot!(1,1, label="   t_elas_y", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightgreen")
# 	plot!(1,1, label="   elas_y", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange")
# 	plot!(1,1, label="   elas_va", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55")
# 	plot!(1,1, label="   t_elas_m", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightblue")
# 	plot!(1,1, label="   elas_m", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue")
# 	plot!(1,1, label="   t_elas_a", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple")
# 	plot!(1,1, label="   elas_a", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red" )
# 	plot!(1,1, label="   elas_dm", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow")
# 	plot!(1,1, label="   d_elas_ra", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="cyan") # Any value other than 1 makes all values 1...
# 	push!(groupplots, plt)
# 	for i in 2:length(pltvarsA[1,:])
# 	# 	# Start with original value from the model
# 		plt =	scatter([pltvarsM[1,1]],[pltvarsM[1,i]], color = "darkgreen", label = false, markersize = 9, 
# 		# set limits from data
# 		ylimits=(min(pltvarsA[1,i],minimum(pltvarsA[2:end,i]),minimum(pltvarsM[2:end,i]),minimum(pltvarsY[2:end,i]),minimum(pltvarsDM[2:end,i]),minimum(pltvarsVA[2:end,i]),minimum(pltvarsT_A[2:end,i]),minimum(pltvarsT_M[2:end,i]),minimum(pltvarsT_Y[2:end,i]),minimum(pltvarsD_RA[2:end,i]))-.1,
# 		max(pltvarsA[1,i],maximum(pltvarsA[2:end,i]),maximum(pltvarsM[2:end,i]),maximum(pltvarsY[2:end,i]),maximum(pltvarsDM[2:end,i]),maximum(pltvarsVA[2:end,i]),maximum(pltvarsT_A[2:end,i]),maximum(pltvarsT_M[2:end,i]),maximum(pltvarsT_Y[2:end,i]),maximum(pltvarsD_RA[2:end,i]))+.1))
# 		# Add all saved Monte Carlo
# 		plt =	plot!(pltvarsT_Y[2:end,i].-1,pltvarsT_Y[2:end,i], seriestype=:scatter, title=(replace(names(pltvarsT_Y)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>"")),
# 				label=false, xaxis = ("Change in each Elasticity of Substitution",font(7)), yaxis=("Diff in Sector activity",font(7)),# ylim=(minimum(Senseoutput[:,i]), maximum(Senseoutput[:,i])),
# 				titlefontsize=10, markersize=3, markerstrokewidth=0, color="lightgreen", alpha=0.3)
# 		plt =	plot!(pltvarsY[2:end,1],pltvarsY[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange", alpha=0.3)
# 		plt =	plot!(pltvarsVA[2:end,1].-1,pltvarsVA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55", alpha=0.3)
# 		plt =	plot!(pltvarsT_M[2:end,1],pltvarsT_M[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="light blue", alpha=0.3)
# 		plt =	plot!(pltvarsM[2:end,1],pltvarsM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue", alpha=0.3)
# 		plt =	plot!(pltvarsT_A[2:end,1].-2,pltvarsT_A[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple", alpha=0.3)
# 		plt =	plot!(pltvarsA[2:end,1],pltvarsA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red", alpha=0.3)
# 		plt =	plot!(pltvarsDM[2:end,1].-2,pltvarsDM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow", alpha=0.3)
# 		plt =	plot!(pltvarsD_RA[2:end,1].-1,pltvarsD_RA[2:end,i], label=false, seriestype=:scatter, markersize=2, markerstrokewidth=0, color="cyan", alpha=0.1) # Any value other than 1 makes all values 1...
# 			#   display(plt)		
# 	# sleep(.5)
# 		push!(groupplots, plt)
# 	end
# 	l = @layout([a{0.1w} [a d e; f g h ; i j k ; l m n; o p q; r s t; u v w]])
# 	plot((groupplots...), layout=l, size=(1000,1200), left_margin=2*Plots.mm, bottom_margin=-1*Plots.mm, top_margin=3*Plots.mm, show=true)
# # display(dspplt)
# 	savefig("Sens_SubSet.png")
# end

# # SensePlotInd(2,13,16,58,73,84,87,129,146,210,204,214,281,275,285,286,287,288,144,145,289)
