# Replication of the WiNDC national MGE model
using MPSGE, JLD2, CSV

using JuMP,PATHSolver


""" 
This function generates a GAMS-like report detailing both the values of the variables,
but also the evaluation of the complementary constraint. Very useful for debugging as the
complementary constraint should always have value 0.
"""

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

	out = "var_name\t value\t\t margin\n"
	for elm in all_variables(m)

		# val = round(value(elm),digits = decimals)
		val = JuMP.is_parameter(elm) ? round(JuMP.parameter_value(elm), digits = decimals) : round(JuMP.value(elm), digits=decimals)

		margin = "."
		try
			margin = round(value(mapping[elm]),digits = mdecimals)
		catch
			margin = "."
		end
		

		out = out*"$elm\t\t $val\t\t $margin\n"
	end

	return(out)
end


# cd(dirname(Base.source_path()))
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./national_ls/DAAData.jld2"))["data"] # load in date from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./national_ls/Indices.jld2"))["data"] # load in date from saved Notebook output Dict, named P
y_ = filter!(x -> x != :oth && x!= :use, S[:i][:]) # These 2 sectors 'use' & 'oth' are in the indices list, but have no data (and therefore cause problems)
a_ = filter!(x -> x != :fbt && x != :mvt && x != :gmt, copy(y_))

# Indexes (set from the data files, via the notebook)
# n = 73   # This is for running with less sectors for quicker troubleshotting etc. Uncomment, set # sectors, and replace 'end' with 'n' in sectorsi  = S[:i][1:end]
n = length(S[:i])
sectorsi  = S[:i][1:n] # "BEA Goods and sectors categories", is "i" in GAMS
sectorsj = copy(sectorsi) # "BEA Goods and sectors categories", is "j" in GAMS, somehow different
xfd = filter!(x -> x != :pce, S[:fd]) # "BEA Final demand categories",
ts = S[:ts] # "BEA Taxes and subsidies categories",
valueadded = filter!(s -> s != :othtax, S[:va]) # "BEA Value added categories excluding othtax", va in GAMS
margin  = S[:m] # "Margins (trade or transport)"; m in GAMS

yr = Symbol(2017)
# PARAMETERS

# Data For a single year, knock out one dimension
y_0 = P[:y_0][yr,:] #	"Gross output",
ys_0 = P[:ys_0][yr,:,:] #	"Sectoral supply",
ty_0 = P[:ty_0][yr,:] #	"Output tax rate"
fs_0 = P[:fs_0][yr,:] #	"Household supply", # All zeros
id_0 = P[:id_0][yr,:,:] #	"Intermediate demand",
fd_0 = P[:fd_0][yr,:,:] #	"Final demand",
va_0 = P[:va_0][yr,:,:] #	"Value added",
# ts_0 = P[:ts_0][yr,:,:] #	"Taxes and subsidies", Not in this model
m_0 = P[:m_0][yr,:] #	"Imports",
x_0 = P[:x_0][yr,:] #	"Exports of goods and services",
# mrg_0 = P[:mrg_0][yr,:] #	"Trade margins", Not in this model
# trn_0 = P[:trn_0][yr,:] #	"Transportation costs",  Not in this model
# duty_0 = P[:ty_0][yr,:] #	"Import duties", Not in this model
# sbd_0 = P[:sbd_0][yr,:] #	"Subsidies on products", Not in this model
# tax_0 = P[:tax_0][yr,:] #	"Taxes on products", Not in this model
ms_0 = P[:ms_0][yr,:,:] #	"Margin supply",
md_0 = P[:md_0][yr,:,:] #	"Margin demand",
s_0 = P[:s_0][yr,:] #	"Aggregate supply",
a_0 = P[:a_0][yr,:][a_]  #	"Armington supply",
bopdef_0 = P[:bopdef_0][yr] #	"Balance of payments deficit",
ta_0 = P[:ta_0][yr,:] #	"Tax net subsidy rate on intermediate demand", Initial, for price
tm_0 = P[:tm_0][yr,:] #	"Import tariff"; Initial, for price 

# ty_0 = add!(WiNnat, Parameter(:ty, indices = (sectorsj,), value=P[:ty_0][year,:].data)) #	"Output tax rate",

# These are filters which are actually set down in lines 269-273 in the gms code  :

# sets	y_(j)	"Sectors with positive production",
# 	a_(i)	"Sectors with absorption",
# 	py_(i)	"Goods with positive supply",
# 	xfd(fd) "Exogenous components of final demand";

# Filters from lines 269-273 in the GAMS version
# 	y_(j) = yes$sum(i,ys0(j,i));
# 	a_(i) = yes$a0(i);
# 	py_(i) = yes$sum(j,ys0(j,i));
# 	xfd(fd) = yes$(not sameas(fd,'pce'));
# *	xfd(fd) = yes$(not sameas('pce', fd));

# function timeWiNnat(n::Int64)
WiNnat = MPSGE.Model()

	# parameters
	ta = add!(WiNnat, MPSGE.Parameter(:ta, indices = (sectorsi,), value=P[:ta_0][yr,sectorsi].data)) #	"Tax net subsidy rate on intermediate demand",
	tm = add!(WiNnat, MPSGE.Parameter(:tm, indices = (sectorsi,), value=P[:tm_0][yr,sectorsi].data)) #	"Import tariff";

	# sectors:
	Y = add!(WiNnat, Sector(:Y, indices=(sectorsj,)))
	A = add!(WiNnat, Sector(:A, indices=(sectorsi,)))

	MS = add!(WiNnat, Sector(:MS, indices=(margin,)))

	# commodities:
	# Should be filtered for sectors in $a0(i)	? Seems to work better to just loop of a_
	PA  = add!(WiNnat, Commodity(:PA, indices=(sectorsi, ))) #	Armington price
	# Should be filtered for sectors in $py_(i)   ? py_ is the same as y_
	PY  = add!(WiNnat, Commodity(:PY, indices=(sectorsi,))) #	Supply
	PVA = add!(WiNnat, Commodity(:PVA, indices=(valueadded,))) #		Value-added
	PM  = add!(WiNnat, Commodity(:PM, indices=(margin,))) #		Margin
	PFX = add!(WiNnat, Commodity(:PFX))	#	Foreign exchnage

	# consumers:
	RA = add!(WiNnat, Consumer(:RA, benchmark = sum(fd_0[:,:pce]) ))

	# production functions
	for j in y_
		@production(WiNnat, Y[j], 0., 0.,
		[	
			Output(PY[i], ys_0[j,i], taxes=[Tax(ty_0[j], RA)]) for i in sectorsi if ys_0[j,i]>0
		], 
		[
			[Input(PA[i], id_0[i,j]) for i in sectorsi if id_0[i,j]>0];
			[Input(Nest(
					Symbol("VA$j"),
					1.,
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
		add!(WiNnat, Production(MS[m], 0., 1., 
			[Output(PM[m], sum(ms_0[:,m]) ) ],
			[Input(PY[i], ms_0[i,m]) for i in sectorsi if ms_0[i,m]>0])) 
	end

	for i in a_  
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
		[Demand(PA[i], fd_0[i,:pce]) for i in sectorsi],
		[
			[Endowment(PY[i], fs_0[i]) for i in sectorsi];
			[Endowment(PA[i], -sum(fd_0[i,x] for x in xfd)) for i in sectorsi];  
			[Endowment(PVA[va], sum(va_0[va,:])) for va in valueadded];
			Endowment(PFX, bopdef_0)
		]
		))

	# MPSGE.build(WiNnat)
	# @time solve!(WiNnat, cumulative_iteration_limit=0)
	# return WiNnat
# end

# WiNnat = timeWiNnat(71)

# Counterfactual solve
# set_value((A[(:gmt)]), 1.0)
# set_value((A[(:mvt)]), 1.0)
# set_value((A[(:fbt)]), 1.0)
# set_fixed!(A[(:gmt)], true)
# set_fixed!(A[(:mvt)], true)
# set_fixed!(A[(:fbt)], true)

set_fixed!(RA, true)
solve!(WiNnat, cumulative_iteration_limit=0);

Report = CSV.File(IOBuffer(generate_report(WiNnat._jump_model, mdecimals = 6)));
CSV.write("FullReportBmrk.csv", Report, missingstring="missing", bom=true)

# Counterfactual
for i in sectorsi
	set_value(ta[i], 0.)
	set_value(tm[i], 0.)
end
set_value(RA,  12453.8963) #So far, this updated default normalization value needs to be set, value from GAMS output. 
solve!(WiNnat, cumulative_iteration_limit=10000);#, convergence_tolerance=1e-0);

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

	Report = CSV.File(IOBuffer(generate_report(WiNnat._jump_model, mdecimals=6)));
	CSV.write("FullReportCounter.csv", (Report), missingstring="missing", bom=true)
	
## For testing with variable numbers of sectors	
	# timeWiNnat()
	# [@elapsed timeWiNnat(t) for t in [2 2 2 8 16]]
	# []@time timeWiNnat(t) for t in [2 2 2 8 16]]
	# [@time timeWiNnat(t) for t in [2 2 8 73]]
	# @profview timeWiNnat(111)
