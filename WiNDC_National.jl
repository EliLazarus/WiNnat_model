# Replication of the WiNDC national MGE model
using MPSGE, JLD2
# cd(dirname(Base.source_path()))
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
P= load("./nationaldata_ls/DAAData.jld2")["data"] # load in date from saved Notebook output Dict, named P
S= load("./nationaldata_ls/Indices.jld2")["data"] # load in date from saved Notebook output Dict, named P
n=8
# function timeWiNnat(n::Int64)
	# Indexes (set from the data files, via the notebook)
	yr = S[:yr] # "Years in WiNDC Database",
	sectorsi  = S[:i][1:n] # "BEA Goods and sectors categories", is "i" in GAMS
	sectorsj = copy(sectorsi) # "BEA Goods and sectors categories", is "j" in GAMS, somehow different
	xfd = S[:fd] # "BEA Final demand categories",
	ts = S[:ts] # "BEA Taxes and subsidies categories",
	valueadded = S[:va] # "BEA Value added categories excluding othtax", va in GAMS
	margin  = S[:m] # "Margins (trade or transport)"; m in GAMS

	WiNnat = Model()

	year = Symbol(2017)
	# PARAMETERS
	# ty = add!(WiNnat, Parameter(:ty, indices = (sectorsj,), value=P[:ty_0][year,:].data)) #	"Output tax rate",
	ta = add!(WiNnat, Parameter(:ta, indices = (sectorsi,), value=P[:ta_0][year,sectorsi].data)) #	"Tax net subsidy rate on intermediate demand",
	tm = add!(WiNnat, Parameter(:tm, indices = (sectorsi,), value=P[:tm_0][year,sectorsi].data)) #	"Import tariff";

	yr = Symbol(2017)

	# Data For a single year, knock out one dimension
		y_0 = P[:y_0][yr,:] #	"Gross output",
		ys_0 = P[:ys_0][yr,:,:] #	"Sectoral supply",
		ty_0 = P[:ty_0][yr,:] #	"Output tax rate"
		fs_0 = P[:fs_0][yr,:] #	"Household supply",
		id_0 = P[:id_0][yr,:,:] #	"Intermediate demand",
		fd_0 = P[:fd_0][yr,:,:] #	"Final demand",
		va_0 = P[:va_0][yr,:,:] #	"Value added",
		ts_0 = P[:ts_0][yr,:,:] #	"Taxes and subsidies",
		m_0 = P[:m_0][yr,:] #	"Imports",
		x_0 = P[:x_0][yr,:] #	"Exports of goods and services",
		mrg_0 = P[:mrg_0][yr,:] #	"Trade margins",
		trn_0 = P[:trn_0][yr,:] #	"Transportation costs",
		duty_0 = P[:ty_0][yr,:] #	"Import duties",
		sbd_0 = P[:sbd_0][yr,:] #	"Subsidies on products",
		tax_0 = P[:tax_0][yr,:] #	"Taxes on products",
		ms_0 = P[:ms_0][yr,:,:] #	"Margin supply",
		md_0 = P[:md_0][yr,:,:] #	"Margin demand",
		s_0 = P[:s_0][yr,:] #	"Aggregate supply",
		#Data Missing
		# d_0 = P[:d_0][yr,:] #	"Sales in the domestic market",
		a_0 = P[:a_0][yr,:] #	"Armington supply",
		bopdef_0 = P[:bopdef_0][yr] #	"Balance of payments deficit",
		ta_0 = P[:ta_0][yr,:] #	"Tax net subsidy rate on intermediate demand",
		tm_0 = P[:tm_0][yr,:] #	"Import tariff";

	# y_0 = add!(WiNnat, Parameter(:y_0, indices = (yr,i))) #	"Gross output",
		# ys_0 = add!(WiNnat, Parameter(:ys_0, indices = (yr,j,i))) #	"Sectoral supply",
		# ty_0 = add!(WiNnat, Parameter(:ty_0, indices = (yr,j))) #	"Output tax rate"
		# fs_0 = add!(WiNnat, Parameter(:fs_0, indices = (yr,i))) #	"Household supply",
		# id_0 = add!(WiNnat, Parameter(:id_0, indices = (yr,i,j))) #	"Intermediate demand",
		# fd_0 = add!(WiNnat, Parameter(:fd_0, indices = (yr,i,fd))) #	"Final demand",
		# va_0 = add!(WiNnat, Parameter(:va_0, indices = (yr,va,j))) #	"Value added",
		# ts_0 = add!(WiNnat, Parameter(:ts_0, indices = (yr,ts,i))) #	"Taxes and subsidies",
		# m_0 = add!(WiNnat, Parameter(:m_0, indices = (yr,i))) #	"Imports",
		# x_0 = add!(WiNnat, Parameter(:x_0, indices = (yr,i))) #	"Exports of goods and services",
		# mrg_0 = add!(WiNnat, Parameter(:mrg_0, indices = (yr,i))) #	"Trade margins",
		# trn_0 = add!(WiNnat, Parameter(:trn_0, indices = (yr,i))) #	"Transportation costs",
		# duty_0 = add!(WiNnat, Parameter(:ty_0, indices = (yr,i))) #	"Import duties",
		# sbd_0 = add!(WiNnat, Parameter(:bd_0, indices = (yr,i))) #	"Subsidies on products",
		# tax_0 = add!(WiNnat, Parameter(:ax_0, indices = (yr,i))) #	"Taxes on products",
		# ms_0 = add!(WiNnat, Parameter(:ms_0, indices = (yr,i,m))) #	"Margin supply",
		# md_0 = add!(WiNnat, Parameter(:md_0, indices = (yr,m,i))) #	"Margin demand",
		# s_0 = add!(WiNnat, Parameter(:s_0, indices = (yr,i))) #	"Aggregate supply",
		# d_0 = add!(WiNnat, Parameter(:d_0, indices = (yr,i))) #	"Sales in the domestic market",
		# a_0 = add!(WiNnat, Parameter(:a_0, indices = (yr,i))) #	"Armington supply",
		# bopdef_0 = add!(WiNnat, Parameter(:ef_0, indices = (yr,))) #	"Balance of payments deficit",
		# ta_0 = add!(WiNnat, Parameter(:ta_0, indices = (yr,i))) #	"Tax net subsidy rate on intermediate demand",
		# tm_0 = add!(WiNnat, Parameter(:tm_0, indices = (yr,i))) #	"Import tariff";

	# PARAMETERS for single year version
	# y0 = add!(WiNnat, Parameter(:y0, indices = (sectorsi,), value=P[:y_0][year,:].data)) #	"Gross output",
	# ys0 = add!(WiNnat, Parameter(:ys0, indices = (sectorsj,sectorsi), value=P[:ys_0][year,:,:].data)) #	"Sectoral supply",
	# ty0 = add!(WiNnat, Parameter(:ty0, indices = (sectorsj,), value=P[:ty_0][year,:].data)) #	"Output tax rate"
	# fs0 = add!(WiNnat, Parameter(:fs0, indices = (sectorsi,), value=P[:fs_0][year,:].data)) #	"Household supply",
	# id0 = add!(WiNnat, Parameter(:id0, indices = (sectorsi,sectorsj), value=P[:id_0][year,:,:].data)) #	"Intermediate demand",
	# fd0 = add!(WiNnat, Parameter(:fd0, indices = (sectorsi,fd), value=P[:fd_0][year,:,:].data)) #	"Final demand",
	# va0 = add!(WiNnat, Parameter(:va0, indices = (valueadded,sectorsj), value=P[:va_0][year,:,:].data)) #	"Value added",
	# ts0 = add!(WiNnat, Parameter(:ts0, indices = (ts,sectorsi), value=P[:ts_0][year,:,:].data)) #	"Taxes and subsidies",
	# m0 = add!(WiNnat, Parameter(:m0, indices = (sectorsi,), value=P[:m_0][year,:].data)) #	"Imports",
	# x0 = add!(WiNnat, Parameter(:x0, indices = (sectorsi,), value=P[:x_0][year,:].data)) #	"Exports of goods and services",
	# mrg0 = add!(WiNnat, Parameter(:rg0, indices = (sectorsi,), value=P[:mrg_0][year,:].data)) #	"Trade margins",
	# trn0 = add!(WiNnat, Parameter(:rn0, indices = (sectorsi,), value=P[:trn_0][year,:].data)) #	"Transportation costs",
	# duty0 = add!(WiNnat, Parameter(:ty0, indices = (sectorsi,), value=P[:duty_0][year,:].data)) #	"Import duties",
	# sbd0 = add!(WiNnat, Parameter(:bd0, indices = (sectorsi,), value=P[:sbd_0][year,:].data)) #	"Subsidies on products",
	# tax0 = add!(WiNnat, Parameter(:ax0, indices = (sectorsi,), value=P[:tax_0][year,:].data)) #	"Taxes on products",
	# ms0 = add!(WiNnat, Parameter(:ms0, indices = (sectorsi,margin), value=P[:ms_0][year,:,:].data)) #	"Margin supply",
	# md0 = add!(WiNnat, Parameter(:md0, indices = (margin,sectorsi), value=P[:md_0][year,:,:].data)) #	"Margin demand",
	# s0 = add!(WiNnat, Parameter(:s0, indices = (sectorsi,), value=P[:s_0][year,:].data)) #	"Aggregate supply",

	# d0 = add!(WiNnat, Parameter(:d0, indices = (sectorsi,), value=P[:d_0][year,:].data)) #	"Sales in the domestic market",

	# a0 = add!(WiNnat, Parameter(:a0, indices = (sectorsi,), value=P[:a_0][year,:].data)) #	"Armington supply",
	# bopdef0 = add!(WiNnat, Parameter(:bopdef0, value=(P[:bopdef_0][year]))) #	"Balance of payments deficit",
	# ta0 = add!(WiNnat, Parameter(:ta0, indices = (sectorsi,), value=P[:ta_0][year,:].data)) #	"Tax net subsidy rate on intermediate demand",
	# tm0 = add!(WiNnat, Parameter(:tm0, indices = (sectorsi,), value=P[:tm_0][year,:].data)) #	"Import tariff";

	# TODO Not sure what this is?
	# Looks like a filters, but I don't see where they're set.

	# sets	y_(j)	"Sectors with positive production",
	# 	a_(i)	"Sectors with absorption",
	# 	py_(i)	"Goods with positive supply",
	# 	xfd(fd) "Exogenous components of final demand";

	# sectors:
	Y = add!(WiNnat, Sector(:Y, indices=(sectorsj,)))
	A = add!(WiNnat, Sector(:A, indices=(sectorsi,)))

	MS = add!(WiNnat, Sector(:MS, indices=(margin,)))

	# commodities:
	# Should be filtered for sectors in $a0(i)	?!
	PA  = add!(WiNnat, Commodity(:PA, indices=(sectorsi, ))) #	Armington price
	# Should be filtered for sectors in $py_(i)   ?!
	PY  = add!(WiNnat, Commodity(:PY, indices=(sectorsi,))) #	!	Supply
	PVA = add!(WiNnat, Commodity(:PVA, indices=(valueadded,))) #		!	Value-added
	PM  = add!(WiNnat, Commodity(:PM, indices=(margin,))) #		!	Margin
	PFX = add!(WiNnat, Commodity(:PFX))	#	!	Foreign exchnage

	# consumers:
	RA = add!(WiNnat, Consumer(:RA, benchmark = sum(fd_0) ))

	# production functions
	for j in sectorsj
		@production(WiNnat, Y[j], 0., 0.,
		[Output(PY[i], ys_0[j,i], taxes=[Tax(ty_0[j], RA)]) for i in sectorsi], 
		[
		 [Input(PA[i], id_0[i,j]) for i in sectorsi];
     	 [Input(Nest(:VA, 1., sum(va_0[va,j] for va in valueadded),
		  [Input(PVA[va], va_0[va,j]) for va in valueadded]),sum(va_0[va,j] for va in valueadded))]
		]
			)
	end

	for m in margin
		add!(WiNnat, Production(MS[m], 0., 1., 
		[Output(PM[m], sum(ms_0[i,m] for i in sectorsi) ) for m in margin],
		[Input(PY[i], ms_0[i,m]) for i in sectorsi]))
	end

	for i in sectorsi 
		# add!(WiNnat, Production(A[i], 2., 0.,
		@production(WiNnat, A[i], 2., 0.,
		[[Output(PA[i], a_0[i], taxes=[Tax(:($(ta[i])*1), RA)], price=(1-ta_0[i]))];
		#  for i in sectorsi]; # Question re ta and ta0
		Output(PFX, x_0[i])],

		#For testing without nesting
		[
			# [Input(PY[i], y_0[i]) for i in sectorsi];
		#  [Input(PFX, m_0[i], taxes=[Tax(:($(tm[i])*1), RA)], price=(1+tm_0[i]))];

	    [Input(Nest(:dm, 2., sum(y_0[i]+m_0[i] for i in sectorsi),
		    [Input(PY[i], y_0[i]),
			 Input(PFX, m_0[i], taxes=[Tax(:($(tm[i])*1), RA)], price=:(1+$(tm[i])*1))
			 ]), sum(y_0[i]+m_0[i] for i in sectorsi)
			 )];
		[Input(PM[m], md_0[m,i]) for m in margin]
		]
		)

	end

	add!(WiNnat, DemandFunction(RA, 1.,
		[Demand(PA[i], fd_0[i,:pce]) for i in sectorsi],
		[[Endowment(PY[i], fs_0[i]) for i in sectorsi];
		 [Endowment(PA[i], -sum(fd_0[i,x] for x in xfd)) for i in sectorsi];
		 [Endowment(PVA[va], sum(va_0[va,:])) for va in valueadded];
		 Endowment(PFX, bopdef_0)
		]))

	# solve!(WiNnat, cumulative_iteration_limit=0)

	solve!(WiNnat)
	# solve!(WiNnat)

# end

# timeWiNnat(2)
# for t in [2 2 4 8]
	# [@elapsed timeWiNnat(t) for t in [2 2 2 8 16]]
	# [@time timeWiNnat(t) for t in [2 2 2 8 16]]
	# [@time timeWiNnat(t) for t in [2 2 8]]

# end
# @profview solve!(WiNnat)
# @profview timeWiNnat(8)