using MPSGE, JLD2, CSV, XLSX, DataFrames

# Not best practise for data, 
cd(dirname(Base.source_path()))
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/national_ls/DAAData.jld2"))["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/national_ls/Indices.jld2"))["data"] # load in data from saved Notebook output Dict, named S
# Alternate, Julia WiNDC generated data
# P= load(joinpath(@__DIR__,"./data/nationaldata_julia/JDAAData.jld2"))["data"] # load in Julia-generated data from saved Notebook output Dict, named P
# S= load(joinpath(@__DIR__,"./data/nationaldata_julia/JIndices.jld2"))["data"] # load in Julia-generated data from saved Notebook output Dict, named S : indexes are a bit different re :use, so use original or update filters

y_ = filter!(x -> x != :oth && x!= :use, S[:i][:]) # These 2 sectors 'use' & 'oth' are in the indices list, but have no data (and therefore cause problems)
a_ = filter!(x -> x != :fbt && x != :mvt && x != :gmt, copy(y_))

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
vam_0 = deepcopy(va_0) #copy for structure
va_0.data[:,:] = va_0.data.-va_0.data/10
vam_0.data[:,:]= vam_0.data/10 # Set second VA nest to all zeros in the benchmark.

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

		# ty_0 = add!(MultiNat, Parameter(:ty, indices = (sectorsj,), value=P[:ty_0][year,:].data)) #	"Output tax rate", Not in this model

"""
 Option to set model build and solve as function for benchmarking tests
"""
# function timeMultiNat(n::Int64)
MultiNat = MPSGE.Model()

	# parameters
	ta = add!(MultiNat, Parameter(:ta, indices = (a_,), value=P[:ta_0][yr,a_].data)) #	"Tax net subsidy rate on intermediate demand",
	tm = add!(MultiNat, Parameter(:tm, indices = (a_,), value=P[:tm_0][yr,a_].data)) #	"Import tariff";
	tch4 = add!(MultiNat, Parameter(:tch4, indices = (y_,), value=zeros(length(y_))))
	tco2 = add!(MultiNat, Parameter(:tco2, indices = (y_,), value=zeros(length(y_))))

	# Elasticity parameters
	t_elas_y =  add!(MultiNat, Parameter(:t_elas_y,  value=0.))
	elas_y =    add!(MultiNat, Parameter(:elas_y,    value=0.))
	elas_va =   add!(MultiNat, Parameter(:elas_va,   value=1.))
	t_elas_m =  add!(MultiNat, Parameter(:t_elas_m,  value=0.))
	elas_m =    add!(MultiNat, Parameter(:elas_m,    value=0.))
	t_elas_a =  add!(MultiNat, Parameter(:t_elas_a,  value=2.))
	elas_a =    add!(MultiNat, Parameter(:elas_a,    value=0.))
	elas_dm =   add!(MultiNat, Parameter(:elas_dm,   value=2.))
	d_elas_ra = add!(MultiNat, Parameter(:d_elas_ra, value=1.))

	# sectors:
	Y = add!(MultiNat, Sector(:Y, indices=(sectorsj,)))
	A = add!(MultiNat, Sector(:A, indices=(a_,)))

	MS = add!(MultiNat, Sector(:MS, indices=(margin,)))

	# commodities:
	PA  = add!(MultiNat, Commodity(:PA, indices=(a_, ))) #	Armington price
	PY  = add!(MultiNat, Commodity(:PY, indices=(sectorsi,))) #	Supply
	PVA = add!(MultiNat, Commodity(:PVA, indices=(valueadded,))) #		Value-added
	PVAM = add!(MultiNat, Commodity(:PVAM, indices=(valueadded,))) #		Value-added
	PM  = add!(MultiNat, Commodity(:PM, indices=(margin,))) #		Margin
	PFX = add!(MultiNat, Commodity(:PFX))	#	Foreign exchnage

	# consumers:
	RA = add!(MultiNat, Consumer(:RA, benchmark = sum(fd_0[:,:pce]) ))

	# production functions
	for j in y_
	# Options to use Floats, or model parameters for elasticities, with syntax :($(t_elas_y)*1), e.g. to test elasticities in sensitivity tests.	
		@production(MultiNat, Y[j], 0., 0.,
		[	
			Output(PY[i], ys_0[j,i], taxes=[Tax(ty_0[j], RA)]) for i in sectorsi if ys_0[j,i]>0
		], 
		[
			[Input(PA[i], id_0[i,j], taxes=[Tax(:($(tco2[j])*1), RA)]) for i in a_ if id_0[i,j]>0];  # filtered to A
			
            [Input(Nest(
                Symbol("VAtop$j"), #
                2.,
                sum(va_0[:,j])+sum(vam_0[:,j]),
                [       
                Input(Nest( #
                        Symbol("VA$j"),
                        1., # or :($(elas_va)*1),
                        sum(va_0[:,j]),
                            [Input(PVA[va], va_0[va,j], taxes=[Tax(:($(tch4[j])*1), RA)]) for va in valueadded if va_0[va,j]>0.] 
                        ),
                        sum(va_0[:,j] )
                        ) #
                        ,
                    Input(Nest(
                        Symbol("VAm$j"),
                        1., # or :($(elas_va)*1),
                        sum(vam_0[:,j]),
                            [Input(PVAM[va], vam_0[va,j], price=4.) for va in valueadded if va_0[va,j]>0.] # Check only for va_0, to match, bc all vam_0 == 0
                        ),
                        sum(vam_0[:,j] )
                        )
                ]
                        ),
            sum(va_0[:,j])+sum(vam_0[:,j]) #
            )
        ]
    ]
)
end

	for m in margin
		add!(MultiNat, Production(MS[m], 0., 0., 
			[Output(PM[m], sum(ms_0[:,m]) ) ],
			[Input(PY[i], ms_0[i,m]) for i in sectorsi if ms_0[i,m]>0])) 
	end

	for i in a_  
		@production(MultiNat, A[i], 2., 0.,
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
						2., # or :($(elas_dm)*1),
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

	add!(MultiNat, DemandFunction(RA, 1.,
		[Demand(PA[i], fd_0[i,:pce]) for i in a_],
		[
			[Endowment(PY[i], fs_0[i]) for i in a_];
			[Endowment(PA[i], -sum(fd_0[i,xfd])) for i in a_];  
			[Endowment(PVA[va], sum(va_0[va,sectorsi])) for va in valueadded];
            [Endowment(PVAM[va], sum(vam_0[va,sectorsi])) for va in valueadded];
			Endowment(PFX, bopdef_0)
		]
		))

set_value(RA, 13138.7573)
# set_fixed!(RA, true)
#Re-set to benchmark
for i in a_
	set_value(ta[i], P[:ta_0][yr,a_][i])
	set_value(tm[i], P[:tm_0][yr,a_][i])
end
for i in y_
    set_value(tch4[i], 0.)
    set_value(tco2[i], 0.)
end
solve!(MultiNat, cumulative_iteration_limit=0)
println("$datayear ","benchmark")

fullvrbnch = var_report(MultiNat, true)
rename!(fullvrbnch, :value => :bnchmrk, :margin => :bmkmarg)
# set_value(RA, 12453.9) # But why...
# set_fixed!(RA, false)
# solve!(MultiNat, cumulative_iteration_limit=10000);
# println("$datayear ","Zero tax counterfactual")

# Counterfactual Value Added for non-methane mitigation is taxed at mitigation cost
set_value(tch4[:agr], 0.4)
set_value(tch4[:min], 0.2)
set_value(tch4[:oil], 0.4)
set_value(tch4[:uti], 0.1)
set_value(tch4[:wst], 0.15)
set_value(tch4[:pip], 0.3)

solve!(MultiNat, cumulative_iteration_limit=10000)
println("$datayear ","Add Methane Tax counterfactual")
fullvrch4 = var_report(MultiNat, true)
rename!(fullvrch4, :value => :ch4, :margin => :ch4marg)

# Counterfactual Fossil fuel extraction is taxed at (VERY) ~ carbon content of combustion
for i in y_
    set_value(tch4[i], 0.0)
end
# test to double-check that we're back to benchmark
solve!(MultiNat, cumulative_iteration_limit=10000);
set_value(tco2[:oil], 0.2)
set_value(tco2[:min], 0.4)
solve!(MultiNat, cumulative_iteration_limit=10000) #;
fullvrco2 = var_report(MultiNat, true)
rename!(fullvrco2, :value => :co2, :margin => :co2marg)

# Counterfactual Carbon Tax AND Value Added for non-methane mitigation is taxed at methane mitigation cost
set_value(tch4[:agr], 0.4)
set_value(tch4[:min], 0.2)
set_value(tch4[:oil], 0.4)
set_value(tch4[:uti], 0.1)
set_value(tch4[:wst], 0.15)
set_value(tch4[:pip], 0.3)

solve!(MultiNat, cumulative_iteration_limit=10000) #;
fullvrboth = var_report(MultiNat, true)
rename!(fullvrboth, :value => :both, :margin => :bothmarg)
FullResults = innerjoin(fullvrbnch, fullvrch4, fullvrco2, fullvrboth, on = [:var], makeunique=true)
CompareFullResults = FullResults[288:end,[1,2,4,6,8]]

Results = innerjoin(vrbnch, vrch4, vrco2, vrboth, on = [:var], makeunique=true)
CompareResults = Results[288:end,[1,2,4,6,8]]
# vr = var_report(MultiNat, true)
# vr2 = deepcopy(vr)
# vr2.var = string.(vr2.var)
# gams_results = XLSX.readxlsx(joinpath(@__DIR__, "Results\\GAMSResults-2024-1-12.xlsx"))
#     a_table = gams_results["Sheet1"][3:490,1:3]  # Generated from JPMGE_MPSGE
#     # WNDCnat = DenseAxisArray(a_table[2:end,2:end],string.(a_table[2:end,1],".",a_table[2:end,2]),a_table[1,2:end])
# Compvr = DataFrame(innerjoin(DataFrame(a_table,[:var, :bnch, :counter]),vr2, on = [:var], makeunique=true))
# Compvr.Counterdiff = Compvr.counter.-Compvr.value  # Comparing set counterfactual only
# print(Compvr)
# print(Compvr[coalesce.(abs.(Compvr.Counterdiff).>0.001, true),:]) # or abs.(Compvr.margin).>
# print(filter(row->any(occursin.("ρ",string(row.var))),CompareFullResults))