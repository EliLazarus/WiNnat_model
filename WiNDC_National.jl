# Replication of the WiNDC national MGE model
using MPSGE, JLD2, CSV

using JuMP,PATHSolver
# For sensitivity tests, timing, and plotting
using DataFrames, Plots, Tables, Dates, Distributions

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

"""
 Option to set model build and solve as function for benchmarking tests
"""
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

	# return WiNnat
# end
# @time solve!(WiNnat, cumulative_iteration_limit=0)
# @time MPSGE.build(WiNnat)
		
"""
 End of optional function version of model used for benchmarking
"""

set_value(RA, 13138.7573)
set_fixed!(RA, true)

for i in a_
	set_value(ta[i], P[:ta_0][yr,a_][i])
	set_value(tm[i], P[:tm_0][yr,a_][i])
end
solve!(WiNnat, cumulative_iteration_limit=0)#, convergence_tolerance=1e-9)
println("$year ","bnchmk")

Report = CSV.File(IOBuffer(generate_report(WiNnat._jump_model, mdecimals = 6)))
CSV.write("./Results/FullReportCounter.csv", Report, missingstring="missing", bom=true)

# Counterfactual
for i in a_
	set_value(ta[i], 0.)
	set_value(tm[i], 0.)
end

set_fixed!(RA, false)
solve!(WiNnat, cumulative_iteration_limit=10000)
println("$year ","counter")

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
