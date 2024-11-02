using MPSGE
import JuMP
# For sensitivity tests, timing, and plotting
using DataFrames, Plots, Tables, Dates, Distributions, CSV
set_value!(CH4_tax, CH4_taxrate)
set_value!(CO2_tax, CO2_taxrate)
# set_value!(CH4_tax, 0)
solve!(MultiNat)
###  SENSITIVITY TESTS
## First, save counterfactual baseline. Run this once, save and then comment out so it doesn't overwrite with values from random sampling
	# rp = generate_report(MultiNat)
	# basecountervalues = [:t_elas_y           :elas_y        :t_elas_a      :elas_a      :elas_dm       :d_elas_ra              transpose(rp[:,:var])   
	# # :elas_va           :t_elas_m           :elas_m           
	# 					value(t_elas_y) value(elas_y) value(t_elas_a) value(elas_a) value(elas_dm) value(d_elas_ra) transpose(rp[:,:value])]
	# # 					# value(elas_va) value(t_elas_m) value(elas_m) 
	# CSV.write("./SensitivityTest_Outputs/MultiNat_Both_CounterVals.csv", Tables.table(basecountervalues), missingstring="missing", bom=true)

	# ## Before 1st run, set and write headers and baseline values into CSV for sensitivity test data output
	# basecountervalues = DataFrame(CSV.File("./SensitivityTest_Outputs/MultiNat_Both_CounterVals.csv", header = 2, skipto=3))
	# CSV.write("./SensitivityTest_Outputs/MultiNat_Both_ExpectedVal.csv", basecountervalues, missingstring="missing", bom=true)

## Visualise Distributions by taking a random sample and plotting
histogram(rand(truncated(Gamma(1,2),upper=2),10^6), label = "Gamma(1,2)")
histogram(rand(truncated(Normal(2,2),lower=0, upper=4),10^6), label = "truncated Normal(2,2)")
histogram(rand(truncated(Normal(1,2),lower=0, upper=2),10^6), label = "truncated Normal(2,2)")

## For this sensetivity test, we alter all elasticities simultaneously for each Monte Carlo run
function ExpectedVal(n, model)
	starttime = now()
	out = Array{Float64}(undef, n, length(rp[:,1])+6) # 9 elasticity values plus the 288 variables I'll report
	for i in 1:n
		t_elas_y_val = rand(truncated(Gamma(1,2),upper=2)) # Yt_elas, 0.
		elas_y_val = rand(truncated(Gamma(1,2),upper=2))  # Ys_elas, 0.
		# elas_va_val = rand(truncated(Normal(1,2),lower=0)) #1  Not updating for MultiNat
		# t_elas_m_val = rand(truncated(Gamma(1,2),upper=2)) # Not updating for MultiNat
		# elas_m_val = rand(truncated(Gamma(1,2),upper=2))  # Not updating for MultiNat
		t_elas_a_val = rand(truncated(Normal(2,2),lower=0, upper=4)) #2 		At_elas, 2.
		elas_a_val = rand(truncated(Gamma(1,2),upper=2)) # 		As_elas, 0.
		elas_dm_val = rand(truncated(Normal(2,2),lower=0, upper=4)) #2 		Adm_elas, 2.
		d_elas_ra_val = rand(truncated(Normal(1,2),lower=0, upper=2)) #1 		dm_elas, 1.

		set_value!(t_elas_y, t_elas_y_val)
		set_value!(elas_y, elas_y_val)
		# set_value!(elas_va, elas_va_val)
		# set_value!(t_elas_m, t_elas_m_val)
		# set_value!(elas_m, elas_m_val)
		set_value!(t_elas_a, t_elas_a_val)
		set_value!(elas_a, elas_a_val)
		set_value!(elas_dm, elas_dm_val)
		set_value!(d_elas_ra, d_elas_ra_val)

		solve!(model, cumulative_iteration_limit=10000);#, convergence_tolerance=1e-0);
		# save sampled elasticity values, and indices for Y through to MS, and PA through PFX in results, variables we're most interested in
		out[i,:] = [value(t_elas_y) value(elas_y) value(t_elas_a) value(elas_a) value(elas_dm) value(d_elas_ra) transpose(generate_report(model)[!,:value])]
		# value(elas_va) value(t_elas_m) value(elas_m) 
	println(i,": So far ",now()-starttime)
	end
	println("It took",now()-starttime)
	return out
end

## Run Monte Carlo n times, and append to csv (can run many small batches) & write the results from the sensitivity test to csv (command on 1 line so alway write after running)
outall = ExpectedVal(5, MultiNat); CSV.write("./SensitivityTest_Outputs/MultiNat_Both_ExpectedVal.csv", Tables.table(outall), missingstring="missing", bom=true, append=true)

## get all results back for plotting
ExVal = DataFrame(CSV.File("./SensitivityTest_Outputs/MultiNat_Both_ExpectedVal.csv", header = 1, skipto=2))

## df data frame; fv first variable index from SensElas; lv last variable index, for loop, 21 variables at a time for plotting in 3x7
function plotgroups(df, fv, lv)
	pltvarsAll = df[:,fv:lv]
	saveplots =[]
	# Loop over 21 variables for plotting
	for i in 1:length(pltvarsAll[1,:])
	# 	# Start with original value from the model, and set the axis ranges from the values (including the originals)
		plt =	scatter([pltvarsAll[1,i]],[pltvarsAll[1,i]], size=(600,1400),color = "green",
		 label = false, markersize = 9,
		ylimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(maximum(maximum(eachcol(pltvarsAll)))/8)),
		xlimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(minimum(minimum(eachcol(pltvarsAll)))/8)), 
		 title=replace(names(pltvarsAll)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>""))
		 plt =	plot!(pltvarsAll[2:end,i],pltvarsAll[2:end,i], seriestype=:scatter,
				label=false, xaxis = ("",font(7)), yaxis=("",font(7)),
				titlefontsize=10, markersize=3, markerstrokewidth=0, color="red", alpha=0.4)
		# display(plt) # If you want to have a quick look at them as they're generated
		# sleep(.5)
		push!(saveplots, plt)
	end
	pltAll = plot((saveplots...), layout=(7,3), size=(900,1000),  left_margin=1*Plots.mm,right_margin=1*Plots.mm,
	 bottom_margin=-1*Plots.mm, top_margin=1*Plots.mm, plot_title="Counterfactual Base Value,\nand the impact on each variable from varying all elasticities",
	  plot_titlefontsize=9, plot_titlefonthalign=:left, fontfamily="Palatino Roman")
	  display(pltAll)
	savefig(replace("EvVal"*names(ExVal)[fv]*"-"*names(ExVal)[lv]*".png",":" => ""))
end

### Run plots in groups
# plotgroups(ExVal, 10,30) # fv = 10; lv=30; # Y
# plotgroups(ExVal, 31,51) # Y
# plotgroups(ExVal, 52,72) # Y
# plotgroups(ExVal, 73,80) # Y
# plotgroups(ExVal, 81,101) # A
# plotgroups(ExVal, 102,122) # A
# plotgroups(ExVal, 123,143) # A
# plotgroups(ExVal, 144,151) # A
# plotgroups(ExVal, 152,153) # MS
# plotgroups(ExVal, 154,174) # PA
# plotgroups(ExVal, 175,195) # PA
# plotgroups(ExVal, 196,216) # PA
# plotgroups(ExVal, 217,221) # PA
# plotgroups(ExVal, 222,242) # PY
# plotgroups(ExVal, 243,263) # PY
# plotgroups(ExVal, 264,284) # PY
# plotgroups(ExVal, 285,292) # PY
# plotgroups(ExVal, 293,297) # PVA -> PFX


## Selection of variables to show as a sample in one collection
# plotgroups(ExVal,10,21,24,66,81,92,95,137,154,218,212,222,289,283,293,294,295,296,152,153,297)
## I didn't figure out to make it work using the above function, so just set it separately
# pltvarsAll = ExVal[:,[10,21,24,66,81,92,95,137,154,218,212,222,289,283,293,294,295,296,152,153,297]]
pltvarsAll_base = ExVal[:,[440+6,  436 + 6, 439 + 6]]
# TotEm is 440+6, CO2TotEm is 436 + 6, CH4TotEm is 434 + 6
pltvarsAll = pltvarsAll_base[:,:].*10^3 # Multiply by 1,000 to put in MMtCO2eq

saveplots =[]
	for i in 1:length(pltvarsAll[1,:])
		# plt =	scatter([pltvarsAll[1,i]],[pltvarsAll[1,i]], size=(600,1400),color = "green",
		plt =	scatter([190],[pltvarsAll[1,i]], size=(600,1400),color = "green",
		 label = false, markersize = 9,
		ylimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(maximum(maximum(eachcol(pltvarsAll)))/8)),
		# xlimits=(minimum(minimum(eachcol(pltvarsAll)))-(minimum(minimum(eachcol(pltvarsAll)))/8),maximum(maximum(eachcol(pltvarsAll)))+(minimum(minimum(eachcol(pltvarsAll)))/8)), 
		xlimits=(180,200), 
		title=replace(names(pltvarsAll)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>""))
		plt =	plot!([fill(190,length(pltvarsAll[:,i]))],pltvarsAll[2:end,i], seriestype=:scatter,
		# plt =	plot!(pltvarsAll[2:end,i],pltvarsAll[2:end,i], seriestype=:scatter,
				label=false, xaxis = ("",font(7)), yaxis=("",font(7)),
				titlefontsize=10, markersize=3, markerstrokewidth=0, color="red", alpha=0.4)
		# display(plt)
				# sleep(.5)
		push!(saveplots, plt)
	end
	pltAll = plot((saveplots...), layout=(7,3), size=(900,1000),  left_margin=1*Plots.mm,right_margin=1*Plots.mm,
	 bottom_margin=-1*Plots.mm, top_margin=1*Plots.mm, plot_title="Counterfactual Base Value,\nand the impact on each variable from varying all elasticities",
	  plot_titlefontsize=9, plot_titlefonthalign=:left, fontfamily="Palatino Roman")
	  display(pltAll)
	savefig("ExValSupSet.png")

## SENSETIVITY TEST FOR INDIVIDUAL ELASTICITY PARAMETERS
## First, function to set header and baseline counterfactual values rows.
function SensitivityTestHead(param, m::MPSGEModel)
	out = Array{Any}(undef, 2, 489)#length(rp[:,1])+1) 
	set_value!(t_elas_y, 0.)
	set_value!(elas_y, 0.05)
	# set_value!(elas_va, 1.)
	# set_value!(t_elas_m, 0.)
	# set_value!(elas_m, 0.)
	set_value!(t_elas_a, 2.)
	set_value!(elas_a, 0.)
	set_value!(elas_dm, 2.)
	set_value!(d_elas_ra, 1.)
	solve!(m, cumulative_iteration_limit=10000);
	# paramname = param.model._parameters[param.index].name
	rp = generate_report(m)
	out[1,:] = [param transpose(rp[!,:var])]
	out[2,:] = [value(param) transpose(rp[!,:value])]
	#  = [paramname        transpose(all_variables(WiNnat._jump_model)[[152:295;435:578],:])]
	#  = [get_value(param)  					 transpose(JuMP.value.(all_variables(WiNnat._jump_model)[[152:295;435:578],:]))]
	CSV.write("./SensitivityTest_Outputs/MultiBoth_Sens_"*String(param.name)*".csv", Tables.table(out), missingstring="missing", bom=true)
end

## Here, we set a function to select a single elasticity and vary that according to the distribution argument, n times in each Monte Carlo run, and append values to csv
function SensitivityTest(param, n, dist, m)
	out = Array{Float64}(undef, n, 489)#length(rp[:,1])+1)
	set_value!(t_elas_y, 0.)
	set_value!(elas_y, 0.05)
	# set_value!(elas_va, 1.)
	# set_value!(t_elas_m, 0.)
	# set_value!(elas_m, 0.)
	set_value!(t_elas_a, 2.)
	set_value!(elas_a, 0.)
	set_value!(elas_dm, 2.)
	set_value!(d_elas_ra, 1.)
	for i in 1:n
		pdr = rand(dist)
		set_value!(param, pdr)
		println(i,"/",n,". ",param.name)
		solve!(m, cumulative_iteration_limit=10000);#, convergence_tolerance=1e-0);
	##indexes for Y through to MS, and PA through PFX in results, variables we're most interested in
		out[i,:] = [value(param) ; generate_report(m)[!,:value]]
	end
	CSV.write("./SensitivityTest_Outputs/MultiBoth_Sens_"*String(param.name)*".csv", Tables.table(out), missingstring="missing", bom=true, append=true)
	# return out
end

### Run heading function once, and then can run ST in batches
# Example: for the first elasticity parameter
SensitivityTestHead(t_elas_y, MultiNat)
SensitivityTestHead(elas_y, MultiNat)
# SensitivityTestHead(elas_va, MultiNat)
# SensitivityTestHead(t_elas_m, MultiNat)
# SensitivityTestHead(elas_m, MultiNat)
SensitivityTestHead(t_elas_a, MultiNat)
SensitivityTestHead(elas_a, MultiNat)
SensitivityTestHead(elas_dm, MultiNat)
SensitivityTestHead(d_elas_ra, MultiNat)

histogram(rand(truncated(Gamma(1,2),upper=2),10^6), label = "Gamma(1,2)")
histogram(rand(truncated(Normal(2,2),lower=0, upper=4),10^6), label = "truncated Normal(2,2)")
histogram(rand(truncated(Normal(1,2),lower=0, upper=2),10^6), label = "truncated Normal(2,2)")

SensitivityTest(t_elas_y, 100, truncated(Gamma(1,2),upper=2), MultiNat) # 0 ->2
SensitivityTest(elas_y, 100, truncated(Gamma(1,2),upper=2), MultiNat) # 0 ->2 
# SensitivityTest(elas_va, 100, truncated(Normal(1,2),lower=0,), MultiNat) #1
# SensitivityTest(t_elas_m, 100, truncated(Gamma(1,2),upper=2), MultiNat) 
# SensitivityTest(elas_m, 100, truncated(Gamma(1,2),upper=2), MultiNat) 
SensitivityTest(t_elas_a, 100, truncated(Normal(2,2),lower=0, upper=4), MultiNat) #2 (0->4)
SensitivityTest(elas_a, 100, truncated(Gamma(1,2),upper=2), MultiNat) # 0 ->2 
SensitivityTest(elas_dm, 100, truncated(Normal(2,2),lower=0, upper=4), MultiNat) #2 (0->4)
SensitivityTest(d_elas_ra, 100, truncated(Normal(1,2),lower=0, upper=2), MultiNat) #1 (0->2)
### Get All Sensitivity results for plotting together
SenseoutputT_Y = DataFrame(CSV.File("./SensitivityTest_Outputs/MultiBoth_Sens_t_elas_y.csv", header = 2, skipto=3))
SenseoutputY = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_elas_y.csv"), header = 2, skipto=3))
# SenseoutputVA = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/SensitivityTest_Outputs/MultiBoth_Sens_elas_va.csv"), header = 2, skipto=3))
# SenseoutputT_M = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_t_elas_m.csv"), header = 2, skipto=3))
# SenseoutputM = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_elas_m.csv"), header = 2, skipto=3))
SenseoutputT_A = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_t_elas_a.csv"), header = 2, skipto=3))
SenseoutputA = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_elas_a.csv"), header = 2, skipto=3))
SenseoutputDM = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_elas_dm.csv"), header = 2, skipto=3))
SenseoutputD_RA = DataFrame(CSV.File(joinpath(@__DIR__,"./SensitivityTest_Outputs/MultiBoth_Sens_d_elas_ra.csv"), header = 2, skipto=3))
# # #fv first variable index from SensElas; lv last variable index, 21 variables, for plotting loop and 3x7 plot collection
function SensePlot(fv, lv)
	pltvarsT_Y = [SenseoutputT_Y[:,:t_elas_y] SenseoutputT_Y[:,fv:lv]]
	pltvarsY = [SenseoutputY[:,:elas_y] SenseoutputY[:,fv:lv]]
	pltvarsVA = [SenseoutputVA[:,:elas_va] SenseoutputVA[:,fv:lv]] #1
	pltvarsT_M = [SenseoutputT_M[:,:t_elas_m] SenseoutputT_M[:,fv:lv]]
	pltvarsM = [SenseoutputM[:,:elas_m] SenseoutputM[:,fv:lv]] 
	pltvarsT_A = [SenseoutputT_A[:,:t_elas_a] SenseoutputT_A[:,fv:lv]] #2
	pltvarsA = [SenseoutputA[:,:elas_a] SenseoutputA[:,fv:lv]] 
	pltvarsDM = [SenseoutputDM[:,:elas_dm] SenseoutputDM[:,fv:lv]] #2
	pltvarsD_RA = [SenseoutputD_RA[:,:d_elas_ra] SenseoutputD_RA[:,fv:lv]] #1
	groupplots =[]
	plt = scatter(1,markersize=12, markerstrokewidth=2, size=(300,400), xaxis=false, yaxis=false,
	label="   Baseline\n   Counterfactual\n   Value", color="green", legendfontsize=9,legend_title="Elasticity\n  Parameters",
	legend_title_font_size=11, fg_legend=:transparent,labelspacing=0.8, framestyle=:none,
	fontfamily="Computer Modern")
	plot!(1,1, label="   t_elas_y", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightgreen")
	plot!(1,1, label="   elas_y", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange")
	plot!(1,1, label="   elas_va", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55")
	plot!(1,1, label="   t_elas_m", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightblue")
	plot!(1,1, label="   elas_m", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue")
	plot!(1,1, label="   t_elas_a", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple")
	plot!(1,1, label="   elas_a", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red" )
	plot!(1,1, label="   elas_dm", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow")
	plot!(1,1, label="   d_elas_ra", seriestype=:scatter, markersize=3, markerstrokewidth=0, color="cyan") # Any value other than 1 makes all values 1...
	push!(groupplots, plt)
	for i in 2:length(pltvarsA[1,:])
	# 	# Start with original value from the model
		plt =	scatter([pltvarsM[1,1]],[pltvarsM[1,i]], color = "darkgreen", label = false, markersize = 9, 
		ylimits=(min(pltvarsA[1,i],minimum(pltvarsA[2:end,i]),minimum(pltvarsM[2:end,i]),minimum(pltvarsY[2:end,i]),minimum(pltvarsDM[2:end,i]),minimum(pltvarsVA[2:end,i]),minimum(pltvarsT_A[2:end,i]),minimum(pltvarsT_M[2:end,i]),minimum(pltvarsT_Y[2:end,i]),minimum(pltvarsD_RA[2:end,i]))-.1,
		max(pltvarsA[1,i],maximum(pltvarsA[2:end,i]),maximum(pltvarsM[2:end,i]),maximum(pltvarsY[2:end,i]),maximum(pltvarsDM[2:end,i]),maximum(pltvarsVA[2:end,i]),maximum(pltvarsT_A[2:end,i]),maximum(pltvarsT_M[2:end,i]),maximum(pltvarsT_Y[2:end,i]),maximum(pltvarsD_RA[2:end,i]))+.1))
		plt =	plot!(pltvarsT_Y[2:end,i].-1,pltvarsT_Y[2:end,i], seriestype=:scatter, title=(replace(names(pltvarsT_Y)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>"")),
				label=false, xaxis = ("Change in each Elasticity of Substitution",font(7)), yaxis=("Diff in Sector activity",font(7)),# ylim=(minimum(Senseoutput[:,i]), maximum(Senseoutput[:,i])),
				titlefontsize=10, markersize=3, markerstrokewidth=0, color="lightgreen", alpha=0.3)
		plt =	plot!(pltvarsY[2:end,1],pltvarsY[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange", alpha=0.3)
		plt =	plot!(pltvarsVA[2:end,1].-1,pltvarsVA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55", alpha=0.3)
		plt =	plot!(pltvarsT_M[2:end,1],pltvarsT_M[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="light blue", alpha=0.3)
		plt =	plot!(pltvarsM[2:end,1],pltvarsM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue", alpha=0.3)
		plt =	plot!(pltvarsT_A[2:end,1].-2,pltvarsT_A[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple", alpha=0.3)
		plt =	plot!(pltvarsA[2:end,1],pltvarsA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red", alpha=0.3)
		plt =	plot!(pltvarsDM[2:end,1].-2,pltvarsDM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow", alpha=0.3)
		plt =	plot!(pltvarsD_RA[2:end,1].-1,pltvarsD_RA[2:end,i], label=false, seriestype=:scatter, markersize=2, markerstrokewidth=0, color="cyan", alpha=0.1) # Any value other than 1 makes all values 1...
			#   display(plt)		
	# sleep(.5)
		push!(groupplots, plt)
	end
	l = @layout([a{0.1w} [a d e; f g h ; i j k ; l m n; o p q; r s t; u v w]])
	plot((groupplots...), layout=l, size=(1000,1200), left_margin=2*Plots.mm, bottom_margin=-1*Plots.mm, top_margin=3*Plots.mm, show=true)
	savefig(replace("Sense_"*names(pltvarsT_Y)[2]*"-"*names(pltvarsT_Y)[end]*".png",":" => ""))
end	

### Make plots in groups of up to 21 (3 x 7), less if the section runs out
# SensePlot(2, 22) # Y
# SensePlot(23, 43) # Y
# SensePlot(44, 64) # Y
# SensePlot(65, 72) # Y
# SensePlot(73, 93) # A
# SensePlot(94, 114) # A
# SensePlot(115, 135) # A
# SensePlot(136, 143) # A
# SensePlot(144, 145) # MS
# SensePlot(146, 166) # PA
# SensePlot(167, 187) # PA
# SensePlot(188, 208) # PA
# SensePlot(209, 213) # PA
# SensePlot(214, 234) # PY
# SensePlot(235, 255) # PY
# SensePlot(256, 276) # PY
# SensePlot(277, 284) # PY
# SensePlot(285, 289) # PM -> PFX

## Function to plot non-contiguous set of 21 variables, to show as a sample
function SensePlotInd(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u, filemane, title="Monte Carlo Sensitivity")
	pltvarsT_Y = [SenseoutputT_Y[:,:t_elas_y] SenseoutputT_Y[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]]
	pltvarsY = [SenseoutputY[:,:elas_y] SenseoutputY[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]]
	# pltvarsVA = [SenseoutputVA[:,:elas_va] SenseoutputVA[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #1
	# pltvarsT_M = [SenseoutputT_M[:,:t_elas_m] SenseoutputT_M[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]]
	# pltvarsM = [SenseoutputM[:,:elas_m] SenseoutputM[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] 
	pltvarsT_A = [SenseoutputT_A[:,:t_elas_a] SenseoutputT_A[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #2
	pltvarsA = [SenseoutputA[:,:elas_a] SenseoutputA[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] 
	pltvarsDM = [SenseoutputDM[:,:elas_dm] SenseoutputDM[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #2
	pltvarsD_RA = [SenseoutputD_RA[:,:d_elas_ra] SenseoutputD_RA[:,[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u]]] #1
	# Generate a legend for the whole group of plots as a separate 'plot'
	groupplots =[]
	plt = scatter(1,markersize=12, markerstrokewidth=2, size=(300,400), xaxis=false, yaxis=false,
	label="   Baseline\n   Counterfactual\n   Value", color="green", legendfontsize=9,legend_title="Elasticity\n  Parameters",
	legend_title_font_size=11, fg_legend=:transparent,labelspacing=0.8, framestyle=:none,
	fontfamily="Computer Modern")
	plot!(1,1, label="   t_elas_y:  "*string(pltvarsT_Y[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightgreen")
	plot!(1,1, label="   elas_y:    "*string(pltvarsY[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange")
	# plot!(1,1, label=" elas_va:   "*string(pltvarsVA[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55")
	# plot!(1,1, label=" t_elas_m:  "*string(pltvarsT_M[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="lightblue")
	# plot!(1,1, label=" elas_m:    "*string(pltvarsM[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue")
	plot!(1,1, label="   t_elas_a:  "*string(pltvarsT_A[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple")
	plot!(1,1, label="   elas_a:    "*string(pltvarsA[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red" )
	plot!(1,1, label="   elas_dm:   "*string(pltvarsDM[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow")
	plot!(1,1, label="   d_elas_ra: "*string(pltvarsD_RA[1,1]), seriestype=:scatter, markersize=3, markerstrokewidth=0, color="cyan") # Any value other than 1 makes all values 1...
	push!(groupplots, plt)
	for i in 2:length(pltvarsA[1,:])
	# 	# Start with original value from the model
		plt =	scatter([pltvarsY[1,1]],[pltvarsY[1,i]], color = "darkgreen", label = false, markersize = 9, 
		# set limits from data
		ylimits=(min(pltvarsA[1,i],minimum(pltvarsA[2:end,i]),minimum(pltvarsY[2:end,i]),minimum(pltvarsDM[2:end,i]),minimum(pltvarsT_A[2:end,i]),minimum(pltvarsT_Y[2:end,i]),minimum(pltvarsD_RA[2:end,i]))-.1,
		# minimum(pltvarsT_M[2:end,i]),minimum(pltvarsVA[2:end,i]),minimum(pltvarsM[2:end,i]),
		max(pltvarsA[1,i],maximum(pltvarsA[2:end,i]),maximum(pltvarsY[2:end,i]),maximum(pltvarsDM[2:end,i]),maximum(pltvarsT_A[2:end,i]),maximum(pltvarsT_Y[2:end,i]),maximum(pltvarsD_RA[2:end,i]))+.1))
		# maximum(pltvarsM[2:end,i]),maximum(pltvarsVA[2:end,i]),maximum(pltvarsT_M[2:end,i]),
		# Add all saved Monte Carlo
		plt =	plot!(pltvarsT_Y[2:end,1],pltvarsT_Y[2:end,i], seriestype=:scatter, title=(replace(names(pltvarsT_Y)[i],"["=>"","]"=>"","("=>"",")"=>"",","=>"")),
				label=false, xaxis = ("Change in each Elasticity of Substitution",font(7)), yaxis=("Diff in Sector activity",font(7)),# ylim=(minimum(Senseoutput[:,i]), maximum(Senseoutput[:,i])),
				titlefontsize=10, markersize=3, markerstrokewidth=0, color="lightgreen", alpha=0.3)
		plt =	plot!(pltvarsY[2:end,1],pltvarsY[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="orange", alpha=0.3)
		# plt =	plot!(pltvarsVA[2:end,1].-1,pltvarsVA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="gray55", alpha=0.3)
		# plt =	plot!(pltvarsT_M[2:end,1],pltvarsT_M[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="light blue", alpha=0.3)
		# plt =	plot!(pltvarsM[2:end,1],pltvarsM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="blue", alpha=0.3)
		plt =	plot!(pltvarsT_A[2:end,1].-2,pltvarsT_A[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="purple", alpha=0.3)
		plt =	plot!(pltvarsA[2:end,1],pltvarsA[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="red", alpha=0.3)
		plt =	plot!(pltvarsDM[2:end,1].-2,pltvarsDM[2:end,i], label=false, seriestype=:scatter, markersize=3, markerstrokewidth=0, color="yellow", alpha=0.3)
		plt =	plot!(pltvarsD_RA[2:end,1].-1,pltvarsD_RA[2:end,i], label=false, seriestype=:scatter, markersize=2, markerstrokewidth=0, color="cyan", alpha=0.1) # Any value other than 1 makes all values 1...
			#   display(plt)		
	# sleep(.5)
		push!(groupplots, plt)
	end
	l = @layout([a{0.1w} [a d e; f g h ; i j k ; l m n; o p q; r s t; u v w]])
	plot((groupplots...), layout=l, size=(1000,1200), left_margin=2*Plots.mm, bottom_margin=-1*Plots.mm, top_margin=3*Plots.mm, show=true, plot_title=title)
# display(dspplt)
	savefig(filemane*"Sens_SubSet.png")
end

### 21 variables chosen as sample, fed in as the 21 arguments
# ,[440+6,  436 + 6, 439 + 6]
SensePlotInd(488,480,489,2,13,16,73,87,129,146,210,204,214,281,275,285,286,287,288,144,145, "BothTax_wAbat_noTariffsElas05","CO2 & CH4 tax, with CH4 abatement & no emission tariffs")
# SensePlotInd(88,44,99,135,136,75,60,114,145,212,206,216,283,277,288,289,143,144,77,37,287)
# SensePlotInd(2,13,16,58,73,84,87,129,146,210,204,214,281,275,285,286,287,288,144,145,289)

# pltvarsY = [SenseoutputY[:,:elas_y] SenseoutputY[:,[441,437,440,2,13,16,73,87,129,146,210,204,214,281,275,285,286,287,288,144,145]]]

# a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u = 441,437,440,2,13,16,73,87,129,146,210,204,214,281,275,285,286,287,288,144,145