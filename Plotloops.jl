function plottaxemisscurve(tax1, tax2, start, interval, finish ,vec, cnst=1)
    """# runs a loop increasing each tax by \$1/t and then plotting Total GHG (CO2 & CH4) **incorporated** emissions 
    # Arguments are: which tax to change, other tax to either change simultaneously OR keep at 0, st=initial \$ tax value, fin= final \$ tax value,
    # and final (optional) argument can be set to 0 to remove other tax, by default """
    MargEVt = DataFrame(EV = Float64[], Totreduced=Float64[], CH4EV_t=Float64[], CO2EV_t=Float64[],TotEV_t=Float64[] ) # Welfare cost per ton reduced
    margemiss = DataFrame(tax1=Float64[], tax2=Float64[], Emissions=Float64[], MevCES=[], EVCES2=Float64[], #EV_pcnt=[], EV_pcntCES=[], 
    # margemiss = DataFrame(tax1=Float64[], tax2=Float64[], Emissions=Float64[], utilityCES=[], Mev=[], MevCES=[], Equiv_Variarion=Float64[], EVCES2=Float64[], #EV_pcnt=[], EV_pcntCES=[], 
        CH4Emissions=Float64[],CO2Emissions=Float64[], CH4perc_red=[], CO2perc_red=[])
    rename!(margemiss,:tax1=>tax1.name) , rename!(margemiss,:tax2=>tax2.name)
    Testvars = DataFrame(taxrt=Float64[], 
    # Yagr=Float64[],Ycoa=Float64[], 
    Ypip=Float64[], Yoil=Float64[],Ygas=Float64[],
    # Ypet=Float64[], Apip=Float64[],Aoil=Float64[],Ywst=Float64[], CompDYPApip=Float64[],CompDApipPApip=Float64[],DemRAPApip=Float64[],
    # PAagr=Float64[],PAcoa=Float64[],PApip=Float64[],PAoil=Float64[],PAwst=Float64[],PAuel=Float64[],
    # compDApipPAoil=Float64[],compdDAoilPAoil=Float64[],
    # VAMagr=Float64[],VAScoa=Float64[],VAMcoa=Float64[],VASpip=Float64[],VAMpip=Float64[], VASwst=Float64[],VAMwst=Float64[],
    # VASagr=Float64[],VAM10agr=Float64[],VAM100agr=Float64[],VAM500agr=Float64[],VAMkagr=Float64[],Apip=[],
    # VAScoa=Float64[],VAM10coa=Float64[],VAM100coa=Float64[],VAM500coa=Float64[],VAMkcoa=Float64[],
    # VASpip=Float64[],VAM10pip=Float64[],
    VAM100pip=Float64[],VAM500pip=Float64[],
    # VAMkpip=Float64[],
    VASoil=Float64[],VAM10oil=Float64[],VAM20oil=Float64[],VAM50oil=Float64[],VAM100oil=Float64[],VAM500oil=Float64[],VAMkoil=Float64[],
    VASgas=Float64[],VAM10gas=Float64[],VAM20gas=Float64[],VAM50gas=Float64[],VAM100gas=Float64[],VAM500gas=Float64[],VAMkgas=Float64[],
    # VASwst=Float64[],VAM10wst=Float64[],VAM100wst=Float64[],VAM500wst=Float64[],VAMkwst=Float64[],
    TotEm=Float64[],CH4TotEm=Float64[],CO2TotEm=Float64[]
    # CH4emoil=Float64[],CH4empip=Float64[],CO2ecoa=Float64[],CO2emoil=Float64[]
    )
    ResultsTroubleshoot = DataFrame(var=[], value=Float64[], margin=Float64[], x1=Float64[]) 
    for (i,j) in zip(start:interval:finish,vec)
        print(i,":",j,", ")
        if i>530 # this should really only happen for CO2 tax, but ...
        # set_value!(tax1, .0)
        # set_value!(tax2, .0)
        # solve!(MultiNat, output="no"); # Solving back to benchmark to help with numerical computational issues. Necessary?
        end
        set_value!(tax1, i)
        set_value!(tax2, cnst*j)
        solve!(MultiNat, output="no");
        Results = generate_report(MultiNat)
        Results[!,:var] = Symbol.(Results[:,:var]);
       
        push!(Testvars, [i,                
        # value(Y[:agr]), value(Y[:coa]),
        value(Y[:pip]),value(Y[:oil]),value(Y[:gas]),
        # value(Y[:pet]), value(A[:pip]),value(A[:oil]), value(Y[:wst]), value(compensated_demand(MultiNat[:Y][:pip],MultiNat[:PA][:pip])),value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:pip])),value(demand(MultiNat[:RA],MultiNat[:PA][:pip])),
        # value(PA[:agr]),value(PA[:coa]),value(PA[:pip]),value(PA[:oil]),value(PA[:wst]),value(PA[:uel]),
        # value(compensated_demand(MultiNat[:A][:pip],MultiNat[:PA][:oil])),value(compensated_demand(MultiNat[:A][:oil],MultiNat[:PA][:oil])),
        # value(VAM[:agr]),value(VAS[:coa]),value(VAM[:coa]),value(VAS[:pip]),value(VAM[:pip]),
        # value(VAS[:wst]),value(VAM[:wst]), value(VAS[:agr]),value(VAM10[:agr]),value(VAM100[:agr]),value(VAM500[:agr]),value(VAM1000[:agr]),value(A[:pip]),
        # value(VAS[:coa]),value(VAM10[:coa]),value(VAM100[:coa]),value(VAM500[:coa]),value(VAM1000[:coa]),
        # value(VAS[:pip]),value(VAM10[:pip]),
        value(VAM100[:pip]),value(VAM500[:pip]),
        # value(VAM1000[:pip]),
        value(VAS[:oil]),value(VAM10[:oil]),value(VAM20[:oil]),value(VAM50[:oil]),value(VAM100[:oil]),value(VAM500[:oil]),value(VAM1000[:oil]),
        value(VAS[:gas]),value(VAM10[:gas]),value(VAM20[:gas]),value(VAM50[:gas]),value(VAM100[:gas]),value(VAM500[:gas]),value(VAM1000[:gas]),# value(VAS[:wst]),value(VAM10[:wst]),value(VAM100[:wst]),value(VAM500[:wst]),value(VAM1000[:wst]),
        value(TotEm),value(CH4TotEm),value(CO2TotEm)
        # value(CH4em[:oil]),value(CH4em[:pip]),value(CO2em[:coa]),value(CO2em[:oil])
        ]  )
        totrevboth  = -(sum([value(MPSGE.tax_revenue(MultiNat[:Y][i],MultiNat[:RA])) for i in Ip])+
            sum([value(MPSGE.tax_revenue(MultiNat[:A][i],MultiNat[:RA])) for i in [i for i in Ip if i∉[:fbt,:mvt,:gmt]]])+
            sum([value(MPSGE.tax_revenue(MultiNat[:VAS][i],MultiNat[:RA])) for i in Ip])+
            sum([sum([value(MPSGE.tax_revenue(vam[c],MultiNat[:RA])) for c in CH4sectors if VAM_costover[MPSGE.name(vam),c]>1]) for vam in VAMcommodSet]))
        income = totrevboth + sum(va_0[[:surplus,:compen],:])+ only(bopdef_0) -sum(fd_0)
        elasRA = MPSGE.elasticity(MultiNat.productions[:FDem].input)
        totdem = sum([value(FDem)*value(compensated_demand(FDem,PA[i])) for i in Ip])
        # util    = prod([(value(FDem)*value(compensated_demand(FDem,PA[i])))^(pce_0[i,:pce]/sum(pce_0)) for i in Ip])
        utilCESf    = sum([(pce_0[i,:pce]/sum(pce_0))*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
        # utilCESf    = sum([(pce_0[i,:pce]/sum(pce_0))^(1/elasRA)*(value(FDem)*value(compensated_demand(FDem,PA[i])))^((elasRA-1)/elasRA) for i in Ip if value(compensated_demand(FDem,PA[i]))>0])^(elasRA/(elasRA-1))
        # Mev = prod([(1/value(PA[i]))^(pce_0[i,:pce]/sum(pce_0[:,:pce])) for i in Ip])*income
        MevCES = prod([(1/value(PA[i]))^(value(FDem)*value(compensated_demand(FDem,PA[i]))/totdem) for i in Ip])*income
        # EV  = Mev - value(RA)
        # EVCES  =  value(RA) - MevCES
        EVCES2 = -((utilCESf-utilCES)/utilCES)*100
        ResultsTroubleshoot =vcat(ResultsTroubleshoot, [Results fill(i,length(Results[:,1]))])
        push!(margemiss, [i j only(filter(:var => ==(:TotEm), Results)[:, :value])  MevCES  EVCES2  ;;
               only(filter(:var => ==(:CH4TotEm), Results)[:, :value])    only(filter(:var => ==(:CO2TotEm), Results)[:, :value])    ;;
              100*(TotCH4bnchmk-only(filter(:var => ==(:CH4TotEm), Results)[:, :value]))/(TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:,:value])) ;;
              100*(TotCO2bnchmk-only(filter(:var => ==(:CO2TotEm), Results)[:, :value]))/(TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:,:value]))     ])

        EV = value(RA)*EVCES2*10^-2
        TotRedt = TotGHGbnchmk-only(filter(:var => ==(:TotEm), Results)[:, :value])
        EVperton_CH4  = EV/(TotCH4bnchmk-only(filter(:var => ==(:CH4TotEm), Results)[:, :value]))
        EVperton_CO2  = EV/(TotCO2bnchmk-only(filter(:var => ==(:CO2TotEm), Results)[:, :value]))
        EVperton_both = EV/TotRedt
        push!(MargEVt, [EV TotRedt EVperton_CH4 EVperton_CO2 EVperton_both])
    end
    if cnst==0
        tax2in = "only"
    else
        tax2in = " & $tax2"
    end
    return margemiss, Testvars, ResultsTroubleshoot, MargEVt#, plch4, plco2, plt3 #,pla, plm, plp, plo, plw, plpa, plpm, plpp, plpo, plpw, 
    # plcdyp, plcdap, plfdrap, plAp, plAo
end

# set_upper_bound(MultiNat[:A][:pip], 8)
# set_upper_bound(MultiNat[:A][:pet], 4) # 30 is fine past $1000/t
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pet]))
# set_upper_bound(MultiNat[:Y][:sle], 12)
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(Y[:sle]))
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pip]))
# MPSGE.JuMP.delete_upper_bound(MPSGE.get_variable(A[:pet]))
set_silent(MultiNat)
#
###############################
## This Section for manaully tweaked vectors of carbon taxes that fit with increasing CH4 taxes to meet reduction targets in each setting
###############################
# co2vec = [42.68,	39.344,	39.247,	39.203,	39.068,	38.932,	38.797,	38.662,	38.527,	38.393,	38.258,	38.124,	37.989,	37.855,	37.721,	37.587,	37.454,	37.32,	37.165,	37.032,	36.899,	36.767,	36.634,	36.502,	36.194,	36.066,	35.938,	35.81,	35.682,	35.555,	35.427,	35.3,	35.173,	35.046,	34.919,	34.792,	34.356,	34.237,	34.117,	33.997,	33.878,	33.758,	33.639,	33.52,	33.401,	33.282,	33.163,	33.045,	32.721,	32.607,	32.493,	32.38,	32.266,	32.152,	32.039,	31.925,	31.812,	31.699,	31.586,	31.473,	31.337,	31.225,	31.113,	31,	30.888,	30.777,	30.665,	30.553,	30.442,	30.331,	30.219,	30.108,	29.997,	29.886,	29.775,	29.665,	29.554,	29.443,	29.333,	29.223,	29.113,	29.003,	28.893,	28.783,	28.673,	28.564,	28.454,	28.344,	28.235,	28.126,	28.017,	27.908,	27.799,	27.69,	27.581,	27.473,	27.364,	27.256,	27.148,	27.04,	26.932,	26.824,	26.716,	26.608,	26.5,	26.393,	26.286,	26.178,	26.071,	25.964,	25.857,	25.75,	25.643,	25.537,	25.43,	25.324,	25.217,	25.111,	25.005,	24.899,	24.596,	24.494,	24.392,	24.291,	24.189,	24.094,	23.991,	23.888,	23.786,	23.683,	23.581,	23.478,	23.376,	23.274,	23.172,	23.07,	22.968,	22.867,	22.765,	22.663,	22.562,	22.461,	22.359,	22.258,	22.157,	22.056,	21.956,	21.855,	21.754,	21.654,	21.553,	21.453,	21.353,	21.253,	21.153,	21.053,	20.953,	20.853,	20.754,	20.655,	20.555,	20.456,	20.357,	20.258,	20.159,	20.06,	19.961,	19.863,	19.764,	19.666,	19.567,	19.469,	19.371,	19.273,	19.175,	19.077,	18.979,	18.882,	18.784,	18.686,	18.589,	18.492,	18.395,	18.298,	18.201,	18.104,	18.007,	17.911,	17.814,	17.718,	17.621,	17.525,	17.429,	17.333,	17.237,	17.141,	17.045,	16.95,	16.854,	16.759,	16.663,	16.568,	16.473,	16.378,	16.283,	16.188,	16.093,	15.999,	15.904,	15.81,	15.715,	15.621,	15.527,	15.433,	15.339,	15.245,	15.151,	15.057,	14.964,	14.87,	14.777,	14.684,	14.591,	14.498,	14.405,	14.312,	14.219,	14.126,	14.033,	13.941,	13.849,	13.756,	13.664,	13.572,	13.48,	13.388,	13.296,	13.204,	13.113,	13.021,	12.93,	12.838,	12.747,	12.656,	12.565,	12.474,	12.383,	12.292,	12.202,	12.111,	12.021,	11.93,	11.84,	11.75,	11.659,	11.569,	11.48,	11.39,	11.298,	11.208,	11.118,	11.028,	10.938,	10.848,	10.758,	10.669,	10.579,	10.49,	10.401,	10.311,	10.222,	10.133,	10.045,	9.956,	9.867,	9.779,	9.69,	9.602,	9.513,	9.425,	9.337,	9.249,	9.161,	9.073,	8.986,	8.898,	8.811,	8.723,	8.636,	8.549,	8.462,	8.375,	8.288,	8.201,	8.115,	8.028,	7.941,	7.855,	7.769,	7.683,	7.596,	7.51,	7.425,	7.339,	7.253,	7.167,	7.082,	6.997,	6.911,	6.826,	6.741,	6.656,	6.571,	6.486,	6.402,	6.317,	6.232,	6.148,	6.064,	5.979,	5.895,	5.811,	5.727,	5.643,	5.56,	5.476,	5.392,	5.309,	5.226,	5.142,	5.059,	4.976,	4.893,	4.81,	4.727,	4.645,	4.562,	4.48,	4.397,	4.315,	4.233,	4.15,	4.069,	3.987,	3.905,	3.823,	3.741,	3.66,	3.578,	3.497,	3.416,	3.334,	3.253,	3.172,	3.091,	3.011,	2.93,	2.849,	2.769,	2.688,	2.608,	2.528,	2.448,	2.368,	2.288,	2.208,	2.128,	2.048,	1.969,	1.889,	1.81,	1.731,	1.651,	1.572,	1.493,	1.414,	1.336,	1.257,	1.178,	1.1,	1.021,	0.943,	0.864,	0.786,	0.708,	0.63,	0.552,	0.474,	0.397,	0.319,	0.242,	0.164,	0.087,	0.009,	0
# ] #For oil, gas, pip only up to 395 (393.134)

# co2vec = [
#     45.38,45.2,45.1,44.95,44.62,44.5,44.31,44.15,44,43.75,43.6,43.35,43.1,42.95,42.8,42.7,42.55,42.39,42.15,41.85,41.7,41.55,41.4,41.25,41.1,40.97,40.8,40.5,40.35,40.2,40.05,39.75,39.6,39.45,39.31,39.18,39.23,39.1,38.95,38.85,38.7,38.1,37.95,37.85,37.7,37.37,37.25,37.1,36.95,36.85,36.7,36.4,36.6,36.44,36.35,35.7,35.6,35.45,35.35,35.05,35.3,35.18,34.56,34.44,34.33,34.21,33.89,33.77,33.65,33.54,33.42,33.11,32.99,32.88,32.76,32.45,32.33,32.22,32.1,31.78,31.66,31.55,31.43,31.31,31.19,30.97,30.85,30.64,30.52,30.3,30.19,30.08,29.96,29.85,29.64,29.53,29.31,29.2,29.08,28.87,28.75,28.54,28.32,28.21,28.09,27.98,27.76,27.65,27.53,27.32,27.2,26.99,26.88,26.66,26.55,26.34,26.13,26.01,25.9,25.7,25.55,25.45,25.35,25.22,25.05,24.95,24.83,24.6,24.52,24.32,24.1,23.98,23.87,23.67,23.56,23.36,23.35,23.15,23.05,22.84,22.64,22.53,22.34,22.25,22.12,21.82,21.72,21.62,21.43,21.23,21.2,21,20.9,20.7,20.6,20.5,20.3,20.2,20,19.9,19.7,19.59,19.48,19.28,19.19,19,18.88,18.79,18.59,18.5,18.3,18.19,17.98,17.88,17.77,17.68,17.5,17.39,17.28,17.09,17,16.8,16.7,16.5,16.4,16.2,16.1,16,15.9,15.7,15.6,15.5,15.3,15.2,15.1,14.9,14.8,14.7,14.53,14.44,14.24,14.12,14.02,13.82,13.73,13.64,13.45,13.26,13.16,13.07,12.98,12.77,12.66,12.57,12.38,12.29,12.2,12.01,11.92,11.81,11.62,11.53,11.44,11.24,11.15,11.05,10.85,10.75,10.55,10.55,10.36,10.28,10.1,10,9.91,9.72,9.62,9.53,9.34,9.24,9.15,9.06,8.86,8.77,8.67,8.48,8.39,8.29,8.1,8,7.91,7.82,7.63,7.54,7.46,7.27,7.18,7.09,6.9,6.81,6.72,6.53,6.44,6.25,6.26,6.06,5.97,5.88,5.79,5.59,5.5,5.41,5.23,5.14,5.05,4.97,4.78,4.69,4.6,4.42,4.33,4.14,4.04,3.95,3.87,3.78,3.7,3.51,3.42,3.23,3.24,3.05,2.96,2.87,2.79,2.6,2.51,2.32,2.23,2.14,2.05,1.97,1.88,1.8,1.61,1.43,1.45,1.26,1.18,1.1,1.01,0.83,0.74,0.56,0.48,0.39,0.31,0.22,0
# ] # For oil/gas/pip only, No abatement. (up to 319) 
# co2vec = [49.8 49 19 17.7 17 16.5 16 15.3 14.9 14.4 14 13.5 13 12.5 12.1 11.6 9.45 9 8.5 8 7.7 7.2 6.8 6.4 6 5.5 5.2 4.7 4.4 3.8 3.5 3.1 2.8 2.5 2.1 1.7 1.3 .9 .6 .21 0
# ]# 41 steps 0:0.25:10 GWP=81.2, CH4 oil x 5, econ-wide CH4 tax, w Abatement
# co2vec = [49.9,49.1,19,18,17.8,17.4,16.9,16.6,16.1,15.7,15.5,15,14.6,14.3,13.8,13.3,12.9,12.6,12.2,11.8,11.5,11.1,10.25,10,9.5,9,8.7,8.5,8.2,7.8,7.6,7.2,6.8,6.3,6.2,6,5.5,5.2,4.8,4.6,4.3,3.9,3.5,3.3,3,2.8,2.6,2.2,1.8,1.6,1.1,0.8,0.5,0.21,0
# ]# 0:0.25:13.5 GWP=81.2, CH4 oil x 5, oil/gas CH4 tax, w Abatement
# co2vec = [49.9,49.1,48.4,47.4,46.8,45.8,45.1913,44.5826,43.7739,43.1652,42.3565,41.7478,41.1391,40.3304,39.7217,39.013,38.2043,37.4957,36.687,36.0783,35.4696,34.7609,34.1522,33.5435,32.8348,32.2261,31.4174,30.8087,30.2,29.5913,28.9826,28.3739,27.5652,26.9565,26.3478,25.7391,25.1304,24.5217,23.913,23.3043,22.6957,22.087,21.4783,20.8696,20.2609,19.6522,19.0435,18.4348,17.8261,17.2174,16.6087,16.2,15.5913,14.9826,14.3739,13.7652,13.3565,12.7478,12.2391,11.6304,11.0217,10.513,9.9043,9.49565,8.98696,8.57826,8.06957,7.46087,7.05217,6.44348,6.03478,5.42609,5.01739,4.4087,3.8,3.4,3,2.6,2.1,1.6,1.1,0.61,0.1
# ]# 0:0.25:13.5 GWP=81.2, CH4 oil x 5, oil/gas CH4 tax, No Abatement
# co2vec = [43.784,	29.916,	29.022,	28.665,	28.131,	27.776,	27.422,	27.068,	25.839,	25.497,	25.156,	24.816,	24.476,	24.138,	23.801,	23.464,	23.129,	22.794,	22.46,	22.127,	21.795,	21.464,	21.134,	20.804,	20.171,	19.845,	19.519,	19.195,	18.872,	18.548,	18.227,	17.906,	17.587,	17.267,	16.949,	16.632,	16.034,	15.718,	15.406,	15.095,	14.785,	14.476,	14.168,	13.861,	13.555,	13.25,	12.946,	12.642,	11.384,	11.093,	10.803,	10.514,	10.225,	9.937,	9.65,	9.365,	9.08,	8.795,	8.512,	8.215,	7.932,	7.651,	7.371,	7.091,	6.813,	6.535,	6.258,	5.983,	5.708,	5.434,	5.16,	4.885,	3.169,	2.917,	2.66,	2.408,	2.158,	1.909,	1.66,	1.413,	1.166,	0.884,	0.637,	0.392,	0.149,	0
# ] # x 5, & standard

# co2vec = [45.361; 32; 31.5; collect(30.8:-((30.8-20.5)/38):20.5) ; collect(20.1:-((20.1-9.1)/38):9.1) ; collect(8.8:-((8.8-.45)/38):.45);.24;.17;.1;.0
# ] # x 5 oil/gas tax only

# co2vec = [45.2715; 31.9; 31; 30.8; 30; 29.7; 29.3; 28.8; 27.9; 27.35; collect(26.96:-((26.96-17.92)/28):17.92) ; collect(17.65:-((17.65-9.2)/28):9.2) ; collect(8.9:-((8.9-1.75)/25):1.75);1.5;1.3;.747;.36;.26;0.1 ;0
# ]  # Just Set up  GWP, 0:1:75
# co2vec = [52.1865,51.9,43.3,41.1,40.9,40.7,40.5,39.8,39.7,39.4,39.5,39.1,38.8,38.8222,38.6444,38.4667,36.0889,35.5111,35.5333,35.3556,35.1778,35.1,34.9222,34.8144,34.6367,34.4689,34.3011,34.1333,33.9556,33.7778,33.7,33.6222,33.4444,33.2667,33.2889,33.1111,32.9333,32.7556,32.6778,32.5,32.4222,32.2444,32.0667,31.9889,31.9111,31.7333,31.6556,31.4778,30.38,30.2122,30.0444,29.8667,29.69,29.7111,29.5333,29.3556,29.3778,29.2,29.0222,28.9344,28.7667,28.6889,28.5111,28.4333,28.3556,28.1778,28,27.9222,27.7444,27.6667,27.4889,27.4111,26.6233,26.5356,26.3778,26.2,26.1222,26.0444,25.8667,25.7889,25.6111,25.5333,25.3556,25.2778,25.2,25.0222,24.9444,24.7667,24.6889,24.6111,24.5333,24.3556,24.2678,24.2,24.0222,23.8444,22.2467,22.0889,22.0011,21.9133,21.7556,21.6678,21.6,21.4,21.2744,21.1489,21.1133,20.9978,20.8722,20.7467,20.6211,20.4956,20.46,20.3444,20.2189,20.0933,19.9678,19.9322,19.8167,19.6911,19.5656,19.44,19.3144,19.2789,19.1633,18.9378,19.0122,18.8867,18.7611,18.6356,18.51,18.3844,18.2589,18.1333,18.0078,17.8822,17.9567,17.8311,17.7056,17.58,17.4544,17.3289,17.2933,17.1778,14.7522,14.8267,14.7011,14.6756,14.55,14.4244,14.3689,14.2433,14.1178,14.0022,13.8867,13.7711,13.6456,13.52,13.3944,13.2689,13.3433,13.2178,13.0922,12.9667,12.9311,12.8156,12.77,12.6544,12.5389,12.4133,12.3778,12.2622,12.2167,12.1011,11.9756,11.86,11.7344,11.6989,11.5833,11.5378,11.4322,11.3067,11.1811,11.1456,11.03,10.9044,10.8689,10.8333,10.7178,10.5922,10.5567,10.3511,9.4256,8.9,8.69,8.59667,8.49333,8.39,8.28667,8.27333,8.17,8.07667,8.05333,7.95,7.85667,7.75333,7.66,7.55667,7.45333,7.35,7.33667,7.23333,7.14,7.11667,7.01333,6.92,6.82667,6.72333,6.62,6.51667,6.41333,6.4,6.30667,6.27333,6.1,5.99667,5.98333,5.89,5.85667,5.68333,5.58,5.56667,5.47333,5.44,5.26667,5.16333,5.15,5.05667,5.02333,4.25,4.14667,4.04333,4.03,3.93667,3.90333,3.81,3.71667,3.62333,3.52,3.41667,3.31333,3.3,3.20667,3.10333,3,2.98667,2.89333,2.87,2.68667,2.67333,2.58,2.54667,2.37333,2.36,2.26667,2.23333,2.06,2.04667,1.95333,1.92,1.74667,1.73333,1.64,1.55,1.43333,1.33,1.22667,1.32333,1.22,1.11667,1.01333,0.91,0.80667,0.703333,0.77,0.6,0.6,0.4,0.4,0.4,0.2,0.2,0.1,0,0,0,0
# ]  # GWP, 0:.25:75
# 5408.
# co2vec = [52.1865;48;47;collect(46:-(46-2.1)/150:2.1); 2;1.9;1.4;;2;1.4;1;.6;.3;.2;0

# ]
# co2vec = [52.1865,51.97,51.76,51.52,51.3,51.11,50.89,50.68,50.458,50.25,50.047,49.835,49.59,49.3965,49.2031,48.9596,48.7661,48.5227,48.3292,48.1357,47.9223,47.7288,47.4853,47.2919,47.0984,46.855,46.6615,46.458,46.2646,46.0711,45.8776,45.6342,45.4407,45.2472,45.0538,44.8603,44.6268,44.4334,44.2199,43.9964,43.803,43.6095,43.406,43.2126,43.0191,42.8256,42.6322,42.4387,42.2452,42.0318,41.8383,41.6249,41.4314,41.2379,41.0445,40.851,40.6575,40.4641,40.2706,40.0771,39.8837,39.6902,39.4967,39.3033,39.1098,38.9163,38.7229,38.5294,38.3359,38.1425,37.949,37.7755,37.5821,37.3886,37.1951,37.0017,36.8282,36.6348,36.4413,36.2478,36.0544,35.8709,35.6774,35.484,35.3405,35.147,34.9536,34.7801,34.5866,34.4132,34.2197,34.0262,33.8428,33.6493,33.4958,33.3024,33.1089,32.9254,32.782,32.5885,32.415,32.2016,32.0581,31.8647,31.6912,31.4977,31.3343,31.1608,30.9673,30.7739,30.6304,30.4369,30.2735,30.1,29.9,29.7436,29.5871,29.3807,29.2243,29.0678,28.8614,28.705,28.5485,28.3421,28.1856,28.0292,27.8528,27.6663,27.5099,27.3335,27.177,26.9706,26.8042,26.6477,26.4813,26.2949,26.1284,25.972,25.7955,25.6391,25.4827,25.3062,25.1498,24.9734,24.7869,24.6205,24.4941,24.3176,24.1312,23.9748,23.8183,23.6619,23.5054,23.329,23.1426,22.9861,22.8297,22.6733,22.5168,22.3604,22.204,22.0475,21.8711,21.7147,21.5282,21.3718,21.2153,21.0589,20.9025,20.746,20.5896,20.4332,20.2767,20.1203,19.9639,19.8074,19.651,19.4946,19.3381,19.1817,19.0252,18.8688,18.7124,18.5559,18.3995,18.2431,18.0866,17.9302,17.7738,17.6173,17.4909,17.3345,17.198,17.0416,16.8851,16.7287,16.5723,16.4158,16.2594,16.103,15.9465,15.8201,15.6837,15.5272,15.3708,15.2144,15.0579,14.9315,14.795,14.6386,14.4822,14.3257,14.1693,14.0429,13.8864,13.75,13.6,13.4463,13.3127,13.179,12.9953,12.8816,12.728,12.5943,12.4406,12.2969,12.1633,12.0296,11.8759,11.7122,11.5786,11.4449,11.2912,11.1576,10.9939,10.8602,10.7265,10.5929,10.4592,10.3255,10.1718,10.0382,9.8745,9.7408,9.6071,9.47347,9.3398,9.20612,9.07245,8.93878,8.7851,8.62143,8.48776,8.35408,8.22041,8.08673,7.95306,7.81939,7.68571,7.55204,7.41837,7.28469,7.15102,7.01735,6.88367,6.75,6.61633,6.48265,6.34898,6.21531,6.08163,5.94796,5.81429,5.68061,5.54694,5.41327,5.27959,5.14592,5.01224,4.87857,4.7449,4.61122,4.47755,4.37388,4.2602,4.12653,3.99286,3.85918,3.72551,3.59184,3.45816,3.35449,3.22082,3.08714,2.95347,2.8198,2.68612,2.57245,2.43878,2.3051,2.20143,2.08776,1.95408,1.82041,1.68673,1.58306,1.44939,1.33571,1.20204,1.06837,0.93469,0.83102,0.697347,0.583673,0.45,0.35,0.23,0.09,0
# ] # GWP, No Abatement
# co2vec = [45.38,45.1,44.9,44.7,44.4,44.2,43.9,43.7,43.4,43.2,43,42.7,42.5,42.3,42.1,42.1,41.6,41.4125,41.2375,41.05,40.8625,40.675,40.4875,40,39.6889,39.4778,39.3667,39.0056,38.8944,38.5833,38.4722,38.1611,37.95,37.72,37.49,37.36,37.13,36.9,36.67,36.44,36.21,35.98,35.65,35.4786,35.3071,35.1357,34.8643,34.5929,34.4214,34.25,34.0182,33.8364,33.5545,33.3727,33.0909,32.9091,32.7273,32.5455,32.3636,32.0318,31.85,31.6208,31.3917,31.2625,31.0333,30.7542,30.625,30.3958,30.2667,30.0375,29.8083,29.5792,29.35,29.1208,28.9917,28.7625,28.6333,28.3042,28.175,27.9458,27.7167,27.4875,27.3583,27.1292,26.9,26.7708,26.5417,26.3125,26.0833,25.9542,25.725,25.5958,25.3667,25.1375,24.9083,24.7792,24.55,24.4208,24.1917,23.9625,23.7333,23.6042,23.325,23.1958,23.0667,22.8375,22.6083,22.3792,22.35,22.05,21.85,21.65,21.45,21.3,21.1,20.9,20.7,20.55,20.6,20.5,20,19.8,19.6,19.4,19.2,19,18.8,18.6,18.5,18.3111,18.0722,17.9333,17.9944,17.5556,17.4167,17.1778,17.0389,17.1,16.9611,16.7222,16.5833,16.1444,16.2056,15.7667,15.5278,15.3889,15.25,15.0625,14.875,14.6875,14.6,14.4125,14.225,14.0375,13.85,13.6625,13.475,13.2875,13.15,12.9625,12.775,12.5875,12.4,12.2722,12.1444,12.0167,11.6889,11.5611,11.3833,11.2556,11.1278,10.9,10.6722,10.5444,10.4167,10.1889,10.0611,9.8333,9.80556,9.57778,9.35,9.22059,8.99118,8.86176,8.73235,8.50294,8.37353,8.24412,8.01471,7.88529,7.70588,7.57647,7.44706,7.21765,7.08824,6.95882,6.72941,6.6,6.38095,6.3619,6.14286,5.92381,5.80476,5.63571,5.51667,5.39762,5.17857,4.95952,4.94048,4.72143,4.50238,4.38333,4.26429,4.04524,3.92619,3.70714,3.5881,3.41905,3.4,3.1,3.09,2.9,2.69,2.55,2.40143,2.25,2.09,1.96714,1.80571,1.65,1.48286,1.36,1.21,1.05,0.9,0.765714,0.604286,0.47,0.32,0.16,.0,0
# ]  # Standard but No CH4 abatement, 0:1:241

# co2vec = [42.68,38.749,38.346,38.16,36.969,36.805,36.642,36.479,36.316,36.153,35.99,35.828,35.292,35.133,34.974,34.815,34.656,34.498,34.083,33.929,33.776,33.622,33.469,33.316,32.569,32.421,32.273,32.125,31.978,31.83,31.683,31.536,31.389,31.242,31.096,30.95,29.93,29.791,29.652,29.513,29.375,29.237,29.099,28.961,28.823,28.685,28.548,28.411,27.683,27.551,27.418,27.286,27.155,27.023,26.891,26.76,26.629,26.497,26.367,26.236,25.832,25.703,25.573,25.444,25.315,25.186,25.057,24.928,24.799,24.671,24.543,24.415,24.287,24.159,24.031,23.904,23.776,23.649,23.522,23.395,23.269,23.142,23.016,22.889,22.763,22.637,22.511,22.386,22.26,22.135,22.01,21.884,21.759,21.635,21.51,21.386,21.261,21.137,21.013,20.889,20.765,20.642,20.518,20.395,20.272,20.149,20.026,19.903,19.781,19.658,19.536,19.414,19.292,19.17,19.048,18.927,18.805,18.684,18.563,18.442,17.168,17.053,16.937,16.822,16.706,16.591,16.476,16.361,16.247,16.132,16.017,15.903,15.789,15.675,15.561,15.447,15.333,15.219,15.106,14.993,14.88,14.766,14.654,14.541,14.428,14.316,14.203,14.091,13.979,13.867,13.755,13.643,13.531,13.42,13.309,13.197,13.086,12.975,12.864,12.754,12.643,12.533,12.422,12.312,12.202,12.092,11.982,11.873,11.763,11.654,11.544,11.435,11.326,11.217,11.108,11,10.891,10.783,10.674,10.566,10.458,10.35,10.242,10.135,10.027,9.92,9.813,9.705,9.598,9.491,9.385,9.278,9.171,9.065,8.959,8.853,8.747,8.641,8.535,8.429,8.324,8.218,8.113,8.008,7.903,7.798,7.693,7.588,7.484,7.379,7.275,7.171,7.067,6.963,6.859,6.755,6.652,6.548,6.445,6.342,6.239,6.136,6.033,5.93,5.828,5.725,5.623,5.521,5.419,5.317,5.215,5.113,5.011,4.91,4.809,4.707,4.606,4.505,4.404,4.304,4.203,4.102,4.002,3.902,3.802,3.702,3.602,3.502,3.402,3.303,3.203,3.104,3.005,2.906,2.807,2.708,2.609,2.51,2.412,2.314,2.215,2.117,2.019,1.921,1.824,1.726,1.628,1.531,1.434,1.337,1.24,1.143,1.046,0.949,0.853,0.756,0.66,0.564,0.467,0.371,0.276,0.18,0.084,0
# ] # Standard, updated 0,1,283

# co2vec = [
    # 35.195,	31.483,	31.104,	30.926,	29.781,	29.624,	29.467,	29.311,	29.155,	28.999,	28.843,	28.688,	28.184,	28.032,	27.879,	27.727,	27.576,	27.424,	27.028,	26.881,	26.735,	26.588,	26.442,	26.295,	25.597,	25.456,	25.314,	25.173,	25.032,	24.891,	24.75,	24.61,	24.469,	24.329,	24.189,	24.05,	23.098,	22.965,	22.832,	22.7,	22.568,	22.436,	22.304,	22.172,	22.04,	21.909,	21.778,	21.646,	20.967,	20.84,	20.714,	20.588,	20.462,	20.336,	20.21,	20.085,	19.96,	19.834,	19.709,	19.584,	19.208,	19.084,	18.961,	18.837,	18.714,	18.59,	18.467,	18.345,	18.222,	18.099,	17.977,	17.854,	17.732,	17.61,	17.488,	17.367,	17.245,	17.124,	17.003,	16.882,	16.761,	16.64,	16.519,	16.399,	16.278,	16.158,	16.038,	15.918,	15.798,	15.679,	15.559,	15.44,	15.321,	15.202,	15.083,	14.964,	14.846,	14.727,	14.609,	14.491,	14.373,	14.255,	14.137,	14.02,	13.902,	13.785,	13.668,	13.551,	13.434,	13.317,	13.201,	13.084,	12.968,	12.852,	12.736,	12.62,	12.505,	12.389,	12.274,	12.158,	10.979,	10.869,	10.758,	10.648,	10.539,	10.429,	10.319,	10.21,	10.101,	9.991,	9.883,	9.773,	9.665,	9.556,	9.447,	9.339,	9.231,	9.123,	9.015,	8.907,	8.799,	8.692,	8.584,	8.477,	8.369,	8.262,	8.155,	8.049,	7.942,	7.835,	7.729,	7.623,	7.516,	7.41,	7.304,	7.199,	7.093,	6.987,	6.882,	6.777,	6.672,	6.567,	6.462,	6.357,	6.252,	6.148,	6.043,	5.939,	5.835,	5.731,	5.627,	5.524,	5.42,	5.316,	5.213,	5.11,	5.006,	4.903,	4.801,	4.698,	4.595,	4.493,	4.39,	4.288,	4.186,	4.084,	3.982,	3.88,	3.779,	3.677,	3.576,	3.475,	3.374,	3.273,	3.172,	3.071,	2.97,	2.87,	2.77,	2.669,	2.569,	2.469,	2.369,	2.27,	2.17,	2.07,	1.971,	1.872,	1.773,	1.674,	1.575,	1.476,	1.377,	1.279,	1.18,	1.082,	0.984,	0.886,	0.788,	0.69,	0.593,	0.495,	0.398,	0.3,	0.204,	0.107,	0.01,	0
# ] # Standard, for reductions to SCC for CH4, reductions=3053.83

# ch4vec = [
#     1356.16,	1340.699,	1327.198,	1313.877,	1300.736,	1287.755,	1274.924,	1262.244,	1249.783,	1237.423,	1225.232,	1213.231,	1201.449,	1193.022,	1181.692,	1170.472,	1159.372,	1148.392,	1137.502,	1126.752,	1116.083,	1105.513,	1095.073,	1084.684,	1074.434,	1064.245,	1054.185,	1044.205,	1034.305,	1024.466,	1014.727,	1005.077,	995.485,	985.983,	976.554,	967.296,	958.107,	948.993,	939.945,	930.964,	922.05,	913.202,	904.418,	895.692,	887.037,	878.436,	869.895,	861.453,	853.063,	844.733,	836.455,	828.236,	820.066,	811.952,	803.891,	795.876,	787.917,	779.999,	772.142,	764.321,	756.553,	748.838,	741.154,	733.528,	725.945,	718.393,	710.905,	703.442,	696.031,	688.659,	681.325,	674.034,	666.784,	659.576,	652.398,	645.268,	638.174,	631.109,	624.091,	617.108,	610.158,	603.243,	596.946,	596.878,	592.082,	585.594,	579.11,	572.616,	566.154,	559.711,	553.298,	546.912,	540.552,	534.226,	527.923,	521.651,	515.409,	509.193,	503.007,	496.85,	490.722,	484.626,	478.558,	472.513,	466.497,	460.509,	454.548,	448.614,	442.706,	436.826,	430.973,	425.145,	419.341,	413.566,	407.815,	402.089,	396.385,	390.71,	385.058,	379.434,	373.83,	368.252,	362.698,	357.167,	351.66,	346.178,	340.718,	335.281,	329.868,	324.477,	319.109,	313.765,	308.443,	303.144,	297.865,	292.611,	287.376,	282.164,	276.976,	271.81,	266.666,	261.542,	256.438,	251.355,	246.296,	241.256,	236.238,	231.241,	226.266,	221.31,	216.376,	211.464,	206.57,	201.698,	196.847,	192.013,	187.202,	182.411,	177.641,	172.889,	168.158,	163.449,	158.758,	154.085,	149.436,	144.805,	140.195,	135.6,	131.03,	126.48,	121.948,	119.351,	119.345,	119.348,	119.222,	115.72,	111.558,	107.438,	103.238,	99.093,	94.899,	90.775,	86.615,	82.522,	78.398,	74.334,	70.264,	66.212,	62.174,	59.689,	57.477,	53.641,	49.741,	47.749,	47.744,	44.901,	41.206,	37.506,	35.822,	35.806,	35.717,	32.526,	29.041,	25.558,	23.942,	23.868,	21.533,	18.19,	16.774,	13.527,	11.937,	10.712,	7.53,	4.354,	3.833,	3.734,	1.519,	0.564,	0.537,	0.517,	0.518,	0.499,	0.497,	0.495,	0.492,	0.355,	0.159,	0
#     ] # Standard, for reductions to SCC for CO2, reductions=3053.83

# checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 1, 400,#284,#Int(ceil(200 * 1.130480652)),#283,
# zeros(401)) 
# # # # # # # # collect(0:1:491))
# # # # reverse(co2vec))
# co2vec)
# # EVboth = checkch4CO2[4]
# # print(EVboth)
# # CH4margemiss=checkch4CO2[1]
# # print(CH4margemiss)
# # EVCH4 = checkch4CO2[4]
# # print(EVCH4)
# # checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, 284, reverse(co2vec), 0) #Breaks at $540
# # CO2margemiss=checkCO2[1]
# # # print(CO2margemiss)
# # EVCO2 = checkCO2[4]
# # print(EVCO2)
# # checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, Int(ceil(200 * 1.130480652)), ch4vec)
# # ###checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, 400, zeros(401), 0) #Breaks at $540

# ### Loop to generate full set of combination tax reductions
# # n=226
# # # overcountdf = DataFrame()
# # # for t in 0:n
# # # loop = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, n+1, fill(t,n+1)) #n+1
# # # append!(overcountdf, loop[1])
# # # end
# # #### CSV.write("overcount226.csv", overcountdf, missingstring="missing", bom=true)

# # overcountdf = CSV.read("overcount226.csv", DataFrame)
# # CListB = reshape( range(colorant"darkblue", stop=colorant"lightblue",length=n), 1, n );
# # # Marginal Reductions of CO₂ with simulateous CH₄ taxes at increasing levels
# # ### Plot Reductions with 0 tax on CH₄
# # p13 = (9,"Palatino Roman")
# # plot(overcountdf[1:n,:CO₂_tax],(TotGHGbnchmk .- overcountdf[1:n,:Emissions] .- (TotGHGbnchmk .-overcountdf[1,:Emissions]))*10^3, color=:darkblue, linewidth=3, 
# # legend=false, colorbar=true, yaxis=("marginal GHG reductions\n(MMt CO₂eq)", font(p13)), guidefont=(15,"Palatino Roman"), xtickfont=p13, ytickfont=p13, xlabel="CO₂ tax level", left_margin=2mm, bottom_margin=3mm)
# # for l in n+1:n:length(overcountdf[:,1])-n # TODO This has to change after the larger run, AND the +99 etc below
# #     plot!(overcountdf[l:l+n-1,:CO₂_tax],(TotGHGbnchmk .- overcountdf[l:l+n-1,:Emissions].- (TotGHGbnchmk .-overcountdf[l,:Emissions]))*10^3, color=CListB[Int(round(l/n,digits=0))])
# # end
# # plot!(overcountdf[length(overcountdf[:,1])-n+1:length(overcountdf[:,1]),:CO₂_tax],(TotGHGbnchmk .- overcountdf[length(overcountdf[:,1])-n+1:length(overcountdf[:,1]),:Emissions].- (TotGHGbnchmk .-overcountdf[length(overcountdf[:,1])-n+1,:Emissions]))*10^3, color=:lightblue)
# # Plots.plot!([NaN], [NaN],line_z=range(1.0, stop=n, length=n),  c=cgrad([:dargreen,:light]), legend=false, colorbarxpad=10,
# # # colorbartitle="pre-exisiting CH₄ tax level",colorbar_titlefontsize=13,colorbartickvals=[n:-1:0], colorbarlabelalias= [0:1:n], colorbar_titlefontfamily="Palatino Roman", colorbartitle_leftmargin=20mm)
# # colorbartitle="",colorbar_titlefontsize=12,colorbartickvals=[n:-1:0], colorbarlabelalias= [0:1:n], colorbar_titlefontfamily="Palatino Roman",rightmargin=12mm)
# # annotate!(290,1500,Plots.text("pre-exisiting CH₄ tax level", 18,"Palatino Roman", rotation=90))
# # 1
# # png(joinpath(@__DIR__,"./Results/GHGreductionsCO2_w_CH4taxbigtext"))
# ### % overcounted (difference with CO2 taxes to reductions from CH4 taxes with 0 CO2 taxes [over reductions from CH4 taxes with 0 CO2 taxes])
# # (overcountdfCH4[24971:25197,:Emissions] .-overcountdfCH4[1:227,:Emissions]) ./overcountdfCH4[1:227,:Emissions]

# # ###############################################################################
# # #### Marginal Reductions of CH₄ with simulateous CO₂ taxes at increasing levels
# # ###############################################################################
# # ### Copy combined results sorted for the marginal CH₄ tax reductions
# # overcountdfCH4 = sort(deepcopy(overcountdf),[:CO₂_tax,:CH₄_tax])
# # # ### Set up color gradient
# # CListG = reshape( range(colorant"darkgreen", stop=colorant"lightgreen",length=n+1), 1, n+1 );

# # ### Plot Reductions with 0 tax on CO₂ 
# # plot(overcountdfCH4[1:n+1,:CH₄_tax],(TotGHGbnchmk .- overcountdfCH4[1:n+1,:Emissions] .- (TotGHGbnchmk .-overcountdfCH4[1,:Emissions]))*10^3, color=:darkgreen, linewidth=3, 
# # legend=false, colorbar=true, ylim=(0,3048), ylabel="marginal GHG reductions\n(MMt CO₂eq)", guidefont=(10,"Palatino Roman"), xtickfont=p13, ytickfont=p13, xlabel="CH₄ tax level", left_margin=2mm, bottom_margin=3mm)
# # # ### Loop to generate lines for each marginal reduction (additional reduction of CH₄ beyond what we already get from a CO₂ tax at each level)
# # for l in n+2:n+1:length(overcountdfCH4[:,1])-n-1 
# #     plot!(overcountdfCH4[l:l+n,:CH₄_tax],(TotGHGbnchmk .- overcountdfCH4[l:l+n,:Emissions].- (TotGHGbnchmk .-overcountdfCH4[l,:Emissions]))*10^3, color=CListG[Int(round(l/(n+1),digits=0))])
# # end
# # plot!(overcountdfCH4[length(overcountdfCH4[:,1])-n:length(overcountdfCH4[:,1]),:CH₄_tax],(TotGHGbnchmk .- overcountdfCH4[length(overcountdfCH4[:,1])-n:length(overcountdf[:,1]),:Emissions].- (TotGHGbnchmk .-overcountdfCH4[length(overcountdf[:,1])-n,:Emissions]))*10^3, color=:honeydew)
# # Plots.plot!(collect(0:.-10),collect(0:.-10),line_z=collect(226:-1:0),  c=cgrad([:darkgreen,:lightgreen]),# colorbartitle="pre-exisiting CO₂ tax level", )#, rev=true]))
# # colorbartitle="",colorbar_titlefontsize=12,colorbartickvals=[n:-1:0], colorbarlabelalias= [0:1:n], rightmargin=14mm)
# # annotate!(278,1500,Plots.text("pre-exisiting CO₂ tax level", 12,"Palatino Roman", rotation=90))

# # png(joinpath(@__DIR__,"./Results/GHGreductions_CH4_0CO2tax_matchaxis"))

# # # # # # # # # EVdf2 = checkch4CO2[1]
# # # # # # # # EVdf_slice2= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget) && x.Emissions >(TotGHGbnchmk*10^3-ReductTarget-1),EVdf2[:,[1,2,3,11]])# 4329.824705730001
# # # # # # # # print(EVdf_slice2)
# # # # # # checkch4CO2[7]
# # # # # # # print(EVdf2)
# # # # # # checkch4CO2[1]
# # # # # # ##### png(checkch4CO2[7], "./Results/Bothtax-Allemiss")
# # # # # ##### checkch4CO2[2]
# # # # # # EVdf = checkch4CO2[1]
# # # # # # ##### print("checkch4CO2",checkch4CO2[3])
# # # # # # EVdf_slice= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget) && x.Emissions >(TotGHGbnchmk*10^3-ReductTarget-1),EVdf)# 4329.824705730001
# # # # # # EVdf_cap= filter(x -> x.Emissions <(TotGHGbnchmk*10^3-ReductTarget),EVdf)

# # resultdf = copy(checkCO2[1])#[1:1300,:]
# resultdf = copy(checkch4CO2[1])
# tax1 = names(resultdf)[1]; tax2 = names(resultdf)[2]
# # # print(resultdf)
# # # 1
# # # # # plt = plot(resultdf[!,tax1], resultdf[!,:Emissions].*10^3,  label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) & $(replace(tax2,"_"=>" "))  \$/t", xlims=(0,resultdf[end,:CH₄_tax]))#title= "RA:\$$(value(RA)) fxd:$isfixed",)
# # # # # plt = plot!([CO2_taxrate], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3))
# # # # # plch4 = plot(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, label=false, title="CH₄ Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(names(resultdf)[1]) $(names(resultdf)[2]) \$/t", xlims=(0,resultdf[end,:CH₄_tax]))
# # # # # plch4 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
# # # # # plco2 = plot(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label=false, title="CO₂ Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(names(resultdf)[1]) $(names(resultdf)[2])  \$/t", xlims=(0,resultdf[end,:CH₄_tax]))
# # # # # plco2 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
# # shapes = [:none, :auto, :circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

# # # for i in 1:length(shapes)
# pltEV = plot(repeat([0],length(resultdf[:,1])), repeat([0],length(resultdf[:,1])), label=false, ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", #title= "Emissions with $(names(resultdf)[1]) $(names(resultdf)[2])",
# xlims=(0,resultdf[end,tax1]), color=:darkgrey, linewidth=2, ylab="MMt CO₂eq", linestyle=:dashdot,
# yguidefontsize=9, legendfont=font(10,"Palatino Roman"),guidefont=font(12,"Palatino Roman"), tickfont=(9,"Palantino Roman"))

# pltEV = plot!(resultdf[!,tax1], repeat([TotGHGbnchmk*10^3],length(resultdf[:,1])), label="Benchmark Emissions",  color=:darkgrey, linewidth=2, linestyle=:dashdot, legend=false)#, yguidefontsize=9, legendfont=font(9,"Palatino Roman"),guidefont=font(10,"Palatino Roman"))
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:Emissions].*10^3, label="Total GHG Emissions", color=:black, linewidth=2, linestyle=:dashdot, legend=false)#:topright)#legend_position=(.09,1), legendfont=font(14,"Palatino Roman"),guidefont=font(14,"Palatino Roman"))
# # pltEV = plot!(resultdf[!,tax1], repeat([TotGHGbnchmk*10^3],length(resultdf[:,1])), label="Benchmark Emissions", linestyle=:dashdot, linewidth=2, color=:darkgrey)
# reduction_target = 1117.47  # -961.344/(TotGHGbnchmk*10^3) for CH4 SCC, # 3053.83/(TotGHGbnchmk*10^3)
# # pltEV = plot!(resultdf[!,tax1],repeat([TotGHGbnchmk*10^3-reduction_target],length(resultdf[:,tax1])), color=:lightgreen, label="Reduction target",linewidth=1.5)
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label="CO₂ Emissions", color=:blue, linewidth= 4)#linestyle=:dash)#, ylim=(0,5000))
# pltEV = plot!(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, label="CH₄ Emissions", linewidth=4, color=:darkgreen)#, linestyle=:dashdotdot)
# pltEV = plot!(resultdf[!,tax1],repeat([TotGHGbnchmk*10^3-ReductTarget],length(resultdf[:,tax1])), color=:yellow, label="Reduction target",linewidth=1.5)
# # pltEV = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3+200), color=:red, linewidth=2)
# pltEV = scatter!([NaN], [NaN], color = :red, label="Equivalent Variation % change (right axis)", markershape= :circle,markersize= 10, markerstrokewidth=10, markercolor="red") ## EV line

# #     # pltEV = plot!(resultdf[!,tax1],resultdf[!,:EVCES2], legend=:right, label="Equivalent Variation:\n% change", ylabel= "Percentage change",
# #     # guidefont=font(9,"Palatino Roman"), linewidth=2, xlim=(0,resultdf[end,tax1]), ylim=(0,2),color=:brown,markershape=:x, markersize=10,
# #     # yticks=([0.0,0.5,1.0,1.5,2.0],["0.0%","0.5%","1.%","1.5%","2.0%"]), legendfont=font(9,"Palatino Roman"))
# ### Supplementary spaced line for EV, to use as markers
# setupEV1 = resultdf[!,:EVCES2]
# setupEV = setupEV1[collect(1:10:length(setupEV1)),:]
# setuptaxes = resultdf[:,tax1][collect(1:10:length(resultdf[:,tax1])),:]
  
# pltEV = plot!(twinx(),resultdf[:,tax1],setupEV1,
#  xlim=(0,resultdf[end,tax1]), legend=false,color=:brown,ylab="Total cost as % of initial income",
#  yguidefontsize=9,guidefont=font(12,"Palatino Roman"), tickfont=(9,"Palantino Roman"),
# #legend=:topright,label="CO2 spillover %\nright axis", legendfont=font(9,"Palatino Roman"),
# ylim=(0,1),yticks=([0,.20,.40,.60,.80,1],["-0.0%","-0.2%","-0.4%","-0.6%","-0.8%","-1.0%"]),ytickfontsize=9,ytickfont="Palantino Roman",legendfont=font(9,"Palatino Roman"))
# pltEV = plot!(twinx(), setuptaxes, setupEV,ylim=(0,1),line=false,yticks=false,#yticks=([0,.20,.40,.60,.80,1],["0%","20%","40%","60%","80%","100%"]), tickfontsize=11,,legendfont=font(11,"Palatino Roman")
# xlim=(0,resultdf[end,tax1]+1),markershape= :circle,markersize= 3, markerstrokewidth=3, markerstrokecolor="red",legend=false)#,yticks=([0,.20,.40,.60,.80,1],["0%","20%","40%","60%","80%","100%"]),tickfontsize=11,legendfont=font(11,"Palatino Roman"))
# # # # display(pltEV)
# # # # sleep(2)
# # # # end

# # pltEV = plot!(twinx(),resultdf[!,tax1],resultdf[!,:EVCES2], legend=:right, label="Excess Burden:\nRight axis", ylabel= "Percentage change",
# #  guidefont=font(9,"Palatino Roman"), linewidth=2., xlim=(0,resultdf[end,tax1]), color=:brown,markershape=:x, markersize=10,
# #  ylim=(0,1), yticks=([0.0,0.2,0.4,0.6,0.8,1.0],["0.0%","0.2%","0.4%","0.6%","0.8%","1.0%"]), legendfont=font(9,"Palatino Roman"))

# pltEV = plot!(twiny(),resultdf[!,tax2],resultdf[!,:Emissions].*10^3, xflip=true, xlim=(resultdf[end,tax2],resultdf[1,tax2]),xtickfont=(10,"Palantino Roman"),xlabel="$(replace(tax2,"_"=>" ")) \$/t", linewidth=0, legend=false, guidefont=font(12,"Palatino Roman"),xticks=([0,maximum(resultdf[!,tax2])],["\$0","\$$(Int(round(maximum(resultdf[!,tax2]);digits=0)))"]))

# png(pltEV, joinpath(@__DIR__,"./Results/Target_noEV_oilgasonly"))

# # # #Not working
# # # # pltEV = Plots.scatter!(twiny(),[resultdf[120,tax2]],[resultdf[120,:CO2Emissions].*10^3],  xflip=true, xlim=(resultdf[end,tax2],resultdf[1,tax2]),xticks=false,legend=false, guidefont=font(10,"Palatino Roman"))

# # # # ### % Reduction version 
# pltEV = plot(repeat([0],length(resultdf[:,1])), repeat([0],length(resultdf[:,1])), label=false, ylim=(-1,0), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", yticks=([-1,-0.8,-0.6,-0.4,-0.2,0],["-100%","-80%","-60%","-40%","-20%","0%"]),#title= "Emissions with $(names(resultdf)[1]) $(names(resultdf)[2])",
# xlims=(0,resultdf[end,tax1]), color=:darkgrey, linewidth=2, ylab="% change in emissions", linestyle=:dashdot,
# yguidefontsize=12, legendfont=font(9,"Palatino Roman"),guidefont=font(11,"Palatino Roman"), tickfont=(12,"Palantino Roman"))
# # pltEV = plot!(resultdf[!,tax1], (resultdf[!,:CO2Emissions].-resultdf[1,:CO2Emissions])./resultdf[1,:CO2Emissions], label="CO₂ Emissions", color=:blue, linewidth= 4)#linestyle=:dash)#, ylim=(0,5000))
# # pltEV = plot!(resultdf[!,tax1], (resultdf[!,:CH4Emissions].-resultdf[1,:CH4Emissions])./resultdf[1,:CH4Emissions], label="CH₄ Emissions", linewidth=4, color=:darkgreen)#, linestyle=:dashdotdot)
# # pltEV = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(-1,0), color=:red, linewidth=2)


# pltEV = plot(resultdf[!,tax1], (resultdf[!,:Emissions].-resultdf[1,:Emissions])./resultdf[1,:Emissions], legend=:bottomleft, label="Total GHG Emissions", ylim=(-1,.1), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq",yticks=([-1,-0.8,-0.6,-0.4,-0.2,0],	["-100%","-80%","-60%","-40%","-20%","0%"]), #title= "Emissions with $(names(resultdf)[1]) $(names(resultdf)[2])", 
# xlims=(0,resultdf[end,tax1]), color=:black, linewidth=1, ylab="% change in emissions", linestyle=:solid,
# yguidefontsize=12, legendfont=font(10,"Palatino Roman"),guidefont=font(11,"Palatino Roman"), tickfont=(12,"Palantino Roman"),
# # yguidefontsize=11, legendfont=font(9,"Palatino Roman"),guidefont=font(11,"Palatino Roman")
# )
# pltEV = plot!(resultdf[!,tax1], (resultdf[!,:CH4Emissions].-resultdf[1,:CH4Emissions])./resultdf[1,:CH4Emissions], label="CH₄ Emissions", color=:darkgreen, linestyle=:dash)
# pltEV = plot!(resultdf[!,tax1], (resultdf[!,:CO2Emissions].-resultdf[1,:CO2Emissions])./resultdf[1,:CO2Emissions], label="CO₂ Emissions", color=:blue, linestyle=:dashdot, linewidth=1)#, ylim=(0,5000))
# # # pltEV = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(-1,.1), color=:red, linewidth=0.2)
# # ## Legend entry for right y-axis series
# pltEV = plot!([NaN], [NaN], label="Spillovers: CO₂ % of reductions (right axis)", color=:orange, linewidth=2,markershape=:star4,markersize= 28, markercolor="orange",markerstrokecolor="orange")

# # # # # # # ### Next line for CH4 spillovers, only valid for CO2-only policy
# # # # pltEV = plot!(twinx(),resultdf[:,tax1],[resultdf[2,:CH4perc_red]; resultdf[2:end,:CH4perc_red]],legend=:right, label="CH₄ spillover %\nright axis", color=:orange, linewidth=2,legendfont=font(9,"Palatino Roman"), 
# # # # ylim=(0,100), xlims=(0,resultdf[end,tax1]),yticks=([0,10,20,30,40,50,60,70,80,90,100],["0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"]))
# # # # pltEV = plot!(twinx(),resultdf[2:end,tax1],(resultdf[1,:CH4Emissions] .-resultdf[2:end,:CH4Emissions])./(resultdf[1,:Emissions] .-resultdf[2:end,:Emissions]), label="CO2 spillover %\nright axis", color=:orange, legend=:left)#, ylim=(0,1),xticks=([0,.10,.20,.30,.40,.50,.60,.70],["0%","10%","20%","30%","40%","50%","60%","70%"])

# # # # # # # ### Next line for CO2 spillovers, only valid for CH4-only policy
# # # # # # #  % OF emission reduction version
# setup1 = [0; (resultdf[1,:CO2Emissions] .-resultdf[2:end,:CO2Emissions])./(resultdf[1,:Emissions] .-resultdf[2:end,:Emissions])] # CO2 spillovers version
# # setup1 = [0; (resultdf[1,:CH4Emissions] .-resultdf[2:end,:CH4Emissions])./(resultdf[1,:Emissions] .-resultdf[2:end,:Emissions])] # CH4 Spillovers version

# # # # # # # Alternate: % OVER individual GHG version
# # # # setup1 = [0; (resultdf[1,:CO2Emissions] .-resultdf[2:end,:CO2Emissions])./(resultdf[1,:CH4Emissions] .-resultdf[2:end,:CH4Emissions])] # CO2 spillovers from CH4 tax version
# # # # setup1 = [0; (resultdf[1,:CH4Emissions] .-resultdf[2:end,:CH4Emissions])./(resultdf[1,:CO2Emissions] .-resultdf[2:end,:CO2Emissions])] # CH4 Spillovers from CO2 tax version

# setup = setup1[collect(2:20:length(setup1)),:]
# setuptaxes = resultdf[:,tax1][collect(2:20:length(resultdf[:,tax1])),:]

# pltEV = plot!(twinx(),resultdf[2:end,tax1],setup1[2:end],
#  xlim=(0,resultdf[end,tax1]), legend=false,color=:orange,#legend=:topright,label="CO2 spillover %\nright axis", legendfont=font(9,"Palatino Roman"),
#  y_foreground_color_text=:black, y_foreground_color_axis=:orange,y_foreground_color_border=:orange,y_guidefontcolor=:orange,
#  ylim=(0,1),yticks=([0,.20,.40,.60,.80,1],["0%","20%","40%","60%","80%","100%"]),ytickfontsize=15,ytickfont="Palantino Roman",legendfont=font(11,"Palatino Roman"))
# # # # # # # # % over individual reductions version (yticks) 
# # # # y_foreground_color_text=:black, y_foreground_color_axis=:orange,y_foreground_color_border=:orange,y_guidefontcolor=:orange,
# # # # ylim=(0,3.1),yticks=([0,.5,1,1.5,2,2.5,3],["0%","50%","100%","150%","200%","250%","300%"]),ytickfontsize=10,ytickfont="Palantino Roman",legendfont=font(11,"Palatino Roman"))

# pltEV = plot!(twinx(),resultdf[2:end,tax1],setup1[2:end],
#  xlim=(0,resultdf[end,tax1]), color=:orange, legend=false,#legend=:topright,label="CO2 spillover %\nright axis", legendfont=font(9,"Palatino Roman"),
#  y_foreground_color_text=:orange, y_foreground_color_axis=:orange,y_foreground_color_border=:orange,y_guidefontcolor=:orange,
#  ylim=(0,1),yticks=([0,.20,.40,.60,.80,1],["0%","20%","40%","60%","80%","100%"]),ytickfontsize=15,ytickfont="Palantino Roman",legendfont=font(11,"Palatino Roman"))
# # # # # # % over individual reductions version (yticks) 
# # y_foreground_color_text=:orange, y_foreground_color_axis=:orange,y_foreground_color_border=:orange,y_guidefontcolor=:orange,
# # ylim=(0,3.1),yticks=([0,.5,1,1.5,2,2.5,3],["0%","50%","100%","150%","200%","250%","300%"]),ytickfontsize=10,ytickfont="Palantino Roman",legendfont=font(11,"Palatino Roman"))

# pltEV = plot!(twinx(), setuptaxes, setup,ylim=(0,1),line=false,yticks=false,#yticks=([0,.20,.40,.60,.80,1],["0%","20%","40%","60%","80%","100%"]), tickfontsize=11,,legendfont=font(11,"Palatino Roman")
# xlim=(0,resultdf[end,tax1]+1),markershape=:star4,markersize= 8, markercolor="orange",legend=false, color=:orange, markerstrokecolor="orange",
# y_foreground_color_text=:orange, y_foreground_color_axis=:orange,y_foreground_color_border=:orange,y_guidefontcolor=:orange,
# )#,yticks=([0,.20,.40,.60,.80,1],["0%","20%","40%","60%","80%","100%"]),tickfontsize=11,legendfont=font(11,"Palatino Roman"))

# # # # ## Legend only
# pltEVleg = plot([NaN], [NaN], color = :darkgrey, linestyle=:dashdot,label="Benchmark Emissions",legendfont=font(14,"Palatino Roman"), linewidth=1.7)
# pltEVleg = plot!([NaN], [NaN], color = :black,label="Total GHG Emissions",linestyle=:dashdot,legendfont=font(14,"Palatino Roman"), linewidth=1.7)
# reduction_target = 1117.47  # -961.344/(TotGHGbnchmk*10^3) for CH4 SCC, # 3053.83/(TotGHGbnchmk*10^3)
# pltEVleg = plot!([NaN], [NaN], color=:yellow, label="Reduction target",linewidth=1.5)
# pltEVleg = plot!([NaN], [NaN], color = :blue, label="CO₂ Emissions",linewidth=4)
# pltEVleg = plot!([NaN], [NaN], color = :darkgreen, label="CH₄ Emissions",linewidth=4)
# # # pltEVleg = plot!([NaN], [NaN], color = :red, label="SCCO₂",seriestype=:vline,linewidth=2)
# #     #    pltEVleg = scatter!([NaN], [NaN], markersize = 10, color = :orange, label="CH4 spillover% (right axis)",markershape=:star4)
# #     #     pltEVleg = scatter!([NaN], [NaN], markersize = 10, color = :orange, label="CO2 spillover% (right axis)",markershape=:rtriangle)
# pltEVleg = scatter!([NaN], [NaN], color = :red, label="Equivalent Variation, % change (right axis)", markershape= :circle,markersize= 10, markerstrokewidth=3, markerstrokecolor="red")
# display(pltEVleg)
# sleep(3)
# end
# pltEV = plot!(twinx(),resultdf[2:end,tax1],resultdf[2:end,:CO2perc_red],legend=:right, label="CO₂ spillover %\nright axis", legendfont=font(9,"Palatino Roman"), color=:orange, xlim=(0,maximum(resultdf[:,tax1])), ylim=(0,100),yticks=([0,10,20,30,40,50,60,70,80,90,100],["0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"]))

# # TestVars = copy(checkch4CO2[2])
# # TestVars = copy(checkCO2[2])
# print(TestVars)
# # # filter(:EV_pcnt2 => ==(minimum(EVdf_slice.EV_pcnt2)),EVdf_slice)
# # # print(sort(EVdf,:EV_pcnt2)[1:80,:])
# # # print(sort(EVdf_cap,:EV_pcnt2)[1:40,:])

# png(pltEV, joinpath(@__DIR__,"./Results/Spillover_CH4_oilgas_reduct_perc"))
# savefig(pltEV, joinpath(@__DIR__,"./Results/Opt.svg"))
# filter(:EV_pcnt => ==(minimum(EVdf_slice.EV_pcnt)),EVdf_slice)
# print(sort(EVdf,:EV_pcnt)[1:120,:])
# print(sort(EVdf_cap,:EV_pcnt)[1:40,:])
# pltEV = plot(resultdf[!,tax1], resultdf[!,:CO2Emissions].*10^3, label="CO₂ Emissions", color=:blue, linestyle=:dashdot, ylim=(0,5000),


# pltEV = plot(resultdf[!,tax1], resultdf[!,:CH4Emissions].*10^3, ylim=(0,700), label="CH₄ Emissions", color=:darkgreen, linestyle=:dash,
# legend=:left, xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", #title= "Emissions with $(names(resultdf)[1]) $(names(resultdf)[2])",
# xlims=(0,resultdf[end,tax1]),
#  linewidth=1, ylab="MMt CO₂eq",
# yguidefontsize=9, legendfont=font(9,"Palatino Roman"),guidefont=font(10,"Palatino Roman"))

##### checkch4CO2[7]
###### Add Quant diff and % of consumption columns for Final Demand report dataframe
# FDemand[:,:ch4Qdelta]=FDemand[:,:ch4tax].-FDemand[:,:bnch]
# FDemand[:,:CO2Qdelta]=FDemand[:,:co2tax].-FDemand[:,:bnch]
# FDemand[:,:bothQdelta]=FDemand[:,:both].-FDemand[:,:bnch]
# FDemand[:,:bncQpc]=FDemand[:,:bnch]./sum(FDemand[:,:bnch])*100
# FDemand[:,:ch4Qpc]=FDemand[:,:ch4tax]./sum(FDemand[:,:ch4tax])*100
# FDemand[:,:CO2Qpc]=FDemand[:,:co2tax]./sum(FDemand[:,:co2tax])*100
# FDemand[:,:bothQpc]=FDemand[:,:both]./sum(FDemand[:,:both])*100

# fix(RA, sum(fd_0[Ip,:pce]))

# checkCO2[7]
# checkCO2[2] # Total emissions
# # # println("PApip up lim =", upper_bound(MultiNat[:A][:pip]))
# png(checkCO2[7], joinpath(@__DIR__,"./Results/CO2_spills"))

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
# checkch4CO2 = plottaxemisscurve(CH₄_tax, CO₂_tax, 0, 20, 1200, round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]))
# checkch4CO2[7]

# checkch4 = plottaxemisscurve(CH₄_tax,CO₂_tax, 0, 1, 500, zeros(501),round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# resultdfSpill = copy(checkch4[1])
# # # # checkCO2 = plottaxemisscurve(CO₂_tax, CH₄_tax, 0, 1, 400, zeros(401), round(value(MultiNat[:RA]),digits=2), is_fixed(MultiNat[:RA]), 0)
# # # # resultdfSpill = copy(checkCO2[1])
# 1
# tax1 = names(resultdfSpill)[1]; tax2 = names(resultdfSpill)[2]
# pltSpill = plot(resultdfSpill[1:end,tax1], resultdfSpill[1:end,:Emissions].*10^3,  label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$(replace(tax1,"_"=>" ")) \$/t CO₂eq", #title= "Emissions with $(names(resultdfSpill)[1]) $(names(resultdfSpill)[2])",
# xlims=(0,resultdfSpill[end,tax1]), color=:black, linewidth=2, ylab="MMt CO₂eq",legend=:bottom, yguidefontsize=9, guidefont=font(10,"Palatino Roman"), legendfont=font(9,"Palatino Roman"))
# pltSpill = plot!([226.1], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
# pltSpill = plot!(resultdfSpill[1:end,tax1], resultdfSpill[1:end,:CO2Emissions].*10^3, label="CO₂ Emissions", color=:blue)
# pltSpill = plot!(resultdfSpill[1:end,tax1], resultdfSpill[1:end,:CH4Emissions].*10^3, label="CH₄ Emissions", color=:green)
# ### Next line for CH4 spillovers, only valid for CO2-only policy
# # pltSpill = plot!(twinx(),resultdfSpill[:,tax1],resultdfSpill[:,:CH4perc_red], ylim = (4.96,5.04),label="Ch4 spillover %\nright axis", color=:orange, xlims=(0,resultdfSpill[end,tax1]),legend=:right,guidefont=font(10,"Palatino Roman"),yticks=([4.96,4.98,5,5.02,5.04,],["4.96%","4.98%","5%","5.02%","5.04%"]))
# ### Next line for CO2 spillovers, only valid for CH4-only policy
# pltSpill = plot!(twinx(),resultdfSpill[2:end,tax1],resultdfSpill[2:end,:CO2perc_red], label="CO2 spillover %\nright axis", legendfont=font(9,"Palatino Roman"),color=:orange, xlims=(0,resultdfSpill[end,:CH₄_tax]),ylim=(0,100),yticks=([0,10,20,30,40,50,60,70,80,90,100],["0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"]), legend=:right,guidefont=font(10,"Palatino Roman"))
# 1
# savefig(pltSpill, joinpath(@__DIR__,"./Results/CH4-oilgas-Spillovers.svg"))

# checkch4[2]
# checkch4[5]
# checkch4[6]
# checkch4[7]
# png(checkch4[7], joinpath(@__DIR__,"./Results/CH4_spills"))
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

##############################################################
## OTHER plots from function
# plt = plot(margemiss[!,tax1.name], margemiss[!,:Emissions].*10^3, title= "RA:\$$RAval fxd:$isfixed", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
# plt = plot!([CO2_taxrate], seriestype=:vline, label="SCCO₂", ylim=(0,TotGHGbnchmk*10^3))
# plch4 = plot(margemiss[!,tax1.name], margemiss[!,:CH4Emissions].*10^3, label=false, title="CH4 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
# plch4 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCH4bnchmk*10^3)) 
# plco2 = plot(margemiss[!,tax1.name], margemiss[!,:CO2Emissions].*10^3, label=false, title="CO2 Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="$tax1 $tax2in \$/t", xlims=(0,finish))
# plco2 = plot!([CO2_taxrate], seriestype=:vline, label=false, ylim=(0,TotCO2bnchmk*10^3))
# plt3 = plot(margemiss[!,tax1.name], margemiss[!,:Emissions].*10^3, title= "Emissions with $tax1 $tax2in", label="Total GHG Emissions", ylim=(0,TotGHGbnchmk*10^3+200), xlabel="\$/t", 
# xlims=(0,finish), color=:black, linewidth=2, ylab="MMt CO2eq",legend=:left, linestyle=:dashdot)
# plt3 = plot!(margemiss[!,tax1.name], margemiss[!,:CH4Emissions].*10^3, label="CH4 Emissions", color=:green, linestyle=:dot)
# plt3 = plot!(margemiss[!,tax1.name], margemiss[!,:CO2Emissions].*10^3, label="CO2 Emissions", color=:blue, linestyle=:dash)
# # plt3 = plot!([CO2_taxrate], seriestype=:vline, label="SCCO2", ylim=(0,TotGHGbnchmk*10^3), color=:red, linewidth=0.4)
# plt3 = plot!(margemiss[!,tax1.name],repeat([TotGHGbnchmk*10^3-ReductTarget],length(margemiss[:,tax1.name])), color=:yellow, label="target",linewidth=1.5)
# plt3 = plot!(twinx(),margemiss[!,tax1.name],margemiss[!,:EV_pcnt], label="EV%", linewidth=2., xlim=(0,finish), ylim=(0,1),legend=:right)

### Next line for CH4 spillovers, only valid for CO2-only policy
# plt3 = plot!(twinx(),margemiss[2:end,tax1.name],margemiss[2:end,:CH4perc_red], label="Ch4 spillover %\nright axis", color=:orange, yticks=([4.99,5,5.01,5.02,5.03],["4.99%","5%","5.01%","5.02%","5.03%"]), legend=:bottomright)
### Next line for CO2 spillovers, only valid for CH4-only policy
# plt3 = plot!(twinx(),margemiss[2:end,tax1.name],margemiss[2:end,:CO2perc_red], label="CO2 spillover %\nright axis", color=:orange, yticks=([0,10,20,30,40,50,60,70],["0%","10%","20%","30%","40%","50%","60%","70%"]), legend=:topright)
# plcdyp = plot(margemiss[!,tax1.name],Testvars[!,:CompDYPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(Y:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDYPApip]),maximum(Testvars[!,:CompDYPApip])), xlabel="$tax1 $tax2in \$/t")
# plcdap = plot(margemiss[!,tax1.name],Testvars[!,:CompDApipPApip], title= "RA:\$$RAval fxd:$isfixed", label="comp_dem(A:pip,PA:pip)", ylim=(minimum(Testvars[!,:CompDApipPApip]),maximum(Testvars[!,:CompDApipPApip])), xlabel="$tax1 $tax2in \$/t")
# plfdrap = plot(margemiss[!,tax1.name],Testvars[!,:DemRAPApip], title= "RA:\$$RAval fxd:$isfixed", label="final_dem(RA,PA:pip)", ylim=(minimum(Testvars[!,:DemRAPApip]),maximum(Testvars[!,:DemRAPApip])), xlabel="$tax1 $tax2in \$/t")
# pla = plot(margemiss[!,tax1.name],Testvars[!,:Yagr], label=false, title="Y:agr", ylim=(minimum(Testvars[!,:Yagr]),maximum(Testvars[!,:Yagr])), xlabel="$tax1 $tax2in \$/t")
# plm = plot(margemiss[!,tax1.name],Testvars[!,:Ymin], title= "RA:\$$RAval fxd:$isfixed", label="Y:coa", ylim=(minimum(Testvars[!,:Ymin]),maximum(Testvars[!,:Ymin])), xlabel="$tax1 $tax2in \$/t")
# plp = plot(margemiss[!,tax1.name],Testvars[!,:Ypip], title= "RA:\$$RAval fxd:$isfixed", label="Y:pip", ylim=(minimum(Testvars[!,:Ypip]),maximum(Testvars[!,:Ypip])), xlabel="$tax1 $tax2in \$/t")
# plo = plot(margemiss[!,tax1.name],Testvars[!,:Yoil], title= "RA:\$$RAval fxd:$isfixed", label="Y:oil", ylim=(minimum(Testvars[!,:Yoil]),maximum(Testvars[!,:Yoil])), xlabel="$tax1 $tax2in \$/t")
# plw = plot(margemiss[!,tax1.name],Testvars[!,:Ywst], title= "RA:\$$RAval fxd:$isfixed", label="Y:wst", ylim=(minimum(Testvars[!,:Ywst]),maximum(Testvars[!,:Ywst])), xlabel="$tax1 $tax2in \$/t")
# plpa = plot(margemiss[!,tax1.name],Testvars[!,:PAagr], title= "RA:\$$RAval fxd:$isfixed", label="PA:agr", ylim=(minimum(Testvars[!,:PAagr]),maximum(Testvars[!,:PAagr])), xlabel="$tax1 $tax2in \$/t")
# plpm = plot(margemiss[!,tax1.name],Testvars[!,:PAmin], title= "RA:\$$RAval fxd:$isfixed", label="PA:coa", ylim=(minimum(Testvars[!,:PAmin]),maximum(Testvars[!,:PAmin])), xlabel="$tax1 $tax2in \$/t")
# plpp = plot(margemiss[!,tax1.name],Testvars[!,:PApip], title= "RA:\$$RAval fxd:$isfixed", label="PA:pip", ylim=(minimum(Testvars[!,:PApip]),maximum(Testvars[!,:PApip])), xlabel="$tax1 $tax2in \$/t")
# plpo = plot(margemiss[!,tax1.name],Testvars[!,:PAoil], title= "RA:\$$RAval fxd:$isfixed", label="PA:oil", ylim=(minimum(Testvars[!,:PAoil]),maximum(Testvars[!,:PAoil])), xlabel="$tax1 $tax2in \$/t")
# plpw = plot(margemiss[!,tax1.name],Testvars[!,:PAwst], title= "RA:\$$RAval fxd:$isfixed", label="PA:wst", ylim=(minimum(Testvars[!,:PAwst]),maximum(Testvars[!,:PAwst])), xlabel="$tax1 $tax2in \$/t")
# plAp = plot(margemiss[!,tax1.name],Testvars[!,:Apip], title= "RA:\$$RAval fxd:$isfixed", label="A:pip", ylim=(minimum(Testvars[!,:Apip]),maximum(Testvars[!,:Apip])), xlabel="$tax1 $tax2in \$/t")
# plAo = plot(margemiss[!,tax1.name],Testvars[!,:Aoil], title= "RA:\$$RAval fxd:$isfixed", label="A:oil", ylim=(minimum(Testvars[!,:Aoil]),maximum(Testvars[!,:Aoil])), xlabel="$tax1 $tax2in \$/t")     # Or label=false, title="price of oil commodity"
##############################################
