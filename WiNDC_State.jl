using MPSGE
using JLD2

cd(dirname(Base.source_path()))
## Load all the data: Data was uploaded and structured into Dicts of DenseAxisArrays with a Julia notebook "national_data.ipynb"
# New data from Mitch Oct 11
P= load(joinpath(@__DIR__,"./data/WiNDC_state_data/parms.jld2"))#["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/WiNDC_state_data/sets.jld2"))#["data"] # load in data from saved Notebook output Dict, named S

year = S[:"YR"]; sectors = S[:"S"]; margins = S[:"M"]; gmargins = S[:"GM"] # Margin related sectors;
regions = S[:"R"]; gsectors = S[:"G"] #alias of s 

cd0 = P[:"cd0"]; x0 = P[:"x0"]; bopdef0 = P[:"bopdef0"]; nd0 = P[:"nd0"]; md0 = P[:"md0"]; 
hhadj0 = P[:"hhadj0"]; ys0 = P[:"ys0"]; gdp0 = P[:"gdp0"]; dd0 = P[:"dd0"]; g0 = P[:"g0"]; 
ld0 = P[:"ld0"]; yh0 = P[:"yh0"]; tm0 = P[:"tm0"]; ty0 = P[:"ty0"]; kd0 = P[:"kd0"]; 
i0 = P[:"i0"]; xn0 = P[:"xn0"]; dm0 = P[:"dm0"]; ta0 = P[:"ta0"]; rx0 = P[:"rx0"]; 
xd0 = P[:"xd0"]; id0 = P[:"id0"]; fe0 = P[:"fe0"]; nm0 = P[:"nm0"]; s0 = P[:"s0"]; c0 = P[:"c0"]; m0 = P[:"m0"]; a0 = P[:"a0"]

yr = year[21] #21 = 2017 ## 1 = 1997 
WState = Model()

ta = add!(WState, Parameter(:ta, indices = (regions,gsectors), value=P[:"ta0"][yr,regions,gsectors].data, description="Tax net subsidy rate on intermediate demand"))
ty = add!(WState, Parameter(:ta, indices = (regions,sectors), value=P[:"ty0"][yr,regions,sectors].data, description="Tax net subsidy rate on intermediate demand"))
tm = add!(WState, Parameter(:ta, indices = (regions,gsectors), value=P[:"tm0"][yr,regions,gsectors].data, description="Tax net subsidy rate on intermediate demand"))

# sectors:
Y = add!(WState, Sector(:Y, indices=(regions,sectors), description = "Production"))
X = add!(WState, Sector(:X, indices=(regions,gsectors), description = "Disposition"))
A = add!(WState, Sector(:A, indices=(regions,gsectors), description = "Absorption"))
C = add!(WState, Sector(:C, indices=(regions,), description = "Aggregate final demand"))
MS = add!(WState, Sector(:MS, indices=(regions,margins), description = "Margin supply"))

# commodities:
PA  = add!(WState, Commodity(:PA, indices=(regions,gsectors), description = "Regional market (input)"))
PY  = add!(WState, Commodity(:PY, indices=(regions,gsectors), description = "Regional market (output)"))
PD  = add!(WState, Commodity(:PY, indices=(regions,gsectors), description = "Local market price"))
PN  = add!(WState, Commodity(:PY, indices=(gsectors,), description = "National market"))
PL  = add!(WState, Commodity(:PY, indices=(regions,), description = "Wage rate"))
PK  = add!(WState, Commodity(:PY, indices=(regions,sectors), description = "Rental rate of capital"))
PM  = add!(WState, Commodity(:PM, indices=(regions,margins), description = "Margin price"))
PC  = add!(WState, Commodity(:PY, indices=(regions,), description = "Consumer price index"))
PFX = add!(WState, Commodity(:PFX, description="Foreign exchnage"))

RA = add!(WState, Consumer(:RA, indices = (regions,), description = "Representative agent")) #benchmark = c0[yr,r],

# for r in regions
#     for s in sectors
#     @production(WState, Y[r,s], 0., 0.,
#         [
#          Output(PY[r,g], sum(ys0[yr,r,s,g]), taxes = [Tax(:($(ty[r,s])*1), RA[r])], price=1-(P[:"ty0"][yr,r,s])) for g in gsectors
#          ],
#         [
#          [Input(PA[r,g], id0[yr,r,g,s]) for g in gsectors];
#             [
#                 Input(Nest(Symbol("va$r$s"), 1.,
#             sum(10.), # TODO actually need the total output
#                 [
#                 Input(PL[r],    sum(ld0[yr,r,sectors]) ),
#                 Input(PK[r,s],  kd0[yr,r,s]) #for s in sectors
#                 ]),
#             sum(10.)) # TODO actually need the total output
#             ]
#         ]
#         )
#     end
# end

for r in regions
    for g in gsectors
    @production(WState, X[r,g], 4, 0, 
    [
        [Output(PFX,     x0[yr,r,g] - rx0[yr,r,g])];
        [Output(PN[g],   xn0[yr,r,g])];
        [Output(PD[r,g], xd0[yr,r,g])]
    ],
    [
        Input(PY[r,g],  s0[yr,r,g])
    ]
    )
    end
end

# for r in regions
#     for g in gsectors
#     @production(WState, A[r,g], 0., 0.,
#         [
#             [Output(PA[r,g], a0[yr,r,g], taxes=[Tax(:($(ta[r,g])*1),RA[r])], price = 1-(P[:"ta0"][yr,r,g]))];
#             [Output(PFX,     rx0[yr,r,g])]
#         ],
#         [
#             [Input(PM[r,m], md0[yr,r,m,g]) for m in margins];
#             [Input(Nest(Symbol("d$r$g"), 2., 
#                 sum(10.), # TODO need to put actual output sum,
#                 [Input(PFX, m0[yr,r,g], taxes = [Tax(:($(tm[r,g])*1), RA[r])], price= 1+(P[:"tm0"][yr,r,g]))] # TODO I think this is meant to be a nest in the dm nest
#         ),
#         sum(10.),
#     )
#             ];
#             [Input(Nest(Symbol("dm$r$g"), 4.,
#                 sum(15.), # TODO need to put actual output sum
#                 [
#                 Input(PN[g],    nd0[yr,r,g]),
#                 Input(PD[r,g],  dd0[yr,r,g])
#                 ]
#         ),
#         sum(15.) # TODO need to put actual output total
#     )
#     ]
#         ]  
#         )
#     end
# end


for r in regions
    for m in margins
    @production(WState, MS[r,m], 0, 0,
    [
      Output(PM[r,m],  sum(md0[yr,r,m,gm] for gm in gmargins))
    ],
      [
          [Input(PN[gm],   nm0[yr,r,gm,m]) for gm in gmargins];
          [Input(PD[r,gm], dm0[yr,r,gm,m]) for gm in gmargins]
      ]
    )
    end
end


for r in regions
    @production(WState, C[r], 0, 1,
            [Output(PC[r],   c0[yr,r])],
            [Input(PA[r,g], cd0[yr,r,g]) for g in gsectors]
    )
end


for r in regions
    add!(WState, DemandFunction(RA[r], 0.,
        [Demand(PC[r], c0[yr,r])],
        [
            [Endowment(PY[r,g], yh0[yr,r,g]) for g in gsectors];
            [Endowment(PFX, bopdef0[yr,r] + hhadj0[yr,r])];
            [Endowment(PA[r,g], -g0[yr,r,g] - i0[yr,r,g]) for g in gsectors];
            [Endowment(PL[r], sum(ld0[yr,r,s] for s in sectors))];
            [Endowment(PK[r,s], kd0[yr,r,s]) for s in sectors]
        ]
        )
    )
end

MPSGE.build(WState)
solve!(WState)