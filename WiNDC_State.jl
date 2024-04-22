using MPSGE
using JLD2

cd(dirname(Base.source_path()))

#parms is too big for Github, so I'm running it from local
P= load(joinpath(@__DIR__,"./data/WiNDC_state_data/parms.jld2"))#["data"] # load in data from saved Notebook output Dict, named P
S= load(joinpath(@__DIR__,"./data/WiNDC_state_data/sets.jld2"))#["data"] # load in data from saved Notebook output Dict, named S

cd0 = P[:"cd0"]; x0 = P[:"x0"]; bopdef0 = P[:"bopdef0"]; nd0 = P[:"nd0"]; md0 = P[:"md0"]; 
hhadj0 = P[:"hhadj0"]; ys0 = P[:"ys0"]; gdp0 = P[:"gdp0"]; dd0 = P[:"dd0"]; g0 = P[:"g0"]; 
ld0 = P[:"ld0"]; yh0 = P[:"yh0"]; tm0 = P[:"tm0"]; ty0 = P[:"ty0"]; kd0 = P[:"kd0"]; 
i0 = P[:"i0"]; xn0 = P[:"xn0"]; dm0 = P[:"dm0"]; ta0 = P[:"ta0"]; rx0 = P[:"rx0"]; 
xd0 = P[:"xd0"]; id0 = P[:"id0"]; fe0 = P[:"fe0"]; nm0 = P[:"nm0"]; s0 = P[:"s0"]; c0 = P[:"c0"]; m0 = P[:"m0"]; a0 = P[:"a0"]

# Indicies of sectors for gmargins to align 
gmind = [1, 5, 8, 9, 11, 13, 15, 19, 24, 26, 29, 32, 33, 39, 40, 42, 44, 47, 49, 50, 51, 52, 54, 57, 58, 59, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71]
# Must use a value in gmind, or manually set the end of gmargins to the index at the position before that number in gmid. e.g. n=12 would be gmargins[1:5]
n = 71

year = S[:"YR"]; sectors = S[:"S"]; margins = S[:"M"]; gmargins = S[:"GM"] # Margin related sectors;
regions = S[:"R"]; gsectors = S[:"G"] #alias of s 

# Filter to the specific year
yr = year[21] #21 = 2017 ## 1 = 1997 
cd0 = cd0[yr,:,:]; x0 = x0[yr,:,:]; bopdef0 = bopdef0[yr,:]; nd0 = nd0[yr,:,:]; md0 = md0[yr,:,:,:];
hhadj0 = hhadj0[yr,:]; ys0 = ys0[yr,:,:,:]; gdp0 = gdp0[yr,:]; dd0 = dd0[yr,:,:]; g0 = g0[yr,:,:]; 
ld0 = ld0[yr,:,:]; yh0 = yh0[yr,:,:]; tm0 = tm0[yr,:,:]; ty0 = ty0[yr,:,:]; kd0 = kd0[yr,:,:]; 
i0 = i0[yr,:,:]; xn0 = xn0[yr,:,:]; dm0 = dm0[yr,:,:,:]; ta0 = ta0[yr,:,:]; rx0 = rx0[yr,:,:]; 
xd0 = xd0[yr,:,:]; id0 = id0[yr,:,:,:]; fe0 = fe0[yr,:]; nm0 = nm0[yr,:,:,:]; s0 = s0[yr,:,:]; c0 = c0[yr,:]; m0 = m0[yr,:,:]; a0 = a0[yr,:,:]; 

sectors = sectors[1:n]
gsectors =gsectors[1:n]
gmargins = gmargins[1:findfirst(x->x==n,gmind)]
WState = Model()

ta = add!(WState, Parameter(:ta, indices = (regions,gsectors), value=P[:"ta0"][yr,regions,gsectors].data, description="Tax net subsidy rate on intermediate demand"))
ty = add!(WState, Parameter(:ty, indices = (regions,sectors), value=P[:"ty0"][yr,regions,sectors].data, description="Tax net subsidy rate on intermediate demand"))
tm = add!(WState, Parameter(:tm, indices = (regions,gsectors), value=P[:"tm0"][yr,regions,gsectors].data, description="Tax net subsidy rate on intermediate demand"))

# sectors:
Y = add!(WState, Sector(:Y, indices=(regions,sectors), description = "Production"))
X = add!(WState, Sector(:X, indices=(regions,gsectors), description = "Disposition"))
A = add!(WState, Sector(:A, indices=(regions,gsectors), description = "Absorption"))
C = add!(WState, Sector(:C, indices=(regions,), description = "Aggregate final demand"))
MS = add!(WState, Sector(:MS, indices=(regions,margins), description = "Margin supply"))

# commodities:
PA  = add!(WState, Commodity(:PA, indices=(regions,gsectors), description = "Regional market (input)"))
PY  = add!(WState, Commodity(:PY, indices=(regions,gsectors), description = "Regional market (output)"))
PD  = add!(WState, Commodity(:PD, indices=(regions,gsectors), description = "Local market price"))
PN  = add!(WState, Commodity(:PN, indices=(gsectors,), description = "National market"))
PL  = add!(WState, Commodity(:PL, indices=(regions,), description = "Wage rate"))
PK  = add!(WState, Commodity(:PK, indices=(regions,sectors), description = "Rental rate of capital"))
PM  = add!(WState, Commodity(:PM, indices=(regions,margins), description = "Margin price"))
PC  = add!(WState, Commodity(:PC, indices=(regions,), description = "Consumer price index"))
PFX = add!(WState, Commodity(:PFX, description="Foreign exchnage"))

RA = add!(WState, Consumer(:RA, indices = (regions,), description = "Representative agent")) #benchmark = c0[yr,r],

for r in regions
    for s in sectors
    @production(WState, Y[r,s], 0., 0.,
        [
         Output(PY[r,g], ys0[r,s,g], taxes = [Tax(:($(ty[r,s])*1), RA[r])], price=1-(ty0[r,s])) for g in gsectors# if ys0[r,s,g]>0
         ],
        [
         [Input(PA[r,g], id0[r,g,s]) for g in gsectors];
            [
                Input(Nest(
                    Symbol("va$r$s"), 
                    1.,
            (ld0[r,s] + kd0[r,s]),
                [
                Input(PL[r],    ld0[r,s]),#;# ,# for r in regions),
                Input(PK[r,s],  kd0[r,s])# for r in regions) for s in sectors
                ]
                ),
                (ld0[r,s] + kd0[r,s])
                )
            ]
        ]
        )
    end
end

for r in regions
    for g in gsectors
    @production(WState, X[r,g], 4, 0, 
    [
        [Output(PFX,     x0[r,g] - rx0[r,g])];
        [Output(PN[g],   xn0[r,g])];
        [Output(PD[r,g], xd0[r,g])]
    ],
    [
        Input(PY[r,g],  s0[r,g])
    ]
    )
    end
end

for r in regions
    for g in gsectors
    @production(WState, A[r,g], 0., 0.,
        [
            [Output(PA[r,g], a0[r,g], taxes=[Tax(:($(ta[r,g])*1),RA[r])], price = 1-(ta0[r,g]))];
            [Output(PFX,     rx0[r,g])]
        ],
        [
            [Input(PM[r,m], md0[r,m,g]) for m in margins];
            [Input(Nest(Symbol("d$r$g"), 2., 
                sum(m0[r,g]), # TODO is this right, re the tax amount (being balanced by the price)?
                [Input(PFX, m0[r,g], taxes = [Tax(:($(tm[r,g])*1), RA[r])], price= 1+(tm0[r,g]))] # TODO I think this is meant to be a nest in the dm nest
        ),
        sum(m0[r,g]),
    )
            ];
            [Input(Nest(Symbol("dm$r$g"), 4.,
                (nd0[r,g] + dd0[r,g]),
                [
                Input(PN[g],    nd0[r,g]),
                Input(PD[r,g],  dd0[r,g])
        ]
        ),
        (nd0[r,g] + dd0[r,g])
    )
    ]
        ]  
        )
    end
end


for r in regions
    for m in margins
    @production(WState, MS[r,m], 0, 0,
    [
      Output(PM[r,m],  sum(md0[r,m,gm] for gm in gmargins))
    ],
      [
          [Input(PN[gm],   nm0[r,gm,m]) for gm in gmargins];
          [Input(PD[r,gm], dm0[r,gm,m]) for gm in gmargins]
      ]
    )
    end
end


for r in regions
    @production(WState, C[r], 0, 1,
            [Output(PC[r],   c0[r])],
            [Input(PA[r,g], cd0[r,g]) for g in gsectors]
    )
end


for r in regions
    add!(WState, DemandFunction(RA[r], 0.,
        [Demand(PC[r], c0[r])],
        [
            [Endowment(PY[r,g], yh0[r,g]) for g in gsectors];
            [Endowment(PFX, bopdef0[r] + hhadj0[r])];
            [Endowment(PA[r,g], -g0[r,g] - i0[r,g]) for g in gsectors];
            [Endowment(PL[r], sum(ld0[r,s] for s in sectors))];
            [Endowment(PK[r,s], kd0[r,s]) for s in sectors]
        ]
        )
    )
end

t=time()
MPSGE.build(WState)
println("n=",n," BuildState: ",(time()-t))

t=time()
solve!(WState, cumulative_iteration_limit=0)
println("n=",n," SolveState0: ",(time()-t))

# keys(DenseAxisArray)