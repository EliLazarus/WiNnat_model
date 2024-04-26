using MPSGE
using JLD2

P= load("C:/Users/Eli/GitFolders/WiNnat_model/data/WiNDC_state_data/parms.jld2")
S= load("C:/Users/Eli/GitFolders/WiNnat_model/data/WiNDC_state_data/sets.jld2")

cd0 = P[:"cd0"]; x0 = P[:"x0"]; bopdef0 = P[:"bopdef0"]; nd0 = P[:"nd0"]; md0 = P[:"md0"]; 
hhadj0 = P[:"hhadj0"]; ys0 = P[:"ys0"]; gdp0 = P[:"gdp0"]; dd0 = P[:"dd0"]; g0 = P[:"g0"]; 
ld0 = P[:"ld0"]; yh0 = P[:"yh0"]; tm0 = P[:"tm0"]; ty0 = P[:"ty0"]; kd0 = P[:"kd0"]; 
i0 = P[:"i0"]; xn0 = P[:"xn0"]; dm0 = P[:"dm0"]; ta0 = P[:"ta0"]; rx0 = P[:"rx0"]; 
xd0 = P[:"xd0"]; id0 = P[:"id0"]; fe0 = P[:"fe0"]; nm0 = P[:"nm0"]; s0 = P[:"s0"]; c0 = P[:"c0"]; m0 = P[:"m0"]; a0 = P[:"a0"]

yr = year[21] #21 = 2017 ## 1 = 1997 

# limit resolution for testing
n = 24
for n in [2, 8, 16, 24, 30, 40, 50]
year = S[:"YR"]; sectors = S[:"S"]; margins = S[:"M"]; gmargins = S[:"GM"] # Margin related sectors;
regions = S[:"R"]; gsectors = S[:"G"] #alias of s 
sectors = sectors[1:n]
gsectors =gsectors[1:n]

DIndTest = Model()

Y = add!(DIndTest, Sector(:Y, indices=(regions,sectors), description = "Production"))
C = add!(DIndTest, Sector(:C, indices=(regions,), description = "Aggregate final demand"))

PY  = add!(DIndTest, Commodity(:PY, indices=(regions,gsectors), description = "Regional market (output)"))
PA  = add!(DIndTest, Commodity(:PA, indices=(regions,gsectors), description = "Regional market (input)"))
PC  = add!(DIndTest, Commodity(:PC, indices=(regions,), description = "Consumer price index"))

RA = add!(DIndTest, Consumer(:RA, indices = (regions,), description = "Representative agent")) #benchmark = c0[yr,r],

for r in regions
    for s in sectors
        # for g in gsectors
        # if ys0[yr,r,s,g] >0 && id0[yr,r,g,s] > 0
    @production(DIndTest, Y[r,s], 0., 0.,
        [Output(PY[r,g], ys0[yr,r,s,g]) for g in gsectors],
         [Input(PA[r,g], id0[yr,r,g,s]) for g in gsectors]
        
    )
        # end
        # end
    end
end

for r in regions
    @production(DIndTest, C[r], 0, 1,
            [Output(PC[r],   c0[yr,r]) ],
            [Input(PY[r,g], cd0[yr,r,g]) for g in gsectors if cd0[yr,r,g]>0]
    )
end


for r in regions
    add!(DIndTest, DemandFunction(RA[r], 0.,
        [Demand(PC[r], c0[yr,r])],
            [Endowment(PA[r,g], -g0[yr,r,g] - i0[yr,r,g]) for g in gsectors]
        )
    )
end

t = time()
MPSGE.build(DIndTest)
println("n=",n," Build: ",time() - t)

t = time()
solve!(DIndTest, cumulative_iteration_limit=0);
println("n=",n," Solve0: ",time() - t)
end