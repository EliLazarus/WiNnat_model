    using MPSGE
    # import JuMP

        m = MPSGEModel()
    # Here parameter values are doubled and input data halved from MPSGE version       
    @parameters(m, begin
        # inputcoeff, 2
        endow, 2
        elascoeff, 2
        outputmult, 2
    end)
    
    @sectors(m, begin
        X
        Y
        U
    end)
    
    @commodities(m, begin
        PX
        PY
        PU
        PL
        PK
    end)
    
    @consumer(m, RA)

    @production(m, X, [t = 2, s = 0], begin#s = .5], begin
        @output(PX,100, t)
        @output(PY,25,t)
        @input(PL, 50 , s)
        @input(PK, 75, s)
    end)

    @production(m, Y, [t = 0, s = .3*elascoeff], begin
        @output(PY,50, t)
        @input(PL, 20, s) 
        @input(PK, 30, s)
    end)

    @production(m, U, [t = 0, s = 1], begin
        @output(PU, 175, t)
        @input(PX, 100, s)
        @input(PY, 75,  s)
    end)

    @demand(m, RA, begin
            @final_demand(PU,175)
        end, begin
            @endowment(PL, 35 * endow)
            @endowment(PK, 105)
    end)

    solve!(m, cumulative_iteration=0)

    set_value!(endow, 2.2)
    fix(RA, 182)
    solve!(m)


#    set_value!(tx,1); set_value!(sx, 0)
#     set_value!(tx,.25); set_value!(sx, 0)
#     set_value!(tx,0); set_value!(sx, .4)
#     set_value!(tx,0); set_value!(sx, 1)
#     set_value!(tx,1); set_value!(sx, 1)
#     set_value!(tx,.25); set_value!(sx, .4)
    

#     # zero profit with s = 0.5
#     a(PL,PK,PX,PY) = (( -25 * PY - 100 * PX) + ((75.0 * (((((0.0 + (0.4 * ((PL) ^ 0.5))) + (0.6 * ((PK) ^ 0.5))) ^ 2.0) / (PK)) ^ 0.5)) * PK) + ((50.0 * (((((0.0 + (0.4 * ((PL) ^ 0.5))) + (0.6 * ((PK) ^ 0.5))) ^ 2.0) / (PL)) ^ 0.5)) * PL)) - 0.0
#     zpx( PK, PX,PY, PL)= 75* PK - 100* PX - 25* PY + 50* PL
#     zpx( value(PK), value(PX),value(PY), value(PL))
#     a(value(PL),value(PK),value(PX),value(PY))
#     a(0.897677,1.07504,1.00227,1.00179)
#     b(PL,PK,PX,PY) = ( -25 * PY - 100 * PX) + ((75.0 * (((((0.0 + (0.4 * ((PL) ^ 0.5))) + (0.6 * ((PK) ^ 0.5))) ^ 2.0) / (PK)) ^ 0.5)) * PK)
#     b(0.897677,1.07504,1.00227,1.00179)
# c(PL,PK,PX,PY) =((50.0 * (((((0.0 + (0.4 * ((PL) ^ 0.5))) + (0.6 * ((PK) ^ 0.5))) ^ 2.0) / (PL)) ^ 0.5)) * PL)
# c(0.897677,1.07504,1.00227,1.00179)
# b(1,1,1,1)
# c(1,1,1,1)
# (c(1,1,1,1)-b(1,1,1,1))/(b(0.897677,1.07504,1.00227,1.00179)-c(0.897677,1.07504,1.00227,1.00179))


# d(PYpip,PMtrn , PYtrn , PYtrk , PYair , PYwtt , PYotr)= 47.89999243090504 * PYpip - 441.38467027850635 * PMtrn + 66.55389291647839 * PYtrn + 299.34733851792424 * PYtrk + 3.260336848676209 * PYair + 21.345366411075176 * PYwtt + 2.977743153447304 * PYotr
# d(0.98129409,         0.9576311509725615,        0.961539262,        0.952990573,        0.952383916,        0.957500146,        0.962834479)
# d(0.9812940905227775,  0.9615392616423939, 0.9529905729611398, 0.9523839157641576, 0.9575001460703759, 0.9628344790212803)
# d(value(WiNnat[:PY][:pip]),value(WiNnat[:PM][:trn]) , value(WiNnat[:PY][:trn]) , value(WiNnat[:PY][:trk]) , value(WiNnat[:PY][:air]) , value(WiNnat[:PY][:wtt]) , value(WiNnat[:PY][:otr]))
# zp(PL,PK,PX,PY)= -25* PY + 50 *PL - 100 *PX + 75 *PK
# zp(0.897677,1.07504,1.00227,1.00179)
# zp(1,1,1,1)
-1 * (0 + ((-100 * (((((0 + (0.8 * ((X) ^ 3))) + (0.2 * ((Y) ^ 3))) ^ 0.333333333333333) / (X)) ^ -2)) * M))
-1 * ((-100 * ((((0.8 * X^3 + (0.2 * Y^3))^0.333333333333333) / X) ^ -2)) * M)