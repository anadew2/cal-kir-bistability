include("../../membrane.jl")
include("../../fixedpoints.jl")

include("../fI_curve/bif_kir.jl")
include("../fI_curve/bif_km.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

###### Step of size (dI_step), and no CaL,no Kir, no KM
    dI_step = 2e-2

    ## Initial conditions
    V_ic = -65. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I1_KirCaLs_fI[1,1]  +dI_step/2 # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = gKir_range_KirCaLs_fI[1]        # [mS/cm²]
    gKM = 0.     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[1]   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_hh = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)



prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
Plots.plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KirCaLs_fI[1,1] +dI_step/2 - dI_step .*pulse(t,4000,10000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step_hh = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_step_hh = solve(prob_step_hh,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_step_hh = I.(sol_step_hh.t)
    plt_V = Plots.plot(sol_step_hh.t,sol_step_hh[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step_hh.t,sol_step_hh[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))

###### Step of size (dI_step),, and CaL,no Kir, no KM
    ## Initial conditions
    V_ic = -55. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I2_KirCaLs_fI[2,1]*pulse(t,0,1000) + (I1_KirCaLs_fI[2,1] + dI_step/2)*pulse(t,1000,80000) # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = gKir_range_KirCaLs_fI[1]        # [mS/cm²]
    gKM = 0.     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[2]   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_CaLs = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)


prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
Plots.plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KirCaLs_fI[2,1] +dI_step/2 - dI_step .*pulse(t,4000,20000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step_CaLs = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_step_CaLs = solve(prob_step_CaLs,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_step_CaLs = I.(sol_step_CaLs.t)
    plt_V = Plots.plot(sol_step_CaLs.t,sol_step_CaLs[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step_CaLs.t,sol_step_CaLs[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))

###### Step of size (dI_step), and Kir
    ## Initial conditions
    V_ic = -55. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I2_KirCaLs_fI[1,3]*pulse(t,0,1000) + (I1_KirCaLs_fI[1,3] + dI_step/2)*pulse(t,1000,80000) # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = gKir_range_KirCaLs_fI[3]        # [mS/cm²]
    gKM = 0.     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[1]   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_Kir = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
Plots.plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KirCaLs_fI[1,3] +dI_step/2 - dI_step .*pulse(t,4000,20000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step_Kir = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_step_Kir = solve(prob_step_Kir,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_step_Kir = I.(sol_step_Kir.t)
    plt_V = Plots.plot(sol_step_Kir.t,sol_step_Kir[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step_Kir.t,sol_step_Kir[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))

###### Step of size (dI_step), and KM
    ## Initial conditions
    V_ic = -55. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I2_KMCaLs_fI[1,3]*pulse(t,0,1000) + (I1_KMCaLs_fI[1,3] + dI_step/2)*pulse(t,1000,80000) # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.      # [mS/cm²]
    gKM = gKM_range_KMCaLs_fI[3]      # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KMCaLs_fI[1]   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_KM = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
Plots.plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KMCaLs_fI[1,3] +dI_step/2 - dI_step .*pulse(t,4000,80000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step_KM = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_step_KM = solve(prob_step_KM,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_step_KM = I.(sol_step_KM.t)
    plt_V = Plots.plot(sol_step_KM.t,sol_step_KM[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step_KM.t,sol_step_KM[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))

###### Step of size (dI_step), Fixedpoints - CaLs and Kir
    ## Initial conditions
    V_ic = -55. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I2_KirCaLs_fI[2,3]*pulse(t,0,1000) + (I1_KirCaLs_fI[2,3] + dI_step/2)*pulse(t,1000,80000) # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = gKir_range_KirCaLs_fI[3]        # [mS/cm²]
    gKM = 0.     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[2]   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_KirCaLs_fI_23 = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
Plots.plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KirCaLs_fI[2,3] +dI_step/2 - dI_step .*pulse(t,4000,80000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step_KirCaLs_fI_23 = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_step_KirCaLs_fI_23 = solve(prob_step_KirCaLs_fI_23,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_step_KirCaLs_fI_23 = I.(sol_step_KirCaLs_fI_23.t)
    plt_V = Plots.plot(sol_step_KirCaLs_fI_23.t,sol_step_KirCaLs_fI_23[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step_KirCaLs_fI_23.t,sol_step_KirCaLs_fI_23[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))


###### Step of size (dI_step), Fixedpoints - CaLs and KM
    ## Initial conditions
    V_ic = -55. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I2_KMCaLs_fI[2,3]*pulse(t,0,1000) + (I1_KMCaLs_fI[2,3] + dI_step/2)*pulse(t,1000,80000) # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.        # [mS/cm²]
    gKM = gKM_range_KMCaLs_fI[3]     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KMCaLs_fI[2]   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_KMCaLs_fI_23 = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
Plots.plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KMCaLs_fI[2,3] +dI_step/2 - dI_step .*pulse(t,4000,80000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step_KMCaLs_fI_23 = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_step_KMCaLs_fI_23 = solve(prob_step_KMCaLs_fI_23,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_step_KMCaLs_fI_23 = I.(sol_step_KMCaLs_fI_23.t)
    plt_V = Plots.plot(sol_step_KMCaLs_fI_23.t,sol_step_KMCaLs_fI_23[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step_KMCaLs_fI_23.t,sol_step_KMCaLs_fI_23[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))


###### Fixedpoints - CaLs max
    ## Initial conditions
    V_ic = -55. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I2_KirCaLs_fI[3,1]*pulse(t,0,1000) + (I1_KirCaLs_fI[3,1] + dI_step/2)*pulse(t,1000,80000) # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = gKir_range_KirCaLs_fI[1]        # [mS/cm²]
    gKM = 0.     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p_var_fp = [0.,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Fixed points 
    fixedpoint_KirCaLs_fI_31 = find_fixed_points([p_var_fp,p_fixed],range(-200.,stop=-20.,step=0.001),-4,0)

