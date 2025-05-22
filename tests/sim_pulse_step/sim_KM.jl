include("../../membrane.jl")
include("../../fixedpoints.jl")

include("../fI_curve/bif_km.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)


###### Pulse of size (dI_pulse), and gKM=0.15
    dI_pulse = 3

    ## Initial conditions
    V_ic = -85. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I1_KMCaLs_fI[3,2] + (I2_KMCaLs_fI[3,2]-I1_KMCaLs_fI[3,2])*2/3  -dI_pulse/2 # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.          # [mS/cm²]
    gKM = gKM_range_KMCaLs_fI[2]      # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]


prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KMCaLs_fI[3,2] + (I2_KMCaLs_fI[3,2]-I1_KMCaLs_fI[3,2])*2/3 -dI_pulse/2 + dI_pulse .*pulse(t,250,350)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_pulse_gKM_fI_2 = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_pulse_gKM_fI_2 = solve(prob_pulse_gKM_fI_2,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_pulse_gKM_fI_2 = I.(sol_pulse_gKM_fI_2.t)
    plt_V = plot(sol_pulse_gKM_fI_2.t,sol_pulse_gKM_fI_2[1,:],ylims=(-155,50))
    plt_Ca = plot(sol_pulse_gKM_fI_2.t,sol_pulse_gKM_fI_2[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(plot(plt_Ca,plt_V,layout=(2,1)))


###### Pulse of size (dI_pulse), and gKM=0.3
    dI_pulse = 3

    ## Initial conditions
    V_ic = -85. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = I1_KMCaLs_fI[3,3] + (I2_KMCaLs_fI[3,3]-I1_KMCaLs_fI[3,3])*2/3  -dI_pulse/2 # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.          # [mS/cm²]
    gKM = gKM_range_KMCaLs_fI[3]      # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]


prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
plot(sol.t,sol[1,:])
ic_end = sol[end]


    I(t) = I1_KMCaLs_fI[3,3] + (I2_KMCaLs_fI[3,3]-I1_KMCaLs_fI[3,3])*2/3 -dI_pulse/2 + dI_pulse .*pulse(t,250,350)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_pulse_gKM_fI_3 = ODEProblem(membrane!,ic_end,(0.,10000.),[p_var,p_fixed])
    sol_pulse_gKM_fI_3 = solve(prob_pulse_gKM_fI_3,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    It_pulse_gKM_fI_3 = I.(sol_pulse_gKM_fI_3.t)
    plt_V = plot(sol_pulse_gKM_fI_3.t,sol_pulse_gKM_fI_3[1,:],ylims=(-155,50))
    plt_Ca = plot(sol_pulse_gKM_fI_3.t,sol_pulse_gKM_fI_3[11,:],ylims=(0,5),title=" I in [$(I(0));$(I(500))]")
    display(plot(plt_Ca,plt_V,layout=(2,1)))


