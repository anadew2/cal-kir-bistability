using Plots,LaTeXStrings,DifferentialEquations
include("gates_compart.jl")
include("currents.jl")
include("solver_fun.jl")

include("compartment.jl")

    ## Fixed parameters
    C = 1 #µF/cm²
    eNa = 50 # [mV]
    eK = -90 # [mV]
    eCa = 120 # [mV]
    eCaAN = 0 # [mV]
    eleak = -60.1 # [mV]
    Ca_o = 2 # Extracellular Ca concentration [mM]

    ## Intracellular calcium dynamics parameters
    k = 1.0e7 # [nm.cm^(-1)] 
    F = 96520 # Faraday default value [C/mol]
    d = 1 # [nm]
    Ca_i_0 = 5.0e-5 # Initial intracellular Ca concentration [mM]
    tau_Ca = 10 # [ms]

    ## Compartment parameters
    Ra = 35.4 # Axial resistivity in [ohm.cm]
    Ld = 850 #[µm]  ; 
    Dd = 1 #[µm]
    Ls = 56.279 #[µm] ; 
    Ds = 56.279 #[µm] ; 

    ## Stimulation 
    I(t) = 0.00014502088607533072 -25/10^6 +( 0.0001464798860753307+200/10^6)*pulse(t,5000,10000) #(25+10-190-0)/10^(6) + 1*(35+190)/10^(6) *pulse(t,5000,10000) + -0 *10/10^(6) *pulse(t,25000,28000)    # [µA]

    p_fixed = [I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra]

    ## Initial conditions
    V_ic = -63 # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic,
        V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic]

    ## Varying parameters
    gNa_s = 30 #mS/cm²
    gKDR_s = 4 #mS/cm²
    gleak_s = 0.03268 #mS/cm²
    gKir_s = 0.02*6.5            # mS/cm² 
    gKM_s = 0.            # mS/cm²
    gKCa_s = 0.03*1         # mS/cm²
    gCaAN_s = 0.            # mS/cm²
    pCaLf_s = 0.15*1 /1000    # [cm/s]
    pCaLs_s = 0. /1000   # [cm/s]

    gNa_d = 0 #mS/cm²
    gKDR_d = 0 #mS/cm²
    gleak_d = 0.03 #mS/cm²
    gKir_d = 0.           # mS/cm²
    gKM_d = 0.            # mS/cm²
    gKCa_d = 0.0185*1  # mS/cm²
    gCaAN_d = 0.125*1 # mS/cm²
    pCaLf_d = 0. /1000    # [cm/s]
    pCaLs_d = 0.02*1. /2 /1000   # [cm/s]
    

    p_var_d = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d,gleak_d]
    p_var_s = [gNa_s,gKDR_s,gKir_s,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]


 

prob = ODEProblem(compartment!,ic,(-20000.,0.),[p_fixed,p_var_d,p_var_s])
sol_eq = solve(prob,dtmax=1,Tsit5(),saveat=-2000:0.1:0) #
ic_eq = sol_eq[end]

    tmax = 100000.
    saveat_range = 0:0.1:tmax

prob = ODEProblem(compartment!,ic_eq,(0.,tmax),[p_fixed,p_var_d,p_var_s])
sol_all =solve(prob,dtmax=1,Tsit5(),callback=cb,save_everystep=false,save_start=false,save_end=false,saveat=saveat_range)	

display(plt_nrn_compare_all(sol_all,saveat_range,p_fixed))

#Plots.savefig("tests/LeFranc-multicomp-jl/figures/LeFranc_compart_jl__nrn_param_hCaLs_notshifted_I_35_pulse.pdf")

