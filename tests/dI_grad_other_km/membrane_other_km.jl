using Plots,LaTeXStrings,DifferentialEquations
include("gates_other_km.jl")
include("currents_other_km.jl")
include("../../solver_fun.jl")

function membrane!(du,u,p,t)
    ## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]
    
    ## Gates steady-state and time-constants 
    mNa_inf=mNa(u[1])
	tau_mNa_=tau_mNa(u[1])
    hNa_inf=hNa(u[1])
	tau_hNa_=tau_hNa(u[1])
    
    mKDR_inf=mKDR(u[1])
	tau_mKDR_=tau_mKDR(u[1])

	mCaLf_inf=mCaLf(u[1])
	tau_mCaLf_=tau_mCaLf(u[1])
    hCaLf_inf=hCaLf(u[1])
	tau_hCaLf_=tau_hCaLf(u[1])

    mCaLs_inf=mCaLs(u[1])
	tau_mCaLs_=tau_mCaLs(u[1])
    hCaLs_inf=hCaLs(u[1])
	tau_hCaLs_=tau_hCaLs(u[1])

    mKir_inf=mKir(u[1])
	tau_mKir_=tau_mKir(u[1])

    mKM_inf=mKM(u[1])
    tau_mKM_= tau_mKM(u[1])

    ## Current computation
    INa_ = INa(gNa,u[2],u[3],u[1],eNa)
    IKDR_ = IKDR(gKDR,u[4],u[1],eK)
    ICaLf_ = ICaLf(pCaLf,u[5],u[6],u[1],u[11],Ca_o)
    ICaLs_ = ICaLs(pCaLs,u[7],u[8],u[1],u[11],Ca_o)
    IKir_ = IKir(gKir,u[9],u[1],eK)
    IKM_ = IKM(gKM,u[10],u[1],eK)
    Ileak_ = Ileak(gleak,u[1],eleak)

    ## Gradients
	du[1] = (-INa_ -IKDR_ -ICaLf_ -ICaLs_- IKir_- IKM_ -Ileak_ +I(t))/C

	du[2] = (mNa_inf-u[2])/tau_mNa_
	du[3] = (hNa_inf-u[3])/tau_hNa_
	du[4] = (mKDR_inf-u[4])/tau_mKDR_
    du[5] = (mCaLf_inf-u[5])/tau_mCaLf_
	du[6] = (hCaLf_inf-u[6])/tau_hCaLf_
	du[7] = (mCaLs_inf-u[7])/tau_mCaLs_
	du[8] = (hCaLs_inf-u[8])/tau_hCaLs_
	du[9] = (mKir_inf-u[9])/tau_mKir_
    du[10] = (mKM_inf-u[10])/tau_mKM_

    ICa_ = (ICaLf_ + ICaLs_)/1000  # [mA/cm^2] instead of [µA/cm^2]
    du[11] = - ICa_ * k /(2*F*d) - ((u[11]-Ca_i_0)/tau_Ca)
    #Ca_o = 2 mM (McCornick)
end

    ## Fixed parameters
    C = 1 #µF/cm²
    eNa = 50 # [mV]
    eK = -90 # [mV]
    eCa = 120 # [mV]
    eleak = -60.1 # [mV]
    Ca_o = 2 # Extracellular Ca concentration [mM]

    ## Intracellular calcium dynamics parameters
    k = 1.0e7 # [nm.cm^(-1)] #TYPO
    F = 96520 # Faraday default value [C/mol]
    d = 1 # [nm]
    Ca_i_0 = 5.0e-5 # Initial intracellular Ca concentration [mM]
    tau_Ca = 10 # [ms]

    p_fixed = [C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]

    ## Initial conditions
    V_ic = -50. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = -0.12     # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.           # [µA/cm^2]
    gKM = 0.            # [µA/cm^2]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = 0.05/1000   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

#==
prob = ODEProblem(membrane!,ic,(0.,20000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7)
sol_cb = solve(prob,dtmax=1,callback=cb,save_everystep=false,save_start=false,save_end=false,abstol=1e-7,reltol=1e-7)	

Plots.plot(sol.t,sol[1,:])
scatter!(sol_cb.t,sol_cb[1,:])
Plots.plot(sol.t,sol[11,:])
scatter!(sol_cb.t,sol_cb[11,:])
Plots.plot(sol.t,ICaLs.(pCaLs,sol[7,:],sol[8,:],sol[1,:],sol[11,:],Ca_o))
Plots.plot(sol.t,INa.(gNa,sol[2,:],sol[3,:],sol[1,:],eNa))
Plots.plot(sol.t,IKDR.(gKDR,sol[4,:],sol[1,:],eK))

## Equivalent conductance (from permeability)  
R       = 8.3134
celsius = 36
T=celsius+273.15
gCaLf = pCaLf*(1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T) #mS/cm^2
gCaLs = pCaLs*(1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T) #mS/cm^2
==#