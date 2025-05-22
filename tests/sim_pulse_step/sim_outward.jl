include("../dI_outw/membrane_outw.jl")
include("../dI_outw/fixedpoints_outw.jl")

include("../dI_outw/local_postpro_outputs.jl")
include("../dI_outw/data_loader_outw.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

    ## Initial conditions
    V_ic = -65. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

###### KIR 
###### Descending step
    ## Varying parameters
    I(t) = I2_Kir_outw[16] # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = gKir_KirCaLs_outw[16]  # [mS/cm²]
    gKM = 0.           # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_KirCaLs_outw   # [cm/s]
    
    p_var_Kir_outw = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    p_var_Kir_outw__ = [I2_Kir_outw[16],gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    find_fixed_points([p_var_Kir_outw__,p_fixed],range(-200.,stop=-20.,step=0.1),-4,0)

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var_Kir_outw,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
ic_end = sol[end]

V_desc_step_Kir_outw = []
t_desc_step_Kir_outw = []
It_desc_step_Kir_outw = []

n_desc = 10
for i_desc in 1:n_desc
    I(t) = I2_Kir_outw[16] .- (I2_Kir_outw[16]-I1_Kir_outw[16])/n_desc *(i_desc) .*pulse(t,250,10000)  # [µA/cm^2]
    println(I(1000))
    p_var_Kir_outw_ = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step = ODEProblem(membrane!,ic_end,(0.,1000.),[p_var_Kir_outw_,p_fixed])
    sol_step = solve(prob_step,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    plt_V = Plots.plot(sol_step.t,sol_step[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step.t,sol_step[11,:],ylims=(0,5),title="i_desc = $i_desc ; I in [$(I(500));$(I(0))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))
    push!(V_desc_step_Kir_outw,sol_step[1,:])
    push!(t_desc_step_Kir_outw,sol_step.t)
    push!(It_desc_step_Kir_outw,I.(sol_step.t))
end

###### KM 
###### Descending step
    ## Varying parameters
    I(t) = I2_KM_outw[16] # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.  # [mS/cm²]
    gKM = gKM_KMCaLs_outw[16]           # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_KMCaLs_outw   # [cm/s]
    
    p_var_KM_outw = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
		
        ## Initial conditions
        p_var_KM_outw__ = [I2_KM_outw[16],gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
        find_fixed_points([p_var_KM_outw__,p_fixed],range(-200.,stop=-20.,step=0.1),-4,0)
        FP_I2_KM_outw = nlsolve((F,x) ->nlsolve_membrane_grad_0_V_!(F,x,[p_var_KM_outw__,p_fixed]), [-65,5.0e-5],ftol=1E-7,xtol=1E-7)

        V_ic = FP_I2_KM_outw.zero[1] # [mV]
        Ca_ic = FP_I2_KM_outw.zero[2] # [mM]
        ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]
    

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var_KM_outw,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
#plot(sol.t,sol[1,:],ylims=(-155,50))
ic_end_KM_outw = sol[end]

V_desc_step_KM_outw = []
t_desc_step_KM_outw = []
It_desc_step_KM_outw = []

n_desc = 10
for i_desc in 1:n_desc
    I(t) = I2_KM_outw[16] .- (I2_KM_outw[16]-I1_KM_outw[16])/n_desc *(i_desc) .*pulse(t,250,10000)  # [µA/cm^2]
    println(I(1000))
    p_var_KM_outw_ = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step = ODEProblem(membrane!,ic_end_KM_outw,(0.,1000.),[p_var_KM_outw_,p_fixed])
    sol_step = solve(prob_step,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    plt_V = Plots.plot(sol_step.t,sol_step[1,:],ylims=(-155,50))
    plt_Ca = Plots.plot(sol_step.t,sol_step[11,:],ylims=(0,5),title="i_desc = $i_desc ; I in [$(I(500));$(I(0))]")
    display(Plots.plot(plt_Ca,plt_V,layout=(2,1)))
    push!(V_desc_step_KM_outw,sol_step[1,:])
    push!(t_desc_step_KM_outw,sol_step.t)
    push!(It_desc_step_KM_outw,I.(sol_step.t))
end

