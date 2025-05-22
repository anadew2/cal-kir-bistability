include("../../membrane.jl")
include("../../fixedpoints.jl")

include("../dI_grad/local_postpro_output.jl")
include("../fI_curve/data_loader.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

    ## Initial conditions
    V_ic = -85. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

###### Ascending step
    ## Varying parameters
    I(t) = -3 # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.          # [mS/cm²]
    gKM = 0.           # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]


prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
ic_end = sol[end]

V_asc_step_CaLs = []
t_asc_step_CaLs = []
I_asc_step_CaLs = []
for i_asc in 1:4
    I(t) = -3  .+ (0.5 +(i_asc-1)) .*pulse(t,250,10000)  # [µA/cm^2]
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step = ODEProblem(membrane!,ic_end,(0.,1000.),[p_var,p_fixed])
    sol_step = solve(prob_step,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    plt_V = plot(sol_step.t,sol_step[1,:],ylims=(-155,50))
    plt_Ca = plot(sol_step.t,sol_step[11,:],ylims=(0,5),title="i_asc = $i_asc ; I in [$(I(0));$(I(500))]")
    display(plot(plt_Ca,plt_V,layout=(2,1)))
    push!(V_asc_step_CaLs,sol_step[1,:])
    push!(t_asc_step_CaLs,sol_step.t)
    push!(I_asc_step_CaLs,I.(sol_step.t))
end

###### Descending step
    ## Varying parameters
    I(t) = 1 # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.          # [mS/cm²]
    gKM = 0.           # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

nspan=4
tspan=(0.,20000)
ic_end_cyc = []
push!(ic_end_cyc,ic)
for ispan = 1:nspan
    println("Current ic = $(ic_end_cyc[end])")
    prob = ODEProblem(membrane!,ic_end_cyc[end],tspan .+(ispan-1)*tspan[2],[p_var,p_fixed])
    sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    push!(ic_end_cyc,sol[end])
    plt_V = plot(sol.t,sol[1,:],ylims=(-155,50))
    plt_hCaLs = plot(sol.t,sol[8,:],ylims=(0,1))
    display(plot(plt_hCaLs,plt_V,layout=(2,1)))
end

V_desc_step_CaLs = []
t_desc_step_CaLs = []
I_desc_step_CaLs = []
for i_desc in 1:4
    I(t) = 1 .- (0.5 +(i_desc-1)) .*pulse(t,250,10000)  # [µA/cm^2]
    println(I(1000))
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    prob_step = ODEProblem(membrane!,ic_end_cyc[end],(0.,1000.),[p_var,p_fixed])
    sol_step = solve(prob_step,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
    plt_V = plot(sol_step.t,sol_step[1,:],ylims=(-155,50))
    plt_Ca = plot(sol_step.t,sol_step[11,:],ylims=(0,5),title="i_desc = $i_desc ; I in [$(I(500));$(I(0))]")
    display(plot(plt_Ca,plt_V,layout=(2,1)))
    push!(V_desc_step_CaLs,sol_step[1,:])
    push!(t_desc_step_CaLs,sol_step.t)
    push!(I_desc_step_CaLs,I.(sol_step.t))
end

