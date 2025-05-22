include("membrane.jl")
include("fixedpoints.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

    ## Fixed parameters
    C = 1 #µF/cm²
    eNa = 50 # [mV]
    eK = -90 # [mV]
    eCa = 120 # [mV]
    eleak = -60.1 # [mV]
    Ca_o = 2 # Extracellular Ca concentration [mM]

    ## Intracellular calcium dynamics parameters
    k = 1.0e7 # [nm.cm^(-1)] 
    F = 96520 # Faraday default value [C/mol]
    d = 1 # [nm]
    Ca_i_0 = 5.0e-5 # Initial intracellular Ca concentration [mM]
    tau_Ca = 10 # [ms]

    p_fixed = [C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]

    ## Initial conditions
    V_ic = -85. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    ## Varying parameters
    I(t) = -0.2516977580264472 + (2.04+1.5) *pulse(t, 5000, 5000+1000)  # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.18          # mS/cm²
    gKM = 0.           # mS/cm²
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = 1.5e-5   # [cm/s]
    
    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

find_fixed_points_ = find_fixed_points([p_var,p_fixed],range(-200.,stop=-20.,step=0.1),-4,0)    

prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
sol = solve(prob,dtmax=1,abstol=1e-7,reltol=1e-7)
sol_cb = solve(prob,dtmax=1,callback=cb,save_everystep=false,save_start=false,save_end=false,abstol=1e-7,reltol=1e-7,Rodas5P())
ic_end = sol[end]

prob2 = ODEProblem(membrane!,ic_end,(0.,240000.),[p_var,p_fixed])
sol2 = solve(prob2,dtmax=1,abstol=1e-7,reltol=1e-7,Rodas5P())
sol_cb2 = solve(prob2,dtmax=1,callback=cb,save_everystep=false,save_start=false,save_end=false,abstol=1e-7,reltol=1e-7,Rodas5P())


plot(sol.t,sol[1,:])
scatter!(sol_cb.t,sol_cb[1,:])
Plots.scatter(sol_cb.t[2:end],1000 ./ (sol_cb.t[2:end]-sol_cb.t[1:end-1]),xlims=(0.,80000.))


plot(sol2.t,sol2[1,:])
scatter!(sol_cb2.t,sol_cb2[1,:])
Plots.scatter(sol_cb2.t[2:end],1000 ./ (sol_cb2.t[2:end]-sol_cb2.t[1:end-1]),xlims=(0.,20000.))


plot(sol.t,sol[11,:])
scatter!(sol_cb.t,sol_cb[11,:])
plot(sol.t,ICaLs.(pCaLs,sol[7,:],sol[8,:],sol[1,:],sol[11,:],Ca_o))
plot(sol.t,INa.(gNa,sol[2,:],sol[3,:],sol[1,:],eNa))
plot(sol.t,IKDR.(gKDR,sol[4,:],sol[1,:],eK))

## Equivalent conductance (from permeability)  
R       = 8.3134
celsius = 36
T=celsius+273.15
gCaLf = pCaLf*(1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T) #mS/cm^2
gCaLs = pCaLs*(1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T) #mS/cm^2