using Plots,LaTeXStrings,DifferentialEquations,JLD

include("../solver_fun.jl")

include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("jld_hdle_fI.jl")
include("data_loader.jl")

    ind_pCaLs_simu = [3,3,16,16]
    ind_gKir_simu = [4,38,4,38]
    I1pulse_simu = [150,50,150,50]
    I0pulse_simu = [-0.7,210.7,-1.2,210.2]

    sol_simu = []
    t_sol_simu = []
    It_sol_simu = []
    f_sol_simu = []
    plt_simu = []

    ind_sol_simu=[1,13,14,26]
    
    ## Initial conditions
    V_ic = -63 # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic,
        V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic]


    for i in eachindex(ind_pCaLs_simu)
        i_pCaLs = ind_pCaLs_simu[i]
        i_gKir = ind_gKir_simu[i]

        ## Stimulation 
        println("------------- i_pCaLs = $(i_pCaLs) & i_gKir = $(i_gKir)-------------")
        #println("I1 = $(I1_Kir_s_CaLs_d_J[i_pCaLs,i_gKir]*10^6) pA")
        #println("I2 = $(I2_Kir_s_CaLs_d_J[i_pCaLs,i_gKir]*10^6) pA")
        
        I(t) = I0pulse_simu[i]/10^6 + (I1pulse_simu[i]/10^6)*pulse(t,500,2000)  # [µA]
        p_fixed = [I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra];

        ## Varying parameters
        pCaLs_d_fI_test = pCaLs_d_range_KirCaLs[i_pCaLs]
        gKir_s_fI_test = gKir_s_range_KirCaLs[i_gKir]
        p_var_d_fI = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_fI_test,gleak_d]
        p_var_s_fI = [gNa_s,gKDR_s,gKir_s_fI_test,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]

        prob = ODEProblem(compartment!,ic,(-20000.,0.),[p_fixed,p_var_d_fI,p_var_s_fI])
        sol_eq = solve(prob,dtmax=1,Tsit5(),saveat=-2000:0.1:0) #
        ic_eq = sol_eq[end]

            tmax = 5000.
            saveat_range = 0:0.1:tmax

        prob = ODEProblem(compartment!,ic_eq,(0.,tmax),[p_fixed,p_var_d_fI,p_var_s_fI])
        sol_all =solve(prob,dtmax=1,Tsit5(),callback=cb,save_everystep=false,save_start=false,save_end=false,saveat=saveat_range)	
        tspk_sol = setdiff(sol_all.t,saveat_range)
        plt_simu_ = plt_nrn_compare_all(sol_all,saveat_range,p_fixed)
        display(plt_simu_)
        #Plots.savefig("tests/LeFranc-multicomp-jl/figures/LeFranc_compart_jl__nrn_param_hCaLs_notshifted_I_35_pulse.pdf")
        
        push!(sol_simu,sol_all[ind_sol_simu,:])
        push!(t_sol_simu,sol_all.t)
        push!(It_sol_simu,I.(sol_all.t) .*10^6) #[pA]
        if length(tspk_sol)>1
            push!(f_sol_simu,1000 ./(tspk_sol[2:end]-tspk_sol[1:end-1]))
        else
            push!(f_sol_simu,[])
        end
        push!(plt_simu,plt_simu_)
    end
#

