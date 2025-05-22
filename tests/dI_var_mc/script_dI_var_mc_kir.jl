using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

include("data_loader.jl")

function save_dI(p_var,find_dI, jld_name)
    save(jld_name, "p_var", p_var, "dI", find_dI[1], "I1", find_dI[2], "I_SN", find_dI[3], "V_SN", find_dI[4],  "I2", find_dI[5])
end
function load_jld(jld_name)
	data = []
	try 
		push!(data,load(jld_name))
	catch err
		push!(data ,NaN)
		println("Not found...")
	end
	return data[1]
end
function load_array(data,array_name)
	return data[array_name]
end
function load_dI(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        p_var = []
        dI = []
        I1 = []
        I_SN = []
        V_SN = []
        I2 = []
        push!(p_var,load_array(jld,"p_var"))
        push!(dI,load_array(jld,"dI"))
        push!(I1,load_array(jld,"I1"))
        push!(I_SN,load_array(jld,"I_SN"))
        push!(V_SN,load_array(jld,"V_SN"))
        push!(I2,load_array(jld,"I2"))
    else
        p_var = NaN
        dI = NaN
        I1 = NaN
        I_SN = NaN
        V_SN = NaN
        I2 = NaN
        println("Not found")
    end
    find_dI = []
    push!(find_dI,dI[1])
    push!(find_dI,I1[1])
    push!(find_dI,I_SN[1])
    push!(find_dI,V_SN[1])
    push!(find_dI,I2[1])
    return p_var,find_dI
end


pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)
heaviside(t)=0*(t<0)+1*(t>=0)



function bif_desc_w_pulse(ode_fun,I_list,p,ic_h,I1,tmin)
    p_var,p_fixed = p


	## Varying parameters
	gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p_var[2:end]


    f_list = []
    rev_I_list = []
    tspan_list = []

    t_pulse = 5000;
    dt_pulse = 1000;
    f_tol = 0.1
    lim_cnt=5

    V_max_c = []
    V_min_c = []

    while_list = []

    for i=1:length(I_list)
        println("$i on $(length(I_list))")
        tspan = (0.0,tmin)
        push!(rev_I_list,I_list[end-i+1])

        I0 = I_list[end-i+1]
        It_fixed(t) = I0 
        p_var_i = [It_fixed,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
        It_pulse(t) = I0 + (I1-I0)*pulse(t, t_pulse, t_pulse+dt_pulse)
        p_var_it = [It_pulse,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

        #Cycle convergence
        freq_test_tol= f_tol+1
        f_I = []
        cnt = 0
        test_tol_list = []

        V_max = 0
        V_min = 0 

        prob_pulse = ODEProblem(ode_fun,ic_h,tspan./2,[p_var_it,p_fixed])
        sol_pulse = solve(prob_pulse,abstol=1e-7,reltol=1e-7,dtmax=1)
        #display(plot(sol_pulse.t,sol_pulse[1,:]))
        f_state=sol_pulse[end]

        while freq_test_tol > f_tol && cnt <lim_cnt
            cnt = cnt +1
            println(f_state)

            cb = ContinuousCallback(condition,affect!,nothing,save_positions=(true,false))
            prob = ODEProblem(ode_fun,f_state,tspan,[p_var,p_fixed])
            sol = solve(prob,callback=cb,save_everystep=false,save_start=false,save_end=false,abstol=1e-7,reltol=1e-7)
            sol_full = solve(prob,abstol=1e-7,reltol=1e-7)

            if length(sol.t) >3 
                t_2 = sol.t[3:end-1]
                t_1 = sol.t[2:end-2]
                freq_isi = 1 ./ (t_2-t_1) .*1000
                freq_test_tol = (maximum(freq_isi) - minimum(freq_isi))#/minimum(freq_isi)

                V_max = maximum(sol_full[1,:])
                V_min = minimum(sol_full[1,:])
                f_state[:]=sol[end-1][:]
            else
                freq_test_tol = 0.
                freq_isi = []
            end

            tspan = tspan .+ tmin

            push!(test_tol_list,sol)
            if length(freq_isi)>1
                mean_f = sum(freq_isi)/length(freq_isi)
            else
                mean_f = 0
            end
            println(mean_f)
            push!(f_I, mean_f)

        end
        push!(while_list,test_tol_list)
        push!(tspan_list,tspan[2]-tmin)
        push!(f_list,f_I[end])

        if f_I[end]>0
            push!(V_max_c,V_max)
            push!(V_min_c,V_min)
        else
            push!(V_max_c,NaN)
            push!(V_min_c,NaN)
        end
    end

    plt = []

    return plt,f_list,tspan_list,while_list,V_max_c,V_min_c
end
function bif_asc_w_pulse_if_bist(ode_fun,I_list,p,ic_h,I1,tmin)
    p_var,p_fixed = p

	## Varying parameters
	gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p_var[2:end]
    
    f_list = []
    rev_I_list = []
    I_tested_list = []

    t_pulse = 5000;
    dt_pulse = 1000;

    f_tol = 0.1
    lim_cnt=5

    for i in eachindex(I_list)

        println("$i on $(length(I_list))")
        println("I[i] = $(I_list[i])")
        push!(rev_I_list,I_list[i])

        tspan = (0.0,tmin)

        I0 = I_list[i]
        It_fixed(t) = I0 
        p_var_i = [It_fixed,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
        It_pulse(t) = I0 + (I1-I0)*pulse(t, t_pulse, t_pulse+dt_pulse)
        p_var_it = [It_pulse,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
        
        push!(I_tested_list,I0)

        stability_FP_I_i = stability_FP([p_var_i,p_fixed],I_list[i],I_list[i],1)

        if length(stability_FP_I_i[6]) >0
            #Cycle convergence
            freq_test_tol= f_tol+1
            f_I = []
            cnt = 0
            test_tol_list = []

            V_max = 0
            V_min = 0 

            prob_pulse = ODEProblem(ode_fun,ic_h,tspan./2,[p_var_it,p_fixed])
            sol_pulse = solve(prob_pulse,abstol=1e-7,reltol=1e-7,dtmax=1)
            #display(plot(sol_pulse.t,sol_pulse[1,:]))
            f_state=sol_pulse[end]

            while freq_test_tol > f_tol && cnt <lim_cnt
                cnt = cnt +1
                println("f_state = $(f_state)")

                cb = ContinuousCallback(condition,affect!,nothing,save_positions=(true,false))
                prob = ODEProblem(ode_fun,f_state,tspan,[p_var_i,p_fixed])
                sol = solve(prob,callback=cb,save_everystep=false,save_start=false,save_end=false,abstol=1e-7,reltol=1e-7)
                sol_full = solve(prob,abstol=1e-7,reltol=1e-7)
                #plt__ = plot(sol_full.t,sol_full[1,:],label="I=$(I_list[end-i+1]) , cnt = $(cnt)")
                #scatter!(sol.t,sol[1,:])
                #display(plt__)
                

                if length(sol.t) >3 
                    t_2 = sol.t[3:end-1]
                    t_1 = sol.t[2:end-2]
                    freq_isi = 1 ./ (t_2-t_1) .*1000
                    freq_test_tol = (maximum(freq_isi) - minimum(freq_isi))#/minimum(freq_isi)

                    V_max = maximum(sol_full[1,:])
                    V_min = minimum(sol_full[1,:])
                    f_state[:]=sol[end-1][:]
                else
                    freq_test_tol = 0.
                    freq_isi = []
                end

                tspan = tspan .+ tmin

                push!(test_tol_list,sol)
                if length(freq_isi)>1
                    mean_f = sum(freq_isi)/length(freq_isi)
                else
                    mean_f = 0
                end
                println("Average frequency : $(mean_f)")
                push!(f_I, mean_f)

            end
            push!(f_list,f_I[end])

            if f_I[end]>0
                println("Cycle found")
                break
            else
                if i==length(I_list)
                    println("No cycle found on the I_list given")
                end
            end
        else
            if length(stability_FP_I_i[3]) >0
                println("is a saddle")
            end
            if length(stability_FP_I_i[9]) >0
                println("is unstable")
            end
        end
    end

    plt = []

    return plt,f_list,I_tested_list
end
function find_bif(p,ic,I_list_0f,I_list_lf,I_list_hf,I1)
	bif_0f = bif_desc_w_pulse(membrane!,I_list_0f,p,ic,I1,50000)
	bif_lf = bif_asc_w_pulse_if_bist(membrane!,I_list_lf,p,ic,I1,120000)
	println(bif_lf[2])
	bif_hf = bif_desc_w_pulse(membrane!,I_list_hf,p,ic,I1,20000)

	I1_bist = bif_lf[3][end] #cycle appears

	return  I1_bist,bif_0f,bif_lf,bif_hf 
end
function find_dI(p,ic,pCaLf_,pCaLs_,gKDR_,gKir_,gKM_,I_X)
    ## Varying parameters
    I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    p_fixed = p[2]

    ## Set of parameters to test 
    p_to_test = []
    for _pCaLf_ in pCaLf_
        p_ = [I,gNa,gKDR,gKir,gKM,_pCaLf_,pCaLs,gleak]
        push!(p_to_test,p_)
    end
    for _pCaLs_ in pCaLs_
        p_ = [I,gNa,gKDR,gKir,gKM,pCaLf,_pCaLs_,gleak]
        push!(p_to_test,p_)
    end
    for _g_ in gKDR_
        p_ = [I,gNa,_g_,gKir,gKM,pCaLf,pCaLs,gleak]
        push!(p_to_test,p_)
    end
    for _g_ in gKir_
        p_ = [I,gNa,gKDR,_g_,gKM,pCaLf,pCaLs,gleak]
        push!(p_to_test,p_)
    end
    for _g_ in gKM_
        p_ = [I,gNa,gKDR,gKir,_g_,pCaLf,pCaLs,gleak]
        push!(p_to_test,p_)
    end
    println(p_to_test)


    stp_V_FP_range=0.01
    V_FP_range = range(-150.,stop=-40.,step=stp_V_FP_range)

    I_SN_ = []
    V_SN_ = []

    I1_ = []
    I2_ = []
    dI_ = []

    println("Nombre de sets de paramètres : $(length(p_to_test))")
    for i in eachindex(p_to_test)
        println("-----$(i) out of $(length(p_to_test))-----")
        I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable= find_fixed_points([p_to_test[i],p_fixed],V_FP_range,I_X[i],0)
        println("I_X = $(I_X[i])  &   I_SN = $(I_SN)")

        I_l_lf = I_X[i] -0.1
        I_h_lf = maximum(I_stable[I_stable.<1000])
    
        if I_h_lf-I_l_lf > 2*10^(0)
            big_Istep = [10^(0),10^(-1),10^(-2),10^(-3)]
        else
            if I_h_lf-I_l_lf > 2*10^(-1)
                big_Istep = [10^(-1),10^(-2),10^(-3)]
            else
                if I_h_lf-I_l_lf > 2*10^(-2)
                    big_Istep = [10^(-2),10^(-3)]
                else
                    big_Istep = [10^(-3)]
                end
            end
        end

        println("big_Istep = $(big_Istep)")
        println(I_h_lf-I_l_lf )
        
        I_list_lf_big = reverse(collect(range(I_h_lf,stop=I_l_lf-big_Istep[1], step=-big_Istep[1])))
        println("$I_l_lf && $(I_h_lf)")

        
        popfirst!(big_Istep)
        I1_est,bif_0f,bif_lf,bif_hf = find_bif([p_to_test[i],p_fixed],ic,[],I_list_lf_big,[],ceil(I_SN)+1.5)
        
        if bif_lf[2][end]>0 
            if length(big_Istep)>0

                I_list_lf_big = vcat(bif_lf[3][end-1],reverse(collect(range(I1_est,stop=bif_lf[3][end-1], step=-big_Istep[1]))))
                popfirst!(big_Istep)                
                I1_est,bif_0f,bif_lf,bif_hf = find_bif([p_to_test[i],p_fixed],ic,[],I_list_lf_big,[],ceil(I_SN)+1.5)
                
                if length(big_Istep)>0
                    
                    I_list_lf_big = vcat(bif_lf[3][end-1],reverse(collect(range(I1_est,stop=bif_lf[3][end-1], step=-big_Istep[1]))))
                    popfirst!(big_Istep)
                    I1_est,bif_0f,bif_lf,bif_hf = find_bif([p_to_test[i],p_fixed],ic,[],I_list_lf_big,[],ceil(I_SN)+1.5)
                    
                    if length(big_Istep)>0
                    
                        I_list_lf_big = vcat(bif_lf[3][end-1],reverse(collect(range(I1_est,stop=bif_lf[3][end-1], step=-big_Istep[1]))))
                        popfirst!(big_Istep)
                        I1_est,bif_0f,bif_lf,bif_hf = find_bif([p_to_test[i],p_fixed],ic,[],I_list_lf_big,[],ceil(I_SN)+1.5)
                   
                    end
                end
            end

            I2 = maximum(I_stable[I_stable.<1000]) 

            push!(I_SN_,I_SN)
            push!(V_SN_,V_SN)

            push!(I1_,I1_est)
            push!(I2_,I2)
            push!(dI_,maximum([0,I2-I1_est]))

            println("I1 = $(I1_est)   &    I2 = $(I2)   &   dI = $(maximum([0,I2-I1_est]))")   
        else
            #par hypothèse sur le diagramme de bifurcation, si pas de cycle au point stable le plus à droite alors I1=I_SN 
            I2 = maximum(I_stable[I_stable.<1000])

            push!(I_SN_,I_SN)
            push!(V_SN_,V_SN)

            push!(I1_,I_SN)
            push!(I2_,I2)
            push!(dI_,maximum([0,I2-I_SN]))

            println("I1 = $(I_SN)   &    I2 = $(I2)   &   dI = $(maximum([0,I2-I_SN]))")    
        end

    end
    return dI_,I1_,I_SN_,V_SN_,I2_
end

### SCRIPT INPUTS ###
i_gNa_var_mat = parse(Int64, ARGS[1]) #1:length(gNa_var_mat)
i_C_var_mat = i_gNa_var_mat 
i_gleak_var_mat = i_gNa_var_mat 

### INPUTS ###
    n_batch_mc_var_mat = 5

    V_ic = -90.     # [mV]
    Ca_ic = 5.0e-5  # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKM(V_ic),mKM(V_ic),Ca_ic]

    jld_name = "data/jld/dI_var$(Int(var_level*100))_mc_KirCaLs"

### TASK ###
    batch_mc_var_mat = i_gNa_var_mat:(minimum([i_gNa_var_mat+n_batch_mc_var_mat-1,length(gNa_var_mat)]))
    for i_batch in batch_mc_var_mat
        jld_name_i = string(jld_name,"_i_C_gleak_gNa_$(i_batch).jld")
        p_fixed_i_batch = [C_var_mat[i_batch],eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]
        p_var_i_batch = [I,gNa_var_mat[i_batch],gKDR,gKir_KirCaLs,gKM_KirCaLs,pCaLf,pCaLs_KirCaLs,gleak_var_mat[i_batch]]
        find_dI_ = find_dI([p_var_i_batch,p_fixed_i_batch],ic,[],[],[],[gKir_KirCaLs],[],[-4])
        save_dI(p_var_i_batch[2:end], find_dI_, jld_name_i)
    end