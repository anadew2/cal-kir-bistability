
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

include("data_loader.jl")

function save_full_bif(p_var, find_full_bif, jld_name)
    save(jld_name, "p_var", p_var, "f_list", find_full_bif[1], "V_max_c", find_full_bif[2], "V_min_c", find_full_bif[3])
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
function load_full_bif(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        p_var = []
        f_list = []
        V_max_c = []
        V_min_c = []
        push!(p_var,load_array(jld,"p_var"))
        push!(f_list,load_array(jld,"f_list"))
        push!(V_max_c,load_array(jld,"V_max_c"))
        push!(V_min_c,load_array(jld,"V_min_c"))
    else
        p_var = NaN
        f_list = NaN
        V_max_c = NaN
        V_min_c = NaN
        println("Not found")
    end
    find_full_bif = []
    push!(find_full_bif,f_list[1])
    push!(find_full_bif,V_max_c[1])
    push!(find_full_bif,V_min_c[1])
    return p_var,find_full_bif
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

    t_pulse = 1000;
    dt_pulse = 1000;

    f_tol = 0.1
    lim_cnt=5

    V_max_c = []
    V_min_c = []

    while_list = []

    for i in eachindex(I_list)

        println("$i on $(length(I_list))")
        println("I[i] = $(I_list[end-i+1])")
        push!(rev_I_list,I_list[end-i+1])

        tspan = (0.0,tmin)

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
        sol_pulse = solve(prob_pulse,abstol=1e-7,reltol=1e-7,dtmax=1,Rodas5P())
        #display(plot(sol_pulse.t,sol_pulse[1,:]))
        f_state=sol_pulse[end]

        while freq_test_tol > f_tol && cnt <lim_cnt
            cnt = cnt +1
            println("f_state = $(f_state)")

            cb = ContinuousCallback(condition,affect!,nothing,save_positions=(true,false))
            prob = ODEProblem(ode_fun,f_state,tspan,[p_var_i,p_fixed])
            sol_full = solve(prob,abstol=1e-7,reltol=1e-7,Rodas5P())

            spikes_from_sol = []
            spik_times = sol_full.t[abs.(sol_full[1,:]) .<=5]
            if length(spik_times)>1
                push!(spikes_from_sol,spik_times[1])
                for i in 2:length(spik_times)
                    if spik_times[i]- spik_times[i-1]> 2 
                        push!(spikes_from_sol,spik_times[i])
                    end
                end
            end

            if length(spikes_from_sol) >3 
                t_2 = spikes_from_sol[3:end-1]
                t_1 = spikes_from_sol[2:end-2]
                freq_isi = 1 ./ (t_2-t_1) .*1000
                freq_test_tol = (maximum(freq_isi) - minimum(freq_isi))#/minimum(freq_isi)

                V_max = maximum(sol_full[1,:])
                V_min = minimum(sol_full[1,:])
                f_state[:]=sol_full[end-1][:]
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
            println("Max t : $(maximum(sol_full.t))")
            println("Average frequency : $(mean_f)")
            push!(f_I, mean_f)

        end
        push!(while_list,test_tol_list)
        push!(tspan_list,tspan[2]-tmin)
        push!(f_list,f_I[end])

        if f_I[end]>0
            println("Cycle found")
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
function find_full_bif(p,ic,I_list)

    I1 = maximum(I_list)+1.5

	bif = bif_desc_w_pulse(membrane!,I_list,p,ic,I1,20000)
	println("f_list = $(bif[2])")

    f_list = bif[2]
    V_max_c = bif[5]
    V_min_c = bif[6]

	return f_list,V_max_c,V_min_c
end


### SCRIPT INPUTS ###
    i_pCaLs_range_KMCaLs= parse(Int64, ARGS[1])     #1:length(pCaLs_range_KMCaLs_fI)
    i_gKM_range_KMCaLs = parse(Int64, ARGS[2])      #1:length(gKM_range_KMCaLs_fI)
    i_batch = parse(Int64, ARGS[3])                   #1:n_batch

### INPUTS ###
    pCaLs = pCaLs_range_KMCaLs_fI[i_pCaLs_range_KMCaLs]
    gKM = gKM_range_KMCaLs_fI[i_gKM_range_KMCaLs]

    n_batch = 1

    stp_I_list_batch_1=0.01
    stp_I_list_batch_2=0.1
    I_list_batch_1 = collect((I1_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]-0.5):stp_I_list_batch_1:(I1_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]+0.5))
    I_list_batch_2 = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I2_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]+0.5))
    I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2)

    ind_batch = Int(ceil(length(I_list_batch_)/n_batch))
    I_list_batch = I_list_batch_[( (Int(1+ ind_batch*(i_batch-1))) : Int(minimum([ind_batch + ind_batch*(i_batch-1),length(I_list_batch_)]))) ]

    V_ic = -50.     # [mV]
    Ca_ic = 5.0e-5  # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

    jld_name = "data/jld/fI_KMCaLs__"

### TASK ###
    jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KMCaLs)_i_gKM_$(i_gKM_range_KMCaLs)_i_batch_$(i_batch).jld")

    p_var_i_batch = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    find_full_bif_ = find_full_bif([p_var_i_batch,p_fixed],ic,I_list_batch)
    save_full_bif(p_var_i_batch[2:end], find_full_bif_, jld_name_i)
   