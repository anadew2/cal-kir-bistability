

function bif_desc_w_pulse(ode_fun,I_list,p,ic_h,I1,tmin)
    ## Fixed parameters
    p_fixed = p[1]

    ## Varying parameters -- dendrite
	p_dend= p[2]

    ## Varying parameters -- soma
	p_soma = p[3]

    f_list = []
    rev_I_list = []
    f_state_list = []

    t_pulse = 1000;
    dt_pulse = 5000;

    f_tol = 0.01
    lim_cnt=10

    Vd_max_c = []
    Vd_min_c = []
    Cad_max_c = []
    Cad_min_c = []
    Vs_max_c = []
    Vs_min_c = []
    Cas_max_c = []
    Cas_min_c = []

    for i in eachindex(I_list)

        println("$i on $(length(I_list))")
        println("I[i] = $(I_list[end-i+1])")

        tspan = (0.0,tmin)
        dt_saveat_range = 0.1
        saveat_range = 0:dt_saveat_range:tmin

        I0 = I_list[end-i+1]
        push!(rev_I_list,I0)

        It_fixed(t) = I0 
        p_fixed_it_cst = vcat(It_fixed,p_fixed[2:end])
        It_pulse(t) = I0 + (I1-I0)*pulse(t, t_pulse, t_pulse+dt_pulse)
        p_fixed_it_pulse = vcat(It_pulse,p_fixed[2:end])

        #Cycle convergence
        freq_test_tol= f_tol+1
        f_I = []
        cnt = 0

        Vd_max = 0
        Vd_min = 0 
        Cad_max = 0
        Cad_min = 0 
        Vs_max = 0
        Vs_min = 0 
        Cas_max = 0
        Cas_min = 0 

        if length(f_list)>0 && f_list[end]<=10.5
            println("Initialized at the f_state of the previous current (last f = $(f_list[end]))")
            f_state = f_state_list[end]
        else
            println("Initialized from a pulse")
            prob_pulse = ODEProblem(ode_fun,ic_h,tspan./2,[p_fixed_it_pulse,p_dend,p_soma])
            sol_pulse = solve(prob_pulse,dtmax=1,Tsit5(),save_everystep=false,save_start=false)
            f_state=sol_pulse[end]
        end

        while freq_test_tol > f_tol && cnt <lim_cnt
            cnt = cnt +1

            prob = ODEProblem(ode_fun,f_state,tspan,[p_fixed_it_cst,p_dend,p_soma])
            sol_all =solve(prob,dtmax=1,Tsit5(),callback=cb,save_everystep=false,save_start=false,save_end=false,saveat=saveat_range)	
            tspk_sol = setdiff(sol_all.t,saveat_range)
            
            if length(sol_all)>0
                println("size(sol_all) = $(size(sol_all))")
                println("Max t : $(maximum(sol_all.t))")
            else
                println("length(sol_all) = $(length(sol_all))")
            end

            if length(tspk_sol) >3 
                freq_isi = 1000 ./(tspk_sol[2:end]-tspk_sol[1:end-1])
                freq_test_tol = (maximum(freq_isi) - minimum(freq_isi))#/minimum(freq_isi)

                T_min = 1000 ./minimum(freq_isi)
                l_T_min = Int(ceil(1.5*T_min/dt_saveat_range))

                if length(sol_all[1,:])>=l_T_min
                    Vd_max = maximum(sol_all[1,end-l_T_min+1:end])
                    Vd_min = minimum(sol_all[1,end-l_T_min+1:end])
                    Cad_max = maximum(sol_all[13,end-l_T_min+1:end])
                    Cad_min = minimum(sol_all[13,end-l_T_min+1:end])
                    Vs_max = maximum(sol_all[14,end-l_T_min+1:end])
                    Vs_min = minimum(sol_all[14,end-l_T_min+1:end])
                    Cas_max = maximum(sol_all[26,end-l_T_min+1:end])
                    Cas_min = minimum(sol_all[26,end-l_T_min+1:end])
                else
                    Vd_max = maximum(sol_all[1,:])
                    Vd_min = minimum(sol_all[1,:])
                    Cad_max = maximum(sol_all[13,:])
                    Cad_min = minimum(sol_all[13,:])
                    Vs_max = maximum(sol_all[14,:])
                    Vs_min = minimum(sol_all[14,:])
                    Cas_max = maximum(sol_all[26,:])
                    Cas_min = minimum(sol_all[26,:])
                end
                if cnt==1
                    tspan=tspan./2
                end
                
                f_state[:]=sol_all[end][:]
            else
                freq_test_tol = 0.
                freq_isi = []
            end


            if length(freq_isi)>1
                mean_f = sum(freq_isi)/length(freq_isi)
            else
                mean_f = 0
            end
            println("Average frequency : $(mean_f)")
            push!(f_I, mean_f)

        end
        push!(f_state_list,f_state)
        push!(f_list,f_I[end])

        if f_I[end]>0
            println("Cycle found")
            push!(Vd_max_c,Vd_max)
            push!(Vd_min_c,Vd_min)
            push!(Cad_max_c,Cad_max)
            push!(Cad_min_c,Cad_min)
            push!(Vs_max_c,Vs_max)
            push!(Vs_min_c,Vs_min)
            push!(Cas_max_c,Cas_max)
            push!(Cas_min_c,Cas_min)
        else
            push!(Vd_max_c,NaN)
            push!(Vd_min_c,NaN)
            push!(Cad_max_c,NaN)
            push!(Cad_min_c,NaN)
            push!(Vs_max_c,NaN)
            push!(Vs_min_c,NaN)
            push!(Cas_max_c,NaN)
            push!(Cas_min_c,NaN)
        end
    end

    return f_list,Vd_max_c,Vd_min_c,Cad_max_c,Cad_min_c,Vs_max_c,Vs_min_c,Cas_max_c,Cas_min_c
end
function find_full_bif(p,ic,I_list)

    I1 = maximum(I_list)+15/10^6

	bif = bif_desc_w_pulse(compartment!,I_list,p,ic,I1,20000)

	return bif
end