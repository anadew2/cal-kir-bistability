#Functions used in script_dI_grad.jl with compartment.jl

pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)
heaviside(t)=0*(t<0)+1*(t>=0)


function bif_asc_w_pulse_if_bist(ode_fun,I_list,p,ic,I1,tmin,fp_V)    
    ## Fixed parameters
    p_fixed = p[1]

    ## Varying parameters -- dendrite
	p_dend= p[2]

    ## Varying parameters -- soma
	p_soma = p[3]
    
    f_list = []
    I_tested_list = []
    f_state_list = []

    t_pulse = 5000;
    dt_pulse = 5000;

    f_tol = 0.01
    lim_cnt=5

    for i in eachindex(I_list)

        println("$i on $(length(I_list))")
        println("I[i] = $(I_list[i])")

        tspan = (0.0,tmin)
        saveat_range = 0:0.1:tmin

        I0 = I_list[i]
        push!(I_tested_list,I0)

        p_fixed_i = vcat(I0,p_fixed[2:end])
        It_fixed(t) = I0 
        p_fixed_it_cst = vcat(It_fixed,p_fixed[2:end])
        It_pulse(t) = I0 + (I1-I0)*pulse(t, t_pulse, t_pulse+dt_pulse)
        p_fixed_it_pulse = vcat(It_pulse,p_fixed[2:end])

        I_stable_fp_V = fp_V[10]
        i_I_stable_fp_V_ = collect(1:length(I_stable_fp_V))[I_stable_fp_V.>=I0]
        if length(i_I_stable_fp_V_)>0
            solver_ic = [fp_V[6][i_I_stable_fp_V_[1]],fp_V[7][i_I_stable_fp_V_[1]],fp_V[8][i_I_stable_fp_V_[1]],fp_V[9][i_I_stable_fp_V_[1]]]
        else
            println("max(I_stable_fp_V) = $(maximum(I_stable_fp_V[I_stable_fp_V.<1000]) )")
            println("I0 = $(I0)")
            if i>1
                solver_ic = [f_state_list[end][1],f_state_list[end][13],f_state_list[end][14],f_state_list[end][26]]
            else
                solver_ic = [-60.,5e-5,-60.,5e-5]
            end
        end        
        println("solver_ic = $(solver_ic)")

        fp_i = find_fixed_points_I([p_fixed_i,p_dend,p_soma],[I0],solver_ic)
        println("fp_i = $(fp_i)")
        
        if length(fp_i[10]) >0
            #Cycle convergence
            freq_test_tol= f_tol+1
            f_I = []
            cnt = 0

            Vd_ic = fp_i[6][1]
            Cad_ic = fp_i[7][1]
            Vs_ic = fp_i[8][1]
            Cas_ic = fp_i[9][1]

            ic = [Vd_ic,mNa(Vd_ic),hNa(Vd_ic),mKDR(Vd_ic),mCaLf(Vd_ic),hCaLf(Vd_ic),mCaLs(Vd_ic),hCaLs(Vd_ic),mKir(Vd_ic),mKM(Vd_ic),mKCa(Cad_ic),mCAN(Cad_ic),Cad_ic,
                Vs_ic,mNa(Vs_ic),hNa(Vs_ic),mKDR(Vs_ic),mCaLf(Vs_ic),hCaLf(Vs_ic),mCaLs(Vs_ic),hCaLs(Vs_ic),mKir(Vs_ic),mKM(Vs_ic),mKCa(Cas_ic),mCAN(Cas_ic),Cas_ic]


            prob_pulse = ODEProblem(ode_fun,ic,tspan./2,[p_fixed_it_pulse,p_dend,p_soma])
            sol_pulse = solve(prob_pulse,dtmax=1,Tsit5(),save_everystep=false,save_start=false)
            f_state=sol_pulse[end]

            while freq_test_tol > f_tol && cnt <lim_cnt
                cnt = cnt +1
                println("f_state = $(f_state)")

                cb = ContinuousCallback(condition,affect!,nothing,save_positions=(true,false))
                prob = ODEProblem(ode_fun,f_state,tspan,[p_fixed_it_cst,p_dend,p_soma])
                sol_all =solve(prob,dtmax=1,Tsit5(),callback=cb,save_everystep=false,save_start=false,save_end=false,saveat=saveat_range)	
                tspk_sol = setdiff(sol_all.t,saveat_range)

                if length(tspk_sol) >3 
                    freq_isi = 1000 ./(tspk_sol[2:end]-tspk_sol[1:end-1])
                    freq_test_tol = (maximum(freq_isi) - minimum(freq_isi))

                    f_state[:]=sol_all[end-1][:]
                else
                    freq_test_tol = 0.
                    freq_isi = []
                end

                tspan = tspan .+ tmin

                if length(freq_isi)>1
                    mean_f = sum(freq_isi)/length(freq_isi)
                else
                    mean_f = 0
                    if length(sol_all)>0
                        println("size(sol_all) = $(size(sol_all))")
                    end
                end
                println("Average frequency : $(mean_f)")
                push!(f_I, mean_f)

            end
            push!(f_state_list,f_state)
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
            if length(fp_i[5]) >0
                println("is a saddle")
            end
            if length(fp_i[15]) >0
                println("is unstable")
            end
        end
    end

    return f_list,I_tested_list
end
function find_bif(p,ic,I_list_lf,I1,fp_V)

	bif_lf = bif_asc_w_pulse_if_bist(compartment!,I_list_lf,p,ic,I1,120000,fp_V)
	println(bif_lf[1])

	I1_bist = bif_lf[2][end] #cycle appears

	return  I1_bist,bif_lf
end
function find_dI_compart(p,ic,I_X)
    ## Fixed parameters
    p_fixed = p[1]
    p_fixed_0 = vcat(0,p_fixed[2:end])

    ## Varying parameters -- dendrite
	p_dend= p[2]

    ## Varying parameters -- soma
	p_soma = p[3]

    println("p_fixed_0 : $(p_fixed_0)")
    println("p_dend : $(p_dend)")
    println("p_soma : $(p_soma)")

    stp_V_FP_range=0.01
    V_FP_range = range(-100.,stop=-40.,step=stp_V_FP_range)

    I1_ = []
    I2_ = []
    dI_ = []
    f_min_ =  []

        println("---------------------------------------")

        solver_ic_v = [25/10^(6),-59. ,Ca_i_0,Ca_i_0]
        fp_V = find_fixed_points([p_fixed_0,p_dend,p_soma],V_FP_range,solver_ic_v)
        I_stable = fp_V[10]       

        i_I_stable = collect(1:length(I_stable))[I_stable.==maximum(I_stable[I_stable.<1000])]
        if length(i_I_stable)>0
            if fp_V[7][i_I_stable[1]] >=0 && fp_V[9][i_I_stable[1]] >=0
                println("Valid I2, Cad = $(fp_V[7][i_I_stable[1]])>=0 and Cas = $(fp_V[9][i_I_stable[1]])")
                I2 = maximum(I_stable[I_stable.<1000])
            else
                println("max(I_stable) is unvalid, searching for a valid one")
                i_max_unvalid = []
                push!(i_max_unvalid,i_I_stable[1])
                i_max_valid = []
                lim_cnt_I_stable_max = 100
                cnt_I_stable_max = 1
                indices_I_stable = collect(1:length(I_stable))
                while length(i_max_valid)==0 && cnt_I_stable_max<=lim_cnt_I_stable_max 
                    cnt_I_stable_max = cnt_I_stable_max +1
                    #println("cnt_I_stable_max = $cnt_I_stable_max")
                    #println("i_max_unvalid = $(i_max_unvalid)")
                    #println("length(indices_I_stable) = $(length(indices_I_stable))")
                    deleteat!(indices_I_stable,sort(i_max_unvalid))
                    pop!(i_max_unvalid)

                    i_I_stable_ = indices_I_stable[I_stable[indices_I_stable].==maximum(I_stable[indices_I_stable][I_stable[indices_I_stable].<1000])]
                    #println("i_I_stable_ = $(i_I_stable_)")
                    #println("maximum(I_stable[indices_I_stable][I_stable[indices_I_stable].<1000]) = $(maximum(I_stable[indices_I_stable][I_stable[indices_I_stable].<1000]))")
                    if length(i_I_stable_)>0 
                        if fp_V[7][i_I_stable_[1]] >=0 && fp_V[9][i_I_stable_[1]] >=0
                            push!(i_max_valid,i_I_stable_[1])
                            println("Valid max(I_stable): (Vd,Cad,Vs,Cas,I)= $((fp_V[6][i_I_stable_[1]],fp_V[7][i_I_stable_[1]],fp_V[8][i_I_stable_[1]],fp_V[9][i_I_stable_[1]],fp_V[10][i_I_stable_[1]]))")
                            I2 = fp_V[10][i_max_valid[1]]
                        else
                            push!(i_max_unvalid,i_I_stable_[1])
                        end
                    else
                        println("Issue in search for a valid I2") 
                    end
                end
                if length(i_max_valid)>0
                    I2 = fp_V[10][i_max_valid[1]]
                else
                    I2 = maximum(I_stable[I_stable.<1000])
                    println("I2 might be unvalid")
                end
            end
        else
            println("i_I_stable = $(i_I_stable)")
            println("--> Could not verify that Cad and Cas are positive at max(I_stable)")
            I2 = maximum(I_stable[I_stable.<1000])
            
        end

        I1_pulse = I2 + 50/10^(6)
        println("I_X = $(I_X), I1_pulse = $(I1_pulse)") 

        I_l_lf = I_X[1] 
        I_h_lf = I2
    
        if I_h_lf-I_l_lf > 2*10^(2)./10^(6)
            big_Istep = [10^(2),10^(1),10^(0),10^(-1),10^(-2),10^(-3)]./10^(6)
        else
            if I_h_lf-I_l_lf > 2*10^(1)./10^(6)
                big_Istep = [10^(1),10^(0),10^(-1),10^(-2),10^(-3)]./10^(6)
            else
                if I_h_lf-I_l_lf > 2*10^(0)./10^(6)
                    big_Istep = [10^(0),10^(-1),10^(-2),10^(-3)]./10^(6)
                else
                    if I_h_lf-I_l_lf > 2*10^(-1)./10^(6)
                        big_Istep = [10^(-1),10^(-2),10^(-3)]./10^(6)
                    else
                        if I_h_lf-I_l_lf > 2*10^(-2)./10^(6)
                            big_Istep = [10^(-2),10^(-3)]./10^(6)
                        else
                            big_Istep = [10^(-3)]./10^(6)
                        end
                    end
                end
            end
        end

        println("big_Istep = $(big_Istep)")
        
        I_list_lf_big = reverse(collect(range(I_h_lf+2*big_Istep[1],stop=I_l_lf-big_Istep[1], step=-big_Istep[1])))
        println("Current test in between: $I_l_lf && $(I_h_lf)")

        
        popfirst!(big_Istep)
        if big_Istep[1]==10^(1)/10^(6)
            I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse + 2*10^2/10^(6),fp_V)
        else
            I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse,fp_V)
        end

        if bif_lf[1][end]>0 
            if length(big_Istep)>0

                I_list_lf_big = vcat(bif_lf[2][end-1],reverse(collect(range(I1_est,stop=bif_lf[2][end-1], step=-big_Istep[1]))))
                popfirst!(big_Istep)                
                I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse,fp_V)
                
                if length(big_Istep)>0
                    if length(bif_lf[2])>1
                        I_list_lf_big = vcat(bif_lf[2][end-1],reverse(collect(range(I1_est,stop=bif_lf[2][end-1], step=-big_Istep[1]))))
                    else
                        I_list_lf_big = vcat(bif_lf[2][end]-5*10*big_Istep[1],reverse(collect(range(I1_est,stop=bif_lf[2][end]-5*10*big_Istep[1], step=-big_Istep[1]))))
                    end
                    popfirst!(big_Istep)
                    I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse,fp_V)
               
                    if length(big_Istep)>0
                    
                        I_list_lf_big = vcat(bif_lf[2][end-1],reverse(collect(range(I1_est,stop=bif_lf[2][end-1], step=-big_Istep[1]))))
                        popfirst!(big_Istep)
                        I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse,fp_V)

                        if length(big_Istep)>0
                    
                            I_list_lf_big = vcat(bif_lf[2][end-1],reverse(collect(range(I1_est,stop=bif_lf[2][end-1], step=-big_Istep[1]))))
                            popfirst!(big_Istep)
                            I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse,fp_V)

                            if length(big_Istep)>0
                    
                                I_list_lf_big = vcat(bif_lf[2][end-1],reverse(collect(range(I1_est,stop=bif_lf[2][end-1], step=-big_Istep[1]))))
                                popfirst!(big_Istep)
                                I1_est,bif_lf = find_bif([p_fixed_0,p_dend,p_soma],ic,I_list_lf_big,I1_pulse,fp_V)
                        
                            end
                        end
                    end
                end
            end

            push!(I1_,I1_est)
            push!(I2_,I2)
            push!(dI_,maximum([0,I2-I1_est]))

            push!(f_min_,bif_lf[1][end])

            println("I1 = $(I1_est)   &    I2 = $(I2)   &   dI = $(maximum([0,I2-I1_est]))   &   f_min = $(bif_lf[1][end]) ")   
        else
            #par hypothèse sur le diagramme de bifurcation, si pas de cycle au point stable le plus à droite alors I1=I_SN 
            push!(I1_,I2)
            push!(I2_,I2)
            push!(dI_,maximum([0,I2-I2]))

            push!(f_min_,bif_lf[1][end])

            println("I1 = $(I2)   &    I2 = $(I2)   &   dI = $(maximum([0,I2-I2]))   &   f_min = $(bif_lf[1][end])  ")
        end

    #end
    return dI_,I1_,I2_,f_min_
end

