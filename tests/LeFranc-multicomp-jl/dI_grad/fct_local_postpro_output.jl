

function save_dI_mat(dI_mat,I1_mat,I2_mat,f_min_mat,jld_name)
    save(jld_name, "dI_mat", dI_mat, "I1_mat", I1_mat, "I2_mat", I2_mat, "f_min_mat", f_min_mat)
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
function load_dI_mat(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        dI = []
        I1 = []
        f_min = []
        I2 = []
        push!(dI,load_array(jld,"dI_mat"))
        push!(I1,load_array(jld,"I1_mat"))
        push!(f_min,load_array(jld,"f_min_mat"))
        push!(I2,load_array(jld,"I2_mat"))
    else
        dI = NaN
        I1 = NaN
        f_min = NaN
        I2 = NaN
        println("Not found")
    end
    return dI[1],I1[1],I2[1],f_min[1]
end

function create_dI_mat(jld_name,ext_name,pCaL_range,gKir_range)
    dI = ones(length(pCaL_range),length(gKir_range)).*NaN
    I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    f_min = ones(length(pCaL_range),length(gKir_range)).*NaN
    for i in 1:length(pCaL_range)
        for j in 1:length(gKir_range)
            jld_name_i = string(jld_name,ext_name[1],"$(i)",ext_name[2],"$(j)",ext_name[3])
            #println(jld_name_i)
            try 
                loaded_find_dI = load_dI(jld_name_i)[2]
                dI[i,j] = loaded_find_dI[1][1]
                I1[i,j] = loaded_find_dI[2][1]
                I2[i,j] = loaded_find_dI[3][1]
                f_min[i,j] = loaded_find_dI[4][1]
            catch err
                println(jld_name_i)
            end
            
        end
    end
  return dI,I1,I2,f_min
end
function create_V_Ix_mat(jld_name,ext_name,pCaL_range,gKir_range,p_fixed,is_gKir_s)
    dI,I1,I2,f_min = load_dI_mat(string(jld_name,".jld"))


    Vd_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cad_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Vs_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cas_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN

    Vd_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cad_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Vs_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cas_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN

    for i in eachindex(pCaL_range)
        for j in eachindex(gKir_range)
            if dI[i,j]!=NaN && I1[i,j]>-200/10^(6)
                jld_name_i = string(jld_name,ext_name[1],"$(i)",ext_name[2],"$(j)",ext_name[3])
                println("---------------- Iteration $i and $j ----------------")

                pCaLs_d_i_batch = pCaL_range[i] 

                if is_gKir_s
                    gKir_s_i_batch = gKir_range[j]

                    p_var_d_i = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
                    p_var_s_i = [gNa_s,gKDR_s,gKir_s_i_batch,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]
                else
                    gKir_d_i_batch = gKir_range[j]

                    p_var_d_i = [gNa_d,gKDR_d,gKir_d_i_batch,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
                    p_var_s_i = [gNa_s,gKDR_s,gKir_s,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]
                end

                if j>1 && Vd_I1[i,j-1] !==NaN 
                    solver_ic_i1 = [Vd_I1[i,j-1]-1,Cad_I1[i,j-1],Vs_I1[i,j-1]-1,Cas_I1[i,j-1]]
                else
                    if i>1 && Vd_I1[i-1,j] !==NaN 
                        solver_ic_i1 = [Vd_I1[i-1,j],Cad_I1[i-1,j],Vs_I1[i-1,j],Cas_I1[i-1,j]]
                    else
                        solver_ic_i1 = [-60.,5e-5,-60.,5e-5]                    
                    end
                end
                if j>1 && Vd_I2[i,j-1] !==NaN 
                    solver_ic_i2 = [Vd_I2[i,j-1]-1,Cad_I2[i,j-1],Vs_I2[i,j-1]-1,Cas_I2[i,j-1]]
                else
                    if i>1 && Vd_I2[i-1,j] !==NaN 
                        solver_ic_i2 = [Vd_I2[i-1,j],Cad_I2[i-1,j],Vs_I2[i-1,j],Cas_I2[i-1,j]]
                    else
                        solver_ic_i2 = [-50.,5e-5,-50.,5e-5]
                    end
                end
                z = find_Vr_I1_I2_comp([p_fixed,p_var_d_i,p_var_s_i],nlsolve_comp_grad_0_V_!,minimum([I1[i,j],I2[i,j]]),I2[i,j],solver_ic_i1,solver_ic_i2)
                println(z)
                
                if z[1][1] !==NaN
                    Vd_I1[i,j] = z[1][1]
                    Cad_I1[i,j] = z[2][1]
                    Vs_I1[i,j] = z[3][1]
                    Cas_I1[i,j] = z[4][1]
                    Vd_I2[i,j] = z[5][1]
                    Cad_I2[i,j] = z[6][1]
                    Vs_I2[i,j] = z[7][1]
                    Cas_I2[i,j] = z[8][1]
                else
                    solver_ic_i1 = [-70.,5e-5,-70.,5e-5]
                    solver_ic_i2 = [-55.,5e-5,-55.,5e-5]
                    z_ = find_Vr_I1_I2_comp([p_fixed,p_var_d_i,p_var_s_i],nlsolve_comp_grad_0_V_!,minimum([I1[i,j],I2[i,j]]),I2[i,j],solver_ic_i1,solver_ic_i2)
                    
                    Vd_I1[i,j] = z_[1][1]
                    Cad_I1[i,j] = z_[2][1]
                    Vs_I1[i,j] = z_[3][1]
                    Cas_I1[i,j] = z_[4][1]
                    println(z_)
                    Vd_I2[i,j] = z_[5][1]
                    Cad_I2[i,j] = z_[6][1]
                    Vs_I2[i,j] = z_[7][1]
                    Cas_I2[i,j] = z_[8][1]
                end
            end
        end
    end
    return Vd_I1,Cad_I1,Vs_I1,Cas_I1,Vd_I2,Cad_I2,Vs_I2,Cas_I2
end
function create_V_Ix_mat(jld_name,ext_name,pCaL_range,gKir_range,p_fixed,is_gKir_s,dI,I1,I2)

    Vd_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cad_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Vs_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cas_I1 = ones(length(pCaL_range),length(gKir_range)).*NaN

    Vd_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cad_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Vs_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN
    Cas_I2 = ones(length(pCaL_range),length(gKir_range)).*NaN

    for i in eachindex(pCaL_range)
        for j in eachindex(gKir_range)
            if dI[i,j]!=NaN && I1[i,j]>-200/10^(6)
                jld_name_i = string(jld_name,ext_name[1],"$(i)",ext_name[2],"$(j)",ext_name[3])
                println("---------------- Iteration $i and $j ----------------")

                pCaLs_d_i_batch = pCaL_range[i] 

                if is_gKir_s
                    gKir_s_i_batch = gKir_range[j]

                    p_var_d_i = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
                    p_var_s_i = [gNa_s,gKDR_s,gKir_s_i_batch,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]
                else
                    gKir_d_i_batch = gKir_range[j]

                    p_var_d_i = [gNa_d,gKDR_d,gKir_d_i_batch,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
                    p_var_s_i = [gNa_s,gKDR_s,gKir_s,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]
                end

                if j>1 && Vd_I1[i,j-1] !==NaN 
                    solver_ic_i1 = [Vd_I1[i,j-1]-1,Cad_I1[i,j-1],Vs_I1[i,j-1]-1,Cas_I1[i,j-1]]
                else
                    if i>1 && Vd_I1[i-1,j] !==NaN 
                        solver_ic_i1 = [Vd_I1[i-1,j],Cad_I1[i-1,j],Vs_I1[i-1,j],Cas_I1[i-1,j]]
                    else
                        solver_ic_i1 = [-60.,5e-5,-60.,5e-5]                    
                    end
                end
                if j>1 && Vd_I2[i,j-1] !==NaN 
                    solver_ic_i2 = [Vd_I2[i,j-1]-1,Cad_I2[i,j-1],Vs_I2[i,j-1]-1,Cas_I2[i,j-1]]
                else
                    if i>1 && Vd_I2[i-1,j] !==NaN 
                        solver_ic_i2 = [Vd_I2[i-1,j],Cad_I2[i-1,j],Vs_I2[i-1,j],Cas_I2[i-1,j]]
                    else
                        solver_ic_i2 = [-50.,5e-5,-50.,5e-5]
                    end
                end
                z = find_Vr_I1_I2_comp([p_fixed,p_var_d_i,p_var_s_i],nlsolve_comp_grad_0_V_!,minimum([I1[i,j],I2[i,j]]),I2[i,j],solver_ic_i1,solver_ic_i2)
                println(z)
                
                if z[1][1] !==NaN
                    Vd_I1[i,j] = z[1][1]
                    Cad_I1[i,j] = z[2][1]
                    Vs_I1[i,j] = z[3][1]
                    Cas_I1[i,j] = z[4][1]
                    Vd_I2[i,j] = z[5][1]
                    Cad_I2[i,j] = z[6][1]
                    Vs_I2[i,j] = z[7][1]
                    Cas_I2[i,j] = z[8][1]
                else
                    solver_ic_i1 = [-70.,5e-5,-70.,5e-5]
                    solver_ic_i2 = [-55.,5e-5,-55.,5e-5]
                    z_ = find_Vr_I1_I2_comp([p_fixed,p_var_d_i,p_var_s_i],nlsolve_comp_grad_0_V_!,minimum([I1[i,j],I2[i,j]]),I2[i,j],solver_ic_i1,solver_ic_i2)
                    
                    Vd_I1[i,j] = z_[1][1]
                    Cad_I1[i,j] = z_[2][1]
                    Vs_I1[i,j] = z_[3][1]
                    Cas_I1[i,j] = z_[4][1]
                    println(z_)
                    Vd_I2[i,j] = z_[5][1]
                    Cad_I2[i,j] = z_[6][1]
                    Vs_I2[i,j] = z_[7][1]
                    Cas_I2[i,j] = z_[8][1]
                end
            end
        end
    end
    return Vd_I1,Cad_I1,Vs_I1,Cas_I1,Vd_I2,Cad_I2,Vs_I2,Cas_I2
end

function find_Vr_I1_I2_comp(p,FP_fun,I1_,I2_,solver_ic_i1,solver_ic_i2)
    # This function computes the value of the stable equilibria associated with a resting state at I1 and I2 
    ## Fixed parameters
    p_fixed = p[1]
    p_fixed_0 = vcat(0,p_fixed[2:end])

    ## Varying parameters -- dendrite
	p_dend= p[2]

    ## Varying parameters -- soma
	p_soma = p[3]

    println("p_dend : $(p_dend)")
    println("p_soma : $(p_soma)")

    #==
    V=-100:0.01:-40
    find_lfp = find_fixed_points([vcat(0,p_fixed[2:end]),p_dend,p_soma],V,[25/10^(6),-59. ,Ca_i_0,Ca_i_0])
    plt_fp = plt_fixed_points(plot(plot(),plot(),plot(),plot(),layout=(4,1)),find_lfp)
    plot!(plt_fp,title="pCaLs_d = $(p_dend[8]) ; gKir_s = $(p_soma[3])")
    display(plt_fp)
    ==#

	Vd_fp_I1 = []
    Cad_fp_I1 = []
    Vs_fp_I1 = []
    Cas_fp_I1 = []
	Vd_fp_I2 = []
    Cad_fp_I2 = []
    Vs_fp_I2 = []
    Cas_fp_I2 = []

    p_fixed_I1 = vcat(I1_,p_fixed[2:end])
    p_fixed_I2 = vcat(I2_,p_fixed[2:end])

    println("p_fixed_I1 : $(p_fixed_I1)")
    println("solver_ic_i1 = $(solver_ic_i1)")
    FP_l = nlsolve((F,x) ->FP_fun(F,x,[p_fixed_I1,p_dend,p_soma]),solver_ic_i1)
    if converged(FP_l)
        println("Fixed point found at Vd=$(FP_l.zero[1]), Cad=$(FP_l.zero[2]), Vs=$(FP_l.zero[3]) & Cas=$(FP_l.zero[4]) for I1")
        lambda = real.(eigvals(J_matrix(FP_l.zero[1],FP_l.zero[2],FP_l.zero[3],FP_l.zero[4],[p_fixed_I1,p_dend,p_soma])) )
        println("lambda = $(lambda)")
        if maximum(lambda) <=0
            println("Valid fixed point for I1")
            push!(Vd_fp_I1,FP_l.zero[1])
            push!(Cad_fp_I1,FP_l.zero[2])
            push!(Vs_fp_I1,FP_l.zero[3])
            push!(Cas_fp_I1,FP_l.zero[4])
        else
            push!(Vd_fp_I1,NaN)
            push!(Cad_fp_I1,NaN)
            push!(Vs_fp_I1,NaN)
            push!(Cas_fp_I1,NaN)
        end
    else
        println("Unvalid fixed point ")
        push!(Vd_fp_I1,NaN)
        push!(Cad_fp_I1,NaN)
        push!(Vs_fp_I1,NaN)
        push!(Cas_fp_I1,NaN)
    end

    println("p_fixed_I2 : $(p_fixed_I2)")
    if Vs_fp_I1[end] !== NaN
        println(Vs_fp_I1[end])
        FP_h = nlsolve((F,x) ->FP_fun(F,x,[p_fixed_I2,p_dend,p_soma]), [Vd_fp_I1[end],Cad_fp_I1[end],Vs_fp_I1[end],Cas_fp_I1[end]])
    else
        FP_h = nlsolve((F,x) ->FP_fun(F,x,[p_fixed_I2,p_dend,p_soma]), solver_ic_i2)
    end
    println("Fixed point found at Vd=$(FP_h.zero[1]), Cad=$(FP_h.zero[2]), Vs=$(FP_h.zero[3]) & Cas=$(FP_h.zero[4]) for I2")
    lambda = real.(eigvals(J_matrix(FP_h.zero[1],FP_h.zero[2],FP_h.zero[3],FP_h.zero[4],[p_fixed_I2,p_dend,p_soma])) )
    println("lambda = $(lambda)")
    if maximum(lambda) <=0 && converged(FP_h)
        println("Valid fixed point for I2")
        push!(Vd_fp_I2,FP_h.zero[1])
        push!(Cad_fp_I2,FP_h.zero[2])
        push!(Vs_fp_I2,FP_h.zero[3])
        push!(Cas_fp_I2,FP_h.zero[4])
    else
        println("Unvalid")
        if Vs_fp_I1[end] !== NaN
            #retry with other x_i
            FP_h = nlsolve((F,x) ->FP_fun(F,x,[p_fixed_I2,p_dend,p_soma]), solver_ic_i2)
            println("Fixed point found at Vd=$(FP_h.zero[1]), Cad=$(FP_h.zero[2]), Vs=$(FP_h.zero[3]) & Cas=$(FP_h.zero[4]) for I2")
            lambda = real.(eigvals(J_matrix(FP_h.zero[1],FP_h.zero[2],FP_h.zero[3],FP_h.zero[4],[p_fixed_I2,p_dend,p_soma])) )
            println("lambda = $(lambda)")
            if maximum(lambda) <=0 && converged(FP_h)
                println("Valid fixed point for I2")
                push!(Vd_fp_I2,FP_h.zero[1])
                push!(Cad_fp_I2,FP_h.zero[2])
                push!(Vs_fp_I2,FP_h.zero[3])
                push!(Cas_fp_I2,FP_h.zero[4])
            else
                println("Unvalid")
                push!(Vd_fp_I2,NaN)
                push!(Cad_fp_I2,NaN)
                push!(Vs_fp_I2,NaN)
                push!(Cas_fp_I2,NaN)
            end
        else
            println("Vs_fp_I1[end] == NaN")
            push!(Vd_fp_I2,NaN)
            push!(Cad_fp_I2,NaN)
            push!(Vs_fp_I2,NaN)
            push!(Cas_fp_I2,NaN)
        end
    end
	return Vd_fp_I1,Cad_fp_I1,Vs_fp_I1,Cas_fp_I1,Vd_fp_I2,Cad_fp_I2,Vs_fp_I2,Cas_fp_I2
end


