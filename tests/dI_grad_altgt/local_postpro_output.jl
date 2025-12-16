using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

include("data_loader.jl")

function save_dI_mat(dI_mat,I1_mat,I2_mat,VSN_mat,jld_name)
    save(jld_name, "dI_mat", dI_mat, "I1_mat", I1_mat, "I2_mat", I2_mat, "VSN_mat", VSN_mat)
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
function load_dI_mat(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        dI = []
        I1 = []
        V_SN = []
        I2 = []
        push!(dI,load_array(jld,"dI_mat"))
        push!(I1,load_array(jld,"I1_mat"))
        push!(V_SN,load_array(jld,"VSN_mat"))
        push!(I2,load_array(jld,"I2_mat"))
    else
        dI = NaN
        I1 = NaN
        V_SN = NaN
        I2 = NaN
        println("Not found")
    end
    return dI[1],I1[1],I2[1],V_SN[1]
end

function find_Vr_I1_I2_clust(p,FP_fun,pCaLf_,pCaLs_,gKDR_,gKir_,gKM_,I1_,I2_,V_est_I1,Ca_est_I1,V_SN_,Ca_SN_)
    # This function computes the value of the stable equilibria associated with a resting state at I1 and I2 
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

	tol_FP = 1e-6
	V_fp_I1 = []
    Ca_fp_I1 = []
	V_fp_I2 = []
    Ca_fp_I2 = []

	for i in eachindex(p_to_test)
        ### print i and length of p_to_test
        ### modify most of print to make them more explicit
        println("------ Parameter vector $i on $(length(p_to_test)) -----")
		p_var_i = p_to_test[i]

        ## Varying parameters
	    gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p_var_i[2:end]


        I_fixed_to_I1 = I1_[i] 
        p_var_i_I1 = [I_fixed_to_I1,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]


		FP_l = nlsolve((F,x) ->FP_fun(F,x,[p_var_i_I1,p_fixed]), [V_est_I1[i],Ca_est_I1[i]],ftol=1E-7,xtol=1E-7)
		res_l = FP_fun([NaN,NaN],FP_l.zero,[p_var_i_I1,p_fixed])
        println("Fixed point found at V=$(FP_l.zero[1]), and Ca=$(FP_l.zero[2]) with residual of $(res_l) for I1")
		if res_l < tol_FP
            lambda = real.(eigvals(J_matrix(FP_l.zero[1],FP_l.zero[2],[p_var_i_I1,p_fixed])) )
            println("lambda = $(lambda)")
            if maximum(lambda) <=0
                println("Valid fixed point for I1")
			    push!(V_fp_I1,FP_l.zero[1])
                push!(Ca_fp_I1,FP_l.zero[2])
            end
		else
            println("Unvalid")
			push!(V_fp_I1,NaN)
            push!(Ca_fp_I1,NaN)
		end

        I_fixed_to_I2 = I2_[i] 
        p_var_i_I2 = [I_fixed_to_I2,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak] 
        println(p_var_i_I2)
        if V_fp_I1[end] !== NaN
            println(V_fp_I1[end])
            FP_h = nlsolve((F,x) ->FP_fun(F,x,[p_var_i_I2,p_fixed]), [V_fp_I1[end],Ca_fp_I1[end]],ftol=1E-7,xtol=1E-7)
        else
            FP_h = nlsolve((F,x) ->FP_fun(F,x,[p_var_i_I2,p_fixed]), [V_SN_[i],Ca_SN_[i]],ftol=1E-7,xtol=1E-7)
        end
        res_h = FP_fun([NaN,NaN],FP_h.zero,[p_var_i_I2,p_fixed])
        println("Fixed point found at V=$(FP_h.zero[1]), and Ca=$(FP_h.zero[2]) with residual of $(res_h) for I2")
        lambda = real.(eigvals(J_matrix(FP_h.zero[1],FP_h.zero[2],[p_var_i_I2,p_fixed])) )
        println("lambda = $(lambda)")
        if maximum(lambda) <=0 && res_h < tol_FP
            println("Valid fixed point for I2")
            push!(V_fp_I2,FP_h.zero[1])
            push!(Ca_fp_I2,FP_h.zero[2])
        else
            println("Unvalid")
            if V_fp_I1[end] !== NaN
                #retry with other x_i
                FP_h = nlsolve((F,x) ->FP_fun(F,x,[p_var_i_I2,p_fixed]), [V_SN_[i],Ca_SN_[i]],ftol=1E-7,xtol=1E-7)
                res_h = FP_fun([NaN,NaN],FP_h.zero,[p_var_i_I2,p_fixed])
                println("Fixed point found at V=$(FP_h.zero[1]), and Ca=$(FP_h.zero[2]) with residual of $(res_h) for I2")
                lambda = real.(eigvals(J_matrix(FP_h.zero[1],FP_h.zero[2],[p_var_i_I2,p_fixed])) )
                println("lambda = $(lambda)")
                if maximum(lambda) <=0 && res_h < tol_FP
                    println("Valid fixed point for I2")
                    push!(V_fp_I2,FP_h.zero[1])
                    push!(Ca_fp_I2,FP_h.zero[2])
                else
                    println("Unvalid")
                    push!(V_fp_I2,NaN)
                    push!(Ca_fp_I2,NaN)
                end
            end
        end
	end
	return V_fp_I1,V_fp_I2,Ca_fp_I1,Ca_fp_I2
end


    jld_name = "tests/dI_grad_altgt/data/jld/dI_KirKMCaLs"
    #
    # This block loads the results from the separated jld files and combine them in a single jld 
        dI_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
        I1_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
        I2_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
        VSN_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
        for i in 1:length(pCaLs_range_KirKMCaLs)
            for j in 1:length(gKir_range_KirKMCaLs)
                for k in 1:length(gKM_range_KirKMCaLs)
                    jld_name_i = string(jld_name,"_i_pCaLs_$(i)_i_Kir_$(j)_i_gKM_$(k).jld")
                    #println(jld_name_i)
                    try 
                        loaded_find_dI = load_dI(jld_name_i)[2]
                        dI_KirKMCaLs[i,j,k] = loaded_find_dI[1][1]
                        I1_KirKMCaLs[i,j,k] = loaded_find_dI[2][1]
                        I2_KirKMCaLs[i,j,k] = loaded_find_dI[5][1]
                        VSN_KirKMCaLs[i,j,k] = loaded_find_dI[4][1]
                    catch err
                        println(jld_name_i)
                    end
                end
            end
        end
        save_dI_mat(dI_KirKMCaLs,I1_KirKMCaLs,I2_KirKMCaLs,VSN_KirKMCaLs,string(jld_name,".jld"))
    ==#

    dI_KirKMCaLs,I1_KirKMCaLs,I2_KirKMCaLs,VSN_KirKMCaLs = load_dI_mat(string(jld_name,".jld"))

    V_I1_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
    V_I2_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
    Ca_I1_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
    Ca_I2_KirKMCaLs = ones(length(pCaLs_range_KirKMCaLs),length(gKir_range_KirKMCaLs),length(gKM_range_KirKMCaLs)).*NaN
    for i in 1:length(pCaLs_range_KirKMCaLs)
        for j in 1:length(gKir_range_KirKMCaLs)
            println("---------------- Iteration $i and $j ----------------")
            for k in 1:length(gKM_range_KirKMCaLs)
                jld_name_i = string(jld_name,"_i_pCaLs_$(i)_i_Kir_$(j)_i_gKM_$(k).jld")
                println("$(k) / $(length(gKM_range_KirKMCaLs))")
                p_var_i_j = [0.,gNa,gKDR,gKir_range_KirKMCaLs[j],gKM_range_KirKMCaLs[k],pCaLf,pCaLs_range_KirKMCaLs[i],gleak]
                if k>1
                    z = find_Vr_I1_I2_clust([p_var_i_j,p_fixed],nlsolve_membrane_grad_0_V_!,[],[],[],[gKir_range_KirKMCaLs[j]],[],[minimum([I1_KirKMCaLs[i,j,k],I2_KirKMCaLs[i,j,k]])],[I2_KirKMCaLs[i,j,k]],[V_I1_KirKMCaLs[i,j,k-1]].-0.02 ,[Ca_I1_KirKMCaLs[i,j,k-1]],[VSN_KirKMCaLs[i,j,k-1]].-0.02,[1e-3])
                else
                    #if i>1
                        #z = find_Vr_I1_I2_clust([p_var_i_j,p_fixed],nlsolve_membrane_grad_0_V_!,[],[],[],[gKir_range_KirKMCaLs[j]],[],[minimum([I1_KirKMCaLs[i,j,k],I2_KirKMCaLs[i,j,k]])],[I2_KirKMCaLs[i,j,k]],[V_I1_KirKMCaLs[i-1,j,k]].-0.02 .-20*(i==2) .+((k-1)*20/31),[Ca_I1_KirKMCaLs[i-1,j,k]],[VSN_KirKMCaLs[i-1,j,k]].-0.02,[1e-3])
                    if j>1
                        z = find_Vr_I1_I2_clust([p_var_i_j,p_fixed],nlsolve_membrane_grad_0_V_!,[],[],[],[gKir_range_KirKMCaLs[j]],[],[minimum([I1_KirKMCaLs[i,j,k],I2_KirKMCaLs[i,j,k]])],[I2_KirKMCaLs[i,j,k]],[V_I1_KirKMCaLs[i,j-1,k]].-1,[0],[VSN_KirKMCaLs[i,j-1,k]].-0.02,[1e-3])
                    else
                        if i>1
                            z = find_Vr_I1_I2_clust([p_var_i_j,p_fixed],nlsolve_membrane_grad_0_V_!,[],[],[],[gKir_range_KirKMCaLs[j]],[],[minimum([I1_KirKMCaLs[i,j,k],I2_KirKMCaLs[i,j,k]])],[I2_KirKMCaLs[i,j,k]],[V_I1_KirKMCaLs[i-1,j,k]].-0.02 .-20*(i-1)/4 ,[Ca_I1_KirKMCaLs[i-1,j,k]],[VSN_KirKMCaLs[i-1,j,k]].-0.02,[1e-3])
                        else
                            z = find_Vr_I1_I2_clust([p_var_i_j,p_fixed],nlsolve_membrane_grad_0_V_!,[],[],[],[gKir_range_KirKMCaLs[j]],[],[minimum([I1_KirKMCaLs[i,j,k],I2_KirKMCaLs[i,j,k]])],[I2_KirKMCaLs[i,j,k]],[-65 -20*(i-1)/2],[1e-5],[-60].-0.02,[1e-3])
                            
                        end
                    end
                end
                V_I1_KirKMCaLs[i,j,k] = z[1][1]
                V_I2_KirKMCaLs[i,j,k] = z[2][1]
                Ca_I1_KirKMCaLs[i,j,k] = z[3][1]
                Ca_I2_KirKMCaLs[i,j,k] = z[4][1]
            end
        end
    end

    pal_V = cgrad([:royalblue2, :white, :springgreen4, :gold, :goldenrod,:hotpink2, :crimson], [0.02, 0.0329,0.1129,0.1229, 0.2057, 0.2929, 0.5, 0.75, 0.9])
  
    i_pCaLs = 5
    clim_ = [0.1,0.2,0.4,0.625,0.75]
    plt = Plots.heatmap(gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,dI_KirKMCaLs[i_pCaLs,:,:],fill=true,levels=100,lw=0,c=:vik,xlabel="gKM",ylabel="gKir",clims=(-1,+1).*clim_[i_pCaLs] .+dI_KirKMCaLs[i_pCaLs,1,1],fontfamily="Computer Modern")
    #savefig(plt,"dI_KirKM_ipCaLs$(i_pCaLs).png")

    plt = Plots.heatmap(gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,V_I1_KirKMCaLs[i_pCaLs,:,:],fill=true,levels=100,lw=0,clims=(-125,-55),c=cgrad(pal_V,rev=true),xlabel="gKM",ylabel="gKir",fontfamily="Computer Modern")
    #savefig(plt,"V_I1_KirKM_ipCaLs$(i_pCaLs).png")
    plt = Plots.heatmap(gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,V_I2_KirKMCaLs[i_pCaLs,:,:],fill=true,levels=100,lw=0,c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),clims=(-66,-56),xlabel="gKM",ylabel="gKir",fontfamily="Computer Modern")
    #savefig(plt,"V_I2_KirKM_ipCaLs$(i_pCaLs).png")


