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


    jld_name = "ChannelUpdate/cluster/dI_grad/data/jld_pCaLf_1.5e-5/dI_KirCaLs_pCaLf_1.5e-5"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KirCaLs = ones(length(pCaLs_range_KirCaLs[1:2:end]),length(gKir_range_KirCaLs[1:3:end])).*NaN
        I1_KirCaLs = ones(length(pCaLs_range_KirCaLs[1:2:end]),length(gKir_range_KirCaLs[1:3:end])).*NaN
        I2_KirCaLs = ones(length(pCaLs_range_KirCaLs[1:2:end]),length(gKir_range_KirCaLs[1:3:end])).*NaN
        VSN_KirCaLs = ones(length(pCaLs_range_KirCaLs[1:2:end]),length(gKir_range_KirCaLs[1:3:end])).*NaN
        for i in 1:length(pCaLs_range_KirCaLs[1:2:end])
            for j in 1:length(gKir_range_KirCaLs[1:3:end])
                jld_name_i = string(jld_name,"_i_pCaLs_$(i)_i_gKir_$(j).jld")
                #println(jld_name_i)
                try 
                    loaded_find_dI = load_dI(jld_name_i)[2]
                    dI_KirCaLs[i,j] = loaded_find_dI[1][1]
                    I1_KirCaLs[i,j] = loaded_find_dI[2][1]
                    I2_KirCaLs[i,j] = loaded_find_dI[5][1]
                    VSN_KirCaLs[i,j] = loaded_find_dI[4][1]
                catch err
                    println(jld_name_i)
                end
            end
        end
        Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],dI_KirCaLs,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
        save_dI_mat(dI_KirCaLs,I1_KirCaLs,I2_KirCaLs,VSN_KirCaLs,string(jld_name,"_1e-2_1e-5.jld"))
    ==#

    dI_KirCaLs_pCaLf_15e_5,I1_KirCaLs_pCaLf_15e_5,I2_KirCaLs_pCaLf_15e_5,VSN_KirCaLs_pCaLf_15e_5 = load_dI_mat(string(jld_name,"_1e-2_1e-5.jld"))
    plt = Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],dI_KirCaLs_pCaLf_15e_5,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(0,3.719))
    plt = Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],dI_KirCaLs[1:2:end,1:3:end],fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(0,3.719))
  
    plt = Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],I1_KirCaLs_pCaLf_15e_5,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(-3.1,3.5))
    plt = Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],I1_KirCaLs[1:2:end,1:3:end],fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(-3.1,3.5))
  
    plt = Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],I2_KirCaLs_pCaLf_15e_5,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    plt = Plots.heatmap(gKir_range_KirCaLs[1:3:end],pCaLs_range_KirCaLs[1:2:end],I2_KirCaLs[1:2:end,1:3:end],fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
  


    jld_name = "ChannelUpdate/cluster/dI_grad/data/jld_pCaLf_1.5e-5/dI_KMCaLs_pCaLf_1.5e-5"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KMCaLs = ones(length(pCaLs_range_KMCaLs[1:2:end]),length(gKM_range_KMCaLs[1:3:end])).*NaN
        I1_KMCaLs = ones(length(pCaLs_range_KMCaLs[1:2:end]),length(gKM_range_KMCaLs[1:3:end])).*NaN
        I2_KMCaLs = ones(length(pCaLs_range_KMCaLs[1:2:end]),length(gKM_range_KMCaLs[1:3:end])).*NaN
        VSN_KMCaLs = ones(length(pCaLs_range_KMCaLs[1:2:end]),length(gKM_range_KMCaLs[1:3:end])).*NaN
        for i in 1:length(pCaLs_range_KMCaLs[1:2:end])
            for j in 1:length(gKM_range_KMCaLs[1:3:end])
                jld_name_i = string(jld_name,"_i_pCaLs_$(i)_i_gKM_$(j).jld")

                try 
                    loaded_find_dI = load_dI(jld_name_i)[2]
                    dI_KMCaLs[i,j] = loaded_find_dI[1][1]
                    I1_KMCaLs[i,j] = loaded_find_dI[2][1]
                    I2_KMCaLs[i,j] = loaded_find_dI[5][1]
                    VSN_KMCaLs[i,j] = loaded_find_dI[4][1]
                catch err
                    println("Not found: ",jld_name_i)
                end
            end
        end
        Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],dI_KMCaLs,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
        save_dI_mat(dI_KMCaLs,I1_KMCaLs,I2_KMCaLs,VSN_KMCaLs,string(jld_name,"_1e-2_1e-5.jld"))
    ==#
    
    dI_KMCaLs_pCaLf_15e_5,I1_KMCaLs_pCaLf_15e_5,I2_KMCaLs_pCaLf_15e_5,VSN_KMCaLs_pCaLf_15e_5 = load_dI_mat(string(jld_name,"_1e-2_1e-5.jld"))
    plt = Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],dI_KMCaLs_pCaLf_15e_5,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(0,3.719))
    plt = Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],dI_KMCaLs[1:2:end,1:3:end],fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(0,3.719))
  
    plt = Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],I1_KMCaLs_pCaLf_15e_5,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(-3.1,6.3))
    plt = Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],I1_KMCaLs[1:2:end,1:3:end],fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(-3.1,6.3))
  
    plt = Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],I2_KMCaLs_pCaLf_15e_5,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    plt = Plots.heatmap(gKM_range_KMCaLs[1:3:end],pCaLs_range_KMCaLs[1:2:end],I2_KMCaLs[1:2:end,1:3:end],fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
  