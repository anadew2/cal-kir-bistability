using JLD

include("membrane_no_ghk_m2.jl")
include("fixedpoints_no_ghk_m2.jl")

include("data_loader_larger_range_no_ghk_m2.jl")

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


    gCaLs_range_KirCaLs_no_ghk_m2 = range(0.,stop=0.5,length=21)
    gKir_range_KirCaLs_no_ghk_m2 = range(0.,stop=1.,length=11)
    gCaLs_range_KMCaLs_no_ghk_m2 = range(0.,stop=0.5,length=21)
    gKM_range_KMCaLs_no_ghk_m2 = range(0.,stop=1.5,length=16)


    jld_name = "ChannelUpdate/cluster/dI_grad_no_ghk_m2/data/jld2/jld2/dI_KirCaLs_larger_range_no_ghk_m2"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KirCaLs = ones(length(gCaLs_range_KirCaLs_no_ghk_m2),length(gKir_range_KirCaLs_no_ghk_m2)).*NaN
        I1_KirCaLs = ones(length(gCaLs_range_KirCaLs_no_ghk_m2),length(gKir_range_KirCaLs_no_ghk_m2)).*NaN
        I2_KirCaLs = ones(length(gCaLs_range_KirCaLs_no_ghk_m2),length(gKir_range_KirCaLs_no_ghk_m2)).*NaN
        VSN_KirCaLs = ones(length(gCaLs_range_KirCaLs_no_ghk_m2),length(gKir_range_KirCaLs_no_ghk_m2)).*NaN
        for i in 1:length(gCaLs_range_KirCaLs_no_ghk_m2)
            for j in 1:length(gKir_range_KirCaLs_no_ghk_m2)
                jld_name_i = string(jld_name,"_i_gCaLs_$(i)_i_gKir_$(j).jld")
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
        save_dI_mat(dI_KirCaLs,I1_KirCaLs,I2_KirCaLs,VSN_KirCaLs,string(jld_name,".jld"))
    ==#

    dI_KirCaLs_no_ghk_m2,I1_KirCaLs_no_ghk_m2,I2_KirCaLs_no_ghk_m2,VSN_KirCaLs_no_ghk_m2 = load_dI_mat(string(jld_name,".jld"))
    plt = Plots.heatmap(gKir_range_KirCaLs_no_ghk_m2,gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    
    plt_kir_no_ghk_m2 = Plots.plot(fontfamily="Computer Modern",xlabel="gCaLs",ylabel="delta I ",title="Original m_inf,  ICaLs =gCaLs.m^2.h.(V-ECa)   ",ylims=(-0.1,3.1),xlims=(0,0.5))
    Plots.plot!(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2[:,1],label="gKir=0",marker=:rect,lw=2,markerstrokewidth=0.,ms=3,fontfamily="Computer Modern",legend=:outertopright,lc=:black,mc=:black,palette=palette(:PuRd_9,rev=true))
    Plots.plot!(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2[:,2],label="gKir=0.1",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2[:,3],label="gKir=0.2",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2[:,4],label="gKir=0.3",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2[:,6],label="gKir=0.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2[:,11],label="gKir=1.",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)


    jld_name = "ChannelUpdate/cluster/dI_grad_no_ghk_m2/data/jld2/jld2/dI_KMCaLs_larger_range_no_ghk_m2"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KMCaLs = ones(length(gCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        I1_KMCaLs = ones(length(gCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        I2_KMCaLs = ones(length(gCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        VSN_KMCaLs = ones(length(gCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        for i in 1:length(gCaLs_range_KMCaLs)
            for j in 1:length(gKM_range_KMCaLs)
                jld_name_i = string(jld_name,"_i_gCaLs_$(i)_i_gKM_$(j).jld")

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
        save_dI_mat(dI_KMCaLs,I1_KMCaLs,I2_KMCaLs,VSN_KMCaLs,string(jld_name,".jld"))
    ==#

    dI_KMCaLs_no_ghk_m2,I1_KMCaLs_no_ghk_m2,I2_KMCaLs_no_ghk_m2,VSN_KMCaLs_no_ghk_m2 = load_dI_mat(string(jld_name,".jld"))
    plt = Plots.heatmap(gKM_range_KMCaLs_no_ghk_m2,gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    
    plt_km_no_ghk_m2 = Plots.plot(fontfamily="Computer Modern",xlabel="gCaLs",ylabel="delta I ",title="Original m_inf,  ICaLs =gCaLs.m^2.h.(V-ECa)   ",ylims=(-0.1,3.1),xlims=(0,0.5))
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,1],label="gKM=0",marker=:rect,lw=2,markerstrokewidth=0.,ms=3,fontfamily="Computer Modern",legend=:outertopright,lc=:black,mc=:black,palette=palette(:roma10,rev=true))
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,2],label="gKM=0.1",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,3],label="gKM=0.2",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,4],label="gKM=0.3",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,6],label="gKM=0.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,11],label="gKM=1.",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2[:,16],label="gKM=1.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    