using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using PlotlyJS, CSV, DataFrames
using JLD

include("data_loader_outw.jl");

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



### KIR GRAD ###
    jld_name = "ChannelUpdate/cluster/dI_outw/data/jld/dI_KirCaLs_outw"

    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KirCaLs_outw = ones(length(gKir_KirCaLs_outw)).*NaN
        I1_KirCaLs_outw = ones(length(gKir_KirCaLs_outw)).*NaN
        I2_KirCaLs_outw = ones(length(gKir_KirCaLs_outw)).*NaN
        VSN_KirCaLs_outw = ones(length(gKir_KirCaLs_outw)).*NaN
        ISN_KirCaLs_outw = ones(length(gKir_KirCaLs_outw)).*NaN


        for i in 1:length(gKir_KirCaLs_outw)
            jld_name_i = string(jld_name,"_i_gKir_$(i).jld")
            println(jld_name_i)

            loaded_find_dI = load_dI(jld_name_i)[2]
            dI_KirCaLs_outw[i] = loaded_find_dI[1][1]
            I1_KirCaLs_outw[i] = loaded_find_dI[2][1]
            I2_KirCaLs_outw[i] = loaded_find_dI[5][1]
            VSN_KirCaLs_outw[i] = loaded_find_dI[4][1]
            ISN_KirCaLs_outw[i] = loaded_find_dI[3][1]
        end

        save_dI_mat(dI_KirCaLs_outw,I1_KirCaLs_outw,I2_KirCaLs_outw,VSN_KirCaLs_outw,string(jld_name,".jld"))
        plt = Plots.scatter(gKir_KirCaLs_outw,dI_KirCaLs_outw,mc=:purple,label=L"\Delta I")
        display(plt) 
    ==#

    #
        jld_name = "ChannelUpdate/cluster/dI_outw/data/jld/dI_KirCaLs_outw"
        dI_Kir_outw,I1_Kir_outw,I2_Kir_outw,VSN_Kir_outw = load_dI_mat(string(jld_name,".jld"))
    ==#
    

### KM GRAD ### 
    jld_name = "ChannelUpdate/cluster/dI_outw/data/jld/dI_KMCaLs_outw"

    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KMCaLs_outw = ones(length(gKM_KMCaLs_outw)).*NaN
        I1_KMCaLs_outw = ones(length(gKM_KMCaLs_outw)).*NaN
        I2_KMCaLs_outw = ones(length(gKM_KMCaLs_outw)).*NaN
        VSN_KMCaLs_outw = ones(length(gKM_KMCaLs_outw)).*NaN
        ISN_KMCaLs_outw = ones(length(gKM_KMCaLs_outw)).*NaN


        for i in 1:length(gKM_KMCaLs_outw)
            jld_name_i = string(jld_name,"_i_gKM_$(i).jld")
            println(jld_name_i)

            loaded_find_dI = load_dI(jld_name_i)[2]
            dI_KMCaLs_outw[i] = loaded_find_dI[1][1]
            I1_KMCaLs_outw[i] = loaded_find_dI[2][1]
            I2_KMCaLs_outw[i] = loaded_find_dI[5][1]
            VSN_KMCaLs_outw[i] = loaded_find_dI[4][1]
            ISN_KMCaLs_outw[i] = loaded_find_dI[3][1]
        end

        save_dI_mat(dI_KMCaLs_outw,I1_KMCaLs_outw,I2_KMCaLs_outw,VSN_KMCaLs_outw,string(jld_name,".jld"))
        plt = Plots.scatter(gKM_KMCaLs_outw,dI_KMCaLs_outw,mc=:purple,label=L"\Delta I")
        display(plt) 
    ==#

    #
        jld_name = "ChannelUpdate/cluster/dI_outw/data/jld/dI_KMCaLs_outw"
        dI_KM_outw,I1_KM_outw,I2_KM_outw,VSN_KM_outw = load_dI_mat(string(jld_name,".jld"))
    ==#


    Plots.plot(gKir_KirCaLs_outw,dI_Kir_inw,marker=:rect,mc=RGB(0.95,0.62,0.75),lw=2,markerstrokewidth=0.,lc=RGB(0.95,0.62,0.75),ms=3,label="Kir, inward",fontfamily="Computer Modern")
    Plots.plot!(gKir_KirCaLs_outw,dI_Kir_outw,marker=:rect,mc=RGB(0.95,0.62,0.75),ls=:dot,lw=3,markerstrokewidth=0.,lc=RGB(0.95,0.62,0.75),ms=3,label="Kir, outward",fontfamily="Computer Modern")
    Plots.plot!(gKM_KMCaLs_inw,dI_KM_inw,marker=:circle,mc=RGB(0.55,0.27,0.57),markerstrokewidth=0.,lc=RGB(0.55,0.27,0.57),ms=3,label="KM, inward",fontfamily="Computer Modern")
    Plots.plot!(gKM_KMCaLs_outw,dI_KM_outw,marker=:circle,mc=RGB(0.55,0.27,0.57),ls=:dot,markerstrokewidth=0.,lc=RGB(0.55,0.27,0.57),ms=3,label="KM, outward",fontfamily="Computer Modern")
    Plots.plot!(gKM_range_KMCaLs,dI_KMCaLs[41,:],label="pCaLs=1.5e-5, KM",marker=:rect,lw=2,markerstrokewidth=0.,ms=3,fontfamily="Computer Modern",legend=:outertopright,palette=palette(:Reds_9,rev=true))
    Plots.plot!(gKir_range_KirCaLs,dI_KirCaLs[41,:],label="pCaLs=1.5e-5, Kir",marker=:rect,lw=2,markerstrokewidth=0.,ms=3,fontfamily="Computer Modern",legend=:outertopright,palette=palette(:Reds_9,rev=true))


    Plots.plot(gKir_KirCaLs_inw,I1_Kir_inw,marker=:rect,mc=palette(:Paired_8)[2],lw=2,markerstrokewidth=0.,lc=palette(:Paired_8)[2],ms=3,label="I1, Kir, inward",fontfamily="Computer Modern")
    Plots.plot!(gKir_KirCaLs_inw,I2_Kir_inw,marker=:rect,mc=palette(:Paired_8)[1],lw=2,markerstrokewidth=0.,lc=palette(:Paired_8)[1],ms=3,label="I2, Kir, inward",fontfamily="Computer Modern")
    Plots.plot!(gKir_KirCaLs_outw,I1_Kir_outw,marker=:rect,mc=palette(:Paired_8)[4],ls=:dot,lw=3,markerstrokewidth=0.,lc=palette(:Paired_8)[4],ms=3,label="I1, Kir, outward",fontfamily="Computer Modern")
    Plots.plot!(gKir_KirCaLs_outw,I2_Kir_outw,marker=:rect,mc=palette(:Paired_8)[3],ls=:dot,lw=3,markerstrokewidth=0.,lc=palette(:Paired_8)[3],ms=3,label="I2, Kir, outward",fontfamily="Computer Modern")
    plot!(xlabel="g",ylabel="I", legend=:topleft,ylims=(-2.2,4.9),grid=:y)

    Plots.plot(gKM_KMCaLs_inw,I1_KM_inw,marker=:circle,mc=palette(:Paired_8)[6],markerstrokewidth=0.,lc=palette(:Paired_8)[6],ms=3,label="I1, KM, inward",fontfamily="Computer Modern")
    Plots.plot!(gKM_KMCaLs_inw,I2_KM_inw,marker=:circle,mc=palette(:Paired_8)[5],markerstrokewidth=0.,lc=palette(:Paired_8)[5],ms=3,label="I2, KM, inward",fontfamily="Computer Modern")
    Plots.plot!(gKM_KMCaLs_outw,I1_KM_outw,marker=:circle,mc=palette(:Paired_8)[8],ls=:dot,markerstrokewidth=0.,lc=palette(:Paired_8)[8],ms=3,label="I1, KM, outward",fontfamily="Computer Modern")
    Plots.plot!(gKM_KMCaLs_outw,I2_KM_outw,marker=:circle,mc=palette(:Paired_8)[7],ls=:dot,markerstrokewidth=0.,lc=palette(:Paired_8)[7],ms=3,label="I2, KM, outward",fontfamily="Computer Modern")
    plot!(xlabel="g",ylabel="I", legend=:topleft,ylims=(-2.2,4.9),grid=:y)


