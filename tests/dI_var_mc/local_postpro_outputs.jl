using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using PlotlyJS, CSV, DataFrames
using JLD

include("data_loader.jl");

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
    var_level =0.1
    jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KirCaLs"

    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KirCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        I1_KirCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        I2_KirCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        VSN_KirCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        ISN_KirCaLs_var_mc = ones(length(gNa_var_mat)).*NaN


        for i in 1:length(gNa_var_mat)
            jld_name_i = string(jld_name,"_i_C_gleak_gNa_$(i).jld")
            println(jld_name_i)

            loaded_find_dI = load_dI(jld_name_i)[2]
            dI_KirCaLs_var_mc[i] = loaded_find_dI[1][1]
            I1_KirCaLs_var_mc[i] = loaded_find_dI[2][1]
            I2_KirCaLs_var_mc[i] = loaded_find_dI[5][1]
            VSN_KirCaLs_var_mc[i] = loaded_find_dI[4][1]
            ISN_KirCaLs_var_mc[i] = loaded_find_dI[3][1]
        end

        save_dI_mat(dI_KirCaLs_var_mc,I1_KirCaLs_var_mc,I2_KirCaLs_var_mc,VSN_KirCaLs_var_mc,string(jld_name,".jld"))
        plt = Plots.scatter(C_var_mat,(dI_KirCaLs_var_mc.-dI_KirCaLs[21,4])./ dI_KirCaLs[21,4],ylims=(-0.3,0.3),mc=:purple,label=L"\Delta I")
        display(plt) 
        plt = Plots.scatter(gNa_var_mat,Veq_I1_KirCaLs_var_mc,mc=:purple,label="Veq(I1)")
        display(plt)
    ==#

    #
    var_level =0.1
        inputs_jld_name = "ChannelUpdate/cluster/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        C_mat_KirCaLs_var10,gleak_mat_KirCaLs_var10,gNa_mat_KirCaLs_var10 = load_inputs_of_script_dI_mc(inputs_jld_name)
        jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KirCaLs"
        dI_KirCaLs_var10,I1_KirCaLs_var10,I2_KirCaLs_var10,VSN_KirCaLs_var10 = load_dI_mat(string(jld_name,".jld"))
    var_level =0.2
        inputs_jld_name = "ChannelUpdate/cluster/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        C_mat_KirCaLs_var20,gleak_mat_KirCaLs_var20,gNa_mat_KirCaLs_var20 = load_inputs_of_script_dI_mc(inputs_jld_name)
        jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KirCaLs"
        dI_KirCaLs_var20,I1_KirCaLs_var20,I2_KirCaLs_var20,VSN_KirCaLs_var20 = load_dI_mat(string(jld_name,".jld"))
    var_level =0.3
        inputs_jld_name = "ChannelUpdate/cluster/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        C_mat_KirCaLs_var30,gleak_mat_KirCaLs_var30,gNa_mat_KirCaLs_var30 = load_inputs_of_script_dI_mc(inputs_jld_name)  
        jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KirCaLs"
        dI_KirCaLs_var30,I1_KirCaLs_var30,I2_KirCaLs_var30,VSN_KirCaLs_var30 = load_dI_mat(string(jld_name,".jld"))
    ==#

    

### KM GRAD ### 
    var_level =0.3
    jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KMCaLs"


    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KMCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        I1_KMCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        I2_KMCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        VSN_KMCaLs_var_mc = ones(length(gNa_var_mat)).*NaN
        ISN_KMCaLs_var_mc = ones(length(gNa_var_mat)).*NaN


        for i in 1:length(gNa_var_mat)
            jld_name_i = string(jld_name,"_i_C_gleak_gNa_$(i).jld")
            println(jld_name_i)

            loaded_find_dI = load_dI(jld_name_i)[2]
            dI_KMCaLs_var_mc[i] = loaded_find_dI[1][1]
            I1_KMCaLs_var_mc[i] = loaded_find_dI[2][1]
            I2_KMCaLs_var_mc[i] = loaded_find_dI[5][1]
            VSN_KMCaLs_var_mc[i] = loaded_find_dI[4][1]
            ISN_KMCaLs_var_mc[i] = loaded_find_dI[3][1]
        end

        save_dI_mat(dI_KMCaLs_var_mc,I1_KMCaLs_var_mc,I2_KMCaLs_var_mc,VSN_KMCaLs_var_mc,string(jld_name,".jld"))
        plt = Plots.scatter(C_var_mat,(dI_KMCaLs_var_mc.-dI_KMCaLs[21,4])./ dI_KMCaLs[21,4],ylims=(-0.3,0.3),mc=:purple,label=L"\Delta I")
        display(plt) 
        plt = Plots.scatter(gNa_var_mat,(dI_KMCaLs_var_mc.-dI_KMCaLs[21,4])./ dI_KMCaLs[21,4],mc=:purple,label="Veq(I1)")
        display(plt)
    ==#

    var_level =0.1
        inputs_jld_name = "ChannelUpdate/cluster/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        C_mat_KMCaLs_var10,gleak_mat_KMCaLs_var10,gNa_mat_KMCaLs_var10 = load_inputs_of_script_dI_mc(inputs_jld_name)
        jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KMCaLs"
        dI_KMCaLs_var10,I1_KMCaLs_var10,I2_KMCaLs_var10,VSN_KMCaLs_var10 = load_dI_mat(string(jld_name,".jld"))
    var_level =0.2
        inputs_jld_name = "ChannelUpdate/cluster/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        C_mat_KMCaLs_var20,gleak_mat_KMCaLs_var20,gNa_mat_KMCaLs_var20 = load_inputs_of_script_dI_mc(inputs_jld_name)
        jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KMCaLs"
        dI_KMCaLs_var20,I1_KMCaLs_var20,I2_KMCaLs_var20,VSN_KMCaLs_var20 = load_dI_mat(string(jld_name,".jld"))
    var_level =0.3
        inputs_jld_name = "ChannelUpdate/cluster/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        C_mat_KMCaLs_var30,gleak_mat_KMCaLs_var30,gNa_mat_KMCaLs_var30 = load_inputs_of_script_dI_mc(inputs_jld_name)  
        jld_name = "ChannelUpdate/cluster/dI_var_mc/data/jld/dI_var$(Int(var_level*100))_mc_KMCaLs"
        dI_KMCaLs_var30,I1_KMCaLs_var30,I2_KMCaLs_var30,VSN_KMCaLs_var30 = load_dI_mat(string(jld_name,".jld"))
