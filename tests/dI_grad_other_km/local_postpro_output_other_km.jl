using JLD

include("../dI_grad/local_postpro_output_larger_range.jl");
include("data_loader_larger_range.jl")

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


    jld_name = "ChannelUpdate/cluster/dI_grad_other_km/data/jld/dI_KMCaLs_other_Kv72"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        dI_KMCaLs = ones(length(pCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        I1_KMCaLs = ones(length(pCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        I2_KMCaLs = ones(length(pCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        VSN_KMCaLs = ones(length(pCaLs_range_KMCaLs),length(gKM_range_KMCaLs)).*NaN
        for i in 1:length(pCaLs_range_KMCaLs)
            for j in 1:length(gKM_range_KMCaLs)
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
        Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,dI_KMCaLs,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
        save_dI_mat(dI_KMCaLs,I1_KMCaLs,I2_KMCaLs,VSN_KMCaLs,string(jld_name,"_1e-2_1e-5.jld"))
    ==#
    
    dI_KMCaLs_other_KM,I1_KMCaLs_other_KM,I2_KMCaLs_other_KM,VSN_KMCaLs_other_KM = load_dI_mat(string(jld_name,"_1e-2_1e-5.jld"))
    plt = Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,dI_KMCaLs_other_KM,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(0,3.719))
    plt = Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,dI_KMCaLs,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(0,3.719))
  
    plt = Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,I1_KMCaLs_other_KM,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(-3.1,6.3))
    plt = Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,I1_KMCaLs,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256,clims=(-3.1,6.3))
  
    plt = Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,I2_KMCaLs_other_KM,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    plt = Plots.heatmap(gKM_range_KMCaLs,pCaLs_range_KMCaLs,I2_KMCaLs,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
  

    Plots.plot(fontfamily="Computer Modern",xlabel=L"p_{\mathrm{CaLs}}",ylabel=L"\Delta I ",title=L"V_{1/2,}_{\mathrm{KM}}="*latexify(26.7(u"mV"))*L"\ ; \ z_{\mathrm{KM}}="*latexify(12.6(u"mV")),ylims=(-0.1,3.1))
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,1],label=L"g_{\mathrm{KM}}=0",marker=:rect,lw=2,markerstrokewidth=0.,ms=3,fontfamily="Computer Modern",legend=:outertopright,lc=:black,mc=:black,palette=palette(:imola10,rev=false))
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,2],label=L"g_{\mathrm{KM}}=0.1",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,3],label=L"g_{\mathrm{KM}}=0.2",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,4],label=L"g_{\mathrm{KM}}=0.3",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,6],label=L"g_{\mathrm{KM}}=0.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,11],label=L"g_{\mathrm{KM}}=1.",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,16],label=L"g_{\mathrm{KM}}=1.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    
