using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

function save_prop_switch(sigma_vec,prop_switch,I0,jld_name)
    save(jld_name, "sigma_vec", sigma_vec, "prop_switch", prop_switch, "I0", I0)
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
function load_prop_switch(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        sigma_vec = []
        prop_switch = []
        I0 = []
        push!(sigma_vec,load_array(jld,"sigma_vec"))
        push!(prop_switch,load_array(jld,"prop_switch"))
        push!(I0,load_array(jld,"I0"))
    else
        sigma_vec = NaN
        prop_switch = NaN
        I0 = NaN
        println("Not found")
    end
    return sigma_vec[1],prop_switch[1],I0[1]
end
function load_howmanyswitch(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        nb_switch = []
        #noisy_sol_t = []
        #noisy_sol_V = []
        #t_noise = []
        #It_noise = []
        switch_noise = []
        std_It_noise = []
        push!(nb_switch,load_array(jld,"nb_switch"))
        #push!(noisy_sol_t,load_array(jld,"noisy_sol_t"))
        #push!(noisy_sol_V,load_array(jld,"noisy_sol_V"))
        #push!(t_noise,load_array(jld,"t_noise"))
        #push!(It_noise,load_array(jld,"It_noise"))
        push!(switch_noise,load_array(jld,"switch_noise"))
        push!(std_It_noise,load_array(jld,"std_It_noise"))
    else
        nb_switch = NaN
        #noisy_sol_t = NaN
        #noisy_sol_V = NaN
        #t_noise = NaN
        #It_noise = NaN
        switch_noise = NaN
        std_It_noise = NaN
        println("Not found")
    end
    return nb_switch[1],switch_noise[1],std_It_noise[1]#nb_switch[1],noisy_sol_t[1],noisy_sol_V[1],t_noise[1],It_noise[1],switch_noise[1]
end
function load_prop_switch(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        sigma_vec = []
        prop_switch = []
        I0 = []
        push!(sigma_vec,load_array(jld,"sigma_vec"))
        push!(prop_switch,load_array(jld,"prop_switch"))
        push!(I0,load_array(jld,"I0"))
    else
        sigma_vec = NaN
        prop_switch = NaN
        I0 = NaN
        println("Not found")
    end
    return sigma_vec[1],prop_switch[1],I0[1]
end

dt = 0.01


## Kir ; from stable state ##
#==
# This block load the results from the separated jld files and combine them in a single jld 
    include("data_loader_larger_range_kir.jl")

    #sigma_vec = 1:0.1:3
    #jld_name = "ChannelUpdate/cluster/dI_noise/data/jld/multi_switch_KirCaLs_N500_n25_"
    sigma_vec = 1.2:0.1:3.5 #24
    jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KirCaLs_N500_n25_"

    n_samples = 25
    nb_switch_std = []
    switch_noise_std = []

    for i_sigma_vec in 1:length(sigma_vec)
        nb_switch_std_i = 0
        switch_noise_std_i = []
        for i_batch_N_samples in 1:20
            jld_name_i = string(jld_name,"sigma_$(i_sigma_vec)_batch_$(i_batch_N_samples).jld")
            println(jld_name_i)
            nb_switch,switch_noise = load_howmanyswitch(jld_name_i)
            nb_switch_std_i = nb_switch_std_i+nb_switch
            switch_noise_std_i = vcat(switch_noise_std_i,switch_noise)
        end
        println(nb_switch_std_i)
        println(switch_noise_std_i)
        push!(nb_switch_std,nb_switch_std_i)
        push!(switch_noise_std_i,switch_noise_std_i)
    end

    prop_switch_std = nb_switch_std ./ (20*n_samples)
    plot!(sigma_vec./sqrt(dt),prop_switch_std,lc=:purple,lw=2,marker=:circle,mc=:purple,ms=5,markerstrokewidth=0,label=:none)

    save_prop_switch(sigma_vec,prop_switch_std,I0,string(jld_name,".jld"))
==#

## Kir ; from spiking ##
#==
# This block load the results from the separated jld files and combine them in a single jld 
    include("data_loader_larger_range_kir.jl")

    #sigma_vec = 1:0.2:3
    #jld_name = "ChannelUpdate/cluster/dI_noise/data/jld/multi_switch_KirCaLs_cyc_N500_n25_"
    sigma_vec = 2:0.2:7 #26
    jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KirCaLs_cyc_N500_n25_"

    n_samples = 25
    nb_switch_std = []
    switch_noise_std = []

    for i_sigma_vec in 1:length(sigma_vec)
        nb_switch_std_i = 0
        switch_noise_std_i = []
        for i_batch_N_samples in 1:20
            jld_name_i = string(jld_name,"sigma_$(i_sigma_vec)_batch_$(i_batch_N_samples).jld")
            println(jld_name_i)
            nb_switch,switch_noise = load_howmanyswitch(jld_name_i)
            nb_switch_std_i = nb_switch_std_i+nb_switch
            switch_noise_std_i = vcat(switch_noise_std_i,switch_noise)
        end
        println(nb_switch_std_i)
        println(switch_noise_std_i)
        push!(nb_switch_std,nb_switch_std_i)
        push!(switch_noise_std_i,switch_noise_std_i)
    end

    prop_switch_std = nb_switch_std ./ (20*n_samples)
    save_prop_switch(sigma_vec,prop_switch_std,I0,string(jld_name,".jld"))
==#

## KM ; from stable state ##
#==
# This block load the results from the separated jld files and combine them in a single jld 
    include("data_loader_larger_range_km.jl")

    sigma_vec = 0.1:0.1:2
    jld_name = "ChannelUpdate/cluster/dI_noise/data/jld/multi_switch_KMCaLs_N500_n25_"
    sigma_vec = 0.4:0.05:1.4 #21
    jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KMCaLs_N500_n25_"

    n_samples = 25
    nb_switch_std = []
    switch_noise_std = []

    for i_sigma_vec in 1:length(sigma_vec)
        nb_switch_std_i = 0
        switch_noise_std_i = []
        for i_batch_N_samples in 1:20
            jld_name_i = string(jld_name,"sigma_$(i_sigma_vec)_batch_$(i_batch_N_samples).jld")
            println(jld_name_i)
            nb_switch,switch_noise = load_howmanyswitch(jld_name_i)
            nb_switch_std_i = nb_switch_std_i+nb_switch
            switch_noise_std_i = vcat(switch_noise_std_i,switch_noise)
        end
        println(nb_switch_std_i)
        println(switch_noise_std_i)
        push!(nb_switch_std,nb_switch_std_i)
        push!(switch_noise_std_i,switch_noise_std_i)
    end

    prop_switch_std = nb_switch_std ./ (20*n_samples)
    plot!(sigma_vec./sqrt(dt),prop_switch_std,lc=:black,lw=2,marker=:circle,mc=:black,ms=5,markerstrokewidth=0,label=:none)

    save_prop_switch(sigma_vec,prop_switch_std,I0,string(jld_name,".jld"))
==#

## KM ; from spiking ##
#==
# This block load the results from the separated jld files and combine them in a single jld 
    include("data_loader_larger_range_km.jl")

    sigma_vec = 0.1:0.2:3
    jld_name = "ChannelUpdate/cluster/dI_noise/data/jld/multi_switch_KMCaLs_cyc_N500_n25_"
    sigma_vec = 2:0.2:7 #26
    jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KMCaLs_cyc_N500_n25_"

    n_samples = 25
    nb_switch_std = []
    switch_noise_std = []

    for i_sigma_vec in 1:length(sigma_vec)
        nb_switch_std_i = 0
        switch_noise_std_i = []
        for i_batch_N_samples in 1:20
            jld_name_i = string(jld_name,"sigma_$(i_sigma_vec)_batch_$(i_batch_N_samples).jld")
            println(jld_name_i)
            nb_switch,switch_noise = load_howmanyswitch(jld_name_i)
            nb_switch_std_i = nb_switch_std_i+nb_switch
            switch_noise_std_i = vcat(switch_noise_std_i,switch_noise)
        end
        println(nb_switch_std_i)
        println(switch_noise_std_i)
        push!(nb_switch_std,nb_switch_std_i)
        push!(switch_noise_std_i,switch_noise_std_i)
    end

    prop_switch_std = nb_switch_std ./ (20*n_samples)
    save_prop_switch(sigma_vec,prop_switch_std,I0,string(jld_name,".jld"))
==#

jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KirCaLs_N500_n25_"
std_vec_Kir,prop_switch_std_Kir,I0_Kir =load_prop_switch(string(jld_name,".jld"))

jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KirCaLs_cyc_N500_n25_"
std_vec_Kir_cyc,prop_switch_std_Kir_cyc,I0_Kir_cyc =load_prop_switch(string(jld_name,".jld"))

jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KMCaLs_N500_n25_"
std_vec_KM,prop_switch_std_KM,I0_KM =load_prop_switch(string(jld_name,".jld"))

jld_name = "ChannelUpdate/cluster/dI_noise/data/jld_burst_both_Ihalf2/switch_KMCaLs_cyc_N500_n25_"
std_vec_KM_cyc,prop_switch_std_KM_cyc,I0_KM_cyc =load_prop_switch(string(jld_name,".jld"))

Plots.plot(fontfamily="Computer Modern")
plot!(std_vec_Kir./sqrt(dt),prop_switch_std_Kir,lc=RGB(0.95,0.62,0.75),lw=2,marker=:circle,mc=RGB(0.95,0.62,0.75),ms=5,markerstrokewidth=0,label=:none)
plot!(std_vec_KM./sqrt(dt),prop_switch_std_KM,lc=:purple,lw=2,marker=:circle,mc=:purple,ms=5,markerstrokewidth=0,label=:none)
plot!(std_vec_Kir_cyc./sqrt(dt),prop_switch_std_Kir_cyc,lc=RGB(0.95,0.62,0.75),lw=2,marker=:circle,mc=RGB(0.95,0.62,0.75),ms=5,markerstrokewidth=0,label=:none,ls=:dot)
plot!(std_vec_KM_cyc./sqrt(dt),prop_switch_std_KM_cyc,lc=:purple,lw=2,marker=:circle,mc=:purple,ms=5,markerstrokewidth=0,label=:none,ls=:dot)

Plots.plot(fontfamily="Computer Modern")
plot!(std_vec_Kir ,prop_switch_std_Kir,lc=RGB(0.95,0.62,0.75),lw=2,marker=:circle,mc=RGB(0.95,0.62,0.75),ms=5,markerstrokewidth=0,label=:none)
plot!(std_vec_KM ,prop_switch_std_KM,lc=:purple,lw=2,marker=:circle,mc=:purple,ms=5,markerstrokewidth=0,label=:none)
plot!(std_vec_Kir_cyc ,prop_switch_std_Kir_cyc,lc=RGB(0.95,0.62,0.75),lw=2,marker=:circle,mc=RGB(0.95,0.62,0.75),ms=5,markerstrokewidth=0,label=:none)
plot!(std_vec_KM_cyc ,prop_switch_std_KM_cyc,lc=:purple,lw=2,marker=:circle,mc=:purple,ms=5,markerstrokewidth=0,label=:none)



