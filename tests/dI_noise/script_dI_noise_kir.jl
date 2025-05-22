using Statistics
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

include("data_loader_larger_range_kir.jl")
#include("data_loader_larger_range_kir_I1.jl")

function save_prop_switch(sigma_vec,prop_switch,I0,jld_name)
    save(jld_name, "sigma_vec", sigma_vec, "prop_switch", prop_switch, "I0", I0)
end
function save_howmanyswitch(nb_switch,switch_noise,std_It_noise,jld_name)
    save(jld_name, "nb_switch", nb_switch,"switch_noise", switch_noise, "std_It_noise", std_It_noise)
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


pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)
heaviside(t)=0*(t<0)+1*(t>=0)


function switch_sim_noisy(state_i,sol_cb_t,f_av,t_av,tspan,f_eq)
    println("________ Switch computation ________")
    T_end = maximum(tspan) .* 1/3
    f_av_th1 = 100
    f_av_th2 = f_av_th1/2
    f_av_th1_neighbours = 1
    if length(sol_cb_t)>1
        sol_cb_dt = sol_cb_t[2:end] .-sol_cb_t[1:end-1]
        sol_cb_f = 1000 ./sol_cb_dt
        mean__sol_cb_f__T_end = mean(sol_cb_f[sol_cb_t[2:end].>=T_end])
        println("Mean f : $(mean__sol_cb_f__T_end)  ;  f_th = $(f_eq .* 2/3)")
    else
        sol_cb_dt = []
    end

    if state_i=="stable"
        if length(sol_cb_t)>1 
            #find points with f_av>100
            #verify if neighbours of these points have f_av>50    ;    if so, switch =1
            ind_f_av_above_th1 = collect(1:length(f_av))[f_av.>=f_av_th1]
            burst_cnt = 0
            if length(ind_f_av_above_th1)>0
                for i_f_av in ind_f_av_above_th1
                    for n in vcat(-f_av_th1_neighbours:-1 , 1:f_av_th1_neighbours)
                        if i_f_av+n<= length(f_av) && i_f_av+n>= 1 && f_av[i_f_av+n] > f_av_th2
                            #increment counter of burst for i_f_av, and maybe break the loop 
                            burst_cnt=burst_cnt+1
                            println("Burst validated at $(t_av[i_f_av]) ms, with $(t_av[i_f_av+n]) above thresh ")
                        end
                    end
                end
                if burst_cnt>0
                     #has switched to cycle
                    switch=1
                    t_switch = sol_cb_t[2]
                else
                    #has not switch to cycle
                    switch=0
                    t_switch = NaN
                end
            else
                #no burst = no switch 
                 #has not switch to cycle
                 switch=0
                 t_switch = NaN
            end
        else
            #has not switch to cycle
            switch=0
            t_switch = NaN
        end
    else
        if length(sol_cb_t)>1 
            #find points with f_av<50
            #verify if neighbours of these points have f_av<50    ;    if so, switch =1
            ind_f_av_under_th2 = collect(1:length(f_av))[f_av.<=f_av_th2]
            quiescent_cnt = 0
            if length(ind_f_av_under_th2)>0
                for i_f_av in ind_f_av_under_th2
                    for n in vcat(-f_av_th1_neighbours:-1 , 1:f_av_th1_neighbours)
                        if i_f_av+n<= length(f_av) && i_f_av+n>= 1 && f_av[i_f_av+n] <= f_av_th2
                            #increment counter of burst for i_f_av, and maybe break the loop 
                            quiescent_cnt=quiescent_cnt+1
                            println("Quasi quiescent window validated at $(t_av[i_f_av]) ms, with $(t_av[i_f_av+n]) under thresh ")
                        end
                    end
                end
                if quiescent_cnt>0
                     #has switched to stable
                    switch=1
                    t_switch = sol_cb_t[2]
                else
                    #has not switch to stable
                    switch=0
                    t_switch = NaN
                end
            else
                #no quiescent period = no switch 
                 #has not switch to stable
                 switch=0
                 t_switch = NaN
            end
        else
            #has not switch to cycle
            switch=0
            t_switch = NaN
        end

    end
    return switch,t_switch
end
function sim_noisy_pulse_cb(p,ic,state_i,tspan,nspan,I1bist,I2bist,f_eq)
    # println("---------------------I0=$(p[1]) & I1=$(It(t_pulse+1)+p[1])---------------------")
    p_sigma = p[3]
    dt_=0.01
    intermediate_state = []
    noisy_sol_cb_t_i = []
    noisy_sol_cb_V_i = []
    It_noise_i = []
    t_noise_i = []

    prob_SDE =  SDEProblem(membrane!,σ_membrane!,ic,tspan,p) 
    sol_SDE_cb = solve(prob_SDE,ImplicitEM(),callback=cb,save_everystep=false,save_start=false,save_end=true,save_noise=true,dt=dt_,adaptive=false)

    It_noise_ = p_var[1].(sol_SDE_cb.W.t[2:end-1]) .+ p_sigma[1].*(sol_SDE_cb.W[1,2:end-1]-sol_SDE_cb.W[1,1:end-2])./(sol_SDE_cb.W.t[2:end-1]-sol_SDE_cb.W.t[1:end-2])  # It_noise = It(t) + σ*dW/dt
    t_noise_ = sol_SDE_cb.W.t[2:end-1]
    push!(It_noise_i,It_noise_)
    push!(t_noise_i,t_noise_)

    if length(sol_SDE_cb.t)>1    
        push!(noisy_sol_cb_t_i, sol_SDE_cb.t[1:end-1]) #removes end point (save_end=true)
        push!(noisy_sol_cb_V_i,sol_SDE_cb[1,1:end-1])
    else
        push!(noisy_sol_cb_t_i,[]) 
        push!(noisy_sol_cb_V_i,[])
    end
     push!(intermediate_state,sol_SDE_cb[end])

     for i=1:nspan-1
        prob_SDE_inter =  SDEProblem(membrane!,σ_membrane!,intermediate_state[i],tspan.+i*tspan[2],[p_var,p_fixed,p_sigma])
        sol_SDE_inter_cb = solve(prob_SDE_inter,ImplicitEM(),callback=cb,save_everystep=false,save_start=false,save_end=true,save_noise=true,dt=dt_,adaptive=false)
        println("---Iter : $(i)")
        println("---Data points : $(length(sol_SDE_inter_cb.t))")
        println("---V points : $(length(sol_SDE_inter_cb[1,:]))")

        It_noise_ = p_var[1].(sol_SDE_inter_cb.W.t[2:end-1]) .+ p_sigma[1].*(sol_SDE_inter_cb.W[1,2:end-1]-sol_SDE_inter_cb.W[1,1:end-2])./(sol_SDE_inter_cb.W.t[2:end-1]-sol_SDE_inter_cb.W.t[1:end-2])  # It_noise = It(t) + σ*dW/dt
        t_noise_ = sol_SDE_inter_cb.W.t[2:end-1]
        push!(It_noise_i,It_noise_)
        push!(t_noise_i,t_noise_)

        if length(sol_SDE_inter_cb.t)>1    
            push!(noisy_sol_cb_t_i,sol_SDE_inter_cb.t[1:end-1]) #removes end point (save_end=true)
            push!(noisy_sol_cb_V_i,sol_SDE_inter_cb[1,1:end-1])
        else
            push!(noisy_sol_cb_t_i,[]) 
            push!(noisy_sol_cb_V_i,[])
        end
        push!(intermediate_state,sol_SDE_inter_cb[end])
     end

    noisy_sol_cb_t = reduce(vcat,noisy_sol_cb_t_i)
    noisy_sol_cb_V = reduce(vcat,noisy_sol_cb_V_i)
    It_noise = reduce(vcat,It_noise_i)
    t_noise = reduce(vcat,t_noise_i)
     
    std_It_noise = std(It_noise)

    T_av = 100
    f_av = []
    t_f_av = []
    f = 1000 ./(noisy_sol_cb_t[2:end]-noisy_sol_cb_t[1:end-1])
    for z_av=0:Int(tspan[2]*nspan/T_av)-1
        ind_av_win_ = collect(1:length(noisy_sol_cb_t[1:end]))[noisy_sol_cb_t[1:end].<=T_av+z_av*T_av]
        ind_av_win = ind_av_win_[noisy_sol_cb_t[ind_av_win_].>=z_av*T_av]
        #println("Number of spikes from $(z_av*T_av) to $(T_av+z_av*T_av) ms : $(length(ind_av_win))")
        if length(ind_av_win)>=1
            push!(f_av, 1000*length(ind_av_win)/T_av)
        else
            push!(f_av,0)   
        end
        push!(t_f_av, T_av/2 + z_av *T_av)
    end
    
    switch,t_switch = switch_sim_noisy(state_i,noisy_sol_cb_t,f_av,t_f_av,(tspan[1],tspan[2]*nspan),f_eq)
    println("------------------------------------------------------------")

    #==

    gr()
    plt_I = plot(t_noise,It_noise,label="Noise")
    plot!(t_noise,I1bist.*ones(size(t_noise)),lc=:red,label="I_1")
    plot!(t_noise,I2bist.*ones(size(t_noise)),lc=:red,label="I_2")
    #plot!(t_switch.*ones(2),[I1bist-0.05,I2bist+0.05],lc=RGB(0.6,0.6,0.6),label="Switch",lw=4)
    #display(plt_I)

    plt_V = scatter(noisy_sol_cb_t,noisy_sol_cb_V,label="V")
    scatter!(noisy_sol_cb_t,zeros(length(noisy_sol_cb_t)),label="spikes",ylims=(-70,50),xlims=(tspan[1],tspan[2]*nspan))
    #plot!(t_switch.*ones(2),[-70,-50],lc=RGB(0.6,0.6,0.6),label="Switch",lw=4)
    #display(plt_V)
    
    plt_f = plot(noisy_sol_cb_t[2:end],f,label="f",markershape=:diamond,xlims=(tspan[1],tspan[2]*nspan))

    plt_f_av = bar(t_f_av,f_av,ylims=(0,150),xlims=(tspan[1],tspan[2]*nspan))

    plt_all = plot(plt_V,plt_f,plt_f_av,layout=(3,1),fontfamily="Computer Modern",title="Switch = $(switch), sigma=$(p[3][1]), std = $(std_It_noise)",size=(400,800))
    savefig(plt_all,"test_kir_st_switch_32.svg")

    for z=0:4
        display(plot(plt_all,xlims=(0,nspan*tspan[2]/5).+z*nspan*tspan[2]/5,xlabel="part $(z+1)",size=(400,800)))
    end
    ==#

 
    return noisy_sol_cb_t,noisy_sol_cb_V,t_noise,It_noise,std_It_noise,switch
end
function howmanyswitchwithstd(p,I1bist,I2bist,f_eq,ic,state_i,tspan,nspan,n_samples)

    noisy_sol_t = []
    noisy_sol_V = []
    t_noise = []
    It_noise = []
    std_It_noise = []
    switch_noise = []
    for i=1:n_samples
        println("----------- Noisy simulation $i over $(n_samples) -----------")
        noisy_sol_t_i,noisy_sol_V_i,t_noise_i,It_noise_i,std_It_noise_i,switch_noise_i = sim_noisy_pulse_cb(p,ic,state_i,tspan,nspan,I1bist,I2bist,f_eq)
        
        push!(noisy_sol_t,noisy_sol_t_i)
        push!(noisy_sol_V,noisy_sol_V_i)
        push!(t_noise,t_noise_i)
        push!(It_noise,It_noise_i)
        push!(std_It_noise,std_It_noise_i)
        push!(switch_noise,switch_noise_i)
        println("            Switch = $(switch_noise_i)            ")
    end
    nb_switch = sum(switch_noise)
    return nb_switch,noisy_sol_t,noisy_sol_V,t_noise,It_noise,std_It_noise,switch_noise
end
function switchwithstd_wsave(p,sigma_vec,I1bist,I2bist,f_eq,ic,state_i,tspan,nspan,n_samples,jld_name)
    p_var = p[1]
    p_fixed = p[2]

    prop_switch = []
    noisy_sol_t = []
    noisy_sol_V = []
    t_noise = []
    It_noise = []
    std_It_noise = []
    switch_noise = []
    for sigma in sigma_vec
        println("--------------------------------- Sigma = $sigma ------------------------------")
        p_sigma = [sigma,0,0,0,0,0,0,0,0,0,0]

        nb_switch_sigma,noisy_sol_t_sigma,noisy_sol_V_sigma,t_noise_sigma,It_noise_sigma,std_It_noise_sigma,switch_noise_sigma = howmanyswitchwithstd([p_var,p_fixed,p_sigma],I1bist,I2bist,f_eq,ic,state_i,tspan,nspan,n_samples)
        push!(prop_switch,nb_switch_sigma/n_samples)
        push!(noisy_sol_t,noisy_sol_t_sigma)
        push!(noisy_sol_V,noisy_sol_V_sigma)
        push!(t_noise,t_noise_sigma)
        push!(It_noise,It_noise_sigma)
        push!(std_It_noise,std_It_noise_sigma)
        push!(switch_noise,switch_noise_sigma)
        
        save_howmanyswitch(nb_switch_sigma,switch_noise_sigma,std_It_noise,jld_name)
        println("Total switch = $nb_switch_sigma")
    end

    return prop_switch
end


#
### SCRIPT INPUTS ###
    i_sigma_vec = parse(Int64, ARGS[1]) #1:length(sigma_vec)
    i_batch_N_samples =parse(Int64, ARGS[2]) #1:20 (for 20*25 n_samples)

### INPUTS ###
    n_samples = 25

    sigma_vec = 1.2:0.1:3.5 #24
    #sigma_vec =2:0.2:7 #26
    tspan=(0,2000) 
    nspan = 25

    jld_name = "data/jld/switch_KirCaLs_N500_n25_"
    #jld_name = "data/jld/switch_KirCaLs_I1_N500_n25_"
    #jld_name = "data/jld/switch_KirCaLs_cyc_N500_n25_"
    #jld_name = "data/jld/switch_KirCaLs_cyc_I1_N500_n25_"

    jld_name_i = string(jld_name,"sigma_$(i_sigma_vec)_batch_$(i_batch_N_samples).jld")

### TASK ###
    sigma_vec_ = [sigma_vec[i_sigma_vec]]
    println("---------------------- I0=$(p_var[1](0)) & sigma_vec=$(sigma_vec_) -----------------------")
    println("...................... with I1=$(I1bist) & I2=$(I2bist) ......................")
    switchwithsigma_noise_i = switchwithstd_wsave([p_var,p_fixed],sigma_vec_,I1bist,I2bist,f_eq,ic_st,"stable",tspan,nspan,n_samples,jld_name_i)
    println("---------------------- prop_switch=$(switchwithsigma_noise_i) -----------------------")
==#   