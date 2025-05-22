using Statistics
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

#include("data_loader_larger_range_km.jl")
include("data_loader_larger_range_km_I1.jl")

function save_asc_sim(std,std_t,V,t,jld_name)
    save(jld_name, "std", std, "std_t", std_t,"V", V, "t", t)
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
function load_asc_sim(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        std = []
        std_t = []
        V = []
        t = []
        push!(std,load_array(jld,"std"))
        push!(std_t,load_array(jld,"std_t"))
        push!(V,load_array(jld,"V"))
        push!(t,load_array(jld,"t"))
    else
        std = NaN
        std_t = NaN
        V = NaN
        t = NaN
        println("Not found")
    end
    return std[1],std_t[1],V[1],t[1]
end

pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)
heaviside(t)=0*(t<0)+1*(t>=0)

function σ_membrane_t!(du,u,p,t)
   
    ## Noise parameters
    sigma_V,sigma_mNa,sigma_hNa,sigma_mKDR,sigma_mCaLf,sigma_hCaLf,sigma_mCaLs,sigma_hCaLs,sigma_mKir,sigma_mKM,sigma_Ca = p[3]

    ## Gradients
    du[1] = sigma_V(t)
    du[2] = sigma_mNa
    du[3] = sigma_hNa
    du[4] = sigma_mKDR
    du[5] = sigma_mCaLf
    du[6] = sigma_hCaLf
    du[7] = sigma_mCaLs
    du[8] = sigma_hCaLs
    du[9] = sigma_mKir
    du[10] = sigma_mKM
    du[11] = sigma_Ca

end
function sim_asc_std(p,ic,tspan,nspan,I1bist,I2bist)
     #println("---------------------I0=$(p[1]) & I1=$(It(t_pulse+1)+p[1])---------------------")
     p_sigma = p[3]
     dt_=0.01
     intermediate_state = []
     noisy_sol_cb_t_i = []
     noisy_sol_cb_V_i = []
     It_noise_i = []
     t_noise_i = []

     prob_SDE =  SDEProblem(membrane!,σ_membrane_t!,ic,tspan,p) 
     sol_SDE = solve(prob_SDE,dtmax=1,ImplicitEM(),save_noise=true,dt=dt_,adaptive=false)

     It_noise_ = p_var[1].(sol_SDE.W.t[2:end-1]) .+ p_sigma[1].(sol_SDE.W.t[2:end-1]) .*(sol_SDE.W[1,2:end-1]-sol_SDE.W[1,1:end-2])./(sol_SDE.W.t[2:end-1]-sol_SDE.W.t[1:end-2])  # It_noise = It(t) + σ*dW/dt
     t_noise_ = sol_SDE.W.t[2:end-1]
     push!(It_noise_i,It_noise_)
     push!(t_noise_i,t_noise_)

     push!(noisy_sol_cb_t_i, sol_SDE.t[1:end])
     push!(noisy_sol_cb_V_i,sol_SDE[1,1:end])


     push!(intermediate_state,sol_SDE[end])

     for i=1:nspan-1
        prob_SDE_inter =  SDEProblem(membrane!,σ_membrane_t!,intermediate_state[i],tspan.+i*tspan[2],[p_var,p_fixed,p_sigma])
        sol_SDE_inter = solve(prob_SDE_inter,ImplicitEM(),save_noise=true,dt=dt_,adaptive=false)

        It_noise_ = p_var[1].(sol_SDE_inter.W.t[2:end-1]) .+ p_sigma[1].(sol_SDE_inter.W.t[2:end-1]).*(sol_SDE_inter.W[1,2:end-1]-sol_SDE_inter.W[1,1:end-2])./(sol_SDE_inter.W.t[2:end-1]-sol_SDE_inter.W.t[1:end-2])  # It_noise = It(t) + σ*dW/dt
        t_noise_ = sol_SDE_inter.W.t[2:end-1]
        push!(It_noise_i,It_noise_)
        push!(t_noise_i,t_noise_)

        push!(noisy_sol_cb_t_i, sol_SDE_inter.t[1:end])
        push!(noisy_sol_cb_V_i,sol_SDE_inter[1,1:end])

        push!(intermediate_state,sol_SDE_inter[end])
     end

     noisy_sol_cb_t = reduce(vcat,noisy_sol_cb_t_i)
     noisy_sol_cb_V = reduce(vcat,noisy_sol_cb_V_i)
     It_noise = reduce(vcat,It_noise_i)
     t_noise = reduce(vcat,t_noise_i)
     
     std_It_noise = sigma_asc.(t_noise)./sqrt(dt_)
     println("------------------------------------------------------------")
 
     #==
     gr()
     plt_I = plot(t_noise,It_noise,label="Noise")
     plot!(t_noise,I1bist.*ones(size(t_noise)),lc=:red,label="I_1")
     plot!(t_noise,I2bist.*ones(size(t_noise)),lc=:red,label="I_2")
     #display(plt_I)
     
     plt_sigma = plot(t_noise,std_It_noise,label="std Noise")
     ==#
 
     #==
     plt_V = scatter(noisy_sol_cb_t,noisy_sol_cb_V,label="V")
     scatter!(noisy_sol_cb_t,zeros(length(noisy_sol_cb_t)),label="spikes")
     display(plt_V)
     display(plot(plt_sigma,plt_I,plt_V,layout=(3,1),fontfamily="Computer Modern"))
     ==#
 
     return noisy_sol_cb_t,noisy_sol_cb_V,t_noise,It_noise,std_It_noise
end

### INPUTS ###
    tspan=(0,10000) 
    nspan = 1

    #jld_name = "data/jld/asc_KMCaLs_10s_"
    jld_name = "data/jld/asc_KMCaLs_I1_10s_"
   
    #jld_name_f = "data/fig/asc_KMCaLs_10s_"
    jld_name_f = "data/fig/asc_KMCaLs_I1_10s_"

### TASK ###
## from stable state ##
    sigma_vec = 0.2:0.2:7
    sigma_asc(t) = sigma_vec[1] + (sigma_vec[end]-sigma_vec[1])*t/(tspan[2]*nspan)
    p_sigma = [sigma_asc ,0,0,0,0,0,0,0,0,0,0]
    sim_asc_std_km_st = sim_asc_std([p_var,p_fixed,p_sigma],ic_st,tspan,nspan,I1bist,I2bist)
    noisy_sol_cb_t_km_st,noisy_sol_cb_V_km_st,t_noise_km_st,It_noise_km_st,std_It_noise_km_st = sim_asc_std_km_st
    println("---------------------- Trying to plot the simunlation ----------------------")
    
    gr()
    plt_sigma = plot(t_noise_km_st[1:10:end],std_It_noise_km_st[1:10:end],label="std Noise")
    plt_V = plot(noisy_sol_cb_t_km_st[1:1:end],noisy_sol_cb_V_km_st[1:1:end],label="V")
    #scatter!(noisy_sol_cb_t_km_st[1:10:end],zeros(length(noisy_sol_cb_t_km_st[1:10:end])),label="spikes")

    plt_sub_s_I_V = plot(plt_sigma,plt_V,layout=(2,1),fontfamily="Computer Modern")
    #savefig(plt_sub_s_I_V,string("ChannelUpdate/cluster/dI_noise/",jld_name_f,"st.svg"))
    #save_asc_sim(std_It_noise_km_st,t_noise_km_st,noisy_sol_cb_V_km_st,noisy_sol_cb_t_km_st,string("ChannelUpdate/cluster/dI_noise/",jld_name,"st.jld"))

    std_It_noise_km_st,t_noise_km_st,noisy_sol_cb_V_km_st,noisy_sol_cb_t_km_st = load_asc_sim(string("ChannelUpdate/cluster/dI_noise/",jld_name,"st.jld"))
    plot(t_noise_km_st,std_It_noise_km_st)
    plot(noisy_sol_cb_t_km_st,noisy_sol_cb_V_km_st)

## from spiking state ##
    sigma_vec_cyc = sigma_vec #0.025:0.025:0.7 #28 
    sigma_asc_cyc(t) = sigma_vec_cyc[1] + (sigma_vec_cyc[end]-sigma_vec_cyc[1])*t/(tspan[2]*nspan)
    p_sigma_cyc = [sigma_asc_cyc ,0,0,0,0,0,0,0,0,0,0]
    sim_asc_std_km_cyc = sim_asc_std([p_var,p_fixed,p_sigma_cyc],ic_cyc,tspan,nspan,I1bist,I2bist)

    noisy_sol_cb_t_km_cyc,noisy_sol_cb_V_km_cyc,t_noise_km_cyc,It_noise_km_cyc,std_It_noise_km_cyc = sim_asc_std_km_cyc
    println("---------------------- Trying to plot the simunlation ----------------------")
    
    gr()   
    plt_sigma = plot(t_noise_km_cyc[1:10:end],std_It_noise_km_cyc[1:10:end],label="std Noise")
    plt_V = plot(noisy_sol_cb_t_km_cyc[1:1:end],noisy_sol_cb_V_km_cyc[1:1:end],label="V")
    #scatter!(noisy_sol_cb_t_km_cyc[1:10:end],zeros(length(noisy_sol_cb_t_km_cyc[1:10:end])),label="spikes")

    plt_sub_s_I_V = plot(plt_sigma,plt_V,layout=(2,1),fontfamily="Computer Modern")  
    #savefig(plt_sub_s_I_V,string("ChannelUpdate/cluster/dI_noise/",jld_name_f,"cyc.svg"))
    #save_asc_sim(std_It_noise_km_cyc,t_noise_km_cyc,noisy_sol_cb_V_km_cyc,noisy_sol_cb_t_km_cyc,string("ChannelUpdate/cluster/dI_noise/",jld_name,"cyc.jld"))

    std_It_noise_km_cyc,t_noise_km_cyc,noisy_sol_cb_V_km_cyc,noisy_sol_cb_t_km_cyc = load_asc_sim(string("ChannelUpdate/cluster/dI_noise/",jld_name,"cyc.jld"))
    plot(t_noise_km_cyc,std_It_noise_km_cyc)
    plot(noisy_sol_cb_t_km_cyc,noisy_sol_cb_V_km_cyc)

    