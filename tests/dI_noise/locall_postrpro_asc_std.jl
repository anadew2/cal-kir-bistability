using Statistics
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")


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


### INPUTS ###

    jld_name_kir_Ihalf = "data/asc_sim/asc_KirCaLs_Ihalf_10s_"
    jld_name_kir_I1 = "data/asc_sim/asc_KirCaLs_I1_10s_"

    jld_name_f_kir_Ihalf = "data/fig/asc_KirCaLs_Ihalf_10s_"
    jld_name_f_kir_I1 = "data/fig/asc_KirCaLs_I1_10s_"

    jld_name_km_Ihalf = "data/jld/asc_KMCaLs_Ihalf_10s_"
    jld_name_km_I1 = "data/jld/asc_KMCaLs_I1_10s_"
   
    jld_name_f_km_Ihalf = "data/fig/asc_KMCaLs_Ihalf_10s_"
    jld_name_f_km_I1 = "data/fig/asc_KMCaLs_I1_10s_"

### TASK ###
## Kir from stable state ##

    std_It_noise_kir_Ihalf_st,t_noise_kir_Ihalf_st,noisy_sol_cb_V_kir_Ihalf_st,noisy_sol_cb_t_kir_Ihalf_st = load_asc_sim(string("tests/dI_noise/",jld_name_kir_Ihalf,"st.jld"))
    plt_sigma_ = Plots.plot(t_noise_kir_Ihalf_st,std_It_noise_kir_Ihalf_st,label="std Noise")
    plt_V_ = Plots.plot(noisy_sol_cb_t_kir_Ihalf_st,noisy_sol_cb_V_kir_Ihalf_st,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")

## Kir from spiking state ##
    
    std_It_noise_kir_Ihalf_cyc,t_noise_kir_Ihalf_cyc,noisy_sol_cb_V_kir_Ihalf_cyc,noisy_sol_cb_t_kir_Ihalf_cyc = load_asc_sim(string("tests/dI_noise/",jld_name_kir_Ihalf,"cyc.jld"))
    plt_sigma_ = Plots.plot(t_noise_kir_Ihalf_cyc,std_It_noise_kir_Ihalf_cyc,label="std Noise",title="Kir")
    plt_V_ = Plots.plot(noisy_sol_cb_t_kir_Ihalf_cyc,noisy_sol_cb_V_kir_Ihalf_cyc,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")

## KM from stable state ##

    std_It_noise_km_Ihalf_st,t_noise_km_Ihalf_st,noisy_sol_cb_V_km_Ihalf_st,noisy_sol_cb_t_km_Ihalf_st = load_asc_sim(string("tests/dI_noise/",jld_name_km_Ihalf,"st.jld"))
    plt_sigma_ = Plots.plot(t_noise_km_Ihalf_st,std_It_noise_km_Ihalf_st,label="std Noise")
    plt_V_ = Plots.plot(noisy_sol_cb_t_km_Ihalf_st,noisy_sol_cb_V_km_Ihalf_st,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")

## KM from spiking state ##

    std_It_noise_km_Ihalf_cyc,t_noise_km_Ihalf_cyc,noisy_sol_cb_V_km_Ihalf_cyc,noisy_sol_cb_t_km_Ihalf_cyc = load_asc_sim(string("tests/dI_noise/",jld_name_km_Ihalf,"cyc.jld"))
    plt_sigma_ = Plots.plot(t_noise_km_Ihalf_cyc,std_It_noise_km_Ihalf_cyc,label="std Noise",title="KM")
    plt_V_ = Plots.plot(noisy_sol_cb_t_km_Ihalf_cyc,noisy_sol_cb_V_km_Ihalf_cyc,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")


### TASK ###
## Kir from stable state ##

    std_It_noise_kir_I1_st,t_noise_kir_I1_st,noisy_sol_cb_V_kir_I1_st,noisy_sol_cb_t_kir_I1_st = load_asc_sim(string("tests/dI_noise/",jld_name_kir_I1,"st.jld"))
    plt_sigma_ = Plots.plot(t_noise_kir_I1_st,std_It_noise_kir_I1_st,label="std Noise")
    plt_V_ = Plots.plot(noisy_sol_cb_t_kir_I1_st,noisy_sol_cb_V_kir_I1_st,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")

## Kir from spiking state ##

    std_It_noise_kir_I1_cyc,t_noise_kir_I1_cyc,noisy_sol_cb_V_kir_I1_cyc,noisy_sol_cb_t_kir_I1_cyc = load_asc_sim(string("tests/dI_noise/",jld_name_kir_I1,"cyc.jld"))
    plt_sigma_ = Plots.plot(t_noise_kir_I1_cyc,std_It_noise_kir_I1_cyc,label="std Noise",title="Kir")
    plt_V_ = Plots.plot(noisy_sol_cb_t_kir_I1_cyc,noisy_sol_cb_V_kir_I1_cyc,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")

## KM from stable state ##

    std_It_noise_km_I1_st,t_noise_km_I1_st,noisy_sol_cb_V_km_I1_st,noisy_sol_cb_t_km_I1_st = load_asc_sim(string("tests/dI_noise/",jld_name_km_I1,"st.jld"))
    plt_sigma_ = Plots.plot(t_noise_km_I1_st,std_It_noise_km_I1_st,label="std Noise")
    plt_V_ = Plots.plot(noisy_sol_cb_t_km_I1_st,noisy_sol_cb_V_km_I1_st,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")

## KM from spiking state ##

    std_It_noise_km_I1_cyc,t_noise_km_I1_cyc,noisy_sol_cb_V_km_I1_cyc,noisy_sol_cb_t_km_I1_cyc = load_asc_sim(string("tests/dI_noise/",jld_name_km_I1,"cyc.jld"))
    plt_sigma_ = Plots.plot(t_noise_km_I1_cyc,std_It_noise_km_I1_cyc,label="std Noise",title="KM")
    plt_V_ = Plots.plot(noisy_sol_cb_t_km_I1_cyc,noisy_sol_cb_V_km_I1_cyc,label="V")
    Plots.plot(plt_sigma_,plt_V_,layout=(2,1),fontfamily="Computer Modern")



    