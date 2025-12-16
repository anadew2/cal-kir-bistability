
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../solver_fun.jl")

include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("data_loader_ext.jl")
include("jld_hdle_dI_grad.jl")
include("fct_dI_grad.jl")


### SCRIPT INPUTS ###
    i_pCaLs_d_range_KirCaLs= parse(Int64, ARGS[1])           #1:length(pCaLs_d_range_KirCaLs)
    i_batch_gKir_s_range_KirCaLs = parse(Int64, ARGS[2])      #1:n_batch_gKir_s_range_KirCaLs:length(gKir_s_range_KirCaLs)

### INPUTS ###
    n_batch_gKir_s_range_KirCaLs = 2

    V_ic = -90. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic,
        V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic]

    jld_name = "data/jld_J_I1down/dI_Kir_s_CaLs_d_compart_ext_"

### TASK ###
    batch_gKir_s_range_KirCaLs = i_batch_gKir_s_range_KirCaLs:(minimum([i_batch_gKir_s_range_KirCaLs+n_batch_gKir_s_range_KirCaLs-1,length(gKir_s_range_KirCaLs)]))
    for i_gKir_s_range_KirCaLs in batch_gKir_s_range_KirCaLs
        jld_name_i = string(jld_name,"i_pCaLs_d_$(i_pCaLs_d_range_KirCaLs)_i_gKir_s_$(i_gKir_s_range_KirCaLs).jld")
        println("------------ Simulation for i_pCaLs_d = $(i_pCaLs_d_range_KirCaLs) & i_gKir_s = $(i_gKir_s_range_KirCaLs) ------------")
        
        pCaLs_d_i_batch = pCaLs_d_range_KirCaLs[i_pCaLs_d_range_KirCaLs] 
        gKir_s_i_batch = gKir_s_range_KirCaLs[i_gKir_s_range_KirCaLs]

        p_var_d_i_batch = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
        p_var_s_i_batch = [gNa_s,gKDR_s,gKir_s_i_batch,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]

        find_dI_ = find_dI_compart([p_fixed,p_var_d_i_batch,p_var_s_i_batch],ic,[IX_KirCaLs[i_pCaLs_d_range_KirCaLs][i_gKir_s_range_KirCaLs]])
        save_dI([p_var_d_i_batch,p_var_s_i_batch], find_dI_, jld_name_i)
        println("--------------------------------------------------------------------------------------------------------------")
    end   
   