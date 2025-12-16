
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../solver_fun.jl")

include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("jld_hdle_fI.jl")
include("data_loader.jl")
include("fct_fI_curve.jl")



### SCRIPT INPUTS ###
    i_pCaLs_d_range_KirCaLs = parse(Int64, ARGS[1])     #1:length(pCaLs_d_range_KirCaLs_fI)
    i_gKir_s_range_KirCaLs = parse(Int64, ARGS[2])      #1:length(gKir_s_range_KirCaLs_fI)
    i_batch = 1 #parse(Int64, ARGS[3])                   #1:n_batch

### INPUTS ###
    pCaLs_d_i_batch = pCaLs_d_range_KirCaLs_fI[i_pCaLs_d_range_KirCaLs]
    gKir_s_i_batch = gKir_s_range_KirCaLs_fI[i_gKir_s_range_KirCaLs]

    n_batch = 1

    stp_I_list_batch_1= 0.1/10^6 #1/10^6 #0.1/10^6 #[µA]
    stp_I_list_batch_2= 5/10^6 #[µA]
    stp_I = 50/10^6 #[µA]
    stp_I_2 = 75/10^6 #[µA]
    
    if I2_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]-I1_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs] >stp_I
        I_list_batch_1 = collect((I1_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]-stp_I):stp_I_list_batch_1:(I1_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]+stp_I))
        I_list_batch_2a = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I2_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]))
        I_list_batch_2b = collect((I2_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]):stp_I_list_batch_2:(minimum([275/10^6,I2_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]+stp_I+stp_I_2])))
        I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2a,I_list_batch_2b)
    else
        I_list_batch_1 = collect((I1_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]-stp_I):stp_I_list_batch_1:(I1_Kir_s_CaLs_d_fI[i_pCaLs_d_range_KirCaLs,i_gKir_s_range_KirCaLs]+stp_I))
        I_list_batch_2 = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(minimum([275/10^6,I_list_batch_1[end]+stp_I_list_batch_2+4*stp_I_2])))
        I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2)
    end
    
    ind_batch = Int(ceil(length(I_list_batch_)/n_batch))
    I_list_batch = I_list_batch_[( (Int(1+ ind_batch*(i_batch-1))) : Int(minimum([ind_batch + ind_batch*(i_batch-1),length(I_list_batch_)]))) ]

    V_ic = -90. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic,
        V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic]

    jld_name = "data/jld/fI_Kir_s_CaLs_d_compart_ext_"

### TASK ###
    jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_d_range_KirCaLs)_i_Kir_$(i_gKir_s_range_KirCaLs)_i_batch_$(i_batch).jld")

    p_var_d_i_batch = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
    p_var_s_i_batch = [gNa_s,gKDR_s,gKir_s_i_batch,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]


    find_full_bif_ = find_full_bif([p_fixed,p_var_d_i_batch,p_var_s_i_batch],ic,I_list_batch)
    save_full_bif([p_fixed,p_var_d_i_batch,p_var_s_i_batch], find_full_bif_, jld_name_i)