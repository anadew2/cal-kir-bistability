using JLD

include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("data_loader.jl")
include("jld_hdle_dI_grad.jl")

include("fct_local_postpro_output.jl")


    Ss = 2*pi*(Ds/2)*Ls *10^(-8) # [cm²]

    jld_name = "tests/LeFranc-multicomp-jl/dI_grad/data/jld/dI_Kir_s_CaLs_d_compart_"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        ext_name = ["i_pCaLs_d_","_i_gKir_s_",".jld"]
        dI_Kir_s_CaLs_d,I1_Kir_s_CaLs_d,I2_Kir_s_CaLs_d,f_min_Kir_s_CaLs_d = create_dI_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_s_range_KirCaLs)   
        save_dI_mat(dI_Kir_s_CaLs_d,I1_Kir_s_CaLs_d,I2_Kir_s_CaLs_d,f_min_Kir_s_CaLs_d,string(jld_name,".jld"))
    ==#

    dI_Kir_s_CaLs_d,I1_Kir_s_CaLs_d,I2_Kir_s_CaLs_d,f_min_Kir_s_CaLs_d = load_dI_mat(string(jld_name,".jld"))
    Vd_I1_Kir_s_CaLs_d,Cad_I1_Kir_s_CaLs_d,Vs_I1_Kir_s_CaLs_d,Cas_I1_Kir_s_CaLs_d,Vd_I2_Kir_s_CaLs_d,Cad_I2_Kir_s_CaLs_d,Vs_I2_Kir_s_CaLs_d,Cas_I2_Kir_s_CaLs_d =create_V_Ix_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_s_range_KirCaLs,p_fixed,true)  
    
    #dI in pA
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,dI_Kir_s_CaLs_d.*(10^6),fill=true,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I1 in pA
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,I1_Kir_s_CaLs_d.*(10^6),fill=true,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I2 in pA
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,I2_Kir_s_CaLs_d.*(10^6),fill=true,c=:linear_wcmr_100_45_c42_n256)
    
    #dI in µA/cm²
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,dI_Kir_s_CaLs_d./Ss,fill=true,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I1 in µA/cm²
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,I1_Kir_s_CaLs_d./Ss,fill=true,lw=0,c=:linear_wcmr_100_45_c42_n256) 
    #I2 in µA/cm²
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,I2_Kir_s_CaLs_d./Ss,fill=true,c=:linear_wcmr_100_45_c42_n256)
    
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,f_min_Kir_s_CaLs_d,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)

   
    #
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I1_Kir_s_CaLs_d,fill=true,levels=100,lw=0,c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true))
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I2_Kir_s_CaLs_d,fill=true,levels=100,lw=0,c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true))
    ==#


    jld_name = "tests/LeFranc-multicomp-jl/dI_grad/data/jld/dI_Kir_d_CaLs_d_compart_"
    #==
    # This block load the results from the separated jld files and combine them in a single jld 
        ext_name = ["i_pCaLs_d_","_i_gKir_d_",".jld"]
        dI_Kir_d_CaLs_d,I1_Kir_d_CaLs_d,I2_Kir_d_CaLs_d,f_min_Kir_d_CaLs_d = create_dI_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_d_range_KirCaLs)   
        save_dI_mat(dI_Kir_d_CaLs_d,I1_Kir_d_CaLs_d,I2_Kir_d_CaLs_d,f_min_Kir_d_CaLs_d,string(jld_name,".jld"))
    ==#

    dI_Kir_d_CaLs_d,I1_Kir_d_CaLs_d,I2_Kir_d_CaLs_d,f_min_Kir_d_CaLs_d = load_dI_mat(string(jld_name,".jld"))
    Vd_I1_Kir_d_CaLs_d,Cad_I1_Kir_d_CaLs_d,Vs_I1_Kir_d_CaLs_d,Cas_I1_Kir_d_CaLs_d,Vd_I2_Kir_d_CaLs_d,Cad_I2_Kir_d_CaLs_d,Vs_I2_Kir_d_CaLs_d,Cas_I2_Kir_d_CaLs_d =create_V_Ix_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_d_range_KirCaLs,p_fixed,false)  



    #dI in pA
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,dI_Kir_d_CaLs_d.*(10^6),fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I1 in pA
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,I1_Kir_d_CaLs_d.*(10^6),fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I2 in pA
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,I2_Kir_d_CaLs_d.*(10^6),fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    
    #dI in µA/cm²
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,dI_Kir_d_CaLs_d./(Ss),fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I1 in µA/cm²
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,I1_Kir_d_CaLs_d./(Ss),fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)
    #I2 in µA/cm²
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,I2_Kir_d_CaLs_d./(Ss),fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)


    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,f_min_Kir_d_CaLs_d,fill=true,levels=100,lw=0,c=:linear_wcmr_100_45_c42_n256)


    #
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I1_Kir_d_CaLs_d,fill=true,clims=(-90,-47),c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),xlabel="gKir_d",ylabel="pCaLs_d",title="V(I1)",fontfamily="Computer Modern")
    plt = Plots.heatmap(gKir_d_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I2_Kir_d_CaLs_d,fill=true,clims=(-61,-47),c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),xlabel="gKir_d",ylabel="pCaLs_d",title="V(I2)",fontfamily="Computer Modern")

    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I1_Kir_s_CaLs_d,fill=true,clims=(-61,-47),c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),xlabel="gKir_s",ylabel="pCaLs_d",title="V(I1)",fontfamily="Computer Modern")
    plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I2_Kir_s_CaLs_d,fill=true,clims=(-61,-47),c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),xlabel="gKir_s",ylabel="pCaLs_d",title="V(I2)",fontfamily="Computer Modern")
