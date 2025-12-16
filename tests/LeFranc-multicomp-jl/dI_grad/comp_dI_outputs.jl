using JLD 

include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("data_loader_ext.jl")
include("jld_hdle_dI_grad.jl")

include("fct_local_postpro_output.jl")

Ss = 2*pi*(Ds/2)*Ls *10^(-8) # [cm²]

### I1pulse = max(I_stable)+100 pA
    jld_name = "tests/LeFranc-multicomp-jl/dI_grad/data/jld_s_J_I1/dI_Kir_s_CaLs_d_compart_ext_"
    #==
    ext_name = ["i_pCaLs_d_","_i_gKir_s_",".jld"]
    dI_Kir_s_CaLs_d_J,I1_Kir_s_CaLs_d_J,I2_Kir_s_CaLs_d_J,f_min_Kir_s_CaLs_d_J = create_dI_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_s_range_KirCaLs)   
    save_dI_mat(dI_Kir_s_CaLs_d_J,I1_Kir_s_CaLs_d_J,I2_Kir_s_CaLs_d_J,f_min_Kir_s_CaLs_d_J,string(jld_name,".jld"))
    ==#
    dI_Kir_s_CaLs_d_J,I1_Kir_s_CaLs_d_J,I2_Kir_s_CaLs_d_J,f_min_Kir_s_CaLs_d_J = load_dI_mat(string(jld_name,".jld")) 

    Plots.heatmap(dI_Kir_s_CaLs_d_J.*(10^6),fill=true,lw=0,clims=(0,200))
    Plots.heatmap(I2_Kir_s_CaLs_d_J.*(10^6),fill=true,lw=0,clims=(0,200))
    Plots.heatmap(I1_Kir_s_CaLs_d_J.*(10^6),fill=true,lw=0,clims=(-100,200))


### I1pulse = max(I_stable)+200 pA
    jld_name = "tests/LeFranc-multicomp-jl/dI_grad/data/jld_s_J_I1up/dI_Kir_s_CaLs_d_compart_ext_"
    #==
    ext_name = ["i_pCaLs_d_","_i_gKir_s_",".jld"]
    dI_Kir_s_CaLs_d_J_I1up,I1_Kir_s_CaLs_d_J_I1up,I2_Kir_s_CaLs_d_J_I1up,f_min_Kir_s_CaLs_d_J_I1up = create_dI_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_s_range_KirCaLs)   
    save_dI_mat(dI_Kir_s_CaLs_d_J_I1up,I1_Kir_s_CaLs_d_J_I1up,I2_Kir_s_CaLs_d_J_I1up,f_min_Kir_s_CaLs_d_J_I1up,string(jld_name,".jld"))
    ==#
    dI_Kir_s_CaLs_d_J_I1up,I1_Kir_s_CaLs_d_J_I1up,I2_Kir_s_CaLs_d_J_I1up,f_min_Kir_s_CaLs_d_J_I1up = load_dI_mat(string(jld_name,".jld")) 

    Plots.heatmap(dI_Kir_s_CaLs_d_J_I1up.*(10^6),fill=true,lw=0,clims=(0,200))
    Plots.heatmap(I2_Kir_s_CaLs_d_J_I1up.*(10^6),fill=true,lw=0,clims=(0,200))
    Plots.heatmap(I1_Kir_s_CaLs_d_J_I1up.*(10^6),fill=true,lw=0,clims=(-100,200))


### I1pulse = max(I_stable)+50 pA
    jld_name = "tests/LeFranc-multicomp-jl/dI_grad/data/jld_s_J_I1down/dI_Kir_s_CaLs_d_compart_ext_"
    #==
    ext_name = ["i_pCaLs_d_","_i_gKir_s_",".jld"]
    dI_Kir_s_CaLs_d_J_I1down,I1_Kir_s_CaLs_d_J_I1down,I2_Kir_s_CaLs_d_J_I1down,f_min_Kir_s_CaLs_d_J_I1down = create_dI_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_s_range_KirCaLs)   
    save_dI_mat(dI_Kir_s_CaLs_d_J_I1down,I1_Kir_s_CaLs_d_J_I1down,I2_Kir_s_CaLs_d_J_I1down,f_min_Kir_s_CaLs_d_J_I1down,string(jld_name,".jld"))
    ==#
    dI_Kir_s_CaLs_d_J_I1down,I1_Kir_s_CaLs_d_J_I1down,I2_Kir_s_CaLs_d_J_I1down,f_min_Kir_s_CaLs_d_J_I1down = load_dI_mat(string(jld_name,".jld")) 

    Plots.heatmap(dI_Kir_s_CaLs_d_J_I1down.*(10^6),fill=true,lw=0,clims=(0,200))
    Plots.heatmap(I2_Kir_s_CaLs_d_J_I1down.*(10^6),fill=true,lw=0,clims=(0,200))
    Plots.heatmap(I1_Kir_s_CaLs_d_J_I1down.*(10^6),fill=true,lw=0,clims=(-100,200))


### Finds the best data point in the results with +50, +100, and +200pA in pulse height
        best_I1pulse_Kir_s_CaLs_d_J = ones(length(pCaLs_d_range_KirCaLs),length(gKir_s_range_KirCaLs)).*NaN
        best_dI_Kir_s_CaLs_d_J = ones(length(pCaLs_d_range_KirCaLs),length(gKir_s_range_KirCaLs)).*NaN
        best_I1_Kir_s_CaLs_d_J = ones(length(pCaLs_d_range_KirCaLs),length(gKir_s_range_KirCaLs)).*NaN
        best_I2_Kir_s_CaLs_d_J = ones(length(pCaLs_d_range_KirCaLs),length(gKir_s_range_KirCaLs)).*NaN
        best_f_min_Kir_s_CaLs_d_J = ones(length(pCaLs_d_range_KirCaLs),length(gKir_s_range_KirCaLs)).*NaN
            for i in 1:31
                for j in 1:length(gKir_s_range_KirCaLs)
                    I1_computed = []
                    for mat in [I1_Kir_s_CaLs_d_J_I1down,I1_Kir_s_CaLs_d_J,I1_Kir_s_CaLs_d_J_I1up]
                        if mat[i,j] !== NaN 
                            push!(I1_computed,mat[i,j])
                        end
                    end
                    if length(I1_computed)>0
                        if minimum(I1_computed)==I1_Kir_s_CaLs_d_J_I1down[i,j] && I1_Kir_s_CaLs_d_J_I1down[i,j] !==NaN && minimum(abs.(diff(I1_computed[1:2])))>0
                            best_I1pulse_Kir_s_CaLs_d_J[i,j] = 50 #[pA]
                            best_dI_Kir_s_CaLs_d_J[i,j] = dI_Kir_s_CaLs_d_J_I1down[i,j]
                            best_I1_Kir_s_CaLs_d_J[i,j] = I1_Kir_s_CaLs_d_J_I1down[i,j]
                            best_I2_Kir_s_CaLs_d_J[i,j] = I2_Kir_s_CaLs_d_J_I1down[i,j]
                            best_f_min_Kir_s_CaLs_d_J[i,j] = f_min_Kir_s_CaLs_d_J_I1down[i,j]
                        else
                            if minimum(I1_computed)==I1_Kir_s_CaLs_d_J[i,j] && I1_Kir_s_CaLs_d_J[i,j] !==NaN
                                best_I1pulse_Kir_s_CaLs_d_J[i,j] = 100 #[pA]
                                best_dI_Kir_s_CaLs_d_J[i,j] = dI_Kir_s_CaLs_d_J[i,j]
                                best_I1_Kir_s_CaLs_d_J[i,j] = I1_Kir_s_CaLs_d_J[i,j]
                                best_I2_Kir_s_CaLs_d_J[i,j] = I2_Kir_s_CaLs_d_J[i,j]
                                best_f_min_Kir_s_CaLs_d_J[i,j] = f_min_Kir_s_CaLs_d_J[i,j]
                            else
                                if minimum(I1_computed)==I1_Kir_s_CaLs_d_J_I1up[i,j] && I1_Kir_s_CaLs_d_J_I1up[i,j] !==NaN 
                                    #the minimum(I1_computed)==I1_Kir_s_CaLs_d_I1up[i,j]
                                    best_I1pulse_Kir_s_CaLs_d_J[i,j] = 200 #[pA]
                                    best_dI_Kir_s_CaLs_d_J[i,j] = dI_Kir_s_CaLs_d_J_I1up[i,j]
                                    best_I1_Kir_s_CaLs_d_J[i,j] = I1_Kir_s_CaLs_d_J_I1up[i,j]
                                    best_I2_Kir_s_CaLs_d_J[i,j] = I2_Kir_s_CaLs_d_J_I1up[i,j]
                                    best_f_min_Kir_s_CaLs_d_J[i,j] = f_min_Kir_s_CaLs_d_J_I1up[i,j]
                                end
                            end
                        end            
                    end
                end
            end
            Plots.heatmap(best_I1pulse_Kir_s_CaLs_d_J,fill=true,lw=0,c=:vik,clims=(0,200),title="Best pulse height from max(I_stable)",xlabel="i_gKir_s",ylabel="i_pCaLs_d",fontfamily="Computer Modern")

            plt = Plots.heatmap(best_I1_Kir_s_CaLs_d_J.*(10^6),fill=true,lw=0,clims=(-100,200),c=:diverging_rainbow_bgymr_45_85_c67_n256)
            plt = Plots.heatmap(best_I2_Kir_s_CaLs_d_J.*(10^6),fill=true,lw=0,clims=(0,200),c=:diverging_rainbow_bgymr_45_85_c67_n256)
            plt = Plots.heatmap(best_dI_Kir_s_CaLs_d_J.*10^6,fill=true,lw=0,clims=(0,200))

            plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,best_f_min_Kir_s_CaLs_d_J,fill=true,lw=0,clims=(0,40))

            ext_name = ["i_pCaLs_d_","_i_gKir_s_",".jld"]
            Vd_I1_Kir_s_CaLs_d,Cad_I1_Kir_s_CaLs_d,Vs_I1_Kir_s_CaLs_d,Cas_I1_Kir_s_CaLs_d,Vd_I2_Kir_s_CaLs_d,Cad_I2_Kir_s_CaLs_d,Vs_I2_Kir_s_CaLs_d,Cas_I2_Kir_s_CaLs_d =create_V_Ix_mat(jld_name,ext_name,pCaLs_d_range_KirCaLs,gKir_s_range_KirCaLs,p_fixed,true,best_dI_Kir_s_CaLs_d_J,best_I1_Kir_s_CaLs_d_J,best_I2_Kir_s_CaLs_d_J)  
            plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I2_Kir_s_CaLs_d,fill=true,clims=(-90,-47),c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),title=L"V(I_1) \ \ [mV]",xlabel="gKir_s",ylabel="pCaLs_d",fontfamily="Computer Modern")
            plt = Plots.heatmap(gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,Vs_I1_Kir_s_CaLs_d,fill=true,clims=(-90,-47),c=cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true),title=L"V(I_1) \ \ [mV]",xlabel="gKir_s",ylabel="pCaLs_d",fontfamily="Computer Modern")

    ==#

    