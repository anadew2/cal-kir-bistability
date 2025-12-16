
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

include("data_loader.jl")
include("../dI_grad_altgt/local_postpro_output.jl")
include("local_postpro_outputs_batch.jl")


V_I1_KirKMCaLs_fI = V_I1_KirKMCaLs[3:2:5,3:26:29,3:26:29]
V_I2_KirKMCaLs_fI = V_I2_KirKMCaLs[3:2:5,3:26:29,3:26:29]

stp_V_FP_range=0.01
V_FP_range = range(-150.,stop=-20.,step=stp_V_FP_range)

I_SN_KirKMCaLs_fI = []
V_SN_KirKMCaLs_fI = []
I_saddle_KirKMCaLs_fI = []
V_saddle_KirKMCaLs_fI = []
I_stable_KirKMCaLs_fI = []
V_stable_KirKMCaLs_fI = []
I_unstable_KirKMCaLs_fI = []
V_unstable_KirKMCaLs_fI = []

for i_pCaLs_range_KirKMCaLs in 1:2

    I_SN_KirKMCaLs_fI_ = []
    V_SN_KirKMCaLs_fI_ = []
    I_saddle_KirKMCaLs_fI_ = []
    V_saddle_KirKMCaLs_fI_ = []
    I_stable_KirKMCaLs_fI_ = []
    V_stable_KirKMCaLs_fI_ = []
    I_unstable_KirKMCaLs_fI_ = []
    V_unstable_KirKMCaLs_fI_ = []

    for i_gKir_range_KirKMCaLs in 1:2
        I_SN_KirKMCaLs_fI__ = []
        V_SN_KirKMCaLs_fI__ = []
        I_saddle_KirKMCaLs_fI__ = []
        V_saddle_KirKMCaLs_fI__ = []
        I_stable_KirKMCaLs_fI__ = []
        V_stable_KirKMCaLs_fI__ = []
        I_unstable_KirKMCaLs_fI__ = []
        V_unstable_KirKMCaLs_fI__ = []

        for i_gKM_range_KirKMCaLs in 1:2

            pCaLs = pCaLs_range_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs]
            gKir = gKir_range_KirKMCaLs_fI[i_gKir_range_KirKMCaLs]
            gKM = gKM_range_KirKMCaLs_fI[i_gKM_range_KirKMCaLs]
            I_list_batch = collect((I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]-0.5):0.01:(I2_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]+0.5))[1:end]

            p_var_i_batch = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

            I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable= find_fixed_points([p_var_i_batch,p_fixed],V_FP_range,(I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]-0.5),0)

            push!(I_SN_KirKMCaLs_fI__,I_SN)
            push!(V_SN_KirKMCaLs_fI__,V_SN)
            push!(I_saddle_KirKMCaLs_fI__,I_saddle)
            push!(V_saddle_KirKMCaLs_fI__,V_saddle)
            push!(I_stable_KirKMCaLs_fI__,I_stable)
            push!(V_stable_KirKMCaLs_fI__,V_stable)
            push!(I_unstable_KirKMCaLs_fI__,I_unstable)
            push!(V_unstable_KirKMCaLs_fI__,V_unstable)

        end
        push!(I_SN_KirKMCaLs_fI_,I_SN_KirKMCaLs_fI__)
        push!(V_SN_KirKMCaLs_fI_,V_SN_KirKMCaLs_fI__)
        push!(I_saddle_KirKMCaLs_fI_,I_saddle_KirKMCaLs_fI__)
        push!(V_saddle_KirKMCaLs_fI_,V_saddle_KirKMCaLs_fI__)
        push!(I_stable_KirKMCaLs_fI_,I_stable_KirKMCaLs_fI__)
        push!(V_stable_KirKMCaLs_fI_,V_stable_KirKMCaLs_fI__)
        push!(I_unstable_KirKMCaLs_fI_,I_unstable_KirKMCaLs_fI__)
        push!(V_unstable_KirKMCaLs_fI_,V_unstable_KirKMCaLs_fI__)
    end

    plt_V = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(-2,5),size=(800,400))
    plot!(I_stable_KirKMCaLs_fI_ ,V_stable_KirKMCaLs_fI_ ,palette=palette(:Greens_6,rev=true),lw=2)
	plot!(I_saddle_KirKMCaLs_fI_ ,V_saddle_KirKMCaLs_fI_ , palette=palette([:tan1,:chocolate3],3), ls=:dot,lw=2)
	plot!(I_unstable_KirKMCaLs_fI_,V_unstable_KirKMCaLs_fI_ , palette=palette([:tomato,:darkred],3), ls=:dash,lw=2)
    scatter!(I_SN_KirKMCaLs_fI_ , V_SN_KirKMCaLs_fI_ ,mc=:black,markerstrokewidth=0,label=:none)
    display(plt_V)

    push!(I_SN_KirKMCaLs_fI,I_SN_KirKMCaLs_fI_)
    push!(V_SN_KirKMCaLs_fI,V_SN_KirKMCaLs_fI_)
    push!(I_saddle_KirKMCaLs_fI,I_saddle_KirKMCaLs_fI_)
    push!(V_saddle_KirKMCaLs_fI,V_saddle_KirKMCaLs_fI_)
    push!(I_stable_KirKMCaLs_fI,I_stable_KirKMCaLs_fI_)
    push!(V_stable_KirKMCaLs_fI,V_stable_KirKMCaLs_fI_)
    push!(I_unstable_KirKMCaLs_fI,I_unstable_KirKMCaLs_fI_)
    push!(V_unstable_KirKMCaLs_fI,V_unstable_KirKMCaLs_fI_)

end

## Display the results 
#


    palette_Kir_fI = palette(:Paired_8)
    for i in 2:-1:2
        for j in 2:-1:1
            for k in 2:-1:1
                fI_fill_ = collect(1:length(I_list_KirKMCaLs_fI[i][j][k]))[I_list_KirKMCaLs_fI[i][j][k] .>=I1_KirKMCaLs_fI[i,j,k]]
                fI_fill = fI_fill_[I_list_KirKMCaLs_fI[i][j][k][fI_fill_] .<= I2_KirKMCaLs_fI[i,j,k]]
                plt_fI = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="f",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(I1_KirKMCaLs_fI[i,j,k]-(I2_KirKMCaLs_fI[i,j,k]+0.5 -I1_KirKMCaLs_fI[i,j,k]),I2_KirKMCaLs_fI[i,j,k]+0.5),size=(800,400))
                plot!(plt_fI,I_list_KirKMCaLs_fI[i][j][k],f_list_KirKMCaLs_fI[i][j][k],lc=palette_Kir_fI[Int((i-1)*4+(k-1)*2+j)],lw=2)
                plot!(I_list_KirKMCaLs_fI[i][j][k][fI_fill],f_list_KirKMCaLs_fI[i][j][k][fI_fill],fillrange=zeros(size(fI_fill)),c=palette_Kir_fI[Int((i-1)*4+(k-1)*2+j)],fillalpha=0.15,lw=0,label=:none)
                plot!(title="CaL = $(i); Kir = $(j); KM = $(k)")
                println("xlims=$((I1_KirKMCaLs_fI[i,j,k]-(I2_KirKMCaLs_fI[i,j,k]+0.5 -I1_KirKMCaLs_fI[i,j,k]),I2_KirKMCaLs_fI[i,j,k]+0.5))")
                
                VI_fill_ = collect(1:length(I_stable_KirKMCaLs_fI[i][j][k]))[I_stable_KirKMCaLs_fI[i][j][k] .>=I1_KirKMCaLs_fI[i,j,k]]
                VI_fill = VI_fill_[I_stable_KirKMCaLs_fI[i][j][k][VI_fill_] .<= I2_KirKMCaLs_fI[i,j,k]]  
                plt_VI = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(I1_KirKMCaLs_fI[i,j,k]-(I2_KirKMCaLs_fI[i,j,k] +0.5-I1_KirKMCaLs_fI[i,j,k]),I2_KirKMCaLs_fI[i,j,k]+0.5),ylims=(-150,-50),size=(800,400))
                plot!(I_stable_KirKMCaLs_fI[i][j][k] ,V_stable_KirKMCaLs_fI[i][j][k],ls=:dash,lw=2,lc=palette_Kir_fI[Int((i-1)*4+(k-1)*2+j)])
                plot!(I_stable_KirKMCaLs_fI[i][j][k][VI_fill],V_stable_KirKMCaLs_fI[i][j][k][VI_fill],fillrange=ones(size(VI_fill)).*-55,c=palette_Kir_fI[Int((i-1)*4+(k-1)*2+j)],fillalpha=0.15,lw=0,label=:none)
                scatter!(I1_KirKMCaLs_fI[i,j,k].*ones(2), V_I1_KirKMCaLs_fI[i,j,k].*ones(2),mc=:chartreuse4,markerstrokewidth=0,label=:none)
                scatter!(I2_KirKMCaLs_fI[i,j,k].*ones(2), V_I2_KirKMCaLs_fI[i,j,k].*ones(2),mc=:darkgreen,markerstrokewidth=0,label=:none)    
                plt=Plots.plot(plt_fI,plt_VI,layout=(2,1),xlims=(-3,10))
                display(plt)
            end
        end
    end
    ==#