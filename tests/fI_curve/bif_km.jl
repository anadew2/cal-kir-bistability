
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

include("data_loader.jl")
include("../dI_grad/local_postpro_output.jl")
include("local_postpro_outputs_batch.jl")


V_I1_KMCaLs_fI = V_I1_KMCaLs[1:20:end,1:25:end]
V_I2_KMCaLs_fI = V_I2_KMCaLs[1:20:end,1:25:end]

stp_V_FP_range=0.01
V_FP_range = range(-150.,stop=-20.,step=stp_V_FP_range)

I_SN_KMCaLs_fI = []
V_SN_KMCaLs_fI = []
I_saddle_KMCaLs_fI = []
V_saddle_KMCaLs_fI = []
I_stable_KMCaLs_fI = []
V_stable_KMCaLs_fI = []
I_unstable_KMCaLs_fI = []
V_unstable_KMCaLs_fI = []

for i_pCaLs_range_KMCaLs in 1:3

    I_SN_KMCaLs_fI_ = []
    V_SN_KMCaLs_fI_ = []
    I_saddle_KMCaLs_fI_ = []
    V_saddle_KMCaLs_fI_ = []
    I_stable_KMCaLs_fI_ = []
    V_stable_KMCaLs_fI_ = []
    I_unstable_KMCaLs_fI_ = []
    V_unstable_KMCaLs_fI_ = []

    for i_gKM_range_KMCaLs in 1:3

        pCaLs = pCaLs_range_KMCaLs_fI[i_pCaLs_range_KMCaLs]
        gKM = gKM_range_KMCaLs_fI[i_gKM_range_KMCaLs]

        I_list_batch = collect((I1_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]-0.5):0.01:(I2_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]+0.5))[1:end]

        p_var_i_batch = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

        I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable= find_fixed_points([p_var_i_batch,p_fixed],V_FP_range,(I1_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]-0.5),0)

        push!(I_SN_KMCaLs_fI_,I_SN)
        push!(V_SN_KMCaLs_fI_,V_SN)
        push!(I_saddle_KMCaLs_fI_,I_saddle)
        push!(V_saddle_KMCaLs_fI_,V_saddle)
        push!(I_stable_KMCaLs_fI_,I_stable)
        push!(V_stable_KMCaLs_fI_,V_stable)
        push!(I_unstable_KMCaLs_fI_,I_unstable)
        push!(V_unstable_KMCaLs_fI_,V_unstable)
    end

    plt_V = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(-2,10),size=(800,400))
    plot!(I_stable_KMCaLs_fI_ ,V_stable_KMCaLs_fI_ ,palette=palette(:Greens_6,rev=true),lw=2)
	plot!(I_saddle_KMCaLs_fI_ ,V_saddle_KMCaLs_fI_ , palette=palette([:tan1,:chocolate3],3), ls=:dot,lw=2)
	plot!(I_unstable_KMCaLs_fI_,V_unstable_KMCaLs_fI_ , palette=palette([:tomato,:darkred],3), ls=:dash,lw=2)
    scatter!(I_SN_KMCaLs_fI_ , V_SN_KMCaLs_fI_ ,mc=:black,markerstrokewidth=0,label=:none)
    display(plt_V)

    push!(I_SN_KMCaLs_fI,I_SN_KMCaLs_fI_)
    push!(V_SN_KMCaLs_fI,V_SN_KMCaLs_fI_)
    push!(I_saddle_KMCaLs_fI,I_saddle_KMCaLs_fI_)
    push!(V_saddle_KMCaLs_fI,V_saddle_KMCaLs_fI_)
    push!(I_stable_KMCaLs_fI,I_stable_KMCaLs_fI_)
    push!(V_stable_KMCaLs_fI,V_stable_KMCaLs_fI_)
    push!(I_unstable_KMCaLs_fI,I_unstable_KMCaLs_fI_)
    push!(V_unstable_KMCaLs_fI,V_unstable_KMCaLs_fI_)

end


for i_cal in 1:1
    plt = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(-3,10),size=(800,400))
    plot!(I_list_KMCaLs_fI[i_cal],V_max_c_KMCaLs_fI[i_cal],palette=palette(:BuPu_3,rev=true),lw=2)
    plot!(I_list_KMCaLs_fI[i_cal],V_min_c_KMCaLs_fI[i_cal],palette=palette(:BuPu_3,rev=true),lw=2)
    plot!(I_stable_KMCaLs_fI[i_cal] ,V_stable_KMCaLs_fI[i_cal],palette=palette(:Greens_6,rev=true),lw=2)
    plot!(I_saddle_KMCaLs_fI[i_cal] ,V_saddle_KMCaLs_fI[i_cal], palette=palette([:tan1,:chocolate3],3), ls=:dot,lw=2)
    plot!(I_unstable_KMCaLs_fI[i_cal],V_unstable_KMCaLs_fI[i_cal], palette=palette([:tomato,:darkred],3), ls=:dash,lw=2)
    scatter!(I_SN_KMCaLs_fI[i_cal], V_SN_KMCaLs_fI[i_cal],mc=:black,markerstrokewidth=0,label=:none)
    scatter!(I1_KMCaLs_fI[i_cal,1:3], V_I1_KMCaLs_fI[i_cal,1:3],mc=:chartreuse4,markerstrokewidth=0,label=:none)
    #scatter!(I2_KMCaLs_fI[i_cal,1:3], V_I2_KMCaLs_fI[i_cal,1:3],mc=:darkgreen,markerstrokewidth=0,label=:none)
    display(plt)
end

## Display the results 
#==
    palette_KM_fI = palette([palette(:Reds_9)[4],palette(:Reds_9)[9]],9)
    for i in 1:1
        for j in 3:3
            fI_fill_ = collect(1:length(I_list_KMCaLs_fI[i][j]))[I_list_KMCaLs_fI[i][j] .>=I1_KMCaLs_fI[i,j]]
            fI_fill = fI_fill_[I_list_KMCaLs_fI[i][j][fI_fill_] .<= I2_KMCaLs_fI[i,j]]
            plt_fI = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="f",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(I1_KMCaLs_fI[i,j]-(I2_KMCaLs_fI[i,j]+0.5 -I1_KMCaLs_fI[i,j]),I2_KMCaLs_fI[i,j]+0.5),size=(800,400))
            plot!(plt_fI,I_list_KMCaLs_fI[i][j],f_list_KMCaLs_fI[i][j],lc=palette_KM_fI[Int((i-1)*3+j)],lw=2)
            plot!(I_list_KMCaLs_fI[i][j][fI_fill],f_list_KMCaLs_fI[i][j][fI_fill],fillrange=zeros(size(fI_fill)),c=palette_KM_fI[Int((i-1)*3+j)],fillalpha=0.15,lw=0,label=:none)
            
            VI_fill_ = collect(1:length(I_stable_KMCaLs_fI[i][j]))[I_stable_KMCaLs_fI[i][j] .>=I1_KMCaLs_fI[i,j]]
            VI_fill = VI_fill_[I_stable_KMCaLs_fI[i][j][VI_fill_] .<= I2_KMCaLs_fI[i,j]]  
            plt_VI = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(I1_KMCaLs_fI[i,j]-(I2_KMCaLs_fI[i,j] +0.5-I1_KMCaLs_fI[i,j]),I2_KMCaLs_fI[i,j]+0.5),ylims=(-150,-50),size=(800,400))
            plot!(I_stable_KMCaLs_fI[i][j] ,V_stable_KMCaLs_fI[i][j],ls=:dash,lw=2,lc=palette_KM_fI[Int((i-1)*3+j)])
            plot!(I_stable_KMCaLs_fI[i][j][VI_fill],V_stable_KMCaLs_fI[i][j][VI_fill],fillrange=ones(size(VI_fill)).*-55,c=palette_KM_fI[Int((i-1)*3+j)],fillalpha=0.15,lw=0,label=:none)
            scatter!(I1_KMCaLs_fI[i,j].*ones(2), V_I1_KMCaLs_fI[i,j].*ones(2),mc=:chartreuse4,markerstrokewidth=0,label=:none)
            scatter!(I2_KMCaLs_fI[i,j].*ones(2), V_I2_KMCaLs_fI[i,j].*ones(2),mc=:darkgreen,markerstrokewidth=0,label=:none)    
            plt=Plots.plot(plt_fI,plt_VI,layout=(2,1),xlims=(-3,8))
            display(plt)
        end
    end
==#

