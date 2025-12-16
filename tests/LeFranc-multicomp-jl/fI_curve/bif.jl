
using Plots,ColorSchemes,LaTeXStrings
using DifferentialEquations,NLsolve
using LinearAlgebra
using JLD


include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("jld_hdle_fI.jl")
include("data_loader.jl")
include("fct_fI_curve.jl")

include("local_postpro_outputs_batch.jl")

ind_CaLs_d = [3,16,29]
ind_Kir_s = [4,21,38]
ind_Kir_d = [4,21,38]

Vs_I1_Kir_s_CaLs_d_fI = Vs_I1_Kir_s_CaLs_d[ind_CaLs_d[:],ind_Kir_d[:]]
Vs_I2_Kir_s_CaLs_d_fI = Vs_I2_Kir_s_CaLs_d[ind_CaLs_d[:],ind_Kir_d[:]]

stp_V_FP_range=0.1
Vs_FP_range = range(-100.,stop=-35.,step=stp_V_FP_range)

I_SN_Kir_s_CaLs_d_fI = []
V_SN_Kir_s_CaLs_d_fI = []
I_saddle_Kir_s_CaLs_d_fI = []
Vd_saddle_Kir_s_CaLs_d_fI = []
Cad_saddle_Kir_s_CaLs_d_fI = []
Vs_saddle_Kir_s_CaLs_d_fI = []
Cas_saddle_Kir_s_CaLs_d_fI = []
I_stable_Kir_s_CaLs_d_fI = []
Vd_stable_Kir_s_CaLs_d_fI = []
Cad_stable_Kir_s_CaLs_d_fI = []
Vs_stable_Kir_s_CaLs_d_fI = []
Cas_stable_Kir_s_CaLs_d_fI = []
I_unstable_Kir_s_CaLs_d_fI = []
Vd_unstable_Kir_s_CaLs_d_fI = []
Cad_unstable_Kir_s_CaLs_d_fI = []
Vs_unstable_Kir_s_CaLs_d_fI = []
Cas_unstable_Kir_s_CaLs_d_fI = []

for i_pCaLs_d in 1:3

    I_SN_Kir_s_CaLs_d_fI_ = []
    V_SN_Kir_s_CaLs_d_fI_ = []
    I_saddle_Kir_s_CaLs_d_fI_ = []
    Vd_saddle_Kir_s_CaLs_d_fI_ = []
    Cad_saddle_Kir_s_CaLs_d_fI_ = []
    Vs_saddle_Kir_s_CaLs_d_fI_ = []
    Cas_saddle_Kir_s_CaLs_d_fI_ = []
    I_stable_Kir_s_CaLs_d_fI_ = []
    Vd_stable_Kir_s_CaLs_d_fI_ = []
    Cad_stable_Kir_s_CaLs_d_fI_ = []
    Vs_stable_Kir_s_CaLs_d_fI_ = []
    Cas_stable_Kir_s_CaLs_d_fI_ = []
    I_unstable_Kir_s_CaLs_d_fI_ = []
    Vd_unstable_Kir_s_CaLs_d_fI_ = []
    Cad_unstable_Kir_s_CaLs_d_fI_ = []
    Vs_unstable_Kir_s_CaLs_d_fI_ = []
    Cas_unstable_Kir_s_CaLs_d_fI_ = []

    for i_gKir_s in 1:3 
        println("----------------- i=$(i_pCaLs_d) ; j=$(i_gKir_s) -----------------")
        pCaLs_d_i_batch = pCaLs_d_range_KirCaLs_fI[i_pCaLs_d]
        gKir_s_i_batch = gKir_s_range_KirCaLs_fI[i_gKir_s]

        I_list_batch = I_list_Kir_s_CaLs_d_fI[i_pCaLs_d][i_gKir_s]

        p_var_d_i_batch = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d_i_batch,gleak_d]
        p_var_s_i_batch = [gNa_s,gKDR_s,gKir_s_i_batch,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]

        Vd_saddle,Cad_saddle,Vs_saddle,Cas_saddle,I_saddle,Vd_stable,Cad_stable,Vs_stable,Cas_stable,I_stable,Vd_unstable,Cad_unstable,Vs_unstable,Cas_unstable,I_unstable = find_fixed_points([vcat(0,p_fixed[2:end]),p_var_d_i_batch,p_var_s_i_batch],Vs_FP_range,[25/10^(6),-59. ,Ca_i_0,Ca_i_0])
        I_SN = NaN
        V_SN = NaN

        push!(I_SN_Kir_s_CaLs_d_fI_,I_SN)
        push!(V_SN_Kir_s_CaLs_d_fI_,V_SN)
        push!(I_saddle_Kir_s_CaLs_d_fI_,I_saddle)
        push!(Vd_saddle_Kir_s_CaLs_d_fI_,Vd_saddle)
        push!(Cad_saddle_Kir_s_CaLs_d_fI_,Cad_saddle)
        push!(Vs_saddle_Kir_s_CaLs_d_fI_,Vs_saddle)
        push!(Cas_saddle_Kir_s_CaLs_d_fI_,Cas_saddle)
        push!(I_stable_Kir_s_CaLs_d_fI_,I_stable)
        push!(Vd_stable_Kir_s_CaLs_d_fI_,Vd_stable)
        push!(Cad_stable_Kir_s_CaLs_d_fI_,Cad_stable)
        push!(Vs_stable_Kir_s_CaLs_d_fI_,Vs_stable)
        push!(Cas_stable_Kir_s_CaLs_d_fI_,Cas_stable)
        push!(I_unstable_Kir_s_CaLs_d_fI_,I_unstable)
        push!(Vd_unstable_Kir_s_CaLs_d_fI_,Vd_unstable)
        push!(Cad_unstable_Kir_s_CaLs_d_fI_,Cad_unstable)
        push!(Vs_unstable_Kir_s_CaLs_d_fI_,Vs_unstable)
        push!(Cas_unstable_Kir_s_CaLs_d_fI_,Cas_unstable)

    end

    plt_V = Plots.plot(fontfamily="Computer Modern", xlabel="I [pA]",ylabel="V [mV]",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(xlims=(-100,250))
    plot!(I_stable_Kir_s_CaLs_d_fI_ .*10^6,Vs_stable_Kir_s_CaLs_d_fI_ ,palette=palette(:Greens_6,rev=true),lw=2)
	plot!(I_saddle_Kir_s_CaLs_d_fI_ *10^6,Vs_saddle_Kir_s_CaLs_d_fI_ , palette=palette([:tan1,:chocolate3],3), ls=:dot,lw=2)
	plot!(I_unstable_Kir_s_CaLs_d_fI_*10^6,Vs_unstable_Kir_s_CaLs_d_fI_ , palette=palette([:tomato,:darkred],3), ls=:dash,lw=2)
    scatter!(I_SN_Kir_s_CaLs_d_fI_*10^6, V_SN_Kir_s_CaLs_d_fI_ ,mc=:black,markerstrokewidth=0,label=:none)
    display(plt_V)

    push!(I_SN_Kir_s_CaLs_d_fI,I_SN_Kir_s_CaLs_d_fI_)
    push!(V_SN_Kir_s_CaLs_d_fI,V_SN_Kir_s_CaLs_d_fI_)
    push!(I_saddle_Kir_s_CaLs_d_fI,I_saddle_Kir_s_CaLs_d_fI_)
    push!(Vd_saddle_Kir_s_CaLs_d_fI,Vd_saddle_Kir_s_CaLs_d_fI_)
    push!(Cad_saddle_Kir_s_CaLs_d_fI,Cad_saddle_Kir_s_CaLs_d_fI_)
    push!(Vs_saddle_Kir_s_CaLs_d_fI,Vs_saddle_Kir_s_CaLs_d_fI_)
    push!(Cas_saddle_Kir_s_CaLs_d_fI,Cas_saddle_Kir_s_CaLs_d_fI_)
    push!(I_stable_Kir_s_CaLs_d_fI,I_stable_Kir_s_CaLs_d_fI_)
    push!(Vd_stable_Kir_s_CaLs_d_fI,Vd_stable_Kir_s_CaLs_d_fI_)
    push!(Cad_stable_Kir_s_CaLs_d_fI,Cad_stable_Kir_s_CaLs_d_fI_)
    push!(Vs_stable_Kir_s_CaLs_d_fI,Vs_stable_Kir_s_CaLs_d_fI_)
    push!(Cas_stable_Kir_s_CaLs_d_fI,Cas_stable_Kir_s_CaLs_d_fI_)
    push!(I_unstable_Kir_s_CaLs_d_fI,I_unstable_Kir_s_CaLs_d_fI_)
    push!(Vd_unstable_Kir_s_CaLs_d_fI,Vd_unstable_Kir_s_CaLs_d_fI_)
    push!(Cad_unstable_Kir_s_CaLs_d_fI,Cad_unstable_Kir_s_CaLs_d_fI_)
    push!(Vs_unstable_Kir_s_CaLs_d_fI,Vs_unstable_Kir_s_CaLs_d_fI_)
    push!(Cas_unstable_Kir_s_CaLs_d_fI,Cas_unstable_Kir_s_CaLs_d_fI_)

end

## Display the results 
#

plt_fI = Plots.plot(fontfamily="Computer Modern", xlabel="I [pA]",ylabel="f [Hz]",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
plt_VI = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,ylims=(-150,-40),size=(800,400))
                                
    palette_Kir_fI = palette(:tab10)
    for i in 2:2
        for j in [1,2,3]#1:3
            
            println("-------------- i = $i  ;  j = $j --------------")
                fI_fill_ = collect(1:length(I_list_Kir_s_CaLs_d_fI[i][j]))[I_list_Kir_s_CaLs_d_fI[i][j] .>=I1_Kir_s_CaLs_d_fI[i,j]]
                fI_fill = fI_fill_[I_list_Kir_s_CaLs_d_fI[i][j][fI_fill_] .<= I2_Kir_s_CaLs_d_fI[i,j]]
                
                println(size(f_list_Kir_s_CaLs_d_fI[i][j]))
                println(size(I_list_Kir_s_CaLs_d_fI[i][j]))

                plot!(plt_fI,I_list_Kir_s_CaLs_d_fI[i][j].*10^6,f_list_Kir_s_CaLs_d_fI[i][j],lc=palette_Kir_fI[Int((i-1)*3+j)],lw=2)
                plot!(plt_fI,I_list_Kir_s_CaLs_d_fI[i][j][fI_fill].*10^6,f_list_Kir_s_CaLs_d_fI[i][j][fI_fill],fillrange=zeros(size(fI_fill)),c=palette_Kir_fI[Int((i-1)*3+j)],fillalpha=0.15,lw=0,label=:none)
                println(plt_fI,"xlims=$((I1_Kir_s_CaLs_d_fI[i,j]-(I2_Kir_s_CaLs_d_fI[i,j]+0.5 -I1_Kir_s_CaLs_d_fI[i,j]),I2_Kir_s_CaLs_d_fI[i,j]+0.5).*10^6)")
                
                VI_fill_ = collect(1:length(I_stable_Kir_s_CaLs_d_fI[i][j]))[I_stable_Kir_s_CaLs_d_fI[i][j] .>=I1_Kir_s_CaLs_d_fI[i,j]]
                VI_fill = VI_fill_[I_stable_Kir_s_CaLs_d_fI[i][j][VI_fill_] .<= I2_Kir_s_CaLs_d_fI[i,j]]  
                plot!(plt_VI,I_stable_Kir_s_CaLs_d_fI[i][j].*10^6 ,Vs_stable_Kir_s_CaLs_d_fI[i][j],ls=:dash,lw=2,lc=palette_Kir_fI[Int((i-1)*3+j)])
                plot!(plt_VI,I_stable_Kir_s_CaLs_d_fI[i][j][VI_fill].*10^6,Vs_stable_Kir_s_CaLs_d_fI[i][j][VI_fill],fillrange=ones(size(VI_fill)).*-55,c=palette_Kir_fI[Int((i-1)*3+j)],fillalpha=0.15,lw=0,label=:none)
                scatter!(plt_VI, I1_Kir_s_CaLs_d_fI[i,j].*10^6 .*ones(2), Vs_I1_Kir_s_CaLs_d_fI[i,j].*ones(2),mc=:chartreuse4,markerstrokewidth=0,label=:none)
                scatter!(plt_VI,I2_Kir_s_CaLs_d_fI[i,j].*10^6 .*ones(2), Vs_I2_Kir_s_CaLs_d_fI[i,j].*ones(2),mc=:darkgreen,markerstrokewidth=0,label=:none)    
                
                
            
        end
    end
    plt=Plots.plot(plt_fI,plt_VI,layout=(2,1),xlims=(-150,250))
    display(plt)



    

    for i in 1:3
        plt_V = Plots.plot(fontfamily="Computer Modern", xlabel="I [pA]",ylabel="Vs [mV]",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
        plot!(xlims=(-100,250))
        plot!(I_stable_Kir_s_CaLs_d_fI[i] .*10^6,Vs_stable_Kir_s_CaLs_d_fI[i]  ,palette=palette(:Greens_6,rev=true),lw=2)
        plot!(I_saddle_Kir_s_CaLs_d_fI[i]  *10^6,Vs_saddle_Kir_s_CaLs_d_fI[i] , palette=palette([:tan1,:chocolate3],3), ls=:dot,lw=2)
        plot!(I_unstable_Kir_s_CaLs_d_fI[i] *10^6,Vs_unstable_Kir_s_CaLs_d_fI[i]  , palette=palette([:tomato,:darkred],3), ls=:dash,lw=2)
        scatter!(I_SN_Kir_s_CaLs_d_fI[i] *10^6, V_SN_Kir_s_CaLs_d_fI[i] ,mc=:black,markerstrokewidth=0,label=:none)
        plot!(I_list_Kir_s_CaLs_d_fI[i] *10^6,Vs_max_c_Kir_s_CaLs_d_fI[i],lc=:midnightblue)
        plot!(I_list_Kir_s_CaLs_d_fI[i] *10^6,Vs_min_c_Kir_s_CaLs_d_fI[i],lc=:cornflowerblue,label=["1","2","3"])
        display(plt_V)
    end

    for i in 1:3
        plt_V = Plots.plot(fontfamily="Computer Modern", xlabel="I [pA]",ylabel="Vd [mV]",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
        plot!(xlims=(-100,250))
        plot!(I_stable_Kir_s_CaLs_d_fI[i] .*10^6,Vd_stable_Kir_s_CaLs_d_fI[i]  ,palette=palette(:Greens_6,rev=true),lw=2)
        plot!(I_saddle_Kir_s_CaLs_d_fI[i]  *10^6,Vd_saddle_Kir_s_CaLs_d_fI[i] , palette=palette([:tan1,:chocolate3],3), ls=:dot,lw=2)
        plot!(I_unstable_Kir_s_CaLs_d_fI[i] *10^6,Vd_unstable_Kir_s_CaLs_d_fI[i]  , palette=palette([:tomato,:darkred],3), ls=:dash,lw=2)
        
        display(plt_V)
    end
    ==#