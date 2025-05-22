using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using Statistics

include("../sim_pulse_step/sim_CaLs.jl")
include("../fI_curve/bif_kir.jl")

function plt_fI_VI(I_list_fI,f_list_fI,I_stable_VI,V_stable_VI,I1,I2,V_I1,V_I2,xlims_,color)
    fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
    fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]
    fontsize_= 25

    plt_fI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legend=:outerbottomright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,size=(800,400))
    plot!(plt_fI,I_list_fI,f_list_fI,lc=color,lw=2,label="Spiking")
    plot!(I_stable_VI ,zeros(size(I_stable_VI)),ls=:dash,lw=2,lc=color,label="Resting")
    plot!(I_list_fI[fI_fill],f_list_fI[fI_fill],fillrange=zeros(size(fI_fill)),c=color,fillalpha=0.15,lw=0,label="Bistable")
    if I1 != I2
        plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]))
    else
        plot!(xticks=([xlims_[1],I1,xlims_[2]],[xlims_[1],L"I_1=I_2",xlims_[2]]))
    end
    plot!([I1,I2],220 .*ones(2),arrow=true,color=:black,linewidth=0.5,label=:none)
    plot!([I2,I1],220 .*ones(2),arrow=true,color=:black,linewidth=0.5,label=:none,ylims=(-10,250))
    annotate!(mean([I2,I1]),250,L"\Delta I",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    plot!(yticks=[0,250])

    VI_fill_ = collect(1:length(I_stable_VI))[I_stable_VI .>=I1]
    VI_fill = VI_fill_[I_stable_VI[VI_fill_] .<= I2]  

    plt_VI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legend=:outerbottomright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-150,-50),yticks=-150:50:-50,xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]),size=(800,400))
    plot!(I_stable_VI ,V_stable_VI,ls=:dash,lw=2,lc=color,label="Resting")
    plot!(I_stable_VI[VI_fill],V_stable_VI[VI_fill],fillrange=ones(size(VI_fill)).*-55,c=color,fillalpha=0.15,lw=0,label="Bistable")
    scatter!(I1.*ones(2), V_I1.*ones(2),mc=color,markerstrokewidth=0,label=:none)
    scatter!(I2.*ones(2), V_I2.*ones(2),mc=color,markerstrokewidth=0,label=:none)    
    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims_)
    return plt
end

function plt_V_mult_steps(t,V,It,lc_V,I_lim,yticks_)
    plt_V_list = []
    fontsize_ = 25
    plt_I = Plots.plot(fontfamily="Computer Modern",grid=:none,legend=:none,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)

    lw_V =2
    write_an = true

    dt_guide = 100
    dV_guide = 50
    dI_guide = 1 #std

    I_max= maximum([It[1][end],It[end][end]])
    I_min= I_lim[1]
    V_min = minimum([V[1][1],V[end][end]])
    V_max = 40

    for i in eachindex(V)
        plt_V = Plots.plot(fontfamily="Computer Modern",grid=:none,legend=:none,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)
        plot!(t[i],V[i],lc=lc_V[i],lw=lw_V) 

            plot!(plt_V,maximum(t[i]) *1.05 * ones(2),[V_min,  V_min+dV_guide ],lc=:black,lw=2,label=:none)
            if write_an == true && i==1
                annotate!(plt_V,maximum(t[i]) *1.21,V_min+dV_guide/2, latexify(dV_guide * u"mV") ,annotationfontsize=fontsize_)
            end

        plot!(xlims=(-maximum(t[i]) *0.03,maximum(t[i]) *1.15) )
        plot!(ylims=(minimum(yticks_),maximum(yticks_)))
        plot!(yticks = yticks_,ylabelfontsize=fontsize_)        
        plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
        plot!(foreground_color_border=:black,showaxis=false,foreground_color_axis=:black)
        push!(plt_V_list,plt_V)

        plot!(plt_I,t[i],It[i],lc=lc_V[i],lw=lw_V,label=L"I(t)")
        if i==1
            plot!(plt_I,maximum(t[i]) *1.05 * ones(2),[I_max,  I_max-dI_guide ],lc=:black,lw=2,label=:none)
            if write_an == true
                annotate!(plt_I,maximum(t[i]) *1.19,I_max -dI_guide/2, latexify(dI_guide * u"µA/cm^2") ,annotationfontsize=fontsize_)
            end
            #plot lower guide
            plot!([(maximum(t[i]))-dt_guide , maximum(t[i])],(I_min) .*ones(2),lc=:black,lw=2)
            if write_an == true
                annotate!(((maximum(t[i]))-dt_guide + maximum(t[i]))/2,(I_min)*1.63, latexify(dt_guide * u"ms"),annotationfontsize=fontsize_)
            end
        end
        #axis attributes 
        plot!(plt_I,xlims=(-maximum(t[i]) *0.03,maximum(t[i]) *1.15) )
        plot!(plt_I,xticks=0:10:maximum(t[i]))
        plot!(plt_I,ylims=I_lim)#(minimum(It),maximum(It)))
        plot!(plt_I,ylabelfontsize=fontsize_) #ytciks=round.([minimum(I),maximum(I)].*100)./100
        plot!(plt_I,grid=false,legend=:none,tickfontfamily="Computer Modern",legendfontfamily="Computer Modern",tickfontsize=fontsize_,legendfontsize=fontsize_)
        plot!(plt_I,foreground_color_border=:black,showaxis=false,foreground_color_axis=:black)
        display(plt_I)
    end
    return plt_V_list,plt_I
end


empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

title_new_CaL_A = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.025,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_new_CaL_B = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.05,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_new_CaL_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.3,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")


plt_fI_VI_CaLs = plt_fI_VI(I_list_KirCaLs_fI[3][1],f_list_KirCaLs_fI[3][1],I_stable_KirCaLs_fI[3][1],V_stable_KirCaLs_fI[3][1],I1_KirCaLs_fI[3,1],I2_KirCaLs_fI[3,1],V_I1_KirCaLs_fI[3,1],V_I2_KirCaLs_fI[3,1],(-3,1),:darkorange2)

plt_V_as_CaL,plt_I_as_CaL = plt_V_mult_steps(reverse(t_asc_step_CaLs),reverse(V_asc_step_CaLs),reverse(I_asc_step_CaLs),[palette(:twelvebitrainbow)[4],palette(:twelvebitrainbow)[3],palette(:twelvebitrainbow)[2],palette(:twelvebitrainbow)[1]],(I_asc_step_CaLs[1][1]-0.8,I_desc_step_CaLs[1][1]),-175:25:50)
plt_VI_as = Plots.plot(title_new_CaL_A,plt_V_as_CaL[1],plt_V_as_CaL[2],plt_V_as_CaL[3],plt_V_as_CaL[4],plt_I_as_CaL,layout=@layout[a{0.001h};b;c;d;e;f{0.05h}],size=(600,1000))
plt_V_des_CaL,plt_I_des_CaL = plt_V_mult_steps(t_desc_step_CaLs,V_desc_step_CaLs,I_desc_step_CaLs,[palette(:twelvebitrainbow)[4],palette(:twelvebitrainbow)[3],palette(:twelvebitrainbow)[2],palette(:twelvebitrainbow)[1]],(I_asc_step_CaLs[1][1]-0.8,I_desc_step_CaLs[1][1]),-175:25:50)
plt_VI_des = Plots.plot(title_new_CaL_B,plt_V_des_CaL[1],plt_V_des_CaL[2],plt_V_des_CaL[3],plt_V_des_CaL[4],plt_I_des_CaL,layout=@layout[a{0.001h};b;c;d;e;f{0.05h}],size=(600,1000))
plt_VI_as_des_CaL = Plots.plot(plt_VI_as,empty_plt,plt_VI_des,layout=@layout[a b{0.001w} c],size=(1200,1000))

plt_fI_VI_CaL = Plots.plot(title_new_CaL_C,plot(plt_fI_VI_CaLs[1]),empty_plt,plot(plt_fI_VI_CaLs[2]),empty_plt,layout=@layout[a{0.001h} ; b ; c{0.03h} ; d ; e{0.005h}])
plt_VI_as_des_fI_VI_CaL = Plots.plot(plt_VI_as_des_CaL,empty_plt,empty_plt,plt_fI_VI_CaL,layout=@layout[a{0.6w} b{0.001w} c{0.001w} d],size=(2300,1200))
plt_new_CaL = Plots.plot(empty_plt,Plots.plot(plt_VI_as_des_fI_VI_CaL,empty_plt,layout=@layout[a b{0.001w}]),empty_plt,layout=@layout[a{0.001h} ; b ; c{0.001h}],size=(2300,1100),fontfamily="Computer Modern")

#savefig(plt_new_CaL,"ChannelUpdate/cluster/figures/pdf/fig-1.pdf")
