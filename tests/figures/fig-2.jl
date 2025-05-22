using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

include("../sim_pulse_step/sim_Kir.jl")
include("../sim_pulse_step/sim_KM.jl")
include("../fI_curve/bif_kir.jl")
include("../fI_curve/bif_km.jl")

function plt_fI_VI(I_list_fI,f_list_fI,I_stable_VI,V_stable_VI,I1,I2,V_I1,V_I2,xlims_,color)
    fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
    fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]
    fontsize_= 25

    plt_fI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legend=:outertopright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,size=(800,400))
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

    plt_VI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legend=:outertopright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-150,-50),yticks=-150:50:-50,xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]),size=(800,400))
    plot!(I_stable_VI ,V_stable_VI,ls=:dash,lw=2,lc=color,label="Resting")
    plot!(I_stable_VI[VI_fill],V_stable_VI[VI_fill],fillrange=ones(size(VI_fill)).*-55,c=color,fillalpha=0.15,lw=0,label="Bistable")
    scatter!(I1.*ones(2), V_I1.*ones(2),mc=color,markerstrokewidth=0,label=:none)
    scatter!(I2.*ones(2), V_I2.*ones(2),mc=color,markerstrokewidth=0,label=:none)    
    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims_)
    return plt
end

function plt_V_I_derj(t,V,I,V_ticks_)
    dt_guide = 1000
    dV_guide = 50
    dI_guide = 1.

    lc_V = RGB(0.2,0.2,0.2)
    lc_I = :black
    lw_V = 2
    lw_I = 2

    fontsize_ = 25
    write_an = true

    V_lim_l = minimum(V_ticks_)

    plt_I = Plots.plot()
    plot!(t,I,lc=lc_I,lw=lw_I)
    #plot left guide
    plot!(maximum(t) *1.015 * ones(2),[minimum(I)+0.1 , minimum(I)+0.1+dI_guide],lc=:black,lw=2)
    if write_an == true
        annotate!(maximum(t) *1.15,minimum(I)+0.1 +dI_guide/2, latexify(Int((dI_guide)) * u"µA/cm^2") ,annotationfontsize=fontsize_)
    end
    ==#
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.15) )
    plot!(xticks=0:10:maximum(t))
    plot!(yticks = ([-maximum(t)*0.03,maximum(t)*1.15], [L"~ ~ ~ ~ ~ ~ ~ ~",L"b"]),ylabel=L"I\,~"  ,ylabelfontsize=fontsize_) #ytciks=round.([minimum(I),maximum(I)].*100)./100
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,showaxis=false,foreground_color_axis=:white)

    plt_V = Plots.plot()
    plot!(t,V,lc=lc_V,lw=lw_V,label=:none)

    #plot lower guide
    plot!([(maximum(t))-dt_guide , maximum(t)],(V_lim_l) .*ones(2),lc=:black,lw=2,label=:none)
    if write_an == true
        annotate!(((maximum(t))-dt_guide + maximum(t))/2,(V_lim_l)*1.18, latexify(Int(dt_guide/1000) * u"s"),annotationfontsize=fontsize_)
    end
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.15) )
    plot!(ylims=(minimum(V_ticks_),maximum(V_ticks_)))
    plot!(xticks=0:dt_guide:maximum(t))
    plot!(yticks = V_ticks_,ylabel = L"V\,~"  *L"(" *latexify(u"mV") *L")",ylabelfontsize=fontsize_)
    
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,showaxis=:y,foreground_color_axis=:black)

    plt = Plots.plot(plt_I,plt_V,layout=@layout[a{0.15h} ; b],size=(1000,800))
    return plt
end

empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

title_kir_km_1_A = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.35-0.,0.5,L"\bar{g}_{\mathrm{KM}}=0.15 \,~"*latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_1_B = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_1_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_1_D = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.35-0.,0.5,L"\bar{g}_{\mathrm{Kir}}=0.15 \,~"*latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_1_E = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{E}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_1_F = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{F}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")


plt_fI_VI_Kir = plt_fI_VI(I_list_KirCaLs_fI[3][2],f_list_KirCaLs_fI[3][2],I_stable_KirCaLs_fI[3][2],V_stable_KirCaLs_fI[3][2],I1_KirCaLs_fI[3,2],I2_KirCaLs_fI[3,2],V_I1_KirCaLs_fI[3,2],V_I2_KirCaLs_fI[3,2],(-2,2),:palevioletred2)
plt_fI_VI_KM = plt_fI_VI(I_list_KMCaLs_fI[3][2],f_list_KMCaLs_fI[3][2],I_stable_KMCaLs_fI[3][2],V_stable_KMCaLs_fI[3][2],I1_KMCaLs_fI[3,2],I2_KMCaLs_fI[3,2],V_I1_KMCaLs_fI[3,2],V_I2_KMCaLs_fI[3,2],(0,4),:mediumpurple3)

plt_pulse_kir = plt_V_I_derj(sol_pulse_gKir_fI_2.t,sol_pulse_gKir_fI_2[1,:],It_pulse_gKir_fI_2,50:-50:-150)
plt_pulse_km = plt_V_I_derj(sol_pulse_gKM_fI_2.t,sol_pulse_gKM_fI_2[1,:],It_pulse_gKM_fI_2,50:-50:-150)

plt_pulse_fI_VI_km = Plots.plot(Plots.plot!(title_kir_km_1_D,plt_pulse_kir,layout=@layout[a{0.0003h};b]),title_kir_km_1_E,plot(plt_fI_VI_Kir[1]),title_kir_km_1_F,plot(plt_fI_VI_Kir[2]),empty_plt,layout=@layout[a ; b{0.01h} ; c ; d{0.01h} ; e ; f{0.01h}],size=(1100,1700))
plt_pulse_fI_VI_kir = Plots.plot(Plots.plot!(title_kir_km_1_A,plt_pulse_km,layout=@layout[a{0.003h};b]),title_kir_km_1_B,plot(plt_fI_VI_KM[1]),title_kir_km_1_C,plot(plt_fI_VI_KM[2]),empty_plt,layout=@layout[a ; b{0.01h} ; c ; d{0.01h} ; e ; f{0.01h}],size=(1100,1700))
plt_kir_km_1 = Plots.plot(empty_plt,Plots.plot(empty_plt,plt_pulse_fI_VI_kir,empty_plt,plt_pulse_fI_VI_km,empty_plt,layout=@layout[a{0.075w} b c{0.1w} d e{0.075w}],size=(2300,1700)),layout=@layout[a{0.001h};b],fontfamily="Computer Modern")

#savefig(plt_kir_km_1,"ChannelUpdate/cluster/figures/pdf/fig-2.pdf")
