using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using JLD
using Statistics

include("../sim_pulse_step/sim_switch.jl")

function plt_V_I_derj_nA2(t,V,I,lc_,Iticks,yticks_)
    dt_guide = 1
    dV_guide = 50
    dI_guide = 10.

    lw_V = 1
    lw_I = 3
    fontsize_= 25
    write_an = true

    plt_I = Plots.plot()
    plot!(t,Iticks[1].*ones(length(t)),ls=:dash,lc=:black,linealpha=0.8,lw=lw_I,ylims=(minimum(I)-0.001,maximum(I)+0.001))
    plot!(t,I,lc=lc_,lw=lw_I)
    plot!(yticks = Iticks,tickfontsize=fontsize_,tickfontcolor=:black)
    #plot left guide
    plot!(maximum(t) *1.015 * ones(2),[minimum(I) , minimum(I)+dI_guide/1000],lc=:black,lw=2)
    if write_an == true
        annotate!(maximum(t) *1.175,minimum(I)+dI_guide/1000/2, latexify(Int((10)) * u"nA/cm^2") ,annotationfontsize=fontsize_)
    end
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.15) )
    plot!(xticks=0:10:maximum(t))
     #ytciks=round.([minimum(I),maximum(I)].*100)./100
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,showaxis=:y,foreground_color_axis=:black)

    plt_V = Plots.plot()
    plot!(t,V,lc=lc_,lw=lw_V)
    #plot lower guide
    plot!([(maximum(t))-dt_guide*1000 , maximum(t)],(yticks_[1]) .*ones(2),lc=:black,lw=2)
    if write_an == true
        annotate!(((maximum(t))-dt_guide*1000 + maximum(t))/2,(yticks_[1])*1.1, latexify(Int(dt_guide) * u"s"),annotationfontsize=fontsize_)
    end
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.15) )
    plot!(ylims=(minimum(yticks_),maximum(yticks_)))
    plot!(xticks=0:dt_guide*1000:maximum(t))
    plot!(yticks = yticks_,ylabel=L"V\,~"  *L"(" *latexify(u"mV") *L")" ,ylabelfontsize=fontsize_)
    
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,showaxis=:y,foreground_color_axis=:black)

    plt = Plots.plot(plt_I,plt_V,layout=@layout[a{0.25h} ; b],size=(1000,800))
    return plt
end
function plt_multiple_SN_eig_plt(find_fixed_points_output,I_cyc_,V_cyc_,I1,I2,xlims_,lc_)

    fontsize_ = 25
    lw_=2
    lc_cyc_ = lc_.*0.8
    lc_saddle = colorant"grey56"
    ylims_ = (-75,45)

    VI_fill_ = collect(1:length(I_cyc_))[I_cyc_ .>=I1]
    VI_fill = VI_fill_[I_cyc_[VI_fill_] .<= I2]  
        
    I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable = find_fixed_points_output

	plt = Plots.plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"V\,~"  *L"(" *latexify(u"mV") *L")" ,legend_position=(0.8,0.75),legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,xlims=xlims_,ylims=ylims_,title=" ",size=(800,400))
    plot!(I_stable ,V_stable ,lc=lc_,lw=lw_,label="Stable node")
	plot!(I_saddle ,V_saddle , lc=lc_saddle, ls=:dash,lw=1,label="Saddle node")
    if length(I_unstable[I_unstable .<= xlims_[2]][I_unstable[I_unstable .<= xlims_[2]] .>= xlims_[1]])>0
        plot!(I_unstable[V_unstable.<=50],V_unstable[V_unstable.<=50], lc=lc_saddle, ls=:dot,lw=1,label="Unstable node")
    end
    scatter!(I_SN .*ones(2), V_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
    plot!(I_cyc_,V_cyc_[1],label="Cycle extrema",lc=lc_cyc_,lw=lw_)
    plot!(I_cyc_,V_cyc_[2],label=:none,lc=lc_cyc_,lw=lw_)
    if length(I_cyc_[VI_fill])>1
        plot!(I_cyc_[VI_fill],convert.(Float64,V_cyc_[1][VI_fill]),fillrange=convert.(Float64,V_cyc_[2][VI_fill]),c=lc_,fillalpha=0.15,lw=0,label=:none)
    end		
    if abs(I2-I1)>=0.1
        plot!(xticks=([round(xlims_[1]*10)/10,I1,I2,round(xlims_[2]*10)/10],[round(xlims_[1]*10)/10,L"I_1",L"I_2",round(xlims_[2]*10)/10]))
    else
        plot!(xticks=([round(xlims_[1]*10)/10,I1,round(xlims_[2]*10)/10],[round(xlims_[1]*10)/10,L"I_1=I_2",round(xlims_[2]*10)/10]))
    end
    plot!(yticks=-90:30:ylims_[2])
	return plt
end

empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

plt_step_hh = plt_V_I_derj_nA2(sol_step_hh.t,sol_step_hh[1,:],It_step_hh,:deepskyblue4,([I1_KirCaLs_fI[1,1]],[L"I_1"]),-90:60:-30)
plt_step_Kir = plt_V_I_derj_nA2(sol_step_Kir.t,sol_step_Kir[1,:],It_step_Kir,:palevioletred2,([I1_KirCaLs_fI[1,3]],[L"I_1"]),-90:60:-30)
plt_step_CaLs = plt_V_I_derj_nA2(sol_step_CaLs.t,sol_step_CaLs[1,:],It_step_CaLs,:tan1,([I1_KirCaLs_fI[2,1]],[L"I_1"]),-90:60:-30)
plt_step_KM = plt_V_I_derj_nA2(sol_step_KM.t,sol_step_KM[1,:],It_step_KM,:mediumpurple4,([I1_KMCaLs_fI[1,3]],[L"I_1"]),-90:60:-30)
plt_step_KirCaLs = plt_V_I_derj_nA2(sol_step_KirCaLs_fI_23.t,sol_step_KirCaLs_fI_23[1,:],It_step_KirCaLs_fI_23,:green,([I1_KirCaLs_fI[2,3]],[L"I_1"]),-90:60:-30)
plt_step_KMCaLs = plt_V_I_derj_nA2(sol_step_KMCaLs_fI_23.t,sol_step_KMCaLs_fI_23[1,:],It_step_KMCaLs_fI_23,:navy,([I1_KMCaLs_fI[2,3]],[L"I_1"]),-90:60:-30)

plt_pulse_lCaL = plt_step_hh
plt_pulse_hCaL =plt_step_CaLs
plt_pulse_lCaL_hCaL = Plots.plot(plt_pulse_lCaL,plt_pulse_hCaL,layout=(1,2))

title_A = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_B = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_D = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_E = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{E}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_F = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{F}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_G = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{G}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_H = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{H}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")

plt_bif_hh = plt_multiple_SN_eig_plt(fixedpoint_hh,I_list_KirCaLs_fI[1][1],[V_min_c_KirCaLs_fI[1][1],V_max_c_KirCaLs_fI[1][1]],I1_KirCaLs_fI[1,1],I2_KirCaLs_fI[1,1],0.3 .*(-1,1).+I2_KirCaLs_fI[1,1],colorant"deepskyblue4")
plt_bif_CaLs = plt_multiple_SN_eig_plt(fixedpoint_CaLs,I_list_KirCaLs_fI[2][1],[V_min_c_KirCaLs_fI[2][1],V_max_c_KirCaLs_fI[2][1]],I1_KirCaLs_fI[2,1],I2_KirCaLs_fI[2,1],0.3 .*(-1,1).+mean([I1_KirCaLs_fI[2,1],I2_KirCaLs_fI[2,1]]),colorant"tan1")
plt_bif_Kir = plt_multiple_SN_eig_plt(fixedpoint_Kir,I_list_KirCaLs_fI[1][3],[V_min_c_KirCaLs_fI[1][3],V_max_c_KirCaLs_fI[1][3]],I1_KirCaLs_fI[1,3],I2_KirCaLs_fI[1,3],0.25 .*(-1,1).+I1_KirCaLs_fI[1,3],colorant"palevioletred2")
plt_bif_KM = plt_multiple_SN_eig_plt(fixedpoint_KM,I_list_KMCaLs_fI[1][3],[V_min_c_KMCaLs_fI[1][3],V_max_c_KMCaLs_fI[1][3]],I1_KMCaLs_fI[1,3],I2_KMCaLs_fI[1,3],0.25 .*(-1,1).+I2_KMCaLs_fI[1,3],colorant"mediumpurple3")
plt_bif_KirCaLs = plt_multiple_SN_eig_plt(fixedpoint_KirCaLs_fI_23,I_list_KirCaLs_fI[2][3],[V_min_c_KirCaLs_fI[2][3],V_max_c_KirCaLs_fI[2][3]],I1_KirCaLs_fI[2,3],I2_KirCaLs_fI[2,3],1 .*(-1,1).+mean([I1_KirCaLs_fI[2,3],I2_KirCaLs_fI[2,3]]),colorant"green")
plt_bif_KMCaLs = plt_multiple_SN_eig_plt(fixedpoint_KMCaLs_fI_23,I_list_KMCaLs_fI[2][3],[V_min_c_KMCaLs_fI[2][3],V_max_c_KMCaLs_fI[2][3]],I1_KMCaLs_fI[2,3],I2_KMCaLs_fI[2,3],1. .*(-1,1).+mean([I1_KMCaLs_fI[2,3],I2_KMCaLs_fI[2,3]]),colorant"navy")
plt_bif_CaLs_max = plt_multiple_SN_eig_plt(fixedpoint_KirCaLs_fI_31,I_list_KirCaLs_fI[3][1],[V_min_c_KirCaLs_fI[3][1],V_max_c_KirCaLs_fI[3][1]],I1_KirCaLs_fI[3,1],I2_KirCaLs_fI[3,1],(-1,1).+I2_KirCaLs_fI[3,1],colorant"grey")

fontsize_ = 25

## Subplot 3x4(2)
plt_steps__= Plots.plot(empty_plt,plt_step_hh,empty_plt,plt_step_CaLs,empty_plt,plt_step_Kir,empty_plt,plt_step_KM,empty_plt,layout=@layout[a{0.003w}  b  c{0.03w} d e{0.03w} f g{0.03w} h i{0.003w}],size=(1700,1700))
plt_bif_1__ = Plots.plot(empty_plt,plt_bif_hh,empty_plt,plt_bif_Kir,empty_plt,layout=@layout[a{0.003w}  b  c{0.03w} d e{0.003w}],size=(1700,1700))
plt_bif_2__ = Plots.plot(empty_plt,plt_bif_CaLs,empty_plt,plt_bif_KM,empty_plt,layout=@layout[a{0.003w}  b  c{0.03w} d e{0.003w}],size=(1700,1700))

plt_full_ = Plots.plot(empty_plt,plt_steps__,empty_plt,plt_bif_1__,empty_plt,plt_bif_2__,empty_plt,layout=@layout[a{0.001h};b{0.2h};c{0.001h};d;e{0.001h};f;g{0.001h}],size=(2300,1700),fontfamily="Computer Modern")
plt_full = Plots.plot(plt_full_,empty_plt,layout=@layout[a; b{0.001h}])

## Subplot 4(2)x3
plt_steps__= Plots.plot(title_A,plt_step_hh,title_B,plt_step_CaLs,title_C,plt_step_Kir,title_D,plt_step_KM,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.03h}; f; g{0.03h}; h; i{0.003h}],size=(1700,1700))
plt_bif_1__ = Plots.plot(title_E,plt_bif_hh,title_F,plt_bif_CaLs,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.003h}],size=(1700,1700))
plt_bif_2__ = Plots.plot(title_G,plt_bif_Kir,title_H,plt_bif_KM,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.003h}],size=(1700,1700))

plt_full_ = Plots.plot(empty_plt,plt_steps__,empty_plt,plt_bif_1__,empty_plt,plt_bif_2__,empty_plt,layout=@layout[a{0.001w} b{0.25w} c{0.001w} d e{0.1w} f g{0.1w}],size=(2300,1700),fontfamily="Computer Modern")
plt_full = Plots.plot(plt_full_,empty_plt,layout=@layout[a; b{0.001h}],size=(2300,1500))
#savefig(plt_full,"ChannelUpdate/cluster/figures/pdf/fig-8.pdf")


## Subplot 2(1)x3
plt_steps__= Plots.plot(title_A,plt_step_KirCaLs,title_B,plt_step_KMCaLs,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.003h}],size=(1700,1700))
plt_bif_1__ = Plots.plot(title_C,Plots.plot(plt_bif_KirCaLs,ylims=(-90,45),title="Kir + CaL",titlefontsize=fontsize_),empty_plt,layout=@layout[a{0.003h};  b;  c{0.003h}],size=(1700,1700))
plt_bif_2__ = Plots.plot(title_D,Plots.plot(plt_bif_KMCaLs,ylims=(-90,45),title="KM + CaL",titlefontsize=fontsize_),empty_plt,layout=@layout[a{0.003h};  b;  c{0.003h}],size=(1700,1700))

plt_full_ = Plots.plot(empty_plt,plt_steps__,empty_plt,plt_bif_1__,empty_plt,plt_bif_2__,empty_plt,layout=@layout[a{0.001w} b{0.25w} c{0.001w} d e{0.1w} f g{0.1w}],size=(2300,1700),fontfamily="Computer Modern")
plt_full = Plots.plot(empty_plt,plt_full_,empty_plt,layout=@layout[a{0.001h}; b; c{0.001h}],size=(2300,750),fontfamily="Computer Modern")

#Plots.savefig(plt_full,"ChannelUpdate/cluster/figures/pdf/fig-9.pdf")
