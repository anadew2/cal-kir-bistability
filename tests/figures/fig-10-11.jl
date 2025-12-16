using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using JLD
using Statistics

include("../sim_pulse_step/sim_switch.jl")

function plt_V_I_derj_nA2(t,V,I,lc_,Iticks,yticks_)
    dt_guide = 1

    lw_V = 1
    lw_I = 3
    fontsize_= 25

    plt_I = Plots.plot()
    #plot!(t,(maximum(I)+minimum(I))/2 .*ones(length(t)),ls=:dash,lc=:black,linealpha=0.8,lw=lw_I,ylims=(minimum(I)-0.001,maximum(I)+0.001))
    plot!(t,I,lc=lc_,lw=lw_I)
    plot!(yticks = Iticks,tickfontsize=fontsize_,tickfontcolor=:black,ylabel=L"I\,~"  *L"[" *latexify(u"µA/cm^2") *L"]" ,ylabelfontsize=fontsize_)
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(xticks=0:round(maximum(t)):round(maximum(t)))
     #ytciks=round.([minimum(I),maximum(I)].*100)./100
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,foreground_color_axis=:black)

    plt_V = Plots.plot()
    plot!(t,V,lc=lc_,lw=lw_V)
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(ylims=(minimum(yticks_),maximum(yticks_)))
    plot!(xticks=0:round(maximum(t)):round(maximum(t)),xlabel=L"t \ \ "  *L"[" *latexify(u"ms") *L"]" ,xlabelfontsize=fontsize_)
    plot!(yticks = yticks_,ylabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]" ,ylabelfontsize=fontsize_)
    
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,foreground_color_axis=:black)

    plt = Plots.plot(plt_I,plt_V,layout=@layout[a{0.35h} ; b],size=(1000,800))
    return plt
end
function plt_V_I_derj_nA2_(t,V,I,lc_,Iticks,yticks_)
    dt_guide = 1

    lw_V = 1
    lw_I = 3
    fontsize_= 25

    plt_I = Plots.plot()
    #plot!(t,(maximum(I)+minimum(I))/2 .*ones(length(t)),ls=:dash,lc=:black,linealpha=0.8,lw=lw_I,ylims=(minimum(I)-0.001,maximum(I)+0.001))
    plot!(t,I,lc=lc_,lw=lw_I)
    plot!(yticks = Iticks,tickfontsize=fontsize_,tickfontcolor=:black,ylabel=L"I\,~"  *L"[" *latexify(u"µA/cm^2") *L"]" ,ylabelfontsize=fontsize_)
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(xticks=0:round(maximum(t)):round(maximum(t)))
     #ytciks=round.([minimum(I),maximum(I)].*100)./100
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,foreground_color_axis=:black)

    plt_V = Plots.plot()
    plot!(t,V,lc=lc_,lw=lw_V)
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(ylims=(minimum(yticks_),maximum(yticks_)))
    plot!(xticks=0:round(maximum(t)):round(maximum(t)),xlabel=L"t \ \ "  *L"[" *latexify(u"ms") *L"]" ,xlabelfontsize=fontsize_)
    plot!(yticks = yticks_,ylabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]" ,ylabelfontsize=fontsize_)
    
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,foreground_color_axis=:black)

    plt = Plots.plot(plt_I,empty_plt,plt_V,layout=@layout[a{0.35h} ; b{0.001h} ; c],size=(1000,800))
    return plt
end
function plt_multiple_SN_eig_plt(find_fixed_points_output,I_cyc_,V_cyc_,I1,I2,xlims_,lc__)

    fontsize_ = 25
    lw_=3

    lc_ = lc__[1]
    lc_cyc_ = lc__[2]
    lc_saddle = colorant"black"
    ylims_ = (-90,45)

    VI_fill_ = collect(1:length(I_cyc_))[I_cyc_ .>=I1]
    VI_fill = VI_fill_[I_cyc_[VI_fill_] .<= I2]  
        
    I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable = find_fixed_points_output

	plt = Plots.plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]" ,legend_position=(0.8,0.75),legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,xlims=xlims_,ylims=ylims_,title=" ",size=(800,400))
    plot!(I_stable ,V_stable ,lc=lc_,lw=lw_,label="Stable node")
	plot!(I_saddle ,V_saddle , lc=lc_saddle, ls=:dash,lw=2,label="Saddle node")
    if length(I_unstable[I_unstable .<= xlims_[2]][I_unstable[I_unstable .<= xlims_[2]] .>= xlims_[1]])>0
        #plot!(I_unstable[V_unstable.<=50],V_unstable[V_unstable.<=50], lc=lc_saddle, ls=:dot,lw=1,label="Unstable node")
    end
    #scatter!(I_SN .*ones(2), V_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
    plot!(I_cyc_,V_cyc_[1],label="Cycle extrema",lc=lc_cyc_,lw=lw_)
    plot!(I_cyc_,V_cyc_[2],label=:none,lc=lc_cyc_,lw=lw_)
    if length(I_cyc_[VI_fill])>1
        plot!(I_cyc_[VI_fill],convert.(Float64,V_cyc_[1][VI_fill]),fillrange=convert.(Float64,V_cyc_[2][VI_fill]),c=lc_cyc_,fillalpha=0.15,lw=0,label=:none)
    end		
    if abs(I2-I1)>=0.1
        plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],round(I1*100)/100,round(I2*100)/100,xlims_[2]]))
        plot!(I1 .*ones(2),[-90,45],lc=:grey75,ls=:solid,label=:none,lw=2)
        annotate!(I1,45+6,L"I_1",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
        plot!(I2 .*ones(2),[-90,45],lc=:grey75,ls=:dot,label=:none,lw=3)
        annotate!(I2,45+6,L"I_2",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    else
        plot!(xticks=([xlims_[1],I1,xlims_[2]],[xlims_[1],round(I1*100)/100,xlims_[2]]))
        plot!(I1.*ones(2),[-90,45],lc=:grey75,ls=:solid,label=:none,lw=2)
        annotate!(I1,45+6,L"I_1=I_2",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    end
    plot!(yticks=-90:30:ylims_[2])
    plot!(grid=false)
	return plt
end

empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

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
title_C_ = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.3,0.05,"Kir + CaL",annotationfontsize=25,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_D_ = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.01-0.22,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.3,0.05,"KM + CaL",annotationfontsize=25,annotationhalign=:left,annotationfontfamily="Computer Modern")


plt_step_hh = plt_V_I_derj_nA2(sol_step_hh.t[1:8913],sol_step_hh[1,:][1:8913],It_step_hh[1:8913],:deepskyblue4,([minimum(It_step_hh[1:8913]),maximum(It_step_hh[1:8913])],[round(minimum(It_step_hh[1:8913])*100)/100,round(maximum(It_step_hh[1:8913])*100)/100]),-90:30:-30)
plt_step_Kir = plt_V_I_derj_nA2(sol_step_Kir.t[1:11483],sol_step_Kir[1,:][1:11483],It_step_Kir[1:11483],:maroon3,([minimum(It_step_Kir[1:11483]),maximum(It_step_Kir[1:11483])],[round(minimum(It_step_Kir[1:11483])*100)/100,round(maximum(It_step_Kir[1:11483])*100)/100]),-90:30:-30)
plt_step_CaLs = plt_V_I_derj_nA2(sol_step_CaLs.t[1:56180],sol_step_CaLs[1,:][1:56180],It_step_CaLs[1:56180],:tan1,([minimum(It_step_CaLs[1:56180]),maximum(It_step_CaLs[1:56180])],[round(minimum(It_step_CaLs[1:56180])*100)/100,round(maximum(It_step_CaLs[1:56180])*100)/100]),-90:30:-30)
plt_step_KM = plt_V_I_derj_nA2(sol_step_KM.t[1:16811],sol_step_KM[1,:][1:16811],It_step_KM[1:16811],:mediumpurple4,([minimum(It_step_KM[1:16811]),maximum(It_step_KM[1:16811])],[round(minimum(It_step_KM[1:16811])*100)/100,round(maximum(It_step_KM[1:16811])*100)/100]),-90:30:-30)
plt_step_KirCaLs = plt_V_I_derj_nA2_(sol_step_KirCaLs_fI_23.t[1:70474],sol_step_KirCaLs_fI_23[1,:][1:70474],It_step_KirCaLs_fI_23[1:70474],:green,([minimum(It_step_KirCaLs_fI_23[1:70474]),maximum(It_step_KirCaLs_fI_23[1:70474])],[round(minimum(It_step_KirCaLs_fI_23[1:70474])*100)/100,round(maximum(It_step_KirCaLs_fI_23[1:70474])*100)/100]),-90:30:-30)
plt_step_KMCaLs = plt_V_I_derj_nA2_(sol_step_KMCaLs_fI_23.t[1:48421],sol_step_KMCaLs_fI_23[1,:][1:48421],It_step_KMCaLs_fI_23[1:48421],:navy,([minimum(It_step_KMCaLs_fI_23[1:48421]),maximum(It_step_KMCaLs_fI_23[1:48421])],[round(minimum(It_step_KMCaLs_fI_23[1:48421])*100)/100,round(maximum(It_step_KMCaLs_fI_23[1:48421])*100)/100]),-90:30:-30)

plt_pulse_lCaL = plt_step_hh
plt_pulse_hCaL =plt_step_CaLs
plt_pulse_lCaL_hCaL = Plots.plot(plt_pulse_lCaL,plt_pulse_hCaL,layout=(1,2))

plt_bif_hh = plt_multiple_SN_eig_plt(fixedpoint_hh,I_list_KMCaLs_fI[1][1],[V_min_c_KMCaLs_fI[1][1],V_max_c_KMCaLs_fI[1][1]],I1_KMCaLs_fI[1,1],I2_KMCaLs_fI[1,1],round.((0.3 .*(-1,1).+I2_KMCaLs_fI[1,1]).*10)./10,[colorant"steelblue1",colorant"dodgerblue4"])
plt_bif_CaLs = plt_multiple_SN_eig_plt(fixedpoint_CaLs,I_list_KMCaLs_fI[2][1],[V_min_c_KMCaLs_fI[2][1],V_max_c_KMCaLs_fI[2][1]],I1_KMCaLs_fI[2,1],I2_KMCaLs_fI[2,1],round.((0.3 .*(-1,1).+mean([I1_KMCaLs_fI[2,1],I2_KMCaLs_fI[2,1]])).*10)./10,[colorant"goldenrod1",colorant"chocolate3"])
plt_bif_Kir = plt_multiple_SN_eig_plt(fixedpoint_Kir,I_list_KirCaLs_fI[1][3],[V_min_c_KirCaLs_fI[1][3],V_max_c_KirCaLs_fI[1][3]],I1_KirCaLs_fI[1,3],I2_KirCaLs_fI[1,3],round.((0.3 .*(-1,1).+I1_KirCaLs_fI[1,3]).*10)./10,[colorant"maroon1",colorant"deeppink4"])
plt_bif_KM = plt_multiple_SN_eig_plt(fixedpoint_KM,I_list_KMCaLs_fI[1][3],[V_min_c_KMCaLs_fI[1][3],V_max_c_KMCaLs_fI[1][3]],I1_KMCaLs_fI[1,3],I2_KMCaLs_fI[1,3],round.((0.3 .*(-1,1).+I2_KMCaLs_fI[1,3]).*10)./10,[colorant"mediumpurple2",colorant"purple4"])
plt_bif_KirCaLs = plt_multiple_SN_eig_plt(fixedpoint_KirCaLs_fI_23,I_list_KirCaLs_fI[2][3],[V_min_c_KirCaLs_fI[2][3],V_max_c_KirCaLs_fI[2][3]],I1_KirCaLs_fI[2,3],I2_KirCaLs_fI[2,3],round.((1 .*(-1,1).+mean([I1_KirCaLs_fI[2,3],I2_KirCaLs_fI[2,3]])).*10)./10,[colorant"lightgreen",colorant"green"])
plt_bif_KMCaLs = plt_multiple_SN_eig_plt(fixedpoint_KMCaLs_fI_23,I_list_KMCaLs_fI[2][3],[V_min_c_KMCaLs_fI[2][3],V_max_c_KMCaLs_fI[2][3]],I1_KMCaLs_fI[2,3],I2_KMCaLs_fI[2,3],round.((1. .*(-1,1).+mean([I1_KMCaLs_fI[2,3],I2_KMCaLs_fI[2,3]])).*10)./10,[colorant"lightskyblue",colorant"navy"])
plt_bif_CaLs_max = plt_multiple_SN_eig_plt(fixedpoint_KirCaLs_fI_31,I_list_KirCaLs_fI[3][1],[V_min_c_KirCaLs_fI[3][1],V_max_c_KirCaLs_fI[3][1]],I1_KirCaLs_fI[3,1],I2_KirCaLs_fI[3,1],round.(((-1,1).+I2_KirCaLs_fI[3,1]).*10)./10,[colorant"grey",colorant"black"])

plt_bif_CaLs = plt_multiple_SN_eig_plt(fixedpoint_CaLs,I_list_KMCaLs_fI[2][1],[V_min_c_KMCaLs_fI[2][1],V_max_c_KMCaLs_fI[2][1]],I1_KMCaLs_fI[2,1],I2_KMCaLs_fI[2,1],round.((0.3 .*(-1,1).+mean([I1_KMCaLs_fI[2,1],I2_KMCaLs_fI[2,1]])).*10)./10,[colorant"goldenrod1",colorant"chocolate3"])


fontsize_ = 25

## Subplot 4(2)x3
plt_steps__= Plots.plot(title_A,plt_step_hh,title_B,plt_step_CaLs,title_C,plt_step_Kir,title_D,plt_step_KM,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.03h}; f; g{0.03h}; h; i{0.003h}],size=(1700,1700))
plt_bif_1__ = Plots.plot(title_E,plt_bif_hh,title_F,plt_bif_CaLs,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.003h}],size=(1700,1700))
plt_bif_2__ = Plots.plot(title_G,plt_bif_Kir,title_H,plt_bif_KM,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.003h}],size=(1700,1700))

plt_full_ = Plots.plot(empty_plt,plt_steps__,empty_plt,plt_bif_1__,empty_plt,plt_bif_2__,empty_plt,layout=@layout[a{0.001w} b{0.25w} c{0.001w} d e{0.1w} f g{0.1w}],size=(2300,1700),fontfamily="Computer Modern")
plt_full = Plots.plot(empty_plt,plt_full_,layout=@layout[a{0.001h}; b],size=(2300,2000),fontfamily="Computer Modern")
#savefig(plt_full,"tests/figures/pdf/fig-8.pdf")


## Subplot 2(1)x3
plt_steps__= Plots.plot(title_A,plt_step_KirCaLs,title_B,plt_step_KMCaLs,empty_plt,layout=@layout[a{0.003h};  b;  c{0.03h}; d; e{0.003h}],size=(1700,1700))
plt_bif_1__ = Plots.plot(title_C_,empty_plt,plt_bif_KirCaLs,empty_plt,layout=@layout[a{0.003h};  b{0.003h};  c ;d{0.003h}],size=(1700,1700))
plt_bif_2__ = Plots.plot(title_D_,empty_plt,plt_bif_KMCaLs,empty_plt,layout=@layout[a{0.003h};   b{0.003h};  c ;d{0.003h}],size=(1700,1700))

plt_full_ = Plots.plot(empty_plt,plt_steps__,empty_plt,plt_bif_1__,empty_plt,plt_bif_2__,empty_plt,layout=@layout[a{0.001w} b{0.25w} c{0.001w} d e{0.1w} f g{0.1w}],size=(2300,1700),fontfamily="Computer Modern")
plt_full = Plots.plot(empty_plt,plt_full_,empty_plt,layout=@layout[a{0.001h}; b; c{0.001h}],size=(2300,1000),fontfamily="Computer Modern")

#Plots.savefig(plt_full,"tests/figures/pdf/fig-9.pdf")
