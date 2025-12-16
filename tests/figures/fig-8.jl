using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using PlotlyJS, CSV, DataFrames

include("../sim_pulse_step/steady_state_currents.jl")


function plt_dIss_(V,dIss,lc_,ylims_,yname)
    fontsize_=25
    lw_ =2

    plt = Plots.plot()
    Plots.plot!(V,dIss,lc=lc_,lw=lw_)
    Plots.plot!(V[dIss.<0],dIss[dIss.<0],fillrange=zeros(size(V)),c=lc_,fillalpha=0.15,lw=0,label="NFeedback")
    Plots.plot!(ylabel=yname,ylabelfontsize=fontsize_,yticks=[0],ylims=ylims_)
    Plots.plot!(xlabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]",xlabelfontsize=fontsize_,xticks=-120:30:0)
    Plots.plot!(tickfont = "Computer Modern",tickfontsize=fontsize_)
    Plots.plot!(grid=:y,legend=:none)

    return plt
end

function plt_Iss_in_out(V,Iss,Iss_in,lc_Iss_in,yname,ylims_)
    fontsize_ = 25
    lw = 2
    lc_Iss = RGB(0.7,0.7,0.7)

    plt = Plots.plot()
    Plots.plot!(V,zeros(size(V)),lw=1,lc=RGB(0.9,0.9,0.9))
    Plots.plot!(V,Iss,lc=lc_Iss,lw=lw,ls=:dot)
    Plots.plot!(V,Iss_in,lc=lc_Iss_in,lw=lw)
    Plots.plot!(ylabel=yname,yticks=[0],ylabelfontsize=fontsize_,ylims=ylims_)
    Plots.plot!(xlabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]" ,xlabelfontsize=fontsize_,xticks=-120:30:0)
    Plots.plot!(tickfont = "Computer Modern",tickfontsize=fontsize_)
    Plots.plot!(grid=:y,legend=:none)
    return plt
end

empty_plt = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

title_f5_new_A = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_f5_new_B = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#annotate!(-0.15,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_f5_new_C = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_f5_new_D = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#annotate!(-0.15,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_f5_new_E = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_f5_new_F = Plots.plot() 
Plots.plot!(grid=false,legend=:none)
Plots.plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#annotate!(-0.15,0.5,L"\textbf{F}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")

plt_dIss_(V_ss,dIKir_ss_dV .*0.5,:maroon3,(minimum(dIKir_ss_dV)-0.1,-minimum(dIKir_ss_dV)+0.1),L"\partial I_{\mathrm{Kir}}/\partial V")


#Ajouter la dépendance au Ca ? Surtout le fait que le feedback + se déplace vers V<

plt_dICaL = plt_dIss_(V_ss,dICaLs_ss_dV .*0.6,:darkorange2,(minimum(dICaLs_ss_dV)-0.1,-minimum(dICaLs_ss_dV)+0.1),L"\partial I_{\mathrm{CaL},\infty}/\partial V")
plt_dIKM_ = plt_dIss_(V_ss,dIKM_ss_dV .*0.08,:mediumpurple3,(minimum(dIKM_ss_dV)-0.1,-minimum(dIKM_ss_dV)+0.1),L"\partial I_{\mathrm{KM},\infty}/\partial V")
plt_dIKir = plt_dIss_(V_ss,dIKir_ss_dV .*0.08,:maroon3,(minimum(dIKM_ss_dV)-0.1,-minimum(dIKM_ss_dV)+0.1),L"\partial I_{\mathrm{Kir},\infty}/\partial V")
plt_ICaL_ = plt_Iss_in_out(V_ss,ICaLs_ss.*NaN,ICaLs_ss.*7e-5,:darkorange2,L"I_{\mathrm{CaL},\infty}",(minimum(IKir_ss),-minimum(IKir_ss)))
plt_IKM = plt_Iss_in_out(V_ss,IKM_ss.*NaN,IKM_ss,:mediumpurple3,L"I_{\mathrm{KM},\infty}",(minimum(IKir_ss)-10,-minimum(IKir_ss)+10))
plt_IKir_ = plt_Iss_in_out(V_ss,IKir_ss.*NaN,IKir_ss,:maroon3,L"I_{\mathrm{Kir},\infty}",(minimum(IKir_ss)-10,-minimum(IKir_ss)+10))

plt_Iion_ = Plots.plot(empty_plt,Plots.plot(title_f5_new_A,plt_ICaL_,layout=@layout[a{0.001h};b]),empty_plt,Plots.plot(title_f5_new_C,plt_IKM,layout=@layout[a{0.001h};b]),empty_plt,Plots.plot(title_f5_new_E,plt_IKir_,layout=@layout[a{0.001h};b]),layout=@layout[a{0.001w} b c{0.001w} d e{0.001w} f])
plt_dIion = Plots.plot(empty_plt,Plots.plot(title_f5_new_B,plt_dICaL,layout=@layout[a{0.001h};b]),empty_plt,Plots.plot(title_f5_new_D,plt_dIKM_,layout=@layout[a{0.001h};b]),empty_plt,Plots.plot(title_f5_new_F,plt_dIKir,layout=@layout[a{0.001h};b]),layout=@layout[a{0.001w} b c{0.001w} d e{0.001w} f])
plt_f5_new = Plots.plot(empty_plt,plt_Iion_,empty_plt,plt_dIion,empty_plt,layout=@layout[a{0.001h} ; b ; c{0.003h} ; d ; e{0.003h}],size=(2300,1100))



#Plots.savefig(plt_f5_new,"tests/figures/pdf/fig-6.pdf")

