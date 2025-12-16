using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using PlotlyJS, CSV, DataFrames
using JLD

include("../sim_pulse_step/sim_inward.jl");
include("../sim_pulse_step/sim_outward.jl");
include("../sim_pulse_step/steady_state_currents.jl")


function plt_kir_in_out_Ibist(g,I1,I2,lc_,xname)
    lim_I = (1,3)
    leg_pos = :none
    lw_ = 2
    fontsize_ = 25

    plt = Plots.plot()
    plot!(g,I1,lc=lc_,lw=lw_,label=:none)
    plot!(g,I1,lc=lc_,lw=2,label=L"I_1")
    #plot!(g,I2,lc=lc_,lw=lw_,ls=:dash,label=:none)
    #plot!(g,NaN.*I2,lc=lc_,lw=1.5,ls=:dash,label=L"I_2")
    plot!(g,I1,fillrange=convert.(Float64,I2),c=lc_,fillalpha=0.05,lw=1.,linealpha=0.,lc=lc_,label=:none)
    #plot!(ylabel=L"I\,~"  *L"[" *latexify(u"µA/cm^2") *L"]",ylabelfontsize=fontsize_,ylims=lim_I,yticks=[ceil(lim_I[1]),0,floor(lim_I[2])])
    plot!(ylabel=L"\Delta I\,~"  *L"[" *latexify(u"µA/cm^2") *L"]",ylabelfontsize=fontsize_,ylims=lim_I,yticks=[ceil(lim_I[1]),0,floor(lim_I[2])])
    plot!(xlabel=xname,xlabelfontsize=fontsize_,xticks=[minimum(g),maximum(g)])
    plot!(tickfont = "Computer Modern",tickfontsize=fontsize_)
    plot!(grid=false,legend=leg_pos,legendfontfamily="Computer Modern",tickfontfamily="Computer Modern",legendfontsize=fontsize_,foreground_color_legend=nothing,background_color_legend=nothing)
    return plt
end

function plt_Iss_in_out(V,Iss,Iss_in,lc_Iss_in,yname,ylims_)
    fontsize_ = 25
    lw = 2
    lc_Iss = RGB(0.7,0.7,0.7)

    plt = Plots.plot()
    plot!(V,zeros(size(V)),lw=1,lc=RGB(0.9,0.9,0.9))
    plot!(V,Iss,lc=lc_Iss,lw=lw,ls=:dot)
    plot!(V,Iss_in,lc=lc_Iss_in,lw=lw)
    plot!(ylabel=yname,yticks=[0],ylabelfontsize=fontsize_,ylims=ylims_)
    plot!(xlabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]" ,xlabelfontsize=fontsize_,xticks=-120:30:0)
    plot!(tickfont = "Computer Modern",tickfontsize=fontsize_)
    plot!(grid=:y,legend=:none)
    return plt
end

function plt_step_response_cgrad_yaxis(t,V,I,yticks_,palette)
    dt_guide = 200

    lw_V = 2
    fontsize_ = 25

    V_lim_h = maximum(yticks_)+1
    V_lim_l = minimum(yticks_)-1

    I0 = []
    for i in eachindex(I)
        push!(I0,I[i][end])
    end
    I_lim_h = maximum(maximum(I))+0.05
    I_lim_l = minimum(I0)-0.05
    
    plt_I = Plots.plot()
    plot!(plt_I,t,I,line_z=I0',color=cgrad(palette),lw=lw_V)
    plot!(plt_I,colorbar_title=L"I")
    plot!(plt_I,colorbar_ticks=false,colorbar_titlefontrotation=90.,colorbar_titlefontsize=fontsize_,colorbar_fontfamily="Computer Modern")
    plot!(colorbar=false)
    plot!(xticks=[t[1][1],t[1][end]])
    plot!(ylims=(I_lim_l,I_lim_h),ylabel=L"I\,~"  *L"[" *latexify(u"µA/cm^2") *L"]" )
    plot!(yticks = ([minimum(I0),maximum(I0)],[round(minimum(I0)*10)/10,round(maximum(I0)*10)/10]),ylabelfontsize=fontsize_)
    plot!(grid=false,legend=:none)
    plot!(foreground_color_border=:black,foreground_color_axis=:black,tickfont = "Computer Modern",tickfontsize=fontsize_)

    plt_V = Plots.plot()
    plot!(plt_V,t,V,line_z=I0',color=cgrad(palette),lw=lw_V)
    plot!(plt_V,colorbar_title=L"I")
    plot!(plt_V,colorbar_ticks=false,colorbar_titlefontrotation=90.,colorbar_titlefontsize=fontsize_,colorbar_fontfamily="Computer Modern")
    plot!(colorbar=false)
    #axis attributes 
    plot!(xlims=(-maximum(maximum(t))*0.03,maximum(maximum(t))*1.0) )
    plot!(ylims=(V_lim_l,V_lim_h))
    plot!(draw_arrow=true)
    plot!(xticks=[t[1][1],t[1][end]],xlabel=L"t\,~"  *L"[" *latexify(u"ms") *L"]",xlabelfontsize=fontsize_)
    plot!(yticks = yticks_,ylabel=L"V\,~"  *L"[" *latexify(u"mV") *L"]" ,ylabelfontsize=fontsize_)
    plot!(grid=false,legend=:none)
    plot!(foreground_color_border=:black,foreground_color_axis=:black,tickfont = "Computer Modern",tickfontsize=fontsize_)

    plt_IV = Plots.plot(plt_I,plt_V,layout=@layout[a{0.25h} ;b])
    return plt_IV
end



palette_step_kir_inw = palette([:black,palette(:Paired_8)[2],palette(:Paired_8)[1]],length(It_desc_step_Kir_inw))
palette_step_kir_outw = palette([:black,palette(:Paired_8)[4],palette(:Paired_8)[3]],length(It_desc_step_Kir_outw))
palette_step_km_inw = palette([:black,palette(:Paired_10)[10],palette(:Paired_10)[9]],length(It_desc_step_KM_inw))
palette_step_km_outw = palette([:black,palette(:Paired_8)[6],palette(:Paired_8)[5]],length(It_desc_step_KM_outw))

empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)


title_KM_Kir_inw_outw_full_hz_1 = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_KM_Kir_inw_outw_full_hz_2 = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_KM_Kir_inw_outw_full_hz_3  = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_KM_Kir_inw_outw_full_hz_4  = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.15,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_step_km_inw  = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(0.5,0.,L"\bar{g}_{\mathrm{KM,inward}}=0.15 "*latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
title_step_km_outw  = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(0.5,0.,L"\bar{g}_{\mathrm{KM,outward}}=0.15 "*latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
title_step_kir_inw  = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(0.5,0.,L"\bar{g}_{\mathrm{Kir,inward}}=0.15 "*latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
title_step_kir_outw  = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(0.5,0.,L"\bar{g}_{\mathrm{Kir,outward}}=0.15 " *latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")


plt_outw_I1I2_ = plt_kir_in_out_Ibist(gKir_KirCaLs_outw[1:31],dI_Kir_outw[1:31],dI_Kir_outw[1:31],palette(:Paired_8)[4],L"\bar{g}_{\mathrm{Kir,outward}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]")
plt_inw_I1I2_ = plt_kir_in_out_Ibist(gKir_KirCaLs_inw[1:31],dI_Kir_inw[1:31],dI_Kir_inw[1:31],palette(:Paired_8)[2],L"\bar{g}_{\mathrm{Kir,inward}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]")
plt_IKir_inw_ = plt_Iss_in_out(V_ss,IKir_ss,IKir_inw_ss,palette(:Paired_8)[2],L"I_{\mathrm{Kir,inward},\infty}",(minimum(IKir_ss),maximum(IKM_ss)-20))
plt_IKir_outw_ = plt_Iss_in_out(V_ss,IKir_ss,IKir_outw_ss,palette(:Paired_8)[4],L"I_{\mathrm{Kir,outward},\infty}",(minimum(IKir_ss),maximum(IKM_ss)-20))
plt_step_kir_in_ = plt_step_response_cgrad_yaxis(t_desc_step_Kir_inw,V_desc_step_Kir_inw,It_desc_step_Kir_inw,-50:-20:-130,palette_step_kir_inw)
plt_step_kir_in = Plots.plot(title_step_km_inw,empty_plt,plt_step_kir_in_,layout=@layout[a{0.001h};b{0.003h};c])
plt_step_kir_out_ = plt_step_response_cgrad_yaxis(t_desc_step_Kir_outw,V_desc_step_Kir_outw,It_desc_step_Kir_outw,-50:-20:-130,palette_step_kir_outw)
plt_step_kir_out = Plots.plot(title_step_km_outw,empty_plt,plt_step_kir_out_,layout=@layout[a{0.001h};b{0.003h};c])


plt_KM_inw_I1I2_ = plt_kir_in_out_Ibist(gKM_KMCaLs_inw[1:31],dI_KM_inw[1:31],dI_KM_inw[1:31],palette(:Paired_10)[10],L"\bar{g}_{\mathrm{KM,inward}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]")
plt_KM_outw_I1I2_ = plt_kir_in_out_Ibist(gKM_KMCaLs_outw[1:31],dI_KM_outw[1:31],dI_KM_outw[1:31],palette(:Paired_8)[6],L"\bar{g}_{\mathrm{KM,outward}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]")
plt_IKM_inw_ = plt_Iss_in_out(V_ss,IKM_ss,IKM_inw_ss,palette(:Paired_10)[10],L"I_{\mathrm{KM,inward},\infty}",(minimum(IKir_ss),maximum(IKM_ss)-20))
plt_IKM_outw_ = plt_Iss_in_out(V_ss,IKM_ss,IKM_outw_ss,palette(:Paired_8)[6],L"I_{\mathrm{KM,outward},\infty}",(minimum(IKir_ss),maximum(IKM_ss)-20))
plt_step_km_in_ = plt_step_response_cgrad_yaxis(t_desc_step_KM_inw,V_desc_step_KM_inw,It_desc_step_KM_inw,-50:-20:-130,palette_step_km_inw)
plt_step_km_in = Plots.plot(title_step_km_inw,empty_plt,plt_step_km_in_,layout=@layout[a{0.001h};b{0.003h};c])
plt_step_km_out_ = plt_step_response_cgrad_yaxis(t_desc_step_KM_outw,V_desc_step_KM_outw,It_desc_step_KM_outw,-50:-20:-130,palette_step_km_outw)
plt_step_km_out = Plots.plot(title_step_km_outw,empty_plt,plt_step_km_out_,layout=@layout[a{0.001h};b{0.003h};c])



plt_inw_I1I2__ = Plots.plot(empty_plt,empty_plt,plt_inw_I1I2_,layout=@layout[a{0.001h};b{0.003h};c])
plt_outw_I1I2__ = Plots.plot(empty_plt,empty_plt,plt_outw_I1I2_,layout=@layout[a{0.001h};b{0.003h};c])
plt_IKir_inw = Plots.plot(title_KM_Kir_inw_outw_full_hz_3,empty_plt,plt_IKir_inw_,layout=@layout[a{0.001h};b{0.003h};c])
plt_IKir_outw = Plots.plot(title_KM_Kir_inw_outw_full_hz_4,empty_plt,plt_IKir_outw_,layout=@layout[a{0.001h};b{0.003h};c])
plt_step_kir_in = Plots.plot(title_step_kir_inw,plt_step_kir_in_,layout=@layout[a{0.001h};b])
plt_step_kir_out = Plots.plot(title_step_km_outw,plt_step_kir_out_,layout=@layout[a{0.001h};b])

plt_KM_inw_I1I2__ = Plots.plot(empty_plt,empty_plt,plt_KM_inw_I1I2_,layout=@layout[a{0.001h};b{0.003h};c])
plt_KM_outw_I1I2__ = Plots.plot(empty_plt,empty_plt,plt_KM_outw_I1I2_,layout=@layout[a{0.001h};b{0.003h};c])
plt_IKM_inw = Plots.plot(title_KM_Kir_inw_outw_full_hz_1,empty_plt,plt_IKM_inw_,layout=@layout[a{0.001h};b{0.003h};c])
plt_IKM_outw = Plots.plot(title_KM_Kir_inw_outw_full_hz_2,empty_plt,plt_IKM_outw_,layout=@layout[a{0.001h};b{0.003h};c])
plt_step_km_in = Plots.plot(title_step_km_inw,plt_step_km_in_,layout=@layout[a{0.001h};b])
plt_step_km_out = Plots.plot(title_step_km_outw,plt_step_km_out_,layout=@layout[a{0.001h};b])


plt_KM_Kir_inw_outw_full_hz_1 = Plots.plot(plt_IKM_inw,empty_plt,plt_KM_inw_I1I2__,empty_plt,plt_step_km_in,empty_plt,layout=@layout[ a b{0.001w} c d{0.001w} e f{0.003w} ],size=(2300,1100/2))
plt_KM_Kir_inw_outw_full_hz_2= Plots.plot(plt_IKM_outw,empty_plt,plt_KM_outw_I1I2__,empty_plt,plt_step_km_out,empty_plt,layout=@layout[ a b{0.001w} c d{0.001w} e f{0.003w} ],size=(2300,1100/2))
plt_KM_Kir_inw_outw_full_hz_3 = Plots.plot(plt_IKir_inw,empty_plt,plt_inw_I1I2__,empty_plt,plt_step_kir_in,empty_plt,layout=@layout[ a b{0.001w} c d{0.001w} e f{0.003w} ],size=(2300,1100/2))
plt_KM_Kir_inw_outw_full_hz_4 = Plots.plot(plt_IKir_outw,empty_plt,plt_outw_I1I2__,empty_plt,plt_step_kir_out,empty_plt,layout=@layout[ a b{0.001w} c d{0.001w} e f{0.003w} ],size=(2300,1100/2))
plt_KM_Kir_inw_outw_full_hz = Plots.plot(empty_plt,plt_KM_Kir_inw_outw_full_hz_1,empty_plt,plt_KM_Kir_inw_outw_full_hz_2,empty_plt,plt_KM_Kir_inw_outw_full_hz_3,empty_plt,plt_KM_Kir_inw_outw_full_hz_4,empty_plt,layout=@layout[a{0.003h};b;c{0.001h};d;e{0.001h};f;g{0.001h};h;i{0.001h}],size=(2300,2300))
plt_KM_Kir_inw_outw_full_hz = Plots.plot(empty_plt,Plots.plot(empty_plt,plt_KM_Kir_inw_outw_full_hz_1,plt_KM_Kir_inw_outw_full_hz_2,plt_KM_Kir_inw_outw_full_hz_3,plt_KM_Kir_inw_outw_full_hz_4,empty_plt,layout=@layout[a{0.003h};b;c;d;e;f{0.001h}],size=(2300,2300)),layout=@layout[a{0.001w} b ])

#savefig(plt_KM_Kir_inw_outw_full_hz,"tests/figures/pdf/fig-7.pdf")
