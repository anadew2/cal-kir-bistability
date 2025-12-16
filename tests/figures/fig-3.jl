using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

include("../dI_grad/local_postpro_output.jl")

function plt_V_I_derj(t,V,I,V_ticks_)
    dt_guide = 500
    dV_guide = 50
    dI_guide = 1.

    lc_V = RGB(0.2,0.2,0.2)
    lc_I = :black
    lw_V = 2
    lw_I = 2

    fontsize_ = 15
    write_an = true

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
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=28)
    plot!(foreground_color_border=:black,showaxis=false,foreground_color_axis=:white)

    plt_V = Plots.plot()
    plot!(t,V,lc=lc_V,lw=lw_V,label=:none)
    #plot lower guide
    plot!([(maximum(t))-dt_guide , maximum(t)],(V_ticks_[1]) .*ones(2),lc=:black,lw=2,label=:none)
    if write_an == true
        annotate!(((maximum(t))-dt_guide + maximum(t))/2,(V_ticks_[1])*1.18, latexify(Int(dt_guide) * u"ms"),annotationfontsize=fontsize_)
    end
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.15) )
    plot!(ylims=(minimum(V_ticks_),maximum(V_ticks_)))
    plot!(xticks=0:dt_guide:maximum(t))
    plot!(yticks = V_ticks_,ylabel = L"V\,~"  *L"(" *latexify(u"mV") *L")",ylabelfontsize=fontsize_)
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=28)
    plot!(foreground_color_border=:black,showaxis=:y,foreground_color_axis=:black)

    plt = Plots.plot(plt_I,plt_V,layout=@layout[a{0.15h} ; b],size=(1000,800))
    return plt
end

function plt_dI_grad(dI_grad_mat,gx,gy,clims_dI,xname,yname,palette)
    fontsize_ =25

    plt_dI = Plots.heatmap(gx,gy,dI_grad_mat,fill=true,lw=0,c=palette,clims=clims_dI)
    #plot!(plt_dI,gx_TC,gy_TC,label=:none,lc=:black,lw=2)
    plot!(colorbar=true,legend=:none)
    yaxis!((minimum(gy),maximum(gy)),yticks=[minimum(gy),maximum(gy)])
    xaxis!((minimum(gx),maximum(gx)),xticks=[minimum(gx),maximum(gx)])
    plot!(ylabel=yname,ylabelfontsize=fontsize_)
    plot!(xlabel=xname,xlabelfontsize=fontsize_)
    plot!(fontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)
    
    plt_annotation = Plots.plot()
    annotate!(0.5,0.57, latexstring("     ",L"\Delta I"),annotationfontsize=fontsize_)
    annotate!(0.5,0.45, latexstring("     ",L"[",latexify(u"µA/cm^2"),L"]"),annotationfontsize=fontsize_)
    plot!(yticks =:none)
    plot!(xticks =:none)
    plot!(grid=false,legend=:none)
    plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

    plt = Plots.plot(plt_dI,plt_annotation,layout=@layout[a b{0.08w}])
    return plt
end
function plt_I1_grad(dI_grad_mat,gx,gy,clims_dI,xname,yname,palette)
    fontsize_ =25

    plt_dI = Plots.heatmap(gx,gy,dI_grad_mat,fill=true,lw=0,c=palette,clims=clims_dI)
    #plot!(plt_dI,gx_TC,gy_TC,label=:none,lc=:black,lw=2)
    plot!(colorbar=true,legend=:none)
    yaxis!((minimum(gy),maximum(gy)),yticks=[minimum(gy),maximum(gy)])
    xaxis!((minimum(gx),maximum(gx)),xticks=[minimum(gx),maximum(gx)])
    plot!(ylabel=yname,ylabelfontsize=fontsize_)
    plot!(xlabel=xname,xlabelfontsize=fontsize_)
    plot!(fontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)
    
    plt_annotation = Plots.plot()
    annotate!(0.5,0.57, latexstring("     ",L"I_1"),annotationfontsize=fontsize_)
    annotate!(0.5,0.45, latexstring("     ",L"[",latexify(u"µA/cm^2"),L"]"),annotationfontsize=fontsize_)
    plot!(yticks =:none)
    plot!(xticks =:none)
    plot!(grid=false,legend=:none)
    plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

    plt = Plots.plot(plt_dI,plt_annotation,layout=@layout[a b{0.08w}])
    return plt
end
function plt_Veq_grad(Veq_grad_mat,gx,gy,clims_V,xname,yname)
    palette = cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true) #cgrad(:BuPu,rev=true)
    fontsize_ = 25

    plt_dI = Plots.heatmap(gx,gy,Veq_grad_mat,fill=true,lw=0,c=palette,clims=clims_V)
    #plot!(plt_dI,gx_TC,gy_TC,label=:none,lc=:black,lw=2)
    plot!(colorbar=true,legend=:none)
    yaxis!((minimum(gy),maximum(gy)),yticks=[minimum(gy),maximum(gy)])
    xaxis!((minimum(gx),maximum(gx)),xticks=[minimum(gx),maximum(gx)])
    plot!(ylabel=yname,ylabelfontsize=fontsize_)
    plot!(xlabel=xname,xlabelfontsize=fontsize_)
    plot!(tickfont = "Computer Modern",tickfontsize=fontsize_)
    plot!(legend=:none)
    
    plt_annotation = Plots.plot()
    annotate!(0.5,0.57, latexstring("     ",L"\overline{V}(I_1)"),annotationfontsize=fontsize_)
    annotate!(0.5,0.45, latexstring("     ",L"[",latexify(u"mV"),L"]"),annotationfontsize=fontsize_)
    plot!(yticks =:none)
    plot!(xticks =:none)
    plot!(grid=false,legend=:none)
    plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

    plt = Plots.plot(plt_dI,plt_annotation,layout=@layout[a b{0.08w}])
    return plt
end


empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

title_kir_km_2_A = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_B = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#annotate!(-0.25,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.1,0.5,L"\mathrm{CaL} + \mathrm{Kir} ",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.1,0.5,L"\mathrm{CaL} + \mathrm{KM} ",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_D = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_E = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#annotate!(-0.25,0.5,L"\textbf{E}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_F = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#annotate!(-0.25,0.5,L"\textbf{F}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")



plt_dI_KMCaLs = plt_dI_grad(dI_KMCaLs,gKM_range_KMCaLs,pCaLs_range_KMCaLs,(0,3.),L"\bar{g}_{\mathrm{KM}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]",cgrad(:seaborn_icefire_gradient ,rev=false))
plt_I1_KMCaLs = plt_I1_grad(I1_KMCaLs,gKM_range_KMCaLs,pCaLs_range_KMCaLs,(-2.5,6.5),L"\bar{g}_{\mathrm{KM}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]",cgrad(:linear_worb_100_25_c53_n256 ,rev=true))

plt_V_I1_KMCaLs = plt_Veq_grad(V_I1_KMCaLs,gKM_range_KMCaLs,pCaLs_range_KMCaLs,(-150,-55),L"\bar{g}_{\mathrm{KM}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]")

plt_dI_KirCaLs = plt_dI_grad(dI_KirCaLs,gKir_range_KirCaLs,pCaLs_range_KirCaLs,(0,3.),L"\bar{g}_{\mathrm{Kir}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]",cgrad(:seaborn_icefire_gradient ,rev=false))
plt_I1_KirCaLs = plt_I1_grad(I1_KirCaLs,gKir_range_KirCaLs,pCaLs_range_KirCaLs,(-2.5,6.5),L"\bar{g}_{\mathrm{Kir}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]",cgrad(:linear_worb_100_25_c53_n256 ,rev=true))

plt_V_I1_KirCaLs = plt_Veq_grad(V_I1_KirCaLs,gKir_range_KirCaLs,pCaLs_range_KirCaLs,(-150,-55),L"\bar{g}_{\mathrm{Kir}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]")

plt_kirgrad = Plots.plot(title_kir_km_2_B,plt_dI_KirCaLs,title_kir_km_2_E,plt_I1_KirCaLs,title_kir_km_2_F,plt_V_I1_KirCaLs,layout=@layout[a{0.001h} ; b ; c{0.001h} ; d; e{0.001h} ; f],size=(1100,1700))
plt_kmgrad = Plots.plot(title_kir_km_2_C,plt_dI_KMCaLs,title_kir_km_2_D,plt_I1_KMCaLs,title_kir_km_2_A,plt_V_I1_KMCaLs,layout=@layout[a{0.001h} ; b ; c{0.001h} ; d; e{0.001h} ; f],size=(1100,1700))
plt_kir_km_2_ = Plots.plot(empty_plt,plt_kmgrad,empty_plt,plt_kirgrad,empty_plt,layout=@layout[a{0.001w} b c{0.07w} d e{0.001w}],size=(2300,1800),fontfamily="Computer Modern")
plt_kir_km_2 = Plots.plot(empty_plt,plt_kir_km_2_,empty_plt,layout=@layout[a{0.001h} ; b ; c{0.001h}],fontfamily="Computer Modern")

#Plots.savefig(plt_kir_km_2,"tests/figures/pdf/fig-3.pdf")
