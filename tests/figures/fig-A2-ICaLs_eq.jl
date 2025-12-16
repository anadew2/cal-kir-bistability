using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

function plt_dI_1g(gx,dI,index_,label_,xlims_,xlabel_,title_,palette_)
    ylims_ = (-0.1,3.1)

    plt = Plots.plot(fontfamily="Computer Modern",xlabel=xlabel_,ylabel=L"\Delta I \ \ "*L"["*latexify(u"µA/cm^2")*L"]",title=title_,ylims=ylims_,xlims=xlims_,legend=:outertopright,palette=palette_)
    for ind in index_
        Plots.plot!(plt,gx,dI[:,ind],label=label_[ind],marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
    end
    return plt
end

# Palettes
pal_kir = []
push!(pal_kir,:black)
for i=1:9
    push!(pal_kir,palette(:PuRd_9,rev=true)[i])
end

pal_km = []
push!(pal_km,:black)
for i=2:9
    push!(pal_km,palette(:roma10,rev=true)[i])
end


## With GHK and m^1
include("../dI_grad/local_postpro_output_larger_range.jl")

    ## Kir
    index_kir_ghk_m = vcat(1:6)
    label_kir_ghk_m = []
    for i in eachindex(gKir_range_KirCaLs)
        push!(label_kir_ghk_m,L"\bar{g}_{\mathrm{Kir}}=%$(round(gKir_range_KirCaLs[i]*1e3)*1e-3) ")
    end
    plt_kir_ghk_m_xp = plt_dI_1g(pCaLs_range_KirCaLs,dI_KirCaLs,index_kir_ghk_m,label_kir_ghk_m,(minimum(pCaLs_range_KirCaLs),maximum(pCaLs_range_KirCaLs)),L"\bar{p}_\mathrm{CaL} \ \ "*L"["*latexify(u"cm/s")*L"]",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_kir))

    #gCaLs_eq_KirCaLs = pCaLs_range_KirCaLs .*((1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T)) #mS/cm^2
    #plt_kir_ghk_m_xg = plt_dI_1g(gCaLs_eq_KirCaLs,dI_KirCaLs,index_kir_ghk_m,label_kir_ghk_m,(0,0.5),L"\bar{g}_\mathrm{CaL,eq}",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_kir))

    ## KM
    index_km_ghk_m = vcat(1:6,11:5:16)
    label_km_ghk_m = []
    for i in eachindex(gKM_range_KMCaLs)
        push!(label_km_ghk_m,L"\bar{g}_{\mathrm{KM}}=%$(round(gKM_range_KMCaLs[i]*1e3)*1e-3) ")
    end
    plt_km_ghk_m_xp = plt_dI_1g(pCaLs_range_KMCaLs,dI_KMCaLs,index_km_ghk_m,label_km_ghk_m,(minimum(pCaLs_range_KirCaLs),maximum(pCaLs_range_KirCaLs)),L"\bar{p}_\mathrm{CaL} \ \ "*L"["*latexify(u"cm/s")*L"]",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_km))

    #gCaLs_eq_KMCaLs = pCaLs_range_KMCaLs .*((1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T)) #mS/cm^2
    #plt_km_ghk_m_xg = plt_dI_1g(gCaLs_eq_KMCaLs,dI_KMCaLs,index_km_ghk_m,label_km_ghk_m,(0,0.5),L"\bar{g}_\mathrm{CaL,eq}",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_km))

## With GHK and m^2
include("../dI_grad_m2/local_postpro_output_m2.jl")

    ## Kir
    index_kir_ghk_m2 = vcat(1:6)
    label_kir_ghk_m2 = []
    for i in eachindex(gKir_range_KirCaLs_m2_ghk)
        push!(label_kir_ghk_m2,L"\bar{g}_{\mathrm{Kir}}=%$(round(gKir_range_KirCaLs_m2_ghk[i]*1e3)*1e-3) ")
    end
    plt_kir_ghk_m2_xp = plt_dI_1g(pCaLs_range_KirCaLs_m2_ghk,dI_KirCaLs_m2_ghk,index_kir_ghk_m2,label_kir_ghk_m2,(minimum(pCaLs_range_KirCaLs_m2_ghk),maximum(pCaLs_range_KirCaLs_m2_ghk)),L"\bar{p}_\mathrm{CaL} \ \ "*L"["*latexify(u"cm/s")*L"]",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}^2.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_kir))
    plot!(plt_kir_ghk_m2_xp,xticks=range(0,stop=0.05/1000,length=3))

    #gCaLs_eq_KirCaLs_m2_ghk = pCaLs_range_KirCaLs_m2_ghk .*((1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T)) #mS/cm^2
    #plt_kir_ghk_m2_xg = plt_dI_1g(gCaLs_eq_KirCaLs_m2_ghk,dI_KirCaLs_m2_ghk,index_kir_ghk_m2,label_kir_ghk_m2,(0,1.5),L"\bar{g}_\mathrm{CaL,eq}",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}^2.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_kir))

    ## KM
    index_km_ghk_m2 = vcat(1:6,11:5:16)
    label_km_ghk_m2 = []
    for i in eachindex(gKM_range_KMCaLs_m2_ghk)
        push!(label_km_ghk_m2,L"\bar{g}_{\mathrm{KM}}=%$(round(gKM_range_KMCaLs_m2_ghk[i]*1e3)*1e-3) ")
    end
    plt_km_ghk_m2_xp = plt_dI_1g(pCaLs_range_KMCaLs_m2_ghk,dI_KMCaLs_m2_ghk,index_km_ghk_m2,label_km_ghk_m2,(minimum(pCaLs_range_KMCaLs_m2_ghk),maximum(pCaLs_range_KMCaLs_m2_ghk)),L"\bar{p}_\mathrm{CaL} \ \ "*L"["*latexify(u"cm/s")*L"]",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}^2.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_km))
    plot!(plt_km_ghk_m2_xp,xticks=range(0,stop=0.05/1000,length=3))

    #gCaLs_eq_KMCaLs_m2_ghk = pCaLs_range_KMCaLs_m2_ghk .*((1e-6) *2*F*2*F*Ca_o *(1e3)/(R*T)) #mS/cm^2
    #plt_km_ghk_m2_xg = plt_dI_1g(gCaLs_eq_KMCaLs_m2_ghk,dI_KMCaLs_m2_ghk,index_km_ghk_m2,label_km_ghk_m2,(0,1.5),L"\bar{g}_\mathrm{CaL,eq}",L"I_{\mathrm{CaL}} =\bar{p}_{\mathrm{CaL}}.m_{\mathrm{CaL}}^2.h_{\mathrm{CaL}}.\mathrm{GHK}(V,\left[Ca^{2+} \right]_i) ",palette(pal_km))

## Without GHK and m^1
include("../dI_grad_no_ghk/local_postpro_output_no_ghk.jl")

    ## Kir
    index_kir_no_ghk_m = vcat(1:6)
    label_kir_no_ghk_m = []
    for i in eachindex(gKir_range_KirCaLs_no_ghk)
        push!(label_kir_no_ghk_m,L"\bar{g}_{\mathrm{Kir}}=%$(round(gKir_range_KirCaLs_no_ghk[i]*1e3)*1e-3) ")
    end
    plt_kir_no_ghk_m_xg = plt_dI_1g(gCaLs_range_KirCaLs_no_ghk,dI_KirCaLs_no_ghk,index_kir_no_ghk_m,label_kir_no_ghk_m,(0,0.21),L"\bar{g}_\mathrm{CaL} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"I_{\mathrm{CaL}} =\bar{g}_{\mathrm{CaL}}.m_{\mathrm{CaL}}.h_{\mathrm{CaL}}.(V-E_{\mathrm{Ca}}) ",palette(pal_kir))

    ## KM
    index_km_no_ghk_m = vcat(1:6,11:5:16)
    label_km_no_ghk_m = []
    for i in eachindex(gKM_range_KMCaLs_no_ghk)
        push!(label_km_no_ghk_m,L"\bar{g}_{\mathrm{KM}}=%$(round(gKM_range_KMCaLs_no_ghk[i]*1e3)*1e-3) ")
    end
    plt_km_no_ghk_m_xg = plt_dI_1g(gCaLs_range_KMCaLs_no_ghk,dI_KMCaLs_no_ghk,index_km_no_ghk_m,label_km_no_ghk_m,(0,0.21),L"\bar{g}_\mathrm{CaL} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"I_{\mathrm{CaL}} =\bar{g}_{\mathrm{CaL}}.m_{\mathrm{CaL}}.h_{\mathrm{CaL}}.(V-E_{\mathrm{Ca}}) ",palette(pal_km))

## Without GHK and m^2
include("../dI_grad_no_ghk_m2/local_postpro_output_no_ghk_m2.jl")

    ## Kir
    index_kir_no_ghk_m2 = vcat(1:6)
    label_kir_no_ghk_m2 = []
    for i in eachindex(gKir_range_KirCaLs_no_ghk_m2)
        push!(label_kir_no_ghk_m2,L"\bar{g}_{\mathrm{Kir}}=%$(round(gKir_range_KirCaLs_no_ghk_m2[i]*1e3)*1e-3) ")
    end
    plt_kir_no_ghk_m2_xg = plt_dI_1g(gCaLs_range_KirCaLs_no_ghk_m2,dI_KirCaLs_no_ghk_m2,index_kir_no_ghk_m2,label_kir_no_ghk_m2,(0,0.5),L"\bar{g}_\mathrm{CaL} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"I_{\mathrm{CaL}} =\bar{g}_{\mathrm{CaL}}.m_{\mathrm{CaL}}^2.h_{\mathrm{CaL}}.(V-E_{\mathrm{Ca}}) ",palette(pal_kir))

    ## KM
    index_km_no_ghk_m2 = vcat(1:6,11:5:16)
    label_km_no_ghk_m2 = []
    for i in eachindex(gKM_range_KMCaLs_no_ghk_m2)
        push!(label_km_no_ghk_m2,L"\bar{g}_{\mathrm{KM}}=%$(round(gKM_range_KMCaLs_no_ghk_m2[i]*1e3)*1e-3) ")
    end
    plt_km_no_ghk_m2_xg = plt_dI_1g(gCaLs_range_KMCaLs_no_ghk_m2,dI_KMCaLs_no_ghk_m2,index_km_no_ghk_m2,label_km_no_ghk_m2,(0,0.5),L"\bar{g}_\mathrm{CaL} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"I_{\mathrm{CaL}} =\bar{g}_{\mathrm{CaL}}.m_{\mathrm{CaL}}^2.h_{\mathrm{CaL}}.(V-E_{\mathrm{Ca}}) ",palette(pal_km))


## Altogether

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

fontsize_=25

plt_A1_km_1 = Plots.plot(title_A,plt_km_ghk_m_xp,title_B,plt_km_no_ghk_m_xg,empty_plt,layout=@layout[a{0.003h} ; b ; c{0.003h} ; d ; e{0.003h} ],size=(1150,1800))
plt_A1_km_2 = Plots.plot(title_C,plt_km_ghk_m2_xp,title_D,plt_km_no_ghk_m2_xg,empty_plt,layout=@layout[a{0.003h} ; b ; c{0.003h} ; d ; e{0.003h} ],size=(1150,1800))
plt_A1_km = Plots.plot(empty_plt,plt_A1_km_1,empty_plt,plt_A1_km_2,empty_plt,layout=@layout[a{0.003w} b c{0.003w} d e{0.003w}],size=(2300,1800),xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,legendfontsize=fontsize_,legend=:topleft,titlefontsize=fontsize_,fontfamily="Computer Modern")

plt_A1_kir_1 = Plots.plot(title_A,plt_kir_ghk_m_xp,title_B,plt_kir_no_ghk_m_xg,empty_plt,layout=@layout[a{0.003h} ; b ; c{0.003h} ; d ; e{0.003h} ],size=(1150,1800))
plt_A1_kir_2 = Plots.plot(title_C,plt_kir_ghk_m2_xp,title_D,plt_kir_no_ghk_m2_xg,empty_plt,layout=@layout[a{0.003h} ; b ; c{0.003h} ; d ; e{0.003h} ],size=(1150,1800))
plt_A1_kir = Plots.plot(empty_plt,plt_A1_kir_1,empty_plt,plt_A1_kir_2,empty_plt,layout=@layout[a{0.003w} b c{0.003w} d e{0.003w}],size=(2300,1800),xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,legendfontsize=fontsize_,legend=:topleft,titlefontsize=fontsize_,fontfamily="Computer Modern")


#Plots.savefig(plt_A1_km,"tests/figures/pdf/fig-A1-ICaLs_eq_KM.pdf")
#Plots.savefig(plt_A1_kir,"tests/figures/pdf/fig-A1-ICaLs_eq_Kir.pdf")

