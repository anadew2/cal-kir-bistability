using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

function plt_dI_1g_(plt,gx,dI,index_,label_,ms_,mshape_,lalpha)
    for ind in index_
        Plots.plot!(plt,gx,dI[:,ind],label=label_[ind],marker=mshape_,lw=2,linealpha=lalpha,ms=ms_,markerstrokewidth=0.)
    end
    return plt
end

# Palettes
pal_kir = []
for i in [1,3]
    push!(pal_kir,palette(:PuRd_9,rev=true)[i])
end
for i in [1,4]
    push!(pal_kir,palette(:beach,rev=false)[i])
end


pal_km = []
for i in [1,3]
    push!(pal_km,palette(:roma10,rev=true)[i])
end
for i in [5,8]
    push!(pal_km,palette(:imola10,rev=false)[i])
end


ylims_ = (-0.1,4.1)
xlabel_=L"\bar{p}_\mathrm{CaL} \ \ "*L"["*latexify(u"cm/s")*L"]"
plt_kir_pCaLf_xp = Plots.plot(fontfamily="Computer Modern",xlabel=xlabel_,ylabel=L"\Delta I \ \ "*L"["*latexify(u"µA/cm^2")*L"]",ylims=ylims_,legend=:outertopright,palette=pal_kir)

ylims_ = (-0.1,4.1)
plt_km_pCaLf_xp = Plots.plot(fontfamily="Computer Modern",xlabel=xlabel_,ylabel=L"\Delta I  \ \ "*L"["*latexify(u"µA/cm^2")*L"]",ylims=ylims_,legend=:outertopright,palette=pal_km)


# todo 
# - reprendre le code de la fonction plt de fig-A1 et la structure 
# - faire le (sub?)plot

include("../dI_grad/local_postpro_output.jl")
    ## Kir
    index_kir_pCaLf_0 = vcat(1,17)
    label_kir_pCaLf_0_ = []
    for i in eachindex(gKir_range_KirCaLs[1:3:end])
        push!(label_kir_pCaLf_0_,L"\bar{p}_\mathrm{CaLf}=0, \ \bar{g}_{\mathrm{Kir}}=%$(round(gKir_range_KirCaLs[1:3:end][i]*10)/10) ")
    end
    plt_dI_1g_(plt_kir_pCaLf_xp,pCaLs_range_KirCaLs,dI_KirCaLs[:,1:3:end],index_kir_pCaLf_0,label_kir_pCaLf_0_,3,:rect,1)

    ## KM
    index_km_pCaLf_0 = vcat(1,17)
    label_km_pCaLf_0 = []
    for i in eachindex(gKM_range_KMCaLs[1:3:end])
        push!(label_km_pCaLf_0,L"\bar{p}_\mathrm{CaLf}=0, \ \bar{g}_{\mathrm{KM}}=%$(round(gKM_range_KMCaLs[1:3:end][i]*10)/10) ")
    end
    plt_dI_1g_(plt_km_pCaLf_xp,pCaLs_range_KMCaLs,dI_KMCaLs[:,1:3:end],index_km_pCaLf_0,label_km_pCaLf_0,3,:rect,1)


include("../dI_grad/local_postpro_output_pCaLf_e-4.jl")

    pCaLs_range_KirCaLs__pCaLf_e_4 = pCaLs_range_KirCaLs[1:2:end]
    gKir_range_KirCaLs__pCaLf_e_4 = gKir_range_KirCaLs[1:3:end]
    pCaLs_range_KMCaLs__pCaLf_e_4 = pCaLs_range_KMCaLs[1:2:end]
    gKM_range_KMCaLs__pCaLf_e_4 = gKM_range_KMCaLs[1:3:end]

    ## Kir
    index_kir_pCaLf_1e_4 = vcat(1,17)
    label_kir_pCaLf_1e_4_ = []
    for i in eachindex(gKir_range_KirCaLs[1:3:end])
        push!(label_kir_pCaLf_1e_4_,L"\bar{p}_\mathrm{CaLf}=1.5e-4, \ \bar{g}_{\mathrm{Kir}}=%$(round(gKir_range_KirCaLs[1:3:end][i]*10)/10) ")
    end
    plt_dI_1g_(plt_kir_pCaLf_xp,pCaLs_range_KirCaLs__pCaLf_e_4,dI_KirCaLs_pCaLf_15e_4,index_kir_pCaLf_1e_4,label_kir_pCaLf_1e_4_,3,:circle,0.8)

    ## KM
    index_km_pCaLf_1e_4 = vcat(1,17)
    label_km_pCaLf_1e_4 = []
    for i in eachindex(gKM_range_KMCaLs__pCaLf_e_4)
        push!(label_km_pCaLf_1e_4,L"\bar{p}_\mathrm{CaLf}=1.5e-4, \ \bar{g}_{\mathrm{KM}}=%$(round(gKM_range_KMCaLs__pCaLf_e_4[i]*10)/10) ")
    end
    plt_dI_1g_(plt_km_pCaLf_xp,pCaLs_range_KMCaLs__pCaLf_e_4,dI_KMCaLs_pCaLf_15e_4,index_km_pCaLf_1e_4,label_km_pCaLf_1e_4,3,:circle,0.8)

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

fontsize_=25
plt_full_1 = Plots.plot(title_A,plt_kir_pCaLf_xp,empty_plt,layout=@layout[a{0.003h}; b; c{0.003h}])
plt_full_2 = Plots.plot(title_B,plt_km_pCaLf_xp,empty_plt,layout=@layout[a{0.003h}; b; c{0.003h}])

plt_full = Plots.plot(empty_plt,plt_full_1,empty_plt,plt_full_2,empty_plt,layout=@layout[a{0.003w} b c{0.003w} d e{0.003w}],size=(2300,900),xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,legendfontsize=fontsize_,legend=:topleft,titlefontsize=fontsize_,fontfamily="Computer Modern")
#Plots.savefig(plt_full,"tests/figures/pdf/fig-A2-pCaLf_increased.pdf")
