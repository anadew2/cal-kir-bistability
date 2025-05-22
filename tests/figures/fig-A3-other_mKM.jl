
using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

include("../dI_grad_other_km/local_postpro_output_other_km.jl")


empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)


fontsize_ = 25
Plots.plot(fontfamily="Computer Modern",xlabel=L"p_{\mathrm{CaL}}",ylabel=L"\Delta I ",ylims=(-0.1,3.1)) #,title=L"V_{1/2,}_{\mathrm{KM}}="*latexify(26.7(u"mV"))*L"\ ; \ z_{\mathrm{KM}}="*latexify(12.6(u"mV"))
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,1],label=L"g_{\mathrm{KM}}=0",marker=:rect,lw=2,markerstrokewidth=0.,ms=3,fontfamily="Computer Modern",legend=:outertopright,lc=:black,mc=:black,palette=palette(:imola10,rev=false))
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,2],label=L"g_{\mathrm{KM}}=0.1",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,3],label=L"g_{\mathrm{KM}}=0.2",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,4],label=L"g_{\mathrm{KM}}=0.3",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,6],label=L"g_{\mathrm{KM}}=0.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,11],label=L"g_{\mathrm{KM}}=1.",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
Plots.plot!(pCaLs_range_KMCaLs,dI_KMCaLs_other_KM[:,16],label=L"g_{\mathrm{KM}}=1.5",marker=:rect,lw=2,markerstrokewidth=0.,ms=3)
plt__ = Plots.plot!(size=(1150,900),xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,legendfontsize=fontsize_,legend=:topleft,titlefontsize=fontsize_,fontfamily="Computer Modern")

plt_ = Plots.plot(empty_plt,plt__,empty_plt,layout=@layout[a{0.003w} b c{0.003w}])
plt = Plots.plot(empty_plt,plt_,empty_plt,layout=@layout[a{0.003h}; b; c{0.003h}],fontfamily="Computer Modern")

#Plots.savefig(plt,"ChannelUpdate/cluster/figures/pdf/fig-A3-other_mKM.pdf")