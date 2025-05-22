using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

include("../dI_noise/local_postpro_outputs.jl")

function plt_V_It_noisy_inc_sigma(t_,V,t_It_,It,t_sigma_,sigma,I1bist,I2bist,mean_It,fill_c,yticks_,sigma_ticks_)
    dt_guide = 1
    dV_guide = 50
    dI_guide = 10^(-1) #std

    fontsize_=25
    write_an = true
    lc_V = RGB(0.2,0.2,0.2)
    lc_I = :black
    lw_V =1

    t = t_ ./1000
    t_It = t_It_ ./1000
    t_sigma = t_sigma_ ./1000


    plt_V = Plots.plot()
    plot!(t,V,lc=lc_V,lw=lw_V)     
    #plot lower guide
    plot!([(maximum(t))-dt_guide , maximum(t)],(yticks_[1]) .*ones(2),lc=:black,lw=2)
    if write_an == 1
        annotate!(((maximum(t))-dt_guide + maximum(t))/2,(yticks_[1])*1.13, latexify(Int(dt_guide) * u"s"),annotationfontsize=fontsize_)
    end
    #axis attributes 
    plot!(xlims=(-maximum(t) *0.03,maximum(t) *1.15) )
    plot!(ylims=(minimum(yticks_),maximum(yticks_)))
    plot!(yticks = yticks_,ylabel=L"V\,~"  *L"(" *latexify(u"mV") *L")",ylabelfontsize=fontsize_)
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,showaxis=:y,foreground_color_axis=:black)
#
    plt_sigma=Plots.plot()
    plot!(t_sigma,sigma,lc=lc_V,lw=2)
    plot!(xlims=(-maximum(t) *0.03,maximum(t) *1.15) )
    plot!(ylabel=L"\sigma\,~"  *L"(" *latexify(u"mV") *L")",ylabelfontsize=fontsize_)
    plot!(grid=false,legend=:none,fontfamily="Computer Modern",tickfontsize=fontsize_,legendfontsize=fontsize_)
    plot!(foreground_color_border=:black,showaxis=:y,foreground_color_axis=:black)

    #==
    l = @layout[
        a{1.0*w,0.3*h}
        b{1.0*w,0.3*h}
        b{1.0*w,0.4*h}
    ]
    plt = Plots.plot(plt_sigma,plt_I,plt_V,layout=l,size=(1000,800))==#

    l = @layout[
        a{1.0*w,0.4*h}
        b{1.0*w,0.6*h}
    ]
    plt = Plots.plot(plt_sigma,plt_V,layout=l,size=(1000,800))
    ==#
    return plt
end

include("../dI_noise/asc_std_kir_Ihalf.jl")
plt_noisy_Kir_asc_st = plt_V_It_noisy_inc_sigma(noisy_sol_cb_t_kir_st,noisy_sol_cb_V_kir_st,NaN,NaN,t_noise_kir_st,std_It_noise_kir_st.*sqrt(0.01),I1bist,I2bist,I0,:mediumpurple3,-100:50:50,0:1:7)
plt_noisy_Kir_asc_cyc = plt_V_It_noisy_inc_sigma(noisy_sol_cb_t_kir_cyc,noisy_sol_cb_V_kir_cyc,NaN,NaN,t_noise_kir_cyc,std_It_noise_kir_cyc.*sqrt(0.01),I1bist,I2bist,I0,:mediumpurple3,-100:50:50,0:1:7)
include("../dI_noise/asc_std_km_Ihalf.jl")
plt_noisy_KM_asc_st = plt_V_It_noisy_inc_sigma(noisy_sol_cb_t_km_st,noisy_sol_cb_V_km_st,NaN,NaN,t_noise_km_st,std_It_noise_km_st.*sqrt(0.01),I1bist,I2bist,I0,:mediumpurple3,-100:50:50,0:10:70)
plt_noisy_KM_asc_cyc = plt_V_It_noisy_inc_sigma(noisy_sol_cb_t_km_cyc,noisy_sol_cb_V_km_cyc,NaN,NaN,t_noise_km_cyc,std_It_noise_km_cyc.*sqrt(0.01),I1bist,I2bist,I0,:mediumpurple3,-100:50:50,0:10:70)


empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

#Horizontal version
plt_noisy_KM_full_inc = Plots.plot(plt_noisy_KM_asc_st,empty_plt,plt_noisy_KM_asc_cyc,layout=@layout[a b{0.05w} c])
plt_noisy_Kir_full_inc = Plots.plot(plt_noisy_Kir_asc_st,empty_plt,plt_noisy_Kir_asc_cyc,layout=@layout[a b{0.05w} c])


title_noisy_KM = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.05,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.465,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_prop_noisy_KM = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.375,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_noisy_Kir = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.05,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.465,0.5,L"\textbf{E}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_prop_noisy_Kir = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.3755,0.5,L"\textbf{F}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")

title_lines_KM_Kir = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
#plot!(0.15*ones(2),0.78 .+ [-0.2,0.2],label=:none,lc=:mediumpurple4)
#plot!(0.15*ones(2),0.2 .+ [-0.2,0.2],label=:none,lc=RGB(0.95,0.62,0.75))
annotate!(0.25,0.78,L"\mathrm{KM}",annotationfontsize=45,annotationhalign=:center,annotationfontfamily="Computer Modern",annotationcolor=:mediumpurple3)
annotate!(0.25,0.18,L"\mathrm{Kir}",annotationfontsize=45,annotationhalign=:center,annotationfontfamily="Computer Modern",annotationcolor=:palevioletred2)


fontsize_=25
plt_switch_Kir = Plots.plot(fontfamily="Computer Modern",xlims=(-0.1,7.1))
plot!(std_vec_Kir,prop_switch_std_Kir,lc=:palevioletred2,lw=2,marker=:circle,mc=:palevioletred2,ms=5,markerstrokewidth=0,label="From resting")
plot!(std_vec_Kir_cyc,prop_switch_std_Kir_cyc,lc=:palevioletred2,lw=2,marker=:circle,mc=:palevioletred2,ms=5,markerstrokewidth=0,label="From spiking",ls=:dash)
plot!(xlabel=L"\sigma\,~"  *L"(" *latexify(u"mV") *L")",ylabel="Proportion of state transition",ylabelfontsize=fontsize_,xlabelfontsize=fontsize_,yticks =0:1,tickfontsize=fontsize_,legendfontsize=fontsize_,legend=:bottomright)

plt_switch_KM = Plots.plot(fontfamily="Computer Modern",xlims=(-0.1,7.1))
plot!(std_vec_KM,prop_switch_std_KM,lc=:mediumpurple3,lw=2,marker=:circle,mc=:mediumpurple3,ms=5,markerstrokewidth=0,label="From resting")
plot!(std_vec_KM_cyc,prop_switch_std_KM_cyc,lc=:mediumpurple3,lw=2,marker=:circle,mc=:mediumpurple3,ms=5,markerstrokewidth=0,label="From spiking",ls=:dash)
plot!(xlabel=L"\sigma\,~"  *L"(" *latexify(u"mV") *L")",ylabel="Proportion of state transition",ylabelfontsize=fontsize_,legend=:topright,xlabelfontsize=fontsize_,yticks =0:1,tickfontsize=fontsize_,legendfontsize=fontsize_)



plt_all_noisy_inc = Plots.plot(title_noisy_KM,plt_noisy_KM_full_inc,empty_plt,title_noisy_Kir,plt_noisy_Kir_full_inc,layout=@layout[a{0.003h}; b; c{0.05h}; d{0.003h} ; e],tickfontsize=22)
plt_switch = Plots.plot(title_prop_noisy_KM,plt_switch_KM,empty_plt,title_prop_noisy_Kir,plt_switch_Kir,layout=@layout[a{0.003h}; b; c{0.003h} ; e{0.003h}; f],tickfontsize=22)
plot!(plt_switch[2],xlabel=L"\sigma\,~"  *L"(" *latexify(u"mV") *L")",ylabel="Proportion of state transition",ylabelfontsize=28,legend=:topright,xlabelfontsize=28)
plot!(plt_switch[5],xlabel=L"\sigma\,~"  *L"(" *latexify(u"mV") *L")",ylabel="Proportion of state transition",ylabelfontsize=28,xlabelfontsize=28)

plt_noisy_inc = Plots.plot(plt_all_noisy_inc,empty_plt,plt_switch,layout=@layout[a{0.7w} b{0.05w} c],fontfamily="Computer Modern",size=(2300,1900))
plt_noisy_inc_ = Plots.plot(empty_plt,plt_noisy_inc,empty_plt,layout=@layout[a{0.001w} b c{0.001w}])
plt_f6_v8 = Plots.plot(empty_plt,plt_noisy_inc_,empty_plt,layout=@layout[a{0.003h} ; b ; c{0.001h}],fontfamily="Computer Modern")

#Plots.savefig(plt_f6_v8,"ChannelUpdate/cluster/figures/pdf/fig-4.pdf")
