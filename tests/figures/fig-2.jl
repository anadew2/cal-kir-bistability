using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

include("../sim_pulse_step/sim_Kir.jl")
include("../sim_pulse_step/sim_KM.jl")
include("../fI_curve/bif_kir.jl")
include("../fI_curve/bif_km.jl")

function plt_fI_VI(I_list_fI,f_list_fI,I_stable_VI,V_stable_VI,I1,I2,V_I1,V_I2,I_scat,ms_scat_fI,ms_scat_VI,c_scat,xlims_,color)
    fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
    fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]
    fontsize_= 25

    ms_=10

    #I_list_fI_hh,f_list_fI_hh,I_stable_VI_hh,V_stable_VI_hh,I1_hh,I2_hh,V_I1_hh,V_I2_hh = I_list_KirCaLs_fI[1][1],f_list_KirCaLs_fI[1][1],I_stable_KirCaLs_fI[1][1],V_stable_KirCaLs_fI[1][1],I1_KirCaLs_fI[1,1],I2_KirCaLs_fI[1,1],V_I1_KirCaLs_fI[1,1],V_I2_KirCaLs_fI[1,1]
    color_hh=:gray

    plt_fI = plot(fontfamily="Computer Modern", ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legend=:outerbottomright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,size=(800,400))
    
    if I1 != I2
        #plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]))
        plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],round(I1*10)/10,round(I2*10)/10,xlims_[2]]),xlims=xlims_)
        plot!(I1.*ones(2),[-5,250],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        plot!(I2.*ones(2),[-5,250],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        annotate!(I1,250+15,L"I_1",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
        annotate!(I2,250+15,L"I_2",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    else
        plot!(xticks=([xlims_[1],I1,xlims_[2]],[xlims_[1],L"I_1=I_2",xlims_[2]]))
    end
    plot!([I1,I2],220 .*ones(2),arrow=true,color=:grey25,linewidth=0.5,label=:none)
    plot!([I2,I1],220 .*ones(2),arrow=true,color=:grey25,linewidth=0.5,label=:none,ylims=(-10,250))
    annotate!(mean([I2,I1]),250,L"\Delta I",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    plot!(yticks=[0,250])
    plot!(plt_fI,I_list_fI,f_list_fI,lc=color,lw=2,label="Spiking")
    plot!(I_stable_VI ,zeros(size(I_stable_VI)),ls=:dash,lw=2,lc=color,label="Resting")
        #plot!(plt_fI,I_list_fI_hh,f_list_fI_hh,lc=color_hh,lw=2,label=:none)
        #plot!(I_stable_VI_hh ,zeros(size(I_stable_VI_hh)),ls=:dash,lw=2,lc=color,label=:none)
    plot!(I_list_fI[fI_fill],f_list_fI[fI_fill],fillrange=zeros(size(fI_fill)),c=color,fillalpha=0.15,lw=0,label="Bistable")
    #Display markers for currents used in step responses 
    for i in eachindex(I_scat)
        ind_I_list_fI = (1:length(I_list_fI))[I_list_fI.>I_scat[i]]
        if length(ind_I_list_fI)>0
            I_scat_i = 1/2*(I_list_fI[ind_I_list_fI[1]] + I_list_fI[ind_I_list_fI[1]-1])
            f_scat_i = 1/2*(f_list_fI[ind_I_list_fI[1]] + f_list_fI[ind_I_list_fI[1]-1])
            scatter!(plt_fI,I_scat_i*ones(2),f_scat_i*ones(2),markercolor=c_scat[i],markerstrokewidth=0.1,markerstrokecolor=:white,markershape=ms_scat_fI[i],ms=ms_,label=:none)
        else
            println("Current is too high to be marked")
        end
    end

    VI_fill_ = collect(1:length(I_stable_VI))[I_stable_VI .>=I1]
    VI_fill = VI_fill_[I_stable_VI[VI_fill_] .<= I2]  

    plt_VI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legend=:outerbottomright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-100,-50),yticks=-100:25:-50,size=(800,400))
    #plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]))
    plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],round(I1*10)/10,round(I2*10)/10,xlims_[2]]),xlims=xlims_)
        #plot!(I_stable_VI_hh ,V_stable_VI_hh,ls=:dash,lw=2,lc=color_hh,label=:none)
        #Display I1 & I2 lines
        plot!(I1.*ones(2),[-150,-50],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        plot!(I2.*ones(2),[-150,-50],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        annotate!(I1,-50+6,L"I_1",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
        annotate!(I2,-50+6,L"I_2",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    plot!(I_stable_VI ,V_stable_VI,ls=:dash,lw=2,lc=color,label="Resting")
    plot!(I_stable_VI[VI_fill],V_stable_VI[VI_fill],fillrange=ones(size(VI_fill)).*-50,c=color,fillalpha=0.15,lw=0,label="Bistable")
    scatter!(I1.*ones(2), V_I1.*ones(2),mc=color,markerstrokewidth=0,label=:none)
    scatter!(I2.*ones(2), V_I2.*ones(2),mc=color,markerstrokewidth=0,label=:none)    
    #Display markers for currents used in step responses 
    for i in eachindex(I_scat)

        ind_I_stable_VI = (1:length(I_stable_VI))[I_stable_VI.>I_scat[i] .&& I_stable_VI.<10]
        if length(ind_I_stable_VI)>0
            I_scat_i = 1/2*(I_stable_VI[ind_I_stable_VI[1]] + I_stable_VI[ind_I_stable_VI[1]-1])
            V_scat_i = 1/2*(V_stable_VI[ind_I_stable_VI[1]] + V_stable_VI[ind_I_stable_VI[1]-1])
            scatter!(plt_VI,I_scat_i*ones(2),V_scat_i*ones(2),markercolor=c_scat[i],markerstrokewidth=0.1,markerstrokecolor=:white,markershape=ms_scat_VI[i],ms=ms_,label=:none)
        else
            println("Current is too high to be marked")
        end
    end

    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims_)
    return plt,plt_fI,plt_VI
end

function plt_V_I_derj(t,V,I,V_ticks_)

    lc_V = RGB(0.2,0.2,0.2)
    lc_I = :black
    lw_V = 1
    lw_I = 2

    fontsize_ = 25
    write_an = true

    V_lim_l = minimum(V_ticks_)

    plt_I = Plots.plot()
    plot!(t,I,lc=lc_I,lw=lw_I)
    ==#
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(xticks=0:round(maximum(t)):round(maximum(t)))
    plot!(yticks = ([minimum(I),maximum(I)],[round(minimum(I)*10)/10,round(maximum(I)*10)/10]),ylabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]"  ,ylabelfontsize=fontsize_) #ytciks=round.([minimum(I),maximum(I)].*100)./100
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)

    plt_V = Plots.plot()
    plot!(t,V,lc=lc_V,lw=lw_V,label=:none)
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(ylims=(minimum(V_ticks_),maximum(V_ticks_)))
    plot!(xticks=0:round(maximum(t)):round(maximum(t)))
    plot!(yticks = V_ticks_,ylabel = L"V \ \ "  *L"[" *latexify(u"mV") *L"]",xlabelfontsize=fontsize_,ylabelfontsize=fontsize_)
    
    plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,foreground_color_axis=:black)

    plt = Plots.plot(plt_I,plt_V,layout=@layout[a{0.3h} ; b],size=(1000,800))
    return plt
end

function plt_Iion(t,ICa,IK,lc_Ca,lc_K,label_Ca,label_K,I_ticks_,Ca_before_K)

    lw_V = 1
    fontsize_ = 25

    plt = Plots.plot()
    if Ca_before_K
        plot!(t,ICa,lc=lc_Ca,lw=lw_V,label=label_Ca)
        plot!(t,IK,lc=lc_K,lw=lw_V,label=label_K)
    else
        plot!(t,IK,lc=lc_K,lw=lw_V,label=label_K)
        plot!(t,ICa,lc=lc_Ca,lw=lw_V,label=label_Ca)
    end
    #axis attributes 
    plot!(xlims=(-maximum(t)*0.03,maximum(t)*1.03) )
    plot!(ylims=(minimum(I_ticks_),maximum(I_ticks_)))
    plot!(xticks=0:round(maximum(t)):round(maximum(t)))
    plot!(yticks = I_ticks_,ylabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",xlabel=L"t \ \ "*L"["*latexify(u"ms")*L"]" ,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_)
    
    plot!(grid=false,legend=:bottomright,legendfontsize=fontsize_,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(foreground_color_border=:black,foreground_color_axis=:black)

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


    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gKir = gKir_range_KirCaLs_fI[2]        # [mS/cm²]
    gKM = 0.     # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]
    
    ICaLs_sol_pulse_gKir_fI_2 = ICaLs.(pCaLs,sol_pulse_gKir_fI_2[7,:],sol_pulse_gKir_fI_2[8,:],sol_pulse_gKir_fI_2[1,:],sol_pulse_gKir_fI_2[11,:],Ca_o)
    IKir_sol_pulse_gKir_fI_2 = IKir.(gKir,sol_pulse_gKir_fI_2[9,:],sol_pulse_gKir_fI_2[1,:],eK)
    IKM_sol_pulse_gKir_fI_2 = IKM.(gKM,sol_pulse_gKir_fI_2[10,:],sol_pulse_gKir_fI_2[1,:],eK)
    INa_sol_pulse_gKir_fI_2  = INa.(gNa,sol_pulse_gKir_fI_2[2,:],sol_pulse_gKir_fI_2[3,:],sol_pulse_gKir_fI_2[1,:],eNa)
    IKDR_sol_pulse_gKir_fI_2  = IKDR.(gKDR,sol_pulse_gKir_fI_2[4,:],sol_pulse_gKir_fI_2[1,:],eK)
    
plt_Iion_Kir = plt_Iion(sol_pulse_gKir_fI_2.t[1:167170],ICaLs_sol_pulse_gKir_fI_2[1:167170],IKir_sol_pulse_gKir_fI_2[1:167170],:darkorange2,:maroon3,L"I_{\mathrm{CaL}}",L"I_{\mathrm{Kir}}",-15:15:15,true)

    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gKir = 0.          # [mS/cm²]
    gKM = gKM_range_KMCaLs_fI[2]      # [mS/cm²]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = pCaLs_range_KirCaLs_fI[3]   # [cm/s]

    ICaLs_sol_pulse_gKM_fI_2 = ICaLs.(pCaLs,sol_pulse_gKM_fI_2[7,:],sol_pulse_gKM_fI_2[8,:],sol_pulse_gKM_fI_2[1,:],sol_pulse_gKM_fI_2[11,:],Ca_o)
    IKir_sol_pulse_gKM_fI_2 = IKir.(gKir,sol_pulse_gKM_fI_2[9,:],sol_pulse_gKM_fI_2[1,:],eK)
    IKM_sol_pulse_gKM_fI_2 = IKM.(gKM,sol_pulse_gKM_fI_2[10,:],sol_pulse_gKM_fI_2[1,:],eK)

plt_Iion_KM = plt_Iion(sol_pulse_gKM_fI_2.t[1:82745],ICaLs_sol_pulse_gKM_fI_2[1:82745],IKM_sol_pulse_gKM_fI_2[1:82745],:darkorange2,:mediumpurple3,L"I_{\mathrm{CaL}}",L"I_{\mathrm{KM}}",-15:15:15,false)

    
plt_fI_VI_Kir,plt_fI_VI_Kir_fi,plt_fI_VI_Kir_vi = plt_fI_VI(I_list_KirCaLs_fI[3][2],f_list_KirCaLs_fI[3][2],I_stable_KirCaLs_fI[3][2],V_stable_KirCaLs_fI[3][2],I1_KirCaLs_fI[3,2],I2_KirCaLs_fI[3,2],V_I1_KirCaLs_fI[3,2],V_I2_KirCaLs_fI[3,2],[],NaN,NaN,NaN,(-2,2).+0.4,:maroon3)
plt_fI_VI_KM,plt_fI_VI_KM_fi,plt_fI_VI_KM_vi = plt_fI_VI(I_list_KMCaLs_fI[3][2],f_list_KMCaLs_fI[3][2],I_stable_KMCaLs_fI[3][2],V_stable_KMCaLs_fI[3][2],I1_KMCaLs_fI[3,2],I2_KMCaLs_fI[3,2],V_I1_KMCaLs_fI[3,2],V_I2_KMCaLs_fI[3,2],[],NaN,NaN,NaN,(0,4),:mediumpurple3)

plt_pulse_kir = plt_V_I_derj(sol_pulse_gKir_fI_2.t[1:167170],sol_pulse_gKir_fI_2[1,1:167170],It_pulse_gKir_fI_2[1:167170],50:-50:-100)
plt_pulse_km = plt_V_I_derj(sol_pulse_gKM_fI_2.t[1:82745],sol_pulse_gKM_fI_2[1,:][1:82745],It_pulse_gKM_fI_2[1:82745],50:-50:-100)

plt_pulse_fI_VI_km = Plots.plot(Plots.plot!(title_kir_km_1_D,empty_plt,plt_pulse_kir,plt_Iion_Kir,layout=@layout[a{0.0003h};b{0.0003h};c;d{0.4h}]),title_kir_km_1_E,plt_fI_VI_Kir_fi,title_kir_km_1_F,plt_fI_VI_Kir_vi,empty_plt,layout=@layout[a{0.6h} ; b{0.01h} ; c ; d{0.01h} ; e ; f{0.01h}],size=(1100,1700))
plt_pulse_fI_VI_kir = Plots.plot(Plots.plot!(title_kir_km_1_A,empty_plt,plt_pulse_km,plt_Iion_KM,layout=@layout[a{0.003h};b{0.0003h};c;d{0.4h}]),title_kir_km_1_B,plt_fI_VI_KM_fi,title_kir_km_1_C,plt_fI_VI_KM_vi,empty_plt,layout=@layout[a{0.6h} ; b{0.01h} ; c ; d{0.01h} ; e ; f{0.01h}],size=(1100,1700))
plt_kir_km_1 = Plots.plot(empty_plt,Plots.plot(empty_plt,plt_pulse_fI_VI_kir,empty_plt,plt_pulse_fI_VI_km,empty_plt,layout=@layout[a{0.075w} b c{0.1w} d e{0.075w}],size=(2300,2000)),layout=@layout[a{0.001h};b],fontfamily="Computer Modern",foreground_color=:gray25)

#savefig(plt_kir_km_1,"tests/figures/pdf/fig-2.pdf")
