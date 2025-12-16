using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using Statistics

include("../dI_grad/comp_dI_outputs.jl")
include("../fI_curve/compartment_simulation.jl")
include("../fI_curve/bif.jl")


function plt_V_I_derj_pA(t,V,I,lc_,Iticks,yticks_)
    dt_guide = 1

    lw_V = 1
    lw_I = 3
    fontsize_= 25

    plt_I = Plots.plot()
    #plot!(t,(maximum(I)+minimum(I))/2 .*ones(length(t)),ls=:dash,lc=:black,linealpha=0.8,lw=lw_I,ylims=(minimum(I)-0.001,maximum(I)+0.001))
    plot!(t,I,lc=lc_,lw=lw_I)
    plot!(yticks = Iticks,tickfontsize=fontsize_,tickfontcolor=:black,ylabel=L"I\,~"  *L"[" *latexify(u"pA") *L"]" ,ylabelfontsize=fontsize_)
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
function plt_fI_VI(I_list_fI_,f_list_fI_,I_stable_VI_,V_stable_VI_,I1_,I2_,I_scat,ms_scat_fI,ms_scat_VI,c_scat,xlims_,color_)

    fontsize_= 25

    ms_=10

    plt_fI = plot(ylims=(-1,75),legend=:topleft,xlims=xlims_,fontfamily="Computer Modern", ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,size=(800,400))
    plt_VI = plot(xlims=xlims_,legend=:bottomright,fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"pA") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-100,-50),yticks=-100:25:-50,size=(800,400))
                
    for j in eachindex(I_list_fI_)
        I_list_fI = I_list_fI_[j]
        f_list_fI = f_list_fI_[j]
        I_stable_VI = I_stable_VI_[j]
        V_stable_VI = V_stable_VI_[j]
        I1 = I1_[j]
        I2 = I2_[j]
        color = color_[j]

        fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
        fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]

        if j==1
            plot!(plt_fI,(xlims_[1]-100) *ones(2),-100 .*ones(2),lc=:grey,lw=1,label="Spiking")
            plot!(plt_fI,(xlims_[1]-100) *ones(2) ,-100 .*ones(2),ls=:dash,lw=1,lc=:grey,label="Resting")
                #plot!(plt_fI,I_list_fI_hh,f_list_fI_hh,lc=color_hh,lw=2,label=:none)
                #plot!(I_stable_VI_hh ,zeros(size(I_stable_VI_hh)),ls=:dash,lw=2,lc=color,label=:none)
            plot!(plt_fI,(xlims_[1]-100) *ones(2) ,-100 .*ones(2),fillrange=-100 .*ones(size(fI_fill)),c=:grey,fillalpha=0.15,lw=0,label="Bistable")
        end

        plot!(plt_fI,I_list_fI,f_list_fI,lc=color,lw=2,label=:none) #"Spiking"
        plot!(plt_fI,I_stable_VI ,zeros(size(I_stable_VI)),ls=:dash,lw=2,lc=color,label=:none) #"Resting"
            #plot!(plt_fI,I_list_fI_hh,f_list_fI_hh,lc=color_hh,lw=2,label=:none)
            #plot!(I_stable_VI_hh ,zeros(size(I_stable_VI_hh)),ls=:dash,lw=2,lc=color,label=:none)
        plot!(plt_fI,I_list_fI[fI_fill],f_list_fI[fI_fill],fillrange=zeros(size(fI_fill)),c=color,fillalpha=0.15,lw=0,label=:none) #"Bistable"
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

        if j==1
            plot!(plt_VI,(xlims_[1]-100) *ones(2) ,-1000 .*ones(2),ls=:dash,lw=1,lc=:grey,label="Resting")
                #plot!(plt_fI,I_list_fI_hh,f_list_fI_hh,lc=color_hh,lw=2,label=:none)
                #plot!(I_stable_VI_hh ,zeros(size(I_stable_VI_hh)),ls=:dash,lw=2,lc=color,label=:none)
            plot!(plt_VI,(xlims_[1]-100) *ones(2) ,-1000 .*ones(2),fillrange=-1000 .*ones(size(fI_fill)),c=:grey,fillalpha=0.15,lw=0,label="Bistable")
        end
        plot!(plt_VI,I_stable_VI ,V_stable_VI,ls=:dash,lw=2,lc=color,label=:none) #"Resting"
        plot!(plt_VI,I_stable_VI[VI_fill],V_stable_VI[VI_fill],fillrange=ones(size(VI_fill)).*-50,c=color,fillalpha=0.15,lw=0,label=:none)  #"Bistable"
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
    end

    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims_)
    return plt,plt_fI,plt_VI
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
    annotate!(0.5,0.45, latexstring("     ",L"[",latexify(u"pA"),L"]"),annotationfontsize=fontsize_)
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
    annotate!(0.5,0.45, latexstring("     ",L"[",latexify(u"pA"),L"]"),annotationfontsize=fontsize_)
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
annotate!(-0.5,0.9,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.5,0.,L"\bar{p}_{CaL} = "*"$(round(minimum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7)"*L" \ \ ; \ \bar{g}_{Kir} = "*"$(round(minimum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")


title_kir_km_2_B = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.5,0.9,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.5,0.,L"\bar{p}_{CaL} = "*"$(round(minimum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7)"*L" \ \ ; \ \bar{g}_{Kir} = "*"$(round(maximum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")

title_kir_km_2_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.5,0.9,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.5,0.,L"\bar{p}_{CaL} = "*"$(round(maximum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7)"*L" \ \ ; \ \bar{g}_{Kir} = "*"$(round(minimum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")

title_kir_km_2_D = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.5,0.9,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
#annotate!(0.5,0.4,L"\bar{p}_{CaL} = "*"$(round(maximum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7)"*latexify(u"cm/s"),annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
#annotate!(0.5,0.,L"\bar{g}_{Kir} = "*"$(round(maximum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)"*latexify(u"mS/cm^2"),annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
#annotate!(0.5,0.,L"\{\bar{p}_{CaL};\bar{g}_{Kir}\} = "*L"\{"*"$(round(maximum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7);$(round(maximum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)"*L"\}",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
#annotate!(0.5,0.4,L"\bar{p}_{CaL} = "*"$(round(maximum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7)",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
#annotate!(0.5,0.,L"\bar{g}_{Kir} = "*"$(round(maximum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")
annotate!(0.5,0.,L"\bar{p}_{CaL} = "*"$(round(maximum(pCaLs_d_range_KirCaLs[ind_pCaLs_simu])*1e7)/10^7)"*L" \ \ ; \ \bar{g}_{Kir} = "*"$(round(maximum(gKir_s_range_KirCaLs[ind_gKir_simu])*1e3)/10^3)",annotationfontsize=25,annotationhalign=:center,annotationfontfamily="Computer Modern")


title_kir_km_2_E = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{E}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_F = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{F}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_G = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{G}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_H = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{H}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_I = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{I}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")



plt_fI_VI_hCaLs,plt_fI_VI_hCaLs_fi,plt_fI_VI_hCaLs_vi = plt_fI_VI([I_list_Kir_s_CaLs_d_fI[2][1],I_list_Kir_s_CaLs_d_fI[2][3]].*10^6,[f_list_Kir_s_CaLs_d_fI[2][1],f_list_Kir_s_CaLs_d_fI[2][3]],[I_stable_Kir_s_CaLs_d_fI[2][1],I_stable_Kir_s_CaLs_d_fI[2][3]].*10^6,[Vs_stable_Kir_s_CaLs_d_fI[2][1],Vs_stable_Kir_s_CaLs_d_fI[2][3]],[I1_Kir_s_CaLs_d_fI[2,1],I_list_Kir_s_CaLs_d_fI[2][3][f_list_Kir_s_CaLs_d_fI[2][3].>0][1]].*10^6,[I2_Kir_s_CaLs_d_fI[2,1],I2_Kir_s_CaLs_d_fI[2,3]]*10^6,[],NaN,NaN,NaN,(-150,250),[:tan1,:green])
plt_fI_VI_lCaLs,plt_fI_VI_lCaLs_fi,plt_fI_VI_lCaLs_vi = plt_fI_VI([I_list_Kir_s_CaLs_d_fI[1][1],I_list_Kir_s_CaLs_d_fI[1][3]].*10^6,[f_list_Kir_s_CaLs_d_fI[1][1],f_list_Kir_s_CaLs_d_fI[1][3]],[I_stable_Kir_s_CaLs_d_fI[1][1],I_stable_Kir_s_CaLs_d_fI[1][3]].*10^6,[Vs_stable_Kir_s_CaLs_d_fI[1][1],Vs_stable_Kir_s_CaLs_d_fI[1][3]],[I1_Kir_s_CaLs_d_fI[1,1],I_list_Kir_s_CaLs_d_fI[1][3][f_list_Kir_s_CaLs_d_fI[1][3].>0][1]].*10^6,[I2_Kir_s_CaLs_d_fI[1,1],I2_Kir_s_CaLs_d_fI[1,3]]*10^6,[],NaN,NaN,NaN,(-150,250),[:deepskyblue4,:maroon3])

plt__lCaLs_lKir =plt_V_I_derj_pA(t_sol_simu[1],sol_simu[1][3,:],It_sol_simu[1],:deepskyblue4,([minimum(It_sol_simu[1]),maximum(It_sol_simu[1])],[round(minimum(It_sol_simu[1])*100)/100,round(maximum(It_sol_simu[1])*100)/100]),-100:50:50)
plt__lCaLs_hKir = plt_V_I_derj_pA(t_sol_simu[2],sol_simu[2][3,:],It_sol_simu[2],:maroon3,([minimum(It_sol_simu[2]),maximum(It_sol_simu[2])],[round(minimum(It_sol_simu[2])*100)/100,round(maximum(It_sol_simu[2])*100)/100]),-100:50:50)
plt__hCaLs_lKir = plt_V_I_derj_pA(t_sol_simu[3],sol_simu[3][3,:],It_sol_simu[3],:tan1,([minimum(It_sol_simu[3]),maximum(It_sol_simu[3])],[round(minimum(It_sol_simu[3])*100)/100,round(maximum(It_sol_simu[3])*100)/100]),-100:50:50)
plt__hCaLs_hKir = plt_V_I_derj_pA(t_sol_simu[4],sol_simu[4][3,:],It_sol_simu[4],:green,([minimum(It_sol_simu[4]),maximum(It_sol_simu[4])],[round(minimum(It_sol_simu[4])*100)/100,round(maximum(It_sol_simu[4])*100)/100]),-100:50:50)


plt_dI_comp = plt_dI_grad(best_dI_Kir_s_CaLs_d_J.*10^6,gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,(0,200.),L"\bar{g}_{\mathrm{Kir}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]",cgrad(:inferno ,rev=false))
plt_I1_comp = plt_I1_grad(best_I1_Kir_s_CaLs_d_J.*10^6,gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,(-100,200),L"\bar{g}_{\mathrm{Kir}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]",cgrad(:linear_worb_100_25_c53_n256 ,rev=true))
plt_Vs_comp = plt_Veq_grad(Vs_I1_Kir_s_CaLs_d,gKir_s_range_KirCaLs,pCaLs_d_range_KirCaLs,(-90,-47),L"\bar{g}_{\mathrm{Kir}} \ \ "*L"["*latexify(u"mS/cm^2")*L"]",L"\bar{p}_{\mathrm{CaL}} \ \ "*L"["*latexify(u"cm/s")*L"]")

palette_fI = palette([:deepskyblue4,:maroon3,:tan1,:green])
for plt in [plt_dI_comp,plt_I1_comp,plt_Vs_comp]
    for i in eachindex(ind_gKir_simu)
        println([gKir_s_range_KirCaLs[ind_gKir_simu[i]]],[pCaLs_d_range_KirCaLs[ind_pCaLs_simu[i]]])
        scatter!(plt,[gKir_s_range_KirCaLs[ind_gKir_simu[i]]],[pCaLs_d_range_KirCaLs[ind_pCaLs_simu[i]]],markershape=:circle,markersize=8,mc=palette_fI[i],markerstrokewidth=0.2,markerstrokecolor=:white)
    end
end



plt_right = plot(empty_plt,title_kir_km_2_G,plt_dI_comp,title_kir_km_2_H,plt_I1_comp,title_kir_km_2_I,plt_Vs_comp,empty_plt,layout=@layout[a{0.001h} ; b{0.02h} ;c ; d{0.02h} ; e ; f{0.02h} ; g ; h{0.001h}],size=(800,1800))
plt_center = plot(empty_plt,title_kir_km_2_E,plt_fI_VI_lCaLs_fi,plt_fI_VI_lCaLs_vi,title_kir_km_2_F,plt_fI_VI_hCaLs_fi,plt_fI_VI_hCaLs_vi,empty_plt,layout=@layout[a{0.001h} ; b{0.001h}; c ; d ;e{0.03h} ; f ; g; h{0.001h}],size=(800,1800))
plt_left = plot(empty_plt,title_kir_km_2_A,plt__lCaLs_lKir,title_kir_km_2_B,plt__lCaLs_hKir,title_kir_km_2_C,plt__hCaLs_lKir,title_kir_km_2_D,plt__hCaLs_hKir,empty_plt,layout=@layout[a{0.001h} ; b{0.05h} ; c ; d{0.05h} ;e ; f{0.05h} ; g; h{0.05h}; i ; j{0.001h}],size=(800,1800))

plot(empty_plt,plt_left,empty_plt,plt_center,empty_plt,plt_right,empty_plt,layout=@layout[a{0.001w} b{0.25w} c{0.001w} d{0.26w} e{0.001w} f g{0.001w}],size=(2300,1800),fontfamily="Computer Modern")

plot(empty_plt,plt_left,empty_plt,plt_center,plt_right,empty_plt,layout=@layout[a{0.001w} b{0.25w} c{0.001w} d{0.27w} e f{0.001w}],size=(2300,2000),fontfamily="Computer Modern")

#Plots.savefig("ChannelUpdate/LeFranc-multicomp-jl/figures/pdf/fig-10.pdf")   