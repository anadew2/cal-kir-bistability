using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using Statistics

include("../sim_pulse_step/sim_CaLs.jl")
include("../fI_curve/bif_kir.jl")

function plt_fI_VI(I_list_fI,f_list_fI,I_stable_VI,V_stable_VI,I1,I2,V_I1,V_I2,I_scat,ms_scat_fI,ms_scat_VI,c_scat,xlims_,color)
    fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
    fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]
    fontsize_= 25

    ms_=13

    I_list_fI_hh,f_list_fI_hh,I_stable_VI_hh,V_stable_VI_hh,I1_hh,I2_hh,V_I1_hh,V_I2_hh = I_list_KirCaLs_fI[1][1],f_list_KirCaLs_fI[1][1],I_stable_KirCaLs_fI[1][1],V_stable_KirCaLs_fI[1][1],I1_KirCaLs_fI[1,1],I2_KirCaLs_fI[1,1],V_I1_KirCaLs_fI[1,1],V_I2_KirCaLs_fI[1,1]
    color_hh=:gray

    plt_fI = plot(fontfamily="Computer Modern",ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legend=:outerbottomright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,size=(800,400))
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
            if ms_scat_fI[i]==:diamond
                scatter!(plt_fI,I_scat_i*ones(2),f_scat_i*ones(2),markercolor=c_scat[i],markerstrokewidth=0.1,markerstrokecolor=:white,markershape=ms_scat_fI[i],ms=ms_+2,label=:none)
            else
                scatter!(plt_fI,I_scat_i*ones(2),f_scat_i*ones(2),markercolor=c_scat[i],markerstrokewidth=0.1,markerstrokecolor=:white,markershape=ms_scat_fI[i],ms=ms_,label=:none)
            end
        else
            #println("Current is too high to be marked")
        end
    end

    VI_fill_ = collect(1:length(I_stable_VI))[I_stable_VI .>=I1]
    VI_fill = VI_fill_[I_stable_VI[VI_fill_] .<= I2]  

    plt_VI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legend=:outerbottomright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-150,-50),yticks=-150:50:-50,size=(800,400))
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
            #println("Current is too high to be marked")
        end
    end

    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims_)
    return plt,plt_fI,plt_VI
end

function plt_V_mult_steps(t,V,It,lc_V,I_lim,ms_scat,yticks_)
    plt_V_list = []
    fontsize_ = 25
    plt_I = Plots.plot(fontfamily="Computer Modern",grid=:none,legend=:none,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)

    lw_V =2
    t_scat = maximum(t[1])-20
    ms_=10
    dV_m = 30

    for i in eachindex(V)
        plt_V = Plots.plot(fontfamily="Computer Modern",grid=:none,legend=:none,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylabel = L"V \ \ "*L"["*latexify(u"mV")*L"]")
        plot!(t[i],V[i],lc=lc_V[i],lw=lw_V) 

        plot!(xlims=(-maximum(t[i]) *0.03,maximum(t[i]) *1.03),xticks=(0:500:maximum(t[i]),[0,"",1000]))
        plot!(ylims=(minimum(yticks_),50))
        plot!(yticks=reverse(50:-100:-175),ylabelfontsize=fontsize_)        
        plot!(grid=false,legend=:none,tickfontfamily="Computer Modern",tickfontsize=fontsize_)
            if V[i][end]>-65
                V_scat = minimum(V[i][end-100:end])-dV_m
            else
                V_scat = V[i][end]+dV_m
            end
            if ms_scat[i]==:diamond
                scatter!(plt_V,t_scat*ones(2),V_scat*ones(2),markercolor=lc_V[i],markerstrokewidth=0.1,markerstrokecolor=:white,markershape=ms_scat[i],ms=ms_+2,label=:none)
            else
                scatter!(plt_V,t_scat*ones(2),V_scat*ones(2),markercolor=lc_V[i],markerstrokewidth=0.1,markerstrokecolor=:white,markershape=ms_scat[i],ms=ms_,label=:none)
            end
        push!(plt_V_list,plt_V)

        plot!(plt_I,t[i],It[i],lc=lc_V[i],lw=lw_V,ylabel = L"I \ \ "*L"["*latexify(u"µA/cm^2")*L"]",xlabel=L"t \ \ "*L"["*latexify(u"ms")*L"]")
        #axis attributes 
        plot!(plt_I,xlims=(-maximum(t[i]) *0.03,maximum(t[i]) *1.03) )
        plot!(plt_I,xticks=(0:500:maximum(t[i]),[0,"",1000]))
        plot!(plt_I,ylims=[I_lim[1]-0.5,I_lim[end]], yticks=[I_lim[1],mean([I_lim[1],I_lim[end]]),I_lim[end]])
        plot!(plt_I,ylabelfontsize=fontsize_) #ytciks=round.([minimum(I),maximum(I)].*100)./100
        plot!(plt_I,grid=false,legend=:none,tickfontfamily="Computer Modern",legendfontfamily="Computer Modern",tickfontsize=fontsize_,legendfontsize=fontsize_)
        #display(plt_I)
    end
    return plt_V_list,plt_I
end


empty_plt = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

title_new_CaL_A = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.025,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_new_CaL_B = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.05,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_new_CaL_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.3,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")

I_marked = []
for i in eachindex(I_asc_step_CaLs)
    push!(I_marked,I_asc_step_CaLs[i][end])
end

ms_spiking = [:circle,:diamond,:diamond,:circle]
ms_resting = [:circle,:rect,:rect,:circle]

plt_fI_VI_CaLs,plt_fI_VI_CaLs_fi,plt_fI_VI_CaLs_vi = plt_fI_VI(I_list_KirCaLs_fI[3][1],f_list_KirCaLs_fI[3][1],I_stable_KirCaLs_fI[3][1],V_stable_KirCaLs_fI[3][1],I1_KirCaLs_fI[3,1],I2_KirCaLs_fI[3,1],V_I1_KirCaLs_fI[3,1],V_I2_KirCaLs_fI[3,1],I_marked,ms_spiking,ms_resting,reverse([palette(:twelvebitrainbow)[4],palette(:twelvebitrainbow)[3],palette(:twelvebitrainbow)[2],palette(:twelvebitrainbow)[1]]),(-3,1),:darkorange2)
plt_fI_VI_hh,plt_fI_VI_hh_fi,plt_fI_VI_hh_vi = plt_fI_VI(I_list_KirCaLs_fI[1][1],f_list_KirCaLs_fI[1][1],I_stable_KirCaLs_fI[1][1],V_stable_KirCaLs_fI[1][1],I1_KirCaLs_fI[1,1],I2_KirCaLs_fI[1,1],V_I1_KirCaLs_fI[1,1],V_I2_KirCaLs_fI[1,1],[],NaN,NaN,NaN,(-3,1),:gray)

plt_V_as_CaL,plt_I_as_CaL = plt_V_mult_steps(reverse(t_asc_step_CaLs),reverse(V_asc_step_CaLs),reverse(I_asc_step_CaLs),[palette(:twelvebitrainbow)[4],palette(:twelvebitrainbow)[3],palette(:twelvebitrainbow)[2],palette(:twelvebitrainbow)[1]],(-3,1),ms_resting,-185:25:50)
plt_VI_as = Plots.plot(title_new_CaL_A,empty_plt,plt_V_as_CaL[1],plt_V_as_CaL[2],plt_V_as_CaL[3],plt_V_as_CaL[4],empty_plt,plt_I_as_CaL,layout=@layout[a{0.001h};b{0.001h};c;d;e;f;g{0.001h};h{0.15h}],size=(600,1000))
plt_V_des_CaL,plt_I_des_CaL = plt_V_mult_steps(t_desc_step_CaLs,V_desc_step_CaLs,I_desc_step_CaLs,[palette(:twelvebitrainbow)[4],palette(:twelvebitrainbow)[3],palette(:twelvebitrainbow)[2],palette(:twelvebitrainbow)[1]],(-3,1),ms_spiking,-185:25:50)
plt_VI_des = Plots.plot(title_new_CaL_B,empty_plt,plt_V_des_CaL[1],plt_V_des_CaL[2],plt_V_des_CaL[3],plt_V_des_CaL[4],empty_plt,plt_I_des_CaL,layout=@layout[a{0.001h};b{0.001h};c;d;e;f;g{0.001h};h{0.15h}],size=(600,1000))
plt_VI_as_des_CaL = Plots.plot(plt_VI_as,empty_plt,plt_VI_des,layout=@layout[a b{0.01w} c],size=(1200,1000))

plt_fI_VI_CaL = Plots.plot(title_new_CaL_C,empty_plt,plt_fI_VI_CaLs_fi,empty_plt,plt_fI_VI_CaLs_vi,layout=@layout[a{0.001h} ; b{0.01h}; c ; d{0.03h} ; e ])
plt_VI_as_des_fI_VI_CaL = Plots.plot(empty_plt,plt_VI_as_des_CaL,empty_plt,empty_plt,plt_fI_VI_CaL,layout=@layout[a{0.01w} b{0.6w} c{0.001w} d{0.001w} e])
plt_new_CaL = Plots.plot(empty_plt,Plots.plot(plt_VI_as_des_fI_VI_CaL,empty_plt,layout=@layout[a b{0.001w}]),empty_plt,layout=@layout[a{0.001h} ; b ; c{0.001h}],size=(2300,1500),fontfamily="Computer Modern",foreground_color=:gray25)

#savefig(plt_new_CaL,"tests/figures/pdf/fig-1.pdf")
