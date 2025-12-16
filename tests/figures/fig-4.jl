using Plots,ColorSchemes,LaTeXStrings
using LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using Statistics

include("../dI_grad_altgt/local_postpro_output.jl")
include("../fI_curve_altgt/bif.jl")

function plt_fI_VI(I_list_fI,f_list_fI,I_stable_VI,V_stable_VI,I1,I2,V_I1,V_I2,xlims_,color)
    fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
    fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]
    fontsize_= 25

    plt_fI = plot(fontfamily="Computer Modern", xlabel=" ",ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legend=:outertopright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,size=(800,400))
    if I1 != I2
        plot!(xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]))
        plot!(I1.*ones(2),[-5,250],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        plot!(I2.*ones(2),[-5,250],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        annotate!(I1,250+15,L"I_1",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
        annotate!(I2,250+15,L"I_2",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")

    else
        plot!(xticks=([xlims_[1],I1,xlims_[2]],[xlims_[1],L"I_1=I_2",xlims_[2]]))
    end
    plot!([I1,I2],220 .*ones(2),arrow=true,color=:black,linewidth=0.5,label=:none)
    plot!([I2,I1],220 .*ones(2),arrow=true,color=:black,linewidth=0.5,label=:none,ylims=(-10,250))
    annotate!(mean([I2,I1]),250,L"\Delta I",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    plot!(yticks=[0,250])
    plot!(plt_fI,I_list_fI,f_list_fI,lc=color,lw=2,label="Spiking")
    plot!(I_stable_VI ,zeros(size(I_stable_VI)),ls=:dash,lw=2,lc=color,label="Resting")
    plot!(I_list_fI[fI_fill],f_list_fI[fI_fill],fillrange=zeros(size(fI_fill)),c=color,fillalpha=0.15,lw=0,label="Bistable")

    VI_fill_ = collect(1:length(I_stable_VI))[I_stable_VI .>=I1]
    VI_fill = VI_fill_[I_stable_VI[VI_fill_] .<= I2]  

    plt_VI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legend=:outertopright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-90,-50),yticks=-90:20:-50,xticks=([xlims_[1],I1,I2,xlims_[2]],[xlims_[1],L"I_1",L"I_2",xlims_[2]]),size=(800,400))
        #Display I1 & I2 lines
        plot!(I1.*ones(2),[-150,-50],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        plot!(I2.*ones(2),[-150,-50],lc=:grey25,linestyle=:dot,label=:none,lw=5,linealpha=0.2)
        annotate!(I1,-50+6,L"I_1",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
        annotate!(I2,-50+6,L"I_2",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
    plot!(I_stable_VI ,V_stable_VI,ls=:dash,lw=2,lc=color,label="Resting")
    plot!(I_stable_VI[VI_fill],V_stable_VI[VI_fill],fillrange=ones(size(VI_fill)).*-55,c=color,fillalpha=0.15,lw=0,label="Bistable")
    scatter!(I1.*ones(2), V_I1.*ones(2),mc=color,markerstrokewidth=0,label=:none)
    scatter!(I2.*ones(2), V_I2.*ones(2),mc=color,markerstrokewidth=0,label=:none)    
    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims_)
    return plt
end
function plt_list_fI_VI(I_list_fI_,f_list_fI_,I_stable_VI_,V_stable_VI_,I1_,I2_,V_I1_,V_I2_,xlims__,color_)

    fontsize_= 25
    plt_fI = plot(fontfamily="Computer Modern", xlabel=" ",ylabel=L"f \ \ "*L"[" *latexify(u"Hz") *L"]",grid=:none,legend=:outertopright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-5,250),yticks=0:50:200,size=(800,400))
    plt_VI = plot(fontfamily="Computer Modern", xlabel=L"I \ \ "*L"[" *latexify(u"µA /cm^2") *L"]",ylabel=L"\overline{V} \ \ "*L"[" *latexify(u"mV") *L"]" ,grid=:none,legend=:outertopright,legendfontsize=fontsize_,xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_,ylims=(-90,-50),yticks=-90:10:-50,size=(800,400))

    for i=1:length(I_list_fI_)
        I_list_fI = I_list_fI_[i]
        f_list_fI = f_list_fI_[i]
        I_stable_VI = I_stable_VI_[i]
        V_stable_VI = V_stable_VI_[i]
        I1 = I1_[i]
        I2 = I2_[i]
        V_I1 = V_I1_[i]
        V_I2 = V_I2_[i]
        xlims_ = xlims__[i]
        color = color_[i]

        fI_fill_ = collect(1:length(I_list_fI))[I_list_fI .>=I1]
        fI_fill = fI_fill_[I_list_fI[fI_fill_] .<= I2]

        if i==1
            plot!(plt_fI,(xlims_[1]-100) *ones(2),-100 .*ones(2),lc=:grey,lw=1,label="Spiking")
            plot!(plt_fI,(xlims_[1]-100) *ones(2) ,-100 .*ones(2),ls=:dash,lw=1,lc=:grey,label="Resting")
            plot!(plt_fI,(xlims_[1]-100) *ones(2) ,-100 .*ones(2),fillrange=-100 .*ones(size(fI_fill)),c=:grey,fillalpha=0.15,lw=0,label="Bistable")
        end
        plot!(plt_fI,I_list_fI,f_list_fI,lc=color,lw=2,label=:none) #"Spiking"
        plot!(plt_fI,I_stable_VI ,zeros(size(I_stable_VI)),ls=:dash,lw=2,lc=color,label=:none) #"Resting"
        plot!(plt_fI,I_list_fI[fI_fill],f_list_fI[fI_fill],fillrange=zeros(size(fI_fill)),c=color,fillalpha=0.15,lw=0,label=:none) #"Bistable"
        #plot!(plt_fI,[I1,I2],220 .*ones(2),arrow=true,color=:black,linewidth=0.5,label=:none)
        #plot!(plt_fI,[I2,I1],220 .*ones(2),arrow=true,color=:black,linewidth=0.5,label=:none,ylims=(-10,250))
        #annotate!(plt_fI,mean([I2,I1]),250,L"\Delta I",annotationfontsize=fontsize_,annotationhalign=:center,annotationfontfamily="Computer Modern")
        #plot!(plt_fI,yticks=[0,250])

        VI_fill_ = collect(1:length(I_stable_VI))[I_stable_VI .>=I1]
        VI_fill = VI_fill_[I_stable_VI[VI_fill_] .<= I2]  

        if i==1
            plot!(plt_VI,(xlims_[1]-100) *ones(2) ,-1000 .*ones(2),ls=:dash,lw=1,lc=:grey,label="Resting")
            plot!(plt_VI,(xlims_[1]-100) *ones(2) ,-1000 .*ones(2),fillrange=-1000 .*ones(size(fI_fill)),c=:grey,fillalpha=0.15,lw=0,label="Bistable")
        end
        plot!(plt_VI,I_stable_VI ,V_stable_VI,ls=:dash,lw=2,lc=color,label=:none) #"Resting"
        plot!(plt_VI,I_stable_VI[VI_fill],V_stable_VI[VI_fill],fillrange=ones(size(VI_fill)).*-50,c=color,fillalpha=0.15,lw=0,label=:none) #"Bistable"
        scatter!(plt_VI,I1.*ones(2), V_I1.*ones(2),mc=color,markerstrokewidth=0,label=:none)
        scatter!(plt_VI,I2.*ones(2), V_I2.*ones(2),mc=color,markerstrokewidth=0,label=:none)  

    end
    plt=plot(plt_fI,plt_VI,layout=(2,1),xlims=xlims__[end])
    return plt
end

function plt_dI_grad(dI_grad_mat,gx,gy,clims_dI,xname,yname)
    palette =  :vik#cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true)
    fontsize_ =25

    plt_dI = Plots.heatmap(gx,gy,dI_grad_mat,fill=true,levels=100,lw=0,c=palette,clims=clims_dI)
    #plot!(plt_dI,gx_TC,gy_TC,label=:none,lc=:black,lw=2)
    plot!(colorbar=true,legend=:none)
    yaxis!((minimum(gy),maximum(gy)),yticks=[minimum(gy),maximum(gy)])
    xaxis!((minimum(gx),maximum(gx)),xticks=[minimum(gx),maximum(gx)])
    plot!(ylabel=yname,ylabelfontsize=fontsize_)
    plot!(xlabel=xname,xlabelfontsize=fontsize_)
    plot!(fontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)
    
    plt_annotation = Plots.plot()
    annotate!(0.5,0.52, latexstring("     ",L"\Delta I"),annotationfontsize=fontsize_)
    plot!(yticks =:none)
    plot!(xticks =:none)
    plot!(grid=false,legend=:none)
    plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)

    plt = Plots.plot(plt_dI,plt_annotation,layout=@layout[a b{0.08w}])
    return plt
end
function plt_dI_grad_(dI_grad_mat,gx,gy,clims_dI,xname,yname)
    palette =  :vik#cgrad(:diverging_rainbow_bgymr_45_85_c67_n256,rev=true)
    fontsize_ =25

    plt_dI = Plots.heatmap(gx,gy,dI_grad_mat,fill=true,levels=100,lw=0,c=palette,clims=clims_dI)
    #plot!(plt_dI,gx_TC,gy_TC,label=:none,lc=:black,lw=2)
    plot!(colorbar=true,legend=:none)
    yaxis!((minimum(gy),maximum(gy)),yticks=[minimum(gy),maximum(gy)])
    xaxis!((minimum(gx),maximum(gx)),xticks=[minimum(gx),maximum(gx)])
    plot!(ylabel=yname,ylabelfontsize=fontsize_)
    plot!(xlabel=xname,xlabelfontsize=fontsize_)
    plot!(fontfamily="Computer Modern",tickfontsize=fontsize_)
    plot!(xlabelfontsize=fontsize_,ylabelfontsize=fontsize_,tickfontsize=fontsize_)
    
    plt_annotation = Plots.plot()
    annotate!(0.5,0.57, latexstring("     ",L"\Delta I-\Delta I_0"),annotationfontsize=fontsize_)
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

    plt_dI = Plots.heatmap(gx,gy,Veq_grad_mat,fill=true,levels=100,lw=0,c=palette,clims=clims_V)
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

title_kir_km_2_AB = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(0.,0.5,L"\textbf{A}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
annotate!(0.5,0.5,L"\textbf{B}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_C = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(0.,0.5,L"\textbf{C}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")
title_kir_km_2_D = Plots.plot() 
plot!(grid=false,legend=:none)
plot!(foreground_color_border=:white,showaxis=false,foreground_color_axis=:white)
annotate!(-0.25,0.5,L"\textbf{D}",annotationfontsize=35,annotationhalign=:left,annotationfontfamily="Computer Modern")


clim_ = [0.1,0.2,0.4,0.625,0.75]


i_pCaLs = 5
#plt_dI_KirKMCaLs_5 = plt_dI_grad(dI_KirKMCaLs[i_pCaLs,:,:],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(0,2.75),L"\bar{g}_{\mathrm{KM}}",L"\bar{g}_{\mathrm{Kir}}")
plt_dI_KirKMCaLs_5 = plt_dI_grad_(dI_KirKMCaLs[i_pCaLs,:,:].-dI_KirKMCaLs[i_pCaLs,1,1],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(-1,+1).*clim_[5] ,L"\bar{g}_{\mathrm{KM}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]",L"\bar{g}_{\mathrm{Kir}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]")
#plt_dI_KirKMCaLs_5 = plt_dI_grad(dI_KirKMCaLs[i_pCaLs,:,:],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(-1,+1).*clim_[i_pCaLs] .+dI_KirKMCaLs[i_pCaLs,1,1],L"\bar{g}_{\mathrm{KM}}",L"\bar{g}_{\mathrm{Kir}}")
plot!(plt_dI_KirKMCaLs_5[1],title=L"\bar{p}_{CaLs} = "*"$(round(pCaLs_range_KirKMCaLs[i_pCaLs]*1e7)/10^7)"*latexify(u"cm/s"),titlefontsize=25)
scatter!(plt_dI_KirKMCaLs_5,[gKM_range_KirKMCaLs[3]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:palevioletred2,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_dI_KirKMCaLs_5,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[3]],markershape=:circle,markersize=18,mc=:mediumpurple1,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_dI_KirKMCaLs_5,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:goldenrod2,markerstrokewidth=0.2,markerstrokecolor=:white)


plt_V_I1_KirKMCaLs_5 = plt_Veq_grad(V_I1_KirKMCaLs[i_pCaLs,:,:],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(-125,-55),L"\bar{g}_{\mathrm{KM}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]",L"\bar{g}_{\mathrm{Kir}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]")
plot!(plt_V_I1_KirKMCaLs_5[1],title=L"\bar{p}_{CaLs} = "*"$(round(pCaLs_range_KirKMCaLs[i_pCaLs]*1e7)/10^7)"*latexify(u"cm/s"),titlefontsize=25)
scatter!(plt_V_I1_KirKMCaLs_5,[gKM_range_KirKMCaLs[3]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:palevioletred2,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_V_I1_KirKMCaLs_5,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[3]],markershape=:circle,markersize=18,mc=:mediumpurple1,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_V_I1_KirKMCaLs_5,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:goldenrod2,markerstrokewidth=0.2,markerstrokecolor=:white)

i_pCaLs_fI =2
I_fI_CaLs_2 = [I_list_KirKMCaLs_fI[i_pCaLs_fI][2][1],I_list_KirKMCaLs_fI[i_pCaLs_fI][1][2],I_list_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
f_fI_CaLs_2 = [f_list_KirKMCaLs_fI[i_pCaLs_fI][2][1],f_list_KirKMCaLs_fI[i_pCaLs_fI][1][2],f_list_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
I_VI_CaLs_2 = [I_stable_KirKMCaLs_fI[i_pCaLs_fI][2][1],I_stable_KirKMCaLs_fI[i_pCaLs_fI][1][2],I_stable_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
V_VI_CaLs_2 = [V_stable_KirKMCaLs_fI[i_pCaLs_fI][2][1],V_stable_KirKMCaLs_fI[i_pCaLs_fI][1][2],V_stable_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
I1_CaLs_2 = [I1_KirKMCaLs_fI[i_pCaLs_fI,2,1],I1_KirKMCaLs_fI[i_pCaLs_fI,1,2],I1_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
I2_CaLs_2 = [I2_KirKMCaLs_fI[i_pCaLs_fI,2,1],I2_KirKMCaLs_fI[i_pCaLs_fI,1,2],I2_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
V_I1_CaLs_2 = [V_I1_KirKMCaLs_fI[i_pCaLs_fI,2,1],V_I1_KirKMCaLs_fI[i_pCaLs_fI,1,2],V_I1_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
V_I2_CaLs_2 = [V_I2_KirKMCaLs_fI[i_pCaLs_fI,2,1],V_I2_KirKMCaLs_fI[i_pCaLs_fI,1,2],V_I2_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
xlims__CaLs_2 = [(0,4.5+5.5),(0,4.5+5.5),(0,4.5+5.5)]
color_CaLs_2 = [:palevioletred2,:mediumpurple1,:goldenrod2]

I_fI_CaLs_2 = reverse(I_fI_CaLs_2)
f_fI_CaLs_2 = reverse(f_fI_CaLs_2)
I_VI_CaLs_2 = reverse(I_VI_CaLs_2)
V_VI_CaLs_2 = reverse(V_VI_CaLs_2)
I1_CaLs_2 = reverse(I1_CaLs_2)
I2_CaLs_2 = reverse(I2_CaLs_2)
V_I1_CaLs_2 = reverse(V_I1_CaLs_2)
V_I2_CaLs_2 = reverse(V_I2_CaLs_2)
xlims__CaLs_2 = reverse(xlims__CaLs_2)
color_CaLs_2 = reverse(color_CaLs_2)

plt_fI_VI_list_CaLs_2 = plt_list_fI_VI(I_fI_CaLs_2,f_fI_CaLs_2,I_VI_CaLs_2,V_VI_CaLs_2,I1_CaLs_2,I2_CaLs_2,V_I1_CaLs_2,V_I2_CaLs_2,xlims__CaLs_2,color_CaLs_2)

## Test 
#==
plt_dI = plot(legend=:outertopright,xlabel="gKir")
for i in 1:1:31
    plot!(plt_dI,dI_KirKMCaLs[i_pCaLs,:,i].-dI_KirKMCaLs[i_pCaLs,1,1],label="i=$i")
end
plot(plt_dI,size=(500,400))==#
#==
plt_dI = plot(legend=:outertopright,xlabel="gKM")
for i in 1:1:31
    plot!(plt_dI,dI_KirKMCaLs[1,i,:].-dI_KirKMCaLs[i_pCaLs,1,1],label="i=$i")
end
plot(plt_dI,size=(500,400))==#


i_pCaLs = 3
#plt_dI_KirKMCaLs_3 = plt_dI_grad(dI_KirKMCaLs[i_pCaLs,:,:],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(0,2.75),L"\bar{g}_{\mathrm{KM}}",L"\bar{g}_{\mathrm{Kir}}")
plt_dI_KirKMCaLs_3 = plt_dI_grad_(dI_KirKMCaLs[i_pCaLs,:,:].-dI_KirKMCaLs[i_pCaLs,1,1],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(-1,+1).*clim_[5] ,L"\bar{g}_{\mathrm{KM}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]",L"\bar{g}_{\mathrm{Kir}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]")
#plt_dI_KirKMCaLs_3 = plt_dI_grad(dI_KirKMCaLs[i_pCaLs,:,:],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(-1,+1).*clim_[i_pCaLs] .+dI_KirKMCaLs[i_pCaLs,1,1],L"\bar{g}_{\mathrm{KM}}",L"\bar{g}_{\mathrm{Kir}}")
plot!(plt_dI_KirKMCaLs_3[1],title=L"\bar{p}_{CaLs} = "*"$(round(pCaLs_range_KirKMCaLs[i_pCaLs]*1e7)/10^7)"*latexify(u"cm/s"),titlefontsize=25)
scatter!(plt_dI_KirKMCaLs_3,[gKM_range_KirKMCaLs[3]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:green,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_dI_KirKMCaLs_3,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[3]],markershape=:circle,markersize=18,mc=:deepskyblue3,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_dI_KirKMCaLs_3,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:gray35,markerstrokewidth=0.2,markerstrokecolor=:white)

plt_V_I1_KirKMCaLs_3 = plt_Veq_grad(V_I1_KirKMCaLs[i_pCaLs,:,:],gKM_range_KirKMCaLs,gKir_range_KirKMCaLs,(-125,-55),L"\bar{g}_{\mathrm{KM}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]",L"\bar{g}_{\mathrm{Kir}} \ \ "*L"[" *latexify(u"mS /cm^2") *L"]")
plot!(plt_V_I1_KirKMCaLs_3[1],title=L"\bar{p}_{CaLs} = "*"$(round(pCaLs_range_KirKMCaLs[i_pCaLs]*1e7)/10^7)"*latexify(u"cm/s"),titlefontsize=25)
scatter!(plt_V_I1_KirKMCaLs_3,[gKM_range_KirKMCaLs[3]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:green,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_V_I1_KirKMCaLs_3,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[3]],markershape=:circle,markersize=18,mc=:deepskyblue3,markerstrokewidth=0.2,markerstrokecolor=:white)
scatter!(plt_V_I1_KirKMCaLs_3,[gKM_range_KirKMCaLs[29]],[gKir_range_KirKMCaLs[29]],markershape=:circle,markersize=18,mc=:gray35,markerstrokewidth=0.2,markerstrokecolor=:white)

i_pCaLs_fI =1
I_fI_CaLs_1 = [I_list_KirKMCaLs_fI[i_pCaLs_fI][2][1],I_list_KirKMCaLs_fI[i_pCaLs_fI][1][2],I_list_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
f_fI_CaLs_1 = [f_list_KirKMCaLs_fI[i_pCaLs_fI][2][1],f_list_KirKMCaLs_fI[i_pCaLs_fI][1][2],f_list_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
I_VI_CaLs_1 = [I_stable_KirKMCaLs_fI[i_pCaLs_fI][2][1],I_stable_KirKMCaLs_fI[i_pCaLs_fI][1][2],I_stable_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
V_VI_CaLs_1 = [V_stable_KirKMCaLs_fI[i_pCaLs_fI][2][1],V_stable_KirKMCaLs_fI[i_pCaLs_fI][1][2],V_stable_KirKMCaLs_fI[i_pCaLs_fI][2][2]]
I1_CaLs_1 = [I1_KirKMCaLs_fI[i_pCaLs_fI,2,1],I1_KirKMCaLs_fI[i_pCaLs_fI,1,2],I1_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
I2_CaLs_1 = [I2_KirKMCaLs_fI[i_pCaLs_fI,2,1],I2_KirKMCaLs_fI[i_pCaLs_fI,1,2],I2_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
V_I1_CaLs_1 = [V_I1_KirKMCaLs_fI[i_pCaLs_fI,2,1],V_I1_KirKMCaLs_fI[i_pCaLs_fI,1,2],V_I1_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
V_I2_CaLs_1 = [V_I2_KirKMCaLs_fI[i_pCaLs_fI,2,1],V_I2_KirKMCaLs_fI[i_pCaLs_fI,1,2],V_I2_KirKMCaLs_fI[i_pCaLs_fI,2,2]]
#xlims__CaLs_1 = [(2.5,4+5.5),(2.5,4+5.5),(2.5,4+5.5)]
xlims__CaLs_1 = [(0,4.5+5.5),(0,4.5+5.5),(0,4.5+5.5)]
color_CaLs_1 = [:green,:deepskyblue3,:gray35]

I_fI_CaLs_1 = reverse(I_fI_CaLs_1)
f_fI_CaLs_1 = reverse(f_fI_CaLs_1)
I_VI_CaLs_1 = reverse(I_VI_CaLs_1)
V_VI_CaLs_1 = reverse(V_VI_CaLs_1)
I1_CaLs_1 = reverse(I1_CaLs_1)
I2_CaLs_1 = reverse(I2_CaLs_1)
V_I1_CaLs_1 = reverse(V_I1_CaLs_1)
V_I2_CaLs_1 = reverse(V_I2_CaLs_1)
xlims__CaLs_1 = reverse(xlims__CaLs_1)
color_CaLs_1 = reverse(color_CaLs_1)

plt_fI_VI_list_CaLs_1 = plt_list_fI_VI(I_fI_CaLs_1,f_fI_CaLs_1,I_VI_CaLs_1,V_VI_CaLs_1,I1_CaLs_1,I2_CaLs_1,V_I1_CaLs_1,V_I2_CaLs_1,xlims__CaLs_1,color_CaLs_1)


plt_fI_VI_list_2 =Plots.plot(empty_plt,plt_fI_VI_list_CaLs_2,empty_plt,layout=@layout[a{0.001w} b c{0.001w}])
plt_fI_VI_list_1 =Plots.plot(empty_plt,plt_fI_VI_list_CaLs_1,empty_plt,layout=@layout[a{0.001w} b c{0.001w}])


plt_dI_V_I1_KirKMCaLs_3 = Plots.plot(empty_plt,plt_dI_KirKMCaLs_3,empty_plt,plt_V_I1_KirKMCaLs_3,empty_plt,layout=@layout[a{0.001w}  b  c{0.03w}  d  e{0.001w} ],size=(2300,700),fontfamily="Computer Modern")
plt_3_B2 = Plots.plot(empty_plt,title_kir_km_2_AB,plt_dI_V_I1_KirKMCaLs_3,title_kir_km_2_C,plt_fI_VI_list_1,empty_plt,layout=@layout[a{0.001h} ;b{0.001h}; c{0.4h} ; d{0.001h} ; e ;f{0.001h} ],size=(2300,1800),fontfamily="Computer Modern",foreground_color=:gray25)
#Plots.savefig(plt_3_B2,"tests/figures/pdf/fig-3-B2.pdf")


plt_dI_V_I1_KirKMCaLs_5 = Plots.plot(empty_plt,plt_dI_KirKMCaLs_5,empty_plt,plt_V_I1_KirKMCaLs_5,empty_plt,layout=@layout[a{0.001w}  b  c{0.03w}  d  e{0.001w} ],size=(2300,700),fontfamily="Computer Modern")
plt_3_B1 = Plots.plot(empty_plt,title_kir_km_2_AB,plt_dI_V_I1_KirKMCaLs_5,title_kir_km_2_C,plt_fI_VI_list_2,empty_plt,layout=@layout[a{0.001h} ;b{0.001h}; c{0.4h} ; d{0.001h} ; e ;f{0.001h} ],size=(2300,1800),fontfamily="Computer Modern",foreground_color=:gray25)
#Plots.savefig(plt_3_B1,"tests/figures/pdf/fig-3-B1.pdf")


