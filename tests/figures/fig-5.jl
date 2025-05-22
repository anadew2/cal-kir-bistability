using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify
using PlotlyJS, CSV, DataFrames
using JLD

include("../dI_grad/local_postpro_output.jl")
include("../dI_var_mc/local_postpro_outputs.jl");

trace_sc_Kir3 = PlotlyJS.scatter(x=(C_mat_KirCaLs_var30.-C)./C,y=(dI_KirCaLs_var30.-dI_KirCaLs[41,51])./ dI_KirCaLs[41,51], mode="markers",marker_color=palette(:PuRd_7,rev=true)[2],legendgroup="groupKir",legendgrouptitle_text="Kir",name="30%",yaxis_tickvals = -0.3:0.1:0.3,yaxis_ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"])
trace_sc_KM3 = PlotlyJS.scatter(x=(C_mat_KMCaLs_var30.-C)./C,y=(dI_KMCaLs_var30.-dI_KMCaLs[41,51])./ dI_KMCaLs[41,51] , mode="markers",marker_color=palette(:Purples_7,rev=true)[2],legendgroup="groupKM",legendgrouptitle_text="KM",name="30%")
trace_sc_Kir2 = PlotlyJS.scatter(x=(C_mat_KirCaLs_var20.-C)./C,y=(dI_KirCaLs_var20.-dI_KirCaLs[41,51])./ dI_KirCaLs[41,51] , mode="markers",marker_color=palette(:PuRd_7,rev=true)[3],legendgroup="groupKir",legendgrouptitle_text="Kir",name="20%")
trace_sc_KM2 = PlotlyJS.scatter(x=(C_mat_KMCaLs_var20.-C)./C,y=(dI_KMCaLs_var20.-dI_KMCaLs[41,51])./ dI_KMCaLs[41,51] , mode="markers",marker_color=palette(:Purples_7,rev=true)[3],legendgroup="groupKM",legendgrouptitle_text="KM",name="20%")
trace_sc_Kir1 = PlotlyJS.scatter(x=(C_mat_KirCaLs_var10.-C)./C,y=(dI_KirCaLs_var10.-dI_KirCaLs[41,51])./ dI_KirCaLs[41,51] , mode="markers",marker_color=palette(:PuRd_7,rev=true)[4],legendgroup="groupKir",legendgrouptitle_text="Kir",name="10%")
trace_sc_KM1 = PlotlyJS.scatter(x=(C_mat_KMCaLs_var10.-C)./C,y=(dI_KMCaLs_var10.-dI_KMCaLs[41,51])./ dI_KMCaLs[41,51], mode="markers",marker_color=palette(:Purples_7,rev=true)[4],legendgroup="groupKM",legendgrouptitle_text="KM",name="10%")

trace_Kirfull = box(
    x=vcat("10%".*replace.(string.(Int.(zeros(size(dI_KirCaLs_var10)))),"0" => ""),"20%".*replace.(string.(Int.(zeros(size(dI_KirCaLs_var20)))),"0" => ""),"30%".*replace.(string.(Int.(zeros(size(dI_KirCaLs_var30)))),"0" => "")),
    y=vcat((dI_KirCaLs_var10.-dI_KirCaLs[41,51])./ dI_KirCaLs[41,51],(dI_KirCaLs_var20.-dI_KirCaLs[41,51])./ dI_KirCaLs[41,51],(dI_KirCaLs_var30.-dI_KirCaLs[41,51])./ dI_KirCaLs[41,51]),
    legendgroup="group",  # this can be any string, not just "group"
    legendgrouptitle_text="Model used",
    name="Kir",
    showlegend=true,
    boxpoints="all",
    marker_color=palette(:PuRd_7,rev=true)[2], #RGB(0.95,0.62,0.75),
    fillcolor=palette(:PuRd_7,rev=true)[4],
    #boxpoints="all",
    whiskerwidth=0.4,
    marker_size = 3,
    line_width = 1.5,
    boxmean=true # represent mean
    )
trace_KMfull = box(
    x=vcat("10%".*replace.(string.(Int.(zeros(size(dI_KMCaLs_var10)))),"0" => ""),"20%".*replace.(string.(Int.(zeros(size(dI_KMCaLs_var20)))),"0" => ""),"30%".*replace.(string.(Int.(zeros(size(dI_KMCaLs_var30)))),"0" => "")),
    y=vcat((dI_KMCaLs_var10.-dI_KMCaLs[41,51])./ dI_KMCaLs[41,51],(dI_KMCaLs_var20.-dI_KMCaLs[41,51])./ dI_KMCaLs[41,51],(dI_KMCaLs_var30.-dI_KMCaLs[41,51])./ dI_KMCaLs[41,51]),
    legendgroup="group",
    name="KM",
    showlegend=true,
    marker_color=palette(:Purples_7,rev=true)[2], #RGB(0.95,0.62,0.75),
    fillcolor=palette(:Purples_7,rev=true)[4],
    boxpoints="all",
    whiskerwidth=0.4,
    marker_size = 3,
    line_width = 1.5,
    boxmean=true # represent mean
    )

    layout = PlotlyJS.Layout(
        title="",
        yaxis=attr(
                title="Change in level of bistability",
                showgrid=true,
                gridcolor="rgb(243,243,243)",
                #tickvals=0:0.01:0.1,
                autorange=true,
                dtick=0.1,
                zeroline=true,
                zerolinecolor="rgb(243,243,243)",
                zerolinewidth=1,
                tickmode = "array",
                tickvals = -0.3:0.1:0.3,
                ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"],
                range=[-0.3,0.3]
            ),
        xaxis=attr(title="Model used"),    
        margin=attr(
                l=40,
                r=30,
                b=80,
                t=80
            ),
        paper_bgcolor="rgb(255,255,255)",
        plot_bgcolor="rgb(255,255,255)",
        font_family="Computer Modern",
        font_size=20,
        showlegend=true, 
        boxmode="group"
    )

    layout_sc = PlotlyJS.Layout(
        title="",
        yaxis=attr(
                showgrid=false,
                gridcolor="rgb(243,243,243)",
                #tickvals=0:0.01:0.1,
                autorange=true,
                dtick=0.1,
                zeroline=true,
                zerolinecolor="rgb(255,255,255)",
                zerolinewidth=0,
                tickmode = "array",
                tickvals = -0.3:0.1:0.3,
                ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"],
                range=[-0.3,0.3]
            ),
        xaxis=attr(
                title="Change in C",
                tickvals = -0.3:0.1:0.3,
                ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"]
            ),       
        margin=attr(
                l=40,
                r=30,
                b=80,
                t=80
            ),
        paper_bgcolor="rgb(255,255,255)",
        plot_bgcolor="rgb(255,255,255)",
        font_family="Computer Modern",
        font_size=20,
        showlegend=true
    )

scplots1 = [trace_sc_Kir3, trace_sc_Kir2, trace_sc_Kir1, trace_sc_KM3, trace_sc_KM2, trace_sc_KM1]
pl1=PlotlyJS.plot(scplots1,layout_sc)
boxplots2=[trace_Kirfull, trace_KMfull]
pl2=Plot(boxplots2, layout)
fig = make_subplots(rows=1, cols=2,
specs=reshape([Spec(kind="scatter")
            Spec(kind="box")], 1,2))
            
    add_trace!(fig, trace_sc_Kir3, row=1, col=1)
    add_trace!(fig, trace_sc_Kir2, row=1, col=1)
    add_trace!(fig, trace_sc_Kir1, row=1, col=1)
    add_trace!(fig, trace_sc_KM3, row=1, col=1)
    add_trace!(fig, trace_sc_KM2, row=1, col=1)
    add_trace!(fig, trace_sc_KM1, row=1, col=1)

    for k=1:2
        add_trace!(fig, pl2.data[k], row=1, col=2)
    end 
relayout!(fig, boxmode="group",
    yaxis_tickvals = -0.3:0.1:0.3,yaxis_ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"],
    yaxis2_tickvals = -0.3:0.1:0.3,yaxis2_ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"],
    yaxis_range=[-0.3,0.3],
    yaxis2_range=[-0.3,0.3],
    xaxis_tickvals = -0.3:0.1:0.3,xaxis_ticktext = ["-30%", "-20%", "-10%", "0%", "10%", "20%", "30%"],
    xaxis_title="Change in C",
    xaxis2_title="Intrinsic variability",
    yaxis_title="Change in level of bistability",
    yaxis2_title="Change in level of bistability",
    paper_bgcolor="rgb(255,255,255)",
    plot_bgcolor="rgb(255,255,255)",
    font_family="Computer Modern",
    width=1200,
    height=600,
    font_size=20
)
display(fig)
#PlotlyJS.savefig(fig,"ChannelUpdate/cluster/figures/pdf/fig-5.pdf",width=1000,height=500) #from local_postpro_outputs_71_94
