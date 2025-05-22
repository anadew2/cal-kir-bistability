using Plots,LaTeXStrings
include("gates.jl")

plotly()

V=-100:0.1:50
plt_mCaL = plot(fontfamily="Computer Modern", xlabel="V",ylabel=L"m_{\infty}",title="CaL",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=15,titlefontsize=15,yticks=0:0.5:1,xticks=-100:25:50,xlims=(-100,50),size=(800,400))
plot!(V,mCaLf_LeFrancLeMasson2010.(V),label="LeFranc & LeMasson2010-fast",lw=3,lc=palette(:Paired_10)[5])
plot!(V,mCaLf.(V),label="LeFranc & LeMasson2010-fast-boltz",lw=3,lc=:hotpink1)
plot!(V,hCaLf_LeFrancLeMasson2010.(V),label="LeFranc & LeMasson2010-fast",lw=3,lc=palette(:Paired_10)[5],ls=:dash)
plot!(V,mCaLs_LeFrancLeMasson2010.(V),label="LeFranc & LeMasson2010-slow",lw=3,lc=palette(:Paired_10)[6])
plot!(V,mCaLs.(V),label="LeFranc & LeMasson2010-slow_BOLTZ3",lw=3,lc=:chocolate3)
plot!(V,hCaLs.(V),label="LeFranc & LeMasson2010-slow",lw=3,lc=palette(:Paired_10)[6],ls=:dash)

plot!(V,mKir.(V))