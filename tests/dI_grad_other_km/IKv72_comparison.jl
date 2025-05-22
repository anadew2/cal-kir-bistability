using Plots,ColorSchemes,LaTeXStrings
using Unitful, Latexify, UnitfulLatexify

## Steady-state currents comparison
V_vec = -150:0.1:50

include("../../currents.jl")

    IKM_used = IKM.(1,mKM.(V_vec),V_vec,-90)
    mKM_used = mKM.(V_vec)

include("currents_other_km.jl")

    IKM_root = IKM.(1,mKM.(V_vec),V_vec,NaN)
    mKM_root = mKM.(V_vec)

plot(V_vec,IKM_used,label="IKM used",fontfamily="Computer Modern",xlabel="V",ylabel=L"I_{\infty}",legend=:bottomright,lc=:purple)
plot!(V_vec,IKM_root,label="IKM root",lc=:red)
plot!(V_vec,IKM_root.*5,label="IKM root *5",lc=:pink,ylims=(-10,150))

plot(V_vec,mKM_used,label="mKM used",fontfamily="Computer Modern",xlabel="V",ylabel=L"m_{\infty}",legend=:bottomright,lc=:purple)
plot!(V_vec,mKM_root,label="mKM root",lc=:red)

## Bifurcation diagrams comparison (fixed points)
    # Parameters
        ## Fixed parameters
        C = 1 #µF/cm²
        eNa = 50 # [mV]
        eK = -90 # [mV]
        eCa = 120 # [mV]
        eleak = -60.1 # [mV]
        Ca_o = 2 # Extracellular Ca concentration [mM]
    
        ## Intracellular calcium dynamics parameters
        k = 1.0e7 # [nm.cm^(-1)] #TYPO
        F = 96520 # Faraday default value [C/mol]
        d = 1 # [nm]
        Ca_i_0 = 5.0e-5 # Initial intracellular Ca concentration [mM]
        tau_Ca = 10 # [ms]
    
        p_fixed = [C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]
    
        ## Varying parameters
        I(t) = 0.     # [µA/cm^2]
        gNa = 30 #mS/cm²
        gKDR = 4 #mS/cm²
        gleak = 0.03268 #mS/cm²
        gKir = 0.           # [µA/cm^2]
        gKM = 0.            # [µA/cm^2]
        pCaLf = 0. /1000    # [cm/s]
        pCaLs = 0.015/1000   # [cm/s]
        
        p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    
include("../../fixedpoints.jl")
find_fixed_points_ = find_fixed_points([p_var,p_fixed],range(-200.,stop=-20.,step=0.1),-4,0)

    p_var = [I,gNa,gKDR,gKir,0.3,pCaLf,pCaLs,gleak]
find_fixed_points_ = find_fixed_points([p_var,p_fixed],range(-200.,stop=-20.,step=0.1),-4,0)

    p_var = [I,gNa,gKDR,gKir,1.5,pCaLf,pCaLs,gleak]
include("fixedpoints_other_km.jl")
find_fixed_points_ = find_fixed_points([p_var,p_fixed],range(-200.,stop=60.,step=0.1),-4,0)
