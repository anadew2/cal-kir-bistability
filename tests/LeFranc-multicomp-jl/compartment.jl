using Plots,LaTeXStrings,DifferentialEquations

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

function membrane(x,p)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[1]

    ## Varying parameters 
	gNa,gKDR,gKir,gKM,gKCa,gCaAN,pCaLf,pCaLs,gleak = p[2]

    ## State variables 
    V_,mNa_,hNa_,mKDR_,mCaLf_,hCaLf_,mCaLs_,hCaLs_,mKir_,mKM_,mKCa_,mCaAN_,Ca_=x

    ## Current computation
    INa_ = INa(gNa,mNa_,hNa_,V_,eNa)
    IKDR_ = IKDR(gKDR,mKDR_,V_,eK)
    ICaLf_ = ICaLf(pCaLf,mCaLf_,hCaLf_,V_,Ca_,Ca_o)
    ICaLs_ = ICaLs(pCaLs,mCaLs_,hCaLs_,V_,Ca_,Ca_o)
    IKir_ = IKir(gKir,mKir_,V_,eK)
    IKCa_ = IKCa(gKCa,mKCa_,V_,eK)
    ICaAN_ = ICaAN(gCaAN,mCaAN_,V_,eCaAN)
    IKM_ = IKM(gKM,mKM_,V_,eK)
    Ileak_ = Ileak(gleak,V_,eleak)

    Iion_ = INa_ +IKDR_ +ICaLf_ +ICaLs_ +IKir_ +IKM_ +IKCa_ +ICaAN_ +Ileak_  # [µA/cm²]

    dV_ = (-Iion_)/C  # [mV/ms]

	dmNa_ = (mNa.(V_)-mNa_)/(tau_mNa.(V_))
	dhNa_ = (hNa.(V_)-hNa_)/(tau_hNa.(V_))
	dmKDR_ = (mKDR.(V_)-mKDR_)/(tau_mKDR.(V_))
    dmCaLf_ = (mCaLf.(V_)-mCaLf_)/(tau_mCaLf.(V_))
	dhCaLf_ = (hCaLf.(V_)-hCaLf_)/(tau_hCaLf.(V_))
	dmCaLs_ = (mCaLs.(V_)-mCaLs_)/(tau_mCaLs.(V_))
	dhCaLs_ = (hCaLs.(V_)-hCaLs_)/(tau_hCaLs.(V_))
	dmKir_ = (mKir.(V_)-mKir_)/(tau_mKir.(V_))
    dmKM_ = (mKM.(V_)-mKM_)/(tau_mKM.(V_))
    dmKCa_ = (mKCa.(Ca_)-mKCa_)/(tau_mKCa.(Ca_))
    dmCaAN_ = (mCAN.(Ca_)-mCaAN_)/(tau_mCAN.(Ca_))

    ICa_ = (ICaLf_ + ICaLs_)/1000  # [mA/cm^2] instead of [µA/cm^2]
    dCa_ = - ICa_ * k /(2*F*d) - ((Ca_-Ca_i_0)/tau_Ca)

    dx = [dV_,dmNa_,dhNa_,dmKDR_,dmCaLf_,dhCaLf_,dmCaLs_,dhCaLs_,dmKir_,dmKM_,dmKCa_,dmCaAN_,dCa_]
    return dx
end

function compartment!(du,u,p,t)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra = p[1]
    Inull(t) = 0
    p_fixed = [Inull,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]

    ## Varying parameters -- dendrite
	p_dend= p[2]

    ## Varying parameters -- soma
	p_soma = p[3]

    ## Longitudinal resistance between compartments
    Rd = 10^(-2) * Ra * 1/2 * (Ld)/(pi*(Dd/2)^2) #[Mohm]
    Rs = 10^(-2) * Ra * 1/2 * (Ls)/(pi*(Ds/2)^2) #[Mohm]
    Rds = Rd +Rs #[Mohm]

    ## Surface areas
    Ss = 2*pi*(Ds/2)*Ls *10^(-8) # [cm²]
    Sd = 2*pi*(Dd/2)*Ld *10^(-8) # [cm²]

    ## ------- Dendrite -------
    ## Gradients 
	dVd_,du[2],du[3],du[4],du[5],du[6],du[7],du[8],du[9],du[10],du[11],du[12],du[13] = membrane(u[1:13],[p_fixed,p_dend])

    ## ------- Soma -------
    ## Gradients 
	dVs_,du[15],du[16],du[17],du[18],du[19],du[20],du[21],du[22],du[23],du[24],du[25],du[26] = membrane(u[14:26],[p_fixed,p_soma])


    ## ------- Compartment -------
    Vd = u[1]
    Vs = u[14]
    ## Current computation
    Ids = (10^(-3) * (Vd -Vs)/Rds) #[µA]
    Ids_dir = Ids

    du[1]=  dVd_ -Ids_dir/(C*Sd)
    du[14]= dVs_ +Ids_dir/(C*Ss) + I(t)/(C*Ss)
end

function plt_comp_all(sol,saveat_range,p_fixed)
    It_ = p_fixed[1].(sol.t) *10^6 #[pA]
    println("Minimum I: $(minimum(It_)) ; Maximum I: $(maximum(It_))")
    plt_I = plot(sol.t,It_,ylabel="I [pA]",lw=2,ylims=(minimum(It_)-0.1,maximum(It_)+0.1),lc=:mediumpurple)
    plt_Vd = plot(sol.t,sol[1,:],ylabel="Vd [mV]",ylims=(-90,50),lc=:sienna1)
    plt_Vs = plot(sol.t,sol[14,:],ylabel="Vs [mV]",ylims=(-90,50))

    tspk_sol = setdiff(sol.t,saveat_range)
    plt_f = plot()
    if length(tspk_sol)>1
        scatter!(plt_Vs,tspk_sol,0 .*tspk_sol,markerstrokewidth=0,markersize=2,mc=:royalblue4)
        scatter!(plt_f,tspk_sol[2:end],1000 ./(tspk_sol[2:end]-tspk_sol[1:end-1]),ylabel="f [Hz]",ylims=(-1,100),xlims=(-500,sol.t[end]),markerstrokewidth=0,markersize=2,mc=:royalblue4)

        return plot(plt_f,plt_Vs,plt_Vd,plt_I,layout=@layout[a{0.1h};b{0.45h};c;d{0.1h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    else
        return plot(plt_Vs,plt_Vd,plt_I,layout=@layout[a{0.5h};b;c{0.15h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    end
end
function plt_nrn_compare_all(sol,saveat_range,p_fixed)
    It_ = p_fixed[1].(sol.t) *10^6 #[pA]
    println("Minimum I: $(minimum(It_)) ; Maximum I: $(maximum(It_))")
    plt_I = plot(sol.t,It_,ylabel="I [pA]",lw=2,ylims=(minimum(It_)-0.1,maximum(It_)+0.1),lc=:mediumpurple)
    plt_Vs = plot(sol.t,sol[14,:],ylabel="Vs [mV]",ylims=(-90,50))
    plot!(plt_Vs,sol.t,sol[1,:],ylabel="Vd [mV]",ylims=(-90,50),lc=:sienna1)

    tspk_sol = setdiff(sol.t,saveat_range)
    plt_f = plot()
    if length(tspk_sol)>1
        scatter!(plt_Vs,tspk_sol,0 .*tspk_sol,markerstrokewidth=0,markersize=2,mc=:royalblue4)
        scatter!(plt_f,tspk_sol[2:end],1000 ./(tspk_sol[2:end]-tspk_sol[1:end-1]),ylabel="f [Hz]",ylims=(-1,100),xlims=(-500,sol.t[end]),markerstrokewidth=0,markersize=2,mc=:royalblue4)
        println("final f: $(1000 ./(tspk_sol[end]-tspk_sol[end-1]))")

        return plot(plt_f,plt_Vs,plt_I,layout=@layout[a{0.1h};b{0.45h};c{0.1h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    else
        return plot(plt_Vs,plt_I,layout=@layout[a{0.5h};b{0.15h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    end
end
function plt_comp(sol,sol_cb,p_fixed)
    It_ = p_fixed[1].(sol.t) *10^6 #[pA]
    println("Minimum I: $(minimum(It_)) ; Maximum I: $(maximum(It_))")
    plt_I = plot(sol.t,It_,ylabel="I [pA]",lw=2,ylims=(minimum(It_)-0.1,maximum(It_)+0.1),lc=:mediumpurple)
    plt_Vd = plot(sol.t,sol[1,:],ylabel="Vd [mV]",ylims=(-90,50),lc=:sienna1)
    plt_Vs = plot(sol.t,sol[14,:],ylabel="Vs [mV]",ylims=(-90,50))
    plt_f = plot()
    if length(sol_cb.t)>1
        scatter!(plt_Vs,sol_cb.t,sol_cb[14,:],markerstrokewidth=0,markersize=2,mc=:royalblue4)
        scatter!(plt_f,sol_cb.t[2:end],1000 ./(sol_cb.t[2:end]-sol_cb.t[1:end-1]),ylabel="f [Hz]",ylims=(-1,100),xlims=(-500,sol.t[end]),markerstrokewidth=0,markersize=2,mc=:royalblue4)

        return plot(plt_f,plt_Vs,plt_Vd,plt_I,layout=@layout[a{0.1h};b{0.45h};c;d{0.1h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    else
        return plot(plt_Vs,plt_Vd,plt_I,layout=@layout[a{0.5h};b;c{0.15h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    end
end
function plt_nrn_compare(sol,sol_cb,p_fixed)
    It_ = p_fixed[1].(sol.t) *10^6 #[pA]
    println("Minimum I: $(minimum(It_)) ; Maximum I: $(maximum(It_))")
    plt_I = plot(sol.t,It_,ylabel="I [pA]",lw=2,ylims=(minimum(It_)-0.1,maximum(It_)+0.1),lc=:mediumpurple)
    plt_Vs = plot(sol.t,sol[14,:],ylabel="Vs [mV]",ylims=(-90,50))
    plot!(sol.t,sol[1,:],ylabel="Vd [mV]",ylims=(-90,50),lc=:sienna1)
    plt_f = plot()
    if length(sol_cb.t)>1
        #scatter!(plt_Vs,sol_cb.t,sol_cb[14,:],markerstrokewidth=0,markersize=2,mc=:royalblue4)
        scatter!(plt_f,sol_cb.t[2:end],1000 ./(sol_cb.t[2:end]-sol_cb.t[1:end-1]),ylabel="f [Hz]",ylims=(-1,60),xlims=(-500,sol.t[end]),markerstrokewidth=0,markersize=2,mc=:royalblue4)

        return plot(plt_f,plt_Vs,plt_I,layout=@layout[a{0.1h};b{0.45h};c{0.1h}],legend=:none,fontfamily="Computer Modern",size=(800,800))
    else
        return plot(plt_Vs,plt_I,layout=@layout[a{0.5h};b{0.15h}],legend=:none,fontfamily="Computer Modern",size=(800,600))
    end
end
    ## Fixed parameters
    C = 1 #µF/cm²
    eNa = 50 # [mV]
    eK = -90 # [mV]
    eCa = 120 # [mV]
    eCaAN = 0 # [mV]
    eleak = -60.1 # [mV]
    Ca_o = 2 # Extracellular Ca concentration [mM]

    ## Intracellular calcium dynamics parameters
    k = 1.0e7 # [nm.cm^(-1)] 
    F = 96520 # Faraday default value [C/mol]
    d = 1 # [nm]
    Ca_i_0 = 5.0e-5 # Initial intracellular Ca concentration [mM]
    tau_Ca = 10 # [ms]

    ## Compartment parameters
    Ra = 35.4 # Axial resistivity in [ohm.cm]
    Ld = 850 #[µm]  ;
    Dd = 1 #[µm]
    Ls = 56.279 #[µm] ; 
    Ds = 56.279 #[µm] ; 

    ## Stimulation 
    I(t) = 25/10^(6) + 1*35/10^(6) *pulse(t,5000,8000)    # [µA]

    p_fixed = [I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra]

    ## Initial conditions
    V_ic = -59. # [mV]
    Ca_ic = 5.0e-5 # [mM]
    ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic,
        V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),mKCa(Ca_ic),mCAN(Ca_ic),Ca_ic]

    ## Varying parameters
    gNa_s = 30 #mS/cm²
    gKDR_s = 4 #mS/cm²
    gleak_s = 0.03268 #mS/cm²
    gKir_s = 0.02           # mS/cm² 
    gKM_s = 0.            # mS/cm²
    gKCa_s = 0.03       # mS/cm²
    gCaAN_s = 0.        # mS/cm²
    pCaLf_s = 0.15 /1000    # [cm/s]
    pCaLs_s = 0. /1000   # [cm/s]

    gNa_d = 0 #mS/cm²
    gKDR_d = 0 #mS/cm²
    gleak_d = 0.03 #mS/cm²
    gKir_d = 0.           # mS/cm²
    gKM_d = 0.            # mS/cm²
    gKCa_d = 0.0185  # mS/cm²
    gCaAN_d = 0.125 # mS/cm²
    pCaLf_d = 0. /1000    # [cm/s]
    pCaLs_d = 0.02 /1000   # [cm/s]
    

    p_var_d = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d,gleak_d]
    p_var_s = [gNa_s,gKDR_s,gKir_s,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]



    