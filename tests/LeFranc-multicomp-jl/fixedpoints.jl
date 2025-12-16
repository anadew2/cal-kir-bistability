using NLsolve,LinearAlgebra
using ForwardDiff

function nlsolve_comp_grad_0_I_!(F,x,p_)
    ## Fixed parameters
    I_,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,Faraday,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra = p_[1]
    Inull = 0
    p_fixed = [Inull,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,Faraday,d,Ca_i_0,tau_Ca]
        ## Longitudinal resistance between compartments
        Rd = 10^(-2) * Ra * 1/2 * (Ld)/(pi*(Dd/2)^2) #[Mohm]
        Rs = 10^(-2) * Ra * 1/2 * (Ls)/(pi*(Ds/2)^2) #[Mohm]
        Rds = Rd +Rs #[Mohm]
        ## Surface areas
        Ss = 2*pi*(Ds/2)*Ls *10^(-8) # [cm²]
        Sd = 2*pi*(Dd/2)*Ld *10^(-8) # [cm²]

        # todo : verify this 
        Cd = p_fixed[2] #same capacitance per surface area for soma & dendrite
        Cs = p_fixed[2] #same capacitance per surface area for soma & dendrite

    ## Varying parameters -- dendrite
	p_dend= p_[2]

    ## Varying parameters -- soma
	p_soma = p_[3]

    ## NLSolver parameter
	Vs = p_[4][1]

    I = x[1] #[µA]
    Vd = x[2]
    Cad = x[3]
    Cas = x[4]

    x_dend = [Vd,mNa(Vd),hNa(Vd),mKDR(Vd),mCaLf(Vd),hCaLf(Vd),mCaLs(Vd),hCaLs(Vd),mKir(Vd),mKM(Vd),mKCa(Cad),mCAN(Cad),Cad]
    x_soma = [Vs,mNa(Vs),hNa(Vs),mKDR(Vs),mCaLf(Vs),hCaLf(Vs),mCaLs(Vs),hCaLs(Vs),mKir(Vs),mKM(Vs),mKCa(Cas),mCAN(Cas),Cas]

    dx_dend = membrane(x_dend,[p_fixed,p_dend]) 
    dx_soma = membrane(x_soma,[p_fixed,p_soma]) 

    Ids = (10^(-3) * (Vd -Vs)/Rds) #[µA]
    Ids_dir = Ids#*(Ids>=0)
    
    F[1]= dx_dend[1] -Ids_dir/(Cd*Sd) # Gradient of Vd
    F[2]= dx_dend[end] # Gradient of Cad
    F[3]= dx_soma[1] +Ids_dir/(Cs*Ss) + I/(Cs*Ss) # Gradient of Vd
    F[4]= dx_soma[end] # Gradient of Cas
end
function nlsolve_comp_grad_0_V_!(F,x,p)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,Faraday,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra = p[1]
    Inull(t) = 0
    p_fixed = [Inull,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,Faraday,d,Ca_i_0,tau_Ca]
        ## Longitudinal resistance between compartments
        Rd = 10^(-2) * Ra * 1/2 * (Ld)/(pi*(Dd/2)^2) #[Mohm]
        Rs = 10^(-2) * Ra * 1/2 * (Ls)/(pi*(Ds/2)^2) #[Mohm]
        Rds = Rd +Rs #[Mohm]
        ## Surface areas
        Ss = 2*pi*(Ds/2)*Ls *10^(-8) # [cm²]
        Sd = 2*pi*(Dd/2)*Ld *10^(-8) # [cm²]

    ## Varying parameters -- dendrite
	p_dend= p[2]

    ## Varying parameters -- soma
	p_soma = p[3]

    # todo : verify this 
    Cd = p_fixed[2] #same capacitance per surface area for soma & dendrite
    Cs = p_fixed[2] #same capacitance per surface area for soma & dendrite


    Vd = x[1]
    Cad = x[2]
    Vs = x[3]
    Cas = x[4]

    x_dend = [Vd,mNa(Vd),hNa(Vd),mKDR(Vd),mCaLf(Vd),hCaLf(Vd),mCaLs(Vd),hCaLs(Vd),mKir(Vd),mKM(Vd),mKCa(Cad),mCAN(Cad),Cad]
    x_soma = [Vs,mNa(Vs),hNa(Vs),mKDR(Vs),mCaLf(Vs),hCaLf(Vs),mCaLs(Vs),hCaLs(Vs),mKir(Vs),mKM(Vs),mKCa(Cas),mCAN(Cas),Cas]

    dx_dend = membrane(x_dend,[p_fixed,p_dend]) 
    dx_soma = membrane(x_soma,[p_fixed,p_soma]) 

    Ids = (10^(-3) * (Vd -Vs)/Rds) #[µA]
    Ids_dir = Ids#*(Ids>=0)

	F[1]= dx_dend[1] -Ids_dir/(Cd*Sd) # Gradient of Vd
    F[2]= dx_dend[end] # Gradient of Cad
    F[3]= dx_soma[1] +Ids_dir/(Cs*Ss) + I/(Cs*Ss) # Gradient of Vd
    F[4]= dx_soma[end] # Gradient of Cas
end

function dVdot_dV(x,Ca,p)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[1]

    ## Varying parameters 
	gNa,gKDR,gKir,gKM,gKCa,gCaAN,pCaLf,pCaLs,gleak = p[2]

    dINa_dV_ = dINa_dV(x,x,gNa)
    dIKDR_dV_ = dIKDR_dV(x,x,gKDR)
    dICaLf_dV_ = dICaLf_dV(x,x,Ca,Ca_o,pCaLf)
    dICaLs_dV_ = dICaLs_dV(x,x,Ca,Ca_o,pCaLs) 
    dIKir_dV_ = dIKir_dV(x,x,gKir)    
    dIKM_dV_ = dIKM_dV(x,x,gKM)
    dIKCa_dV_ =  dIKCa_dV(x,Ca,gKCa)
    dICaAN_dV_ = dICaAN_dV(x,Ca,gCaAN)
    dIleak_dV_ = dIleak_dV(x,gleak) 

    dIcell_dV = dINa_dV_ + dIKDR_dV_ + dICaLf_dV_ + dICaLs_dV_ + dIKir_dV_ + dIKM_dV_ + dIKCa_dV_ + dICaAN_dV_ + dIleak_dV_ # [µA/cm²]
    dVdot_dV_ = -(dIcell_dV)/C
    return dVdot_dV_
end
function dVdot_dCa(x,V,p)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[1]

    ## Varying parameters 
	gNa,gKDR,gKir,gKM,gKCa,gCaAN,pCaLf,pCaLs,gleak = p[2]

    dICaLf_dCa_ = dICaLf_dCa(x,V,Ca_o,pCaLf)
    dICaLs_dCa_ = dICaLs_dCa(x,V,Ca_o,pCaLs)

    dIcell_dCa = dICaLf_dCa_ + dICaLs_dCa_ 
    dVdot_dCa_ = -(dIcell_dCa)/C
    return dVdot_dCa_
end
function dCadot_dCa(x,V,p)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[1]

    ## Varying parameters 
	gNa,gKDR,gKir,gKM,gKCa,gCaAN,pCaLf,pCaLs,gleak = p[2]

    dICaLf_dCa_ = dICaLf_dCa(x,V,Ca_o,pCaLf) 
    dICaLs_dCa_ = dICaLs_dCa(x,V,Ca_o,pCaLs)

    dICa_dCa = (dICaLf_dCa_ + dICaLs_dCa_)/1000 # [mA/cm^2] instead of [µA/cm^2]
    dCadot_dCa_ = -(dICa_dCa) *( k /(2*F*d))  -1/tau_Ca
    return dCadot_dCa_
end
function dCadot_dV(x,Ca,p)
    ## Fixed parameters
    I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[1]

    ## Varying parameters 
	gNa,gKDR,gKir,gKM,gKCa,gCaAN,pCaLf,pCaLs,gleak = p[2]

    dICaLf_dV_ = dICaLf_dV(x,x,Ca,Ca_o,pCaLf) 
    dICaLs_dV_ = dICaLs_dV(x,x,Ca,Ca_o,pCaLs) 

    dICa_dV = (dICaLf_dV_ + dICaLs_dV_ )/1000 # [mA/cm^2] instead of [µA/cm^2]
    dCadot_dV_ = -(dICa_dV) *( k /(2*F*d)  )
    return dCadot_dV_
end

dIds_dVd(x,Vs,Rds)  = ForwardDiff.derivative(x -> (10^(-3) * (x -Vs)/Rds), x)
dIds_dVs(x,Vd,Rds)  = ForwardDiff.derivative(x -> (10^(-3) * (Vd -x)/Rds), x)

dVddot_dVd(Vd,Cad,Vs,Rds,Cd,Sd,p) = dVdot_dV(Vd,Cad,p) - dIds_dVd(Vd,Vs,Rds)/(Cd*Sd)
dVsdot_dVs(Vs,Cas,Vd,Rds,Cs,Ss,p) = dVdot_dV(Vs,Cas,p) + dIds_dVs(Vs,Vd,Rds)/(Cs*Ss) 
dVddot_dVs(Vd,Vs,Rds,Cd,Sd) = - dIds_dVs(Vs,Vd,Rds)/(Cd*Sd)
dVsdot_dVd(Vs,Vd,Rds,Cs,Ss) = dIds_dVd(Vd,Vs,Rds)/(Cs*Ss) 
dVddot_dCad(Cad,Vd,p) = dVdot_dCa(Cad,Vd,p)
dVsdot_dCas(Cas,Vs,p) = dVdot_dCa(Cas,Vs,p)
dCaddot_dVd(Vd,Cad,p) = dCadot_dV(Vd,Cad,p)     
dCasdot_dVs(Vs,Cas,p) = dCadot_dV(Vs,Cas,p)    
dCaddot_dCad(Cad,Vd,p) = dCadot_dCa(Cad,Vd,p)
dCasdot_dCas(Cas,Vs,p) = dCadot_dCa(Cas,Vs,p)


epsilon=1e-6

dINa_dV(x,V,gNa)  = ForwardDiff.derivative(x -> INa(gNa,mNa(V),hNa(V),x,eNa), x)
dIKDR_dV(x,V,gKDR)  = ForwardDiff.derivative(x -> IKDR(gKDR,mKDR(V),x,eK), x)
dICaLf_dV(x,V,Ca,Ca_o,pCaLf) = ForwardDiff.derivative(x -> ICaLf(pCaLf,mCaLf(V),hCaLf(V),x,Ca,Ca_o), x)
dICaLs_dV(x,V,Ca,Ca_o,pCaLs)  = ForwardDiff.derivative(x -> ICaLs(pCaLs,mCaLs(V),hCaLs(V),x,Ca,Ca_o), x)
dIKir_dV(x,V,gKir)  = ForwardDiff.derivative(x -> IKir(gKir,mKir(V),x,eK), x)    
dIKM_dV(x,V,gKM)  = ForwardDiff.derivative(x -> IKM(gKM,mKM(V),x,eK), x)  
dIKCa_dV(x,Ca,gKCa) = ForwardDiff.derivative(x -> IKCa(gKCa,mKCa(Ca),x,eK), x)
dICaAN_dV(x,Ca,gCaAN)  = ForwardDiff.derivative(x -> ICaAN(gCaAN,mCAN(Ca),x,eCaAN), x)
dIleak_dV(x,gleak)  = ForwardDiff.derivative(x -> Ileak(gleak,x,eleak) , x)

dICaAN_dCa(x,V,gCaAN)  = 0
dICaLf_dCa(x,V,Ca_o,pCaLf) = ForwardDiff.derivative(x -> ICaLf(pCaLf,mCaLf(V),hCaLf(V),V,x,Ca_o), x)
dICaLs_dCa(x,V,Ca_o,pCaLs)  = ForwardDiff.derivative(x -> ICaLs(pCaLs,mCaLs(V),hCaLs(V),V,x,Ca_o), x)

dINa_dmNa(x,V,gNa) = ForwardDiff.derivative(x -> INa(gNa,x,hNa(V),V,eNa), x)
dINa_dhNa(x,V,gNa) = ForwardDiff.derivative(x -> INa(gNa,mNa(V),x,V,eNa), x)
dIKDR_dmKDR(x,V,gKDR) = ForwardDiff.derivative(x -> IKDR(gKDR,x,V,eK), x)
dICaLf_dmCaLf(x,V,Ca,Ca_o,pCaLf) = ForwardDiff.derivative(x -> ICaLf(pCaLf,x,hCaLf(V),V,Ca,Ca_o), x)
dICaLf_dhCaLf(x,V,Ca,Ca_o,pCaLf) = ForwardDiff.derivative(x -> ICaLf(pCaLf,mCaLf(V),x,V,Ca,Ca_o), x)
dICaLs_dmCaLs(x,V,Ca,Ca_o,pCaLs) = ForwardDiff.derivative(x -> ICaLs(pCaLs,x,hCaLs(V),V,Ca,Ca_o), x)
dICaLs_dhCaLs(x,V,Ca,Ca_o,pCaLs) = ForwardDiff.derivative(x -> ICaLs(pCaLs,mCaLs(V),x,V,Ca,Ca_o), x)
dIKir_dmKir(x,V,gKir) = ForwardDiff.derivative(x -> IKir(gKir,x,V,eK), x)              
dIKM_dmKM(x,V,gKM) = ForwardDiff.derivative(x -> IKM(gKM,x,V,eK), x)
dIKCa_dmKCa(x,V,gKCa) = ForwardDiff.derivative(x -> IKCa(gKCa,x,V,eK), x)
dICaAN_dmCaAN(x,V,gCaAN) = ForwardDiff.derivative(x -> ICaAN(gCaAN,x,V,eCaAN), x)

dmNadot_dV(x) = ForwardDiff.derivative(x -> (mNa(x)), x)/tau_mNa(x)
dhNadot_dV(x) = ForwardDiff.derivative(x -> (hNa(x)), x)/tau_hNa(x) 
dmKDRdot_dV(x) = ForwardDiff.derivative(x -> (mKDR(x)), x)/tau_mKDR(x)
dmCaLfdot_dV(x) = ForwardDiff.derivative(x -> (mCaLf(x)), x)/tau_mCaLf(x)
dhCaLfdot_dV(x) = ForwardDiff.derivative(x -> (hCaLf(x)), x)/tau_hCaLf(x)
dmCaLsdot_dV(x) = ForwardDiff.derivative(x -> (mCaLs(x)), x)/tau_mCaLs(x)
dhCaLsdot_dV(x) = ForwardDiff.derivative(x -> (hCaLs(x)), x)/tau_hCaLs(x)   
dmKirdot_dV(x) = ForwardDiff.derivative(x -> (mKir(x)), x)/tau_mKir(x)
dmKMdot_dV(x) = ForwardDiff.derivative(x -> (mKM(x)), x)/tau_mKM(x)
dmKCadot_dCa(x) = ForwardDiff.derivative(x -> (mKCa(x)), x)/tau_mKCa(x)
dmCaANdot_dCa(x) = ForwardDiff.derivative(x -> (mCAN(x)), x)/tau_mCAN(x)


function J_matrix(Vd,Cad,Vs,Cas,p)
	## Fixed parameters
	I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra = p[1]
    Inull = 0
    p_fixed = [Inull,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]
        ## Longitudinal resistance between compartments
        Rd = 10^(-2) * Ra * 1/2 * (Ld)/(pi*(Dd/2)^2) #[Mohm]
        Rs = 10^(-2) * Ra * 1/2 * (Ls)/(pi*(Ds/2)^2) #[Mohm]
        Rds = Rd +Rs #[Mohm]
        ## Surface areas
        Ss = 2*pi*(Ds/2)*Ls *10^(-8) # [cm²]
        Sd = 2*pi*(Dd/2)*Ld *10^(-8) # [cm²]

    ## Varying parameters -- dendrite
	p_dend= p[2]
    gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d,gleak_d = p_dend

    ## Varying parameters -- soma
	p_soma = p[3]
    gNa_s,gKDR_s,gKir_s,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s = p_soma


    # todo : verify this 
    Cd = p_fixed[2] #same capacitance per surface area for soma & dendrite
    Cs = p_fixed[2] #same capacitance per surface area for soma & dendrite

	J_mat = zeros(26,26)
	J_mat[1,1] = dVddot_dVd(Vd,Cad,Vs,Rds,Cd,Sd,[p_fixed,p_dend])                             
	J_mat[1,2] = -1/Cd * dINa_dmNa(mNa(Vd),Vd,gNa_d)
	J_mat[1,3] = -1/Cd * dINa_dhNa(hNa(Vd),Vd,gNa_d)
	J_mat[1,4] = -1/Cd * dIKDR_dmKDR(mKDR(Vd),Vd,gKDR_d)
    J_mat[1,5] = -1/Cd * dICaLf_dmCaLf(mCaLf(Vd),Vd,Cad,Ca_o,pCaLf_d)
	J_mat[1,6] = -1/Cd * dICaLf_dhCaLf(hCaLf(Vd),Vd,Cad,Ca_o,pCaLf_d)
    J_mat[1,7] = -1/Cd * dICaLs_dmCaLs(mCaLs(Vd),Vd,Cad,Ca_o,pCaLs_d)
	J_mat[1,8] = -1/Cd * dICaLs_dhCaLs(hCaLs(Vd),Vd,Cad,Ca_o,pCaLs_d)
    J_mat[1,9] = -1/Cd * dIKir_dmKir(mKir(Vd),Vd,gKir_d)
	J_mat[1,10] = -1/Cd * dIKM_dmKM(mKM(Vd),Vd,gKM_d)
	J_mat[1,11] = -1/Cd * dIKCa_dmKCa(mKCa(Cad),Vd,gKCa_d)
	J_mat[1,12] = -1/Cd * dICaAN_dmCaAN(mCAN(Cad),Vd,gCaAN_d)
	J_mat[1,13] = dVddot_dCad(Cad,Vd,[p_fixed,p_dend])                          
	J_mat[1,14] = dVddot_dVs(Vd,Vs,Rds,Cd,Sd)
	J_mat[1,15:25] .= 0
    J_mat[1,26] = 0 #dVddot_dCas(Ca,V,p)   
    
    J_mat[14,1] = dVsdot_dVd(Vs,Vd,Rds,Cs,Ss)                             
	J_mat[14,2:12] .= 0
	J_mat[14,13] = 0 #dVsdot_dCad(Ca,V,p)                          
	J_mat[14,14] = dVsdot_dVs(Vs,Cas,Vd,Rds,Cs,Ss,[p_fixed,p_soma])
	J_mat[14,15] = -1/Cs * dINa_dmNa(mNa(Vs),Vs,gNa_s)
    J_mat[14,16] = -1/Cs * dINa_dhNa(hNa(Vs),Vs,gNa_s)
	J_mat[14,17] = -1/Cs * dIKDR_dmKDR(mKDR(Vs),Vs,gKDR_s)
    J_mat[14,18] = -1/Cs * dICaLf_dmCaLf(mCaLf(Vs),Vs,Cas,Ca_o,pCaLf_s)
	J_mat[14,19] = -1/Cs * dICaLf_dhCaLf(hCaLf(Vs),Vs,Cas,Ca_o,pCaLf_s)
    J_mat[14,20] = -1/Cs * dICaLs_dmCaLs(mCaLs(Vs),Vs,Cas,Ca_o,pCaLs_s)
    J_mat[14,21] = -1/Cs * dICaLs_dhCaLs(hCaLs(Vs),Vs,Cas,Ca_o,pCaLs_s)
	J_mat[14,22] = -1/Cs * dIKir_dmKir(mKir(Vs),Vs,gKir_s)
    J_mat[14,23] = -1/Cs * dIKM_dmKM(mKM(Vs),Vs,gKM_s)
	J_mat[14,24] = -1/Cs * dIKCa_dmKCa(mKCa(Cas),Vs,gKCa_s)
	J_mat[14,25] = -1/Cs * dICaAN_dmCaAN(mCAN(Cas),Vs,gCaAN_s)
    J_mat[14,26] = dVsdot_dCas(Cas,Vs,[p_fixed,p_soma])   

	J_mat[2,1] = dmNadot_dV(Vd)
	J_mat[3,1] = dhNadot_dV(Vd)
	J_mat[4,1] = dmKDRdot_dV(Vd)
    J_mat[5,1] = dmCaLfdot_dV(Vd)
	J_mat[6,1] = dhCaLfdot_dV(Vd)
    J_mat[7,1] = dmCaLsdot_dV(Vd)
	J_mat[8,1] = dhCaLsdot_dV(Vd)
    J_mat[9,1] = dmKirdot_dV(Vd)
    J_mat[10,1] = dmKMdot_dV(Vd)
	J_mat[11,1] = 0                             #!! Ca-dep
	J_mat[12,1] = 0                             #!! Ca-dep
	J_mat[13,1] = dCaddot_dVd(Vd,Cad,[p_fixed,p_dend])                    


	J_mat[15,14] = dmNadot_dV(Vs)
    J_mat[16,14] = dhNadot_dV(Vs)
	J_mat[17,14] = dmKDRdot_dV(Vs)
    J_mat[18,14] = dmCaLfdot_dV(Vs)
	J_mat[19,14] = dhCaLfdot_dV(Vs)
    J_mat[20,14] = dmCaLsdot_dV(Vs)
    J_mat[21,14] = dhCaLsdot_dV(Vs)
    J_mat[22,14] = dmKirdot_dV(Vs)
	J_mat[23,14] = dmKMdot_dV(Vs)
    J_mat[24,14] = 0                             #!! Ca-dep
    J_mat[25,14] = 0                             #!! Ca-dep
	J_mat[26,14] = dCasdot_dVs(Vs,Cas,[p_fixed,p_soma])                 

    J_mat[11,13] = dmKCadot_dCa(Cad)      #!! Ca-dep
	J_mat[12,13] = dmCaANdot_dCa(Cad)     #!! Ca-dep

    J_mat[24,26] = dmKCadot_dCa(Cas)      #!! Ca-dep
	J_mat[25,26] = dmCaANdot_dCa(Cas)     #!! Ca-dep

	J_mat[2,2] = -1/tau_mNa(Vd)
	J_mat[3,3] = -1/tau_hNa(Vd)
	J_mat[4,4] = -1/tau_mKDR(Vd)
    J_mat[5,5] = -1/tau_mCaLf(Vd)
	J_mat[6,6] = -1/tau_hCaLf(Vd)
    J_mat[7,7] = -1/tau_mCaLs(Vd)
	J_mat[8,8] = -1/tau_hCaLs(Vd)
    J_mat[9,9] = -1/tau_mKir(Vd)
    J_mat[10,10] = -1/tau_mKM(Vd)
	J_mat[11,11] = -1/tau_mKCa(Cad)             #!! Ca-dep
    J_mat[12,12] = -1/tau_mCAN(Cad)             #!! Ca-dep
    J_mat[13,13] = dCaddot_dCad(Cad,Vd,[p_fixed,p_dend])             
	J_mat[15,15] = -1/tau_mNa(Vs)
	J_mat[16,16] = -1/tau_hNa(Vs)
	J_mat[17,17] = -1/tau_mKDR(Vs)
    J_mat[18,18] = -1/tau_mCaLf(Vs)
	J_mat[19,19] = -1/tau_hCaLf(Vs)
    J_mat[20,20] = -1/tau_mCaLs(Vs)
	J_mat[21,21] = -1/tau_hCaLs(Vs)
    J_mat[22,22] = -1/tau_mKir(Vs)
    J_mat[23,23] = -1/tau_mKM(Vs)
	J_mat[24,24] = -1/tau_mKCa(Cas)             #!! Ca-dep
    J_mat[25,25] = -1/tau_mCAN(Cas)             #!! Ca-dep
    J_mat[26,26] = dCasdot_dCas(Cas,Vs,[p_fixed,p_soma])             

	return J_mat
end

function plt_fixed_points(plt,fixed_points)
    Vd_saddle,Cad_saddle,Vs_saddle,Cas_saddle,I_saddle,Vd_stable,Cad_stable,Vs_stable,Cas_stable,I_stable,Vd_unstable,Cad_unstable,Vs_unstable,Cas_unstable,I_unstable = fixed_points
    
    println(length(I_stable))
    println(length(Vs_stable))

    I_stable = I_stable *10^6 #[pA]
    I_saddle = I_saddle *10^6 #[pA]
    I_unstable = I_unstable *10^6 #[pA]

    println(length(I_stable))
    println(length(Vs_stable))

    plt_Vs = Plots.plot(plt[1],fontfamily="Computer Modern", xlabel=" ",ylabel="Vs",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(I_stable ,Vs_stable ,lc=:green,lw=2)
	plot!(I_saddle ,Vs_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(I_unstable,Vs_unstable , lc=:red, ls=:dash,lw=2)
    #scatter!(I_SN .*ones(2), V_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)
    display(plt_Vs)

    plt_Cas = Plots.plot(plt[2],fontfamily="Computer Modern", xlabel="I",ylabel="Cas",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(I_stable ,Cas_stable ,lc=:green,lw=2)
	plot!(I_saddle ,Cas_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(I_unstable,Cas_unstable , lc=:red, ls=:dash,lw=2)
    #scatter!(I_SN .*ones(2), Ca_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)

    plt_Vd = Plots.plot(plt[1],fontfamily="Computer Modern", xlabel=" ",ylabel="Vd",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(I_stable ,Vd_stable ,lc=:green,lw=2)
	plot!(I_saddle ,Vd_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(I_unstable,Vd_unstable , lc=:red, ls=:dash,lw=2)
    #scatter!(I_SN .*ones(2), V_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)

    plt_Cad = Plots.plot(plt[2],fontfamily="Computer Modern", xlabel="I",ylabel="Cad",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(I_stable ,Cad_stable ,lc=:green,lw=2)
	plot!(I_saddle ,Cad_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(I_unstable,Cad_unstable , lc=:red, ls=:dash,lw=2)
    #scatter!(I_SN .*ones(2), Ca_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)
	
    plt_ = Plots.plot(plt_Vs,plt_Cas,plt_Vd,plt_Cad,layout=(4,1),size=(600,1000))
    return plt_
end
function plt_fixed_points_tr(plt,fixed_points)
    Vd_saddle,Cad_saddle,Vs_saddle,Cas_saddle,I_saddle,Vd_stable,Cad_stable,Vs_stable,Cas_stable,I_stable,Vd_unstable,Cad_unstable,Vs_unstable,Cas_unstable,I_unstable = fixed_points

    plt_Vs = Plots.plot(plt[1],fontfamily="Computer Modern", xlabel=" ",ylabel="Vs",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(Vs_stable,I_stable,lc=:green,lw=2)
	plot!(Vs_saddle ,I_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(Vs_unstable ,I_unstable, lc=:red, ls=:dash,lw=2)
    #scatter!( V_SN .*ones(2),I_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)

    plt_Cas = Plots.plot(plt[2],fontfamily="Computer Modern", xlabel="I",ylabel="Cas",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,ylims=(-0.005,(Ca_SN*10)),size=(800,400))
    plot!(Cas_stable ,I_stable ,lc=:green,lw=2)
	plot!(Cas_saddle ,I_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(Cas_unstable ,I_unstable, lc=:red, ls=:dash,lw=2)
    #scatter!( Ca_SN .*ones(2),I_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)

    plt_Vd = Plots.plot(plt[4],fontfamily="Computer Modern", xlabel=" ",ylabel="Vd",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,size=(800,400))
    plot!(Vd_stable,I_stable,lc=:green,lw=2)
	plot!(Vd_saddle ,I_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(Vd_unstable ,I_unstable, lc=:red, ls=:dash,lw=2)
    #scatter!( V_SN .*ones(2),I_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)

    plt_Cad = Plots.plot(plt[4],fontfamily="Computer Modern", xlabel="I",ylabel="Cad",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,ylims=(-0.005,(Ca_SN*10)),size=(800,400))
    plot!(Cad_stable ,I_stable ,lc=:green,lw=2)
	plot!(Cad_saddle ,I_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(Cad_unstable ,I_unstable, lc=:red, ls=:dash,lw=2)
    #scatter!( Ca_SN .*ones(2),I_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)
	
    plt_ = Plots.plot(plt_Vs,plt_Cas,plt_Vd,plt_Cad,layout=(4,1),size=(600,1000))
    return plt_
end
function find_fixed_points(p,Vs_FP_range,solver_ic)
	## Computation of the SN from maximum I_FP observed when setting all gradients to 0 
	## in the conductance-based model for Vs_FP_range 
	## Outputs may contain NaN for a better and easy result display in a plot

    p_fixed,p_dend,p_soma = p


	I_FP_range = []
    Vd_FP_range = []
    Cad_FP_range = []
    Cas_FP_range = []
	Vd_saddle = []
	Cad_saddle = []
    Vs_saddle = []
	Cas_saddle = []
    I_saddle = []
	Vd_stable = []
	Cad_stable = []
    Vs_stable = []
	Cas_stable = []
    I_stable = []
	Vd_unstable = []
	Cad_unstable = []
    Vs_unstable = []
	Cas_unstable = []
    I_unstable = []

	for vs in Vs_FP_range
        if length(I_FP_range)>0
            solve_comp_grad_0_I_ = nlsolve((F,x) ->nlsolve_comp_grad_0_I_!(F,x,[p_fixed,p_dend,p_soma,vs]), [I_FP_range[end],Vd_FP_range[end],Cad_FP_range[end],Cas_FP_range[end]])
        else
            solve_comp_grad_0_I_ = nlsolve((F,x) ->nlsolve_comp_grad_0_I_!(F,x,[p_fixed,p_dend,p_soma,vs]), solver_ic)
		end
        push!(I_FP_range,solve_comp_grad_0_I_.zero[1])
        push!(Vd_FP_range,solve_comp_grad_0_I_.zero[2]) 
        push!(Cad_FP_range,solve_comp_grad_0_I_.zero[3]) 
        push!(Cas_FP_range,solve_comp_grad_0_I_.zero[4]) 
        #println("Vs = $(vs) ; zeros = $(solve_comp_grad_0_I_.zero)") 
		lambda = real.(eigvals(J_matrix(Vd_FP_range[end],Cad_FP_range[end],vs,Cas_FP_range[end],[p_fixed,p_dend,p_soma])) )

        if maximum(lambda) <=0
			push!(Vs_stable,vs)
			push!(I_stable,solve_comp_grad_0_I_.zero[1])
            push!(Vd_stable,solve_comp_grad_0_I_.zero[2])
            push!(Cad_stable,solve_comp_grad_0_I_.zero[3])
            push!(Cas_stable,solve_comp_grad_0_I_.zero[4])
			push!(Vs_saddle,vs)
			push!(I_saddle,NaN)
            push!(Vd_saddle,NaN)
            push!(Cad_saddle,NaN)
            push!(Cas_saddle,NaN)
			push!(Vs_unstable,vs)
			push!(I_unstable,NaN)
            push!(Vd_unstable,NaN)
            push!(Cad_unstable,NaN)
            push!(Cas_unstable,NaN)
		else
			if maximum(lambda[lambda.!=maximum(lambda)]) <=0
				push!(Vs_saddle,vs)
				push!(I_saddle,solve_comp_grad_0_I_.zero[1])
                push!(Vd_saddle,solve_comp_grad_0_I_.zero[2])
                push!(Cad_saddle,solve_comp_grad_0_I_.zero[3])
                push!(Cas_saddle,solve_comp_grad_0_I_.zero[4])
				push!(Vs_stable,vs)
				push!(I_stable,NaN)
                push!(Vd_stable,NaN)
                push!(Cad_stable,NaN)
                push!(Cas_stable,NaN)
				push!(Vs_unstable,vs)
				push!(I_unstable,NaN)
                push!(Vd_unstable,NaN)
                push!(Cad_unstable,NaN)
                push!(Cas_unstable,NaN)
			else
				push!(Vs_unstable,vs)
				push!(I_unstable,solve_comp_grad_0_I_.zero[1])
                push!(Vd_unstable,solve_comp_grad_0_I_.zero[2])
                push!(Cad_unstable,solve_comp_grad_0_I_.zero[3])
                push!(Cas_unstable,solve_comp_grad_0_I_.zero[4])
				push!(Vs_saddle,vs)
				push!(I_saddle,NaN)
                push!(Vd_saddle,NaN)
                push!(Cad_saddle,NaN)
                push!(Cas_saddle,NaN)
				push!(Vs_stable,vs)
				push!(I_stable,NaN)
                push!(Vd_stable,NaN)
                push!(Cad_stable,NaN)
                push!(Cas_stable,NaN)
			end
		end
	end

    #==
	I_SN = NaN#maximum(I_FP_range[Vs_FP_range.<=-40])
	i_I_SN = NaN#findfirst(x -> x == I_SN, I_FP_range)
	V_SN = NaN#Vs_FP_range[i_I_SN]
    Cad_SN = NaN#Cad_FP_range[i_I_SN]==#

    return Vd_saddle,Cad_saddle,Vs_saddle,Cas_saddle,I_saddle,Vd_stable,Cad_stable,Vs_stable,Cas_stable,I_stable,Vd_unstable,Cad_unstable,Vs_unstable,Cas_unstable,I_unstable
end
function find_fixed_points_I(p,I_FP_range,solver_ic)

    p_fixed,p_dend,p_soma = p

    Vd_FP_range = []
    Cad_FP_range = []
    Vs_FP_range = []
    Cas_FP_range = []

    Vd_saddle = []
	Cad_saddle = []
    Vs_saddle = []
	Cas_saddle = []
    I_saddle = []
	Vd_stable = []
	Cad_stable = []
    Vs_stable = []
	Cas_stable = []
    I_stable = []
	Vd_unstable = []
	Cad_unstable = []
    Vs_unstable = []
	Cas_unstable = []
    I_unstable = []

    for i in I_FP_range
        p_fixed_i = vcat(i,p_fixed[2:end])
        if length(Vs_FP_range)>0
            solve_comp_grad_0_V_ = nlsolve((F,x) ->nlsolve_comp_grad_0_V_!(F,x,[p_fixed_i,p_dend,p_soma]), [Vd_FP_range[end],Cad_FP_range[end],Vs_FP_range[end],Cas_FP_range[end]])
        else
            solve_comp_grad_0_V_ = nlsolve((F,x) ->nlsolve_comp_grad_0_V_!(F,x,[p_fixed_i,p_dend,p_soma]), solver_ic)
		end
        push!(Vd_FP_range,solve_comp_grad_0_V_.zero[1]) 
        push!(Cad_FP_range,solve_comp_grad_0_V_.zero[2]) 
        push!(Vs_FP_range,solve_comp_grad_0_V_.zero[3])
        push!(Cas_FP_range,solve_comp_grad_0_V_.zero[4]) 
        if converged(solve_comp_grad_0_V_)
            lambda = real.(eigvals(J_matrix(Vd_FP_range[end],Cad_FP_range[end],Vs_FP_range[end],Cas_FP_range[end],[p_fixed_i,p_dend,p_soma])) )
            if maximum(lambda) <=0
                push!(I_stable,i)
                push!(Vd_stable,solve_comp_grad_0_V_.zero[1])
                push!(Cad_stable,solve_comp_grad_0_V_.zero[2])
                push!(Vs_stable,solve_comp_grad_0_V_.zero[3])
                push!(Cas_stable,solve_comp_grad_0_V_.zero[4])
                push!(I_saddle,i)
                push!(Vd_saddle,NaN)
                push!(Cad_saddle,NaN)
                push!(Vs_saddle,NaN)
                push!(Cas_saddle,NaN)
                push!(I_unstable,i)
                push!(Vd_unstable,NaN)
                push!(Cad_unstable,NaN)
                push!(Vs_unstable,NaN)
                push!(Cas_unstable,NaN)
            else
                if maximum(lambda[lambda.!=maximum(lambda)]) <=0
                    push!(I_saddle,i)
                    push!(Vd_saddle,solve_comp_grad_0_V_.zero[1])
                    push!(Cad_saddle,solve_comp_grad_0_V_.zero[2])
                    push!(Vs_saddle,solve_comp_grad_0_V_.zero[3])
                    push!(Cas_saddle,solve_comp_grad_0_V_.zero[4])
                    push!(I_stable,i)
                    push!(Vd_stable,NaN)
                    push!(Cad_stable,NaN)
                    push!(Vs_stable,NaN)
                    push!(Cas_stable,NaN)
                    push!(I_unstable,i)
                    push!(Vd_unstable,NaN)
                    push!(Cad_unstable,NaN)
                    push!(Vs_unstable,NaN)
                    push!(Cas_unstable,NaN)
                else
                    push!(I_unstable,i)
                    push!(Vd_unstable,solve_comp_grad_0_V_.zero[1])
                    push!(Cad_unstable,solve_comp_grad_0_V_.zero[2])
                    push!(Vs_unstable,solve_comp_grad_0_V_.zero[3])
                    push!(Cas_unstable,solve_comp_grad_0_V_.zero[4])
                    push!(I_saddle,i)
                    push!(Vd_saddle,NaN)
                    push!(Cad_saddle,NaN)
                    push!(Vs_saddle,NaN)
                    push!(Cas_saddle,NaN)
                    push!(I_stable,i)
                    push!(Vd_stable,NaN)
                    push!(Cad_stable,NaN)
                    push!(Vs_stable,NaN)
                    push!(Cas_stable,NaN)
                end
            end
        end

    end

    return Vd_saddle,Cad_saddle,Vs_saddle,Cas_saddle,I_saddle,Vd_stable,Cad_stable,Vs_stable,Cas_stable,I_stable,Vd_unstable,Cad_unstable,Vs_unstable,Cas_unstable,I_unstable
end


#==
    V=-100:0.1:-40
    find_lfp = find_fixed_points([vcat(0,p_fixed[2:end]),p_var_d,p_var_s],V,[25/10^(6),-59. ,Ca_i_0,Ca_i_0])
    plt_fp = plt_fixed_points(plot(plot(),plot(),plot(),plot(),layout=(4,1)),find_lfp)
==#
#==
    V=-58.9:0.1:-30
    find_fp = find_fixed_points([vcat(0,p_fixed[2:end]),p_var_d,p_var_s],V,[-155/10^(6),-24.3,3.25,Ca_i_0])
    plt_fp = plt_fixed_points(plot(plot(),plot(),plot(),plot(),layout=(4,1)),find_fp)
==#

