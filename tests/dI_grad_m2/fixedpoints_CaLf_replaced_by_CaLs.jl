using Plots,LaTeXStrings
using DifferentialEquations,NLsolve,LinearAlgebra

include("gates_CaLf_replaced_by_CaLs.jl")
include("currents_m2.jl")

function derivative(fct,step,x)
	# This function computes the derivative of fct at x for a given step 
	derivative = (fct(x+step)-fct(x))/step
	return derivative
end
function membrane_grad_0(V,Ca,p)
	# This function computes the gradients of V and Ca, with the gates at steady-state
	## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]
    
    ## Gates steady-state 
    mNa_inf=mNa(V)
    hNa_inf=hNa(V)    
    mKDR_inf=mKDR(V)
	mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
    mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
    mKir_inf=mKir(V)
    mKM_inf=mKM(V)

    ## Current computation
    INa_ = INa(gNa,mNa_inf,hNa_inf,V,eNa)
    IKDR_ = IKDR(gKDR,mKDR_inf,V,eK)
    ICaLf_ = ICaLf(pCaLf,mCaLf_inf,hCaLf_inf,V,Ca,Ca_o)
    ICaLs_ = ICaLs(pCaLs,mCaLs_inf,hCaLs_inf,V,Ca,Ca_o)
    IKir_ = IKir(gKir,mKir_inf,V,eK)
    IKM_ = IKM(gKM,mKM_inf,V,eK)
    Ileak_ = Ileak(gleak,V,eleak)
    
    V_dot = (-INa_ -IKDR_ -ICaLf_ -ICaLs_- IKir_- IKM_ -Ileak_ +I)/C

    ICa_ = (ICaLf_ + ICaLs_)/1000  # [mA/cm^2] instead of [µA/cm^2]
    Ca_dot = - ICa_ * k /(2*F*d) - ((Ca-Ca_i_0)/tau_Ca)

    grad = [V_dot,Ca_dot]	
	return grad
end
function nlsolve_membrane_grad_0_I_!(F,x,p_)
    ## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p_[1]

    ## Fixed parameters
    p_fixed = p_[2]

    ## NLSolver parameter
	V = p_[3][1]

	p_var=[x[1],gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
    p=[p_var,p_fixed]
    
    membrane_grad = membrane_grad_0(V,x[2],p) 
	F[1]= membrane_grad[1] # Gradient of V
    F[2]= membrane_grad[2] # Gradient of Ca
end
function J_dINa_dmNa(V)
	mNa_inf=mNa(V)
    hNa_inf=hNa(V)    
	dVdot_dmNa = 3*mNa_inf^2*hNa_inf*(V-eNa)
	return dVdot_dmNa
end
function J_dINa_dhNa(V)
	mNa_inf=mNa(V)
    hNa_inf=hNa(V)    
	dVdot_dhNa = mNa_inf^3*(V-eNa)
	return dVdot_dhNa
end
function J_dIKDR_dmKDR(V)
	mKDR_inf=mKDR(V)    
	dVdot_dmKDR = 4*mKDR_inf^3*(V-eK)
	return dVdot_dmKDR
end
function J_dIKir_dmKir(V)
	mKir_inf=mKir(V)
	dVdot_dmKir = (V-eK)
	return dVdot_dmKir
end
function J_dICaLf_dmCaLf(V,Ca,Ca_o)
	mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
	dICaLf_dmCaLf = 2*mCaLf_inf*hCaLf_inf*ghk_LeFranc(V, Ca, Ca_o) *1000
	return dICaLf_dmCaLf
end
function J_dICaLf_dhCaLf(V,Ca,Ca_o)
	mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
	dICaLf_dhCaLf = mCaLf_inf^2*ghk_LeFranc(V, Ca, Ca_o) *1000 
	return dICaLf_dhCaLf
end
function J_dICaLs_dmCaLs(V,Ca,Ca_o)
	mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
	dICaLs_dmCaLs = hCaLs_inf*ghk_LeFranc(V, Ca, Ca_o) *1000
	return dICaLs_dmCaLs
end
function J_dICaLs_dhCaLs(V,Ca,Ca_o)
	mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
	dICaLs_dhCaLs = mCaLs_inf*ghk_LeFranc(V, Ca, Ca_o) *1000 
	return dICaLs_dhCaLs
end
function J_dIKM_dmKM(V)
	mKM_inf=mKM(V)
	dIKM_dmKM = (V-eK)
	return dIKM_dmKM
end
function compute_dghk_dV(V,Ca,Ca_o,step)
	dghk_dV = (ghk_LeFranc(V+step, Ca, Ca_o)-ghk_LeFranc(V, Ca, Ca_o))/step
	return dghk_dV
end
function compute_dVdot_dV(V,Ca,p,step)
    ## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]

    mNa_inf=mNa(V)
    hNa_inf=hNa(V)    
    dINa_dV  = gNa*mNa_inf^3*hNa_inf

    mKDR_inf=mKDR(V)    
    dIKDR_dV  = gKDR * mKDR_inf^4

    mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
    dICaLf_dV  = pCaLf*mCaLf_inf^2*hCaLf_inf*compute_dghk_dV(V,Ca,Ca_o,step) *1000 

    mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
    dICaLs_dV  = pCaLs*mCaLs_inf*hCaLs_inf*compute_dghk_dV(V,Ca,Ca_o,step) *1000 

    mKir_inf=mKir(V)
    dIKir_dV = gKir * mKir_inf

	mKM_inf=mKM(V)
    dIKM_dV = gKM * mKM_inf

    dIleak_dV = gleak

    g_static = -(dINa_dV .+ dIKDR_dV .+ dICaLf_dV .+ dICaLs_dV .+ dIKir_dV .+ dIKM_dV .+ dIleak_dV)

	return g_static
end
function J_dVdot_dV(V,Ca,p,step)
	J = compute_dVdot_dV(V,Ca,p,step)
	return J 
end
function J_dmNadot_dV(V,step)
	tau_mNa_=tau_mNa(V)
	dm_inf_dV = derivative(mNa,step,V)
	J = dm_inf_dV/tau_mNa_
	return J
end
function J_dhNadot_dV(V,step)
	tau_hNa_=tau_hNa(V)
	dh_inf_dV = derivative(hNa,step,V)
	J = dh_inf_dV/tau_hNa_
	return J
end
function J_dmKDRdot_dV(V,step)
	tau_mKDR_=tau_mKDR(V)
	dm_inf_dV =  derivative(mKDR,step,V)
	J = dm_inf_dV/tau_mKDR_
	return J
end
function J_dmKirdot_dV(V,step)
	tau_mKir_=tau_mKir(V)
	dm_inf_dV = derivative(mKir,step,V)
	J = dm_inf_dV/tau_mKir_
	return J
end
function J_dmCaLfdot_dV(V,step)
	tau_mCaLf_=tau_mCaLs(V)
	dm_inf_dV = derivative(mCaLs,step,V)
	J = dm_inf_dV/tau_mCaLf_
	return J
end
function J_dhCaLfdot_dV(V,step)
	tau_hCaLf_=tau_hCaLs(V)
	dh_inf_dV = derivative(hCaLs,step,V)
	J = dh_inf_dV/tau_hCaLf_
	return J
end
function J_dmCaLsdot_dV(V,step)
	tau_mCaLs_=tau_mCaLs(V)
	dm_inf_dV = derivative(mCaLs,step,V)
	J = dm_inf_dV/tau_mCaLs_
	return J
end
function J_dhCaLsdot_dV(V,step)
	tau_hCaLs_=tau_hCaLs(V)
	dh_inf_dV = derivative(hCaLs,step,V)
	J = dh_inf_dV/tau_hCaLs_
	return J
end
function J_dmKMdot_dV(V,step)
	tau_mKM_= tau_mKM(V)
	dm_inf_dV = derivative(mKM,step,V)
	J = dm_inf_dV/tau_mKM_
	return J
end
function compute_dghk_dCa(V,Ca,Ca_o,step)
	dghk_dV = (ghk_LeFranc(V, Ca+step, Ca_o)-ghk_LeFranc(V, Ca, Ca_o))/step
	return dghk_dV
end
function compute_dVdot_dCa(V,Ca,p,step)
    ## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]

    mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
    dICaLf_dCa  = pCaLf*mCaLf_inf^2*hCaLf_inf*compute_dghk_dCa(V,Ca,Ca_o,step) *1000 

    mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
    dICaLs_dCa  = pCaLs*mCaLs_inf*hCaLs_inf*compute_dghk_dCa(V,Ca,Ca_o,step) *1000 

    dVdot_dCa = -( dICaLf_dCa .+ dICaLs_dCa)

	return dVdot_dCa
end
function J_dVdot_dCa(V,Ca,p,step)
	J = compute_dVdot_dCa(V,Ca,p,step)
	return J 
end
function compute_dCadot_dV(V,Ca,p,step)
    ## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]

    mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
    dICaLf_dV  = pCaLf*mCaLf_inf^2*hCaLf_inf*compute_dghk_dV(V,Ca,Ca_o,step) *1000 

    mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
    dICaLs_dV  = pCaLs*mCaLs_inf*hCaLs_inf*compute_dghk_dV(V,Ca,Ca_o,step) *1000 


    dCadot_dV = - (dICaLf_dV .+ dICaLs_dV)/1000 * k /(2*F*d)
	return dCadot_dV
end
function J_dCadot_dV(V,Ca,p,step)
	J = compute_dCadot_dV(V,Ca,p,step)
	return J 
end
function compute_dCadot_dCa(V,Ca,p,step)
    ## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]

    mCaLf_inf=mCaLs(V)
    hCaLf_inf=hCaLs(V)
    dICaLf_dCa  = pCaLf*mCaLf_inf^2*hCaLf_inf*compute_dghk_dCa(V,Ca,Ca_o,step) *1000 

    mCaLs_inf=mCaLs(V)
    hCaLs_inf=hCaLs(V)
    dICaLs_dCa  = pCaLs*mCaLs_inf*hCaLs_inf*compute_dghk_dCa(V,Ca,Ca_o,step) *1000 

    dCadot_dCa = - (dICaLf_dCa .+ dICaLs_dCa)/1000 * k /(2*F*d) -1/tau_Ca
	return dCadot_dCa
end
function J_dCadot_dCa(V,Ca,p,step)
	J = compute_dCadot_dCa(V,Ca,p,step)
	return J 
end
function J_matrix(V,Ca,p)
	## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p[1]

    ## Fixed parameters
    C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca = p[2]

	step=0.0000001

	J_mat = zeros(11,11)
	J_mat[1,1] = 1/C * J_dVdot_dV(V,Ca,p,step)
	J_mat[1,2] = -gNa * J_dINa_dmNa(V)
	J_mat[1,3] = -gNa * J_dINa_dhNa(V)
	J_mat[1,4] = -gKDR * J_dIKDR_dmKDR(V)
	J_mat[1,5] = -pCaLf * J_dICaLf_dmCaLf(V,Ca,Ca_o)
	J_mat[1,6] = -pCaLf * J_dICaLf_dhCaLf(V,Ca,Ca_o)
    J_mat[1,7] = -pCaLs * J_dICaLs_dmCaLs(V,Ca,Ca_o)
	J_mat[1,8] = -pCaLs * J_dICaLs_dhCaLs(V,Ca,Ca_o)
    J_mat[1,9] = -gKir * J_dIKir_dmKir(V)
	J_mat[1,10] = -gKM * J_dIKM_dmKM(V)
    J_mat[1,11] = 1/C *J_dVdot_dCa(V,Ca,p,step)

	J_mat[2,1] = J_dmNadot_dV(V,step)
	J_mat[3,1] = J_dhNadot_dV(V,step)
	J_mat[4,1] = J_dmKDRdot_dV(V,step)
	J_mat[5,1] = J_dmCaLfdot_dV(V,step)
	J_mat[6,1] = J_dhCaLfdot_dV(V,step)
    J_mat[7,1] = J_dmCaLsdot_dV(V,step)
	J_mat[8,1] = J_dhCaLsdot_dV(V,step)
    J_mat[9,1] = J_dmKirdot_dV(V,step)
	J_mat[10,1] = J_dmKMdot_dV(V,step)
    J_mat[11,1] = J_dCadot_dV(V,Ca,p,step)

	J_mat[2,2] = -1/tau_mNa(V)
	J_mat[3,3] = -1/tau_hNa(V)
	J_mat[4,4] = -1/tau_mKDR(V)
	J_mat[5,5] = -1/tau_mCaLs(V)
	J_mat[6,6] = -1/tau_hCaLs(V)
    J_mat[7,7] = -1/tau_mCaLs(V)
	J_mat[8,8] = -1/tau_hCaLs(V)
    J_mat[9,9] = -1/tau_mKir(V)
	J_mat[10,10] = -1/tau_mKM(V)
    J_mat[11,11] = J_dCadot_dCa(V,Ca,p,step)

	return J_mat
end
function nlsolve_membrane_grad_0_V_!(F,x,p)
    V = x[1]
    Ca = x[2]

    membrane_grad = membrane_grad_0(V,Ca,p) 
	F[1]= membrane_grad[1] # Gradient of V
    F[2]= membrane_grad[2] # Gradient of Ca
end

function stability_FP(p,I_l_lf,I_h_lf,step_Illf_to_Ihlf)
    p_var,p_fixed = p

	## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p_var

    ## Fixed parameters
    Ca_i_0 = p_fixed[10]

    I_FP_range = []
    Ca_FP_range = []
    V_FP_range = range(-150.,stop=-40.,step=0.01)

    for v in V_FP_range
        solve_membrane_grad_0_I_ = nlsolve((F,x) ->nlsolve_membrane_grad_0_I_!(F,x,[p_var,p_fixed,v]), [-1.,Ca_i_0],ftol=1E-7,xtol=1E-7)
        push!(I_FP_range,solve_membrane_grad_0_I_.zero[1])
        push!(Ca_FP_range,solve_membrane_grad_0_I_.zero[2])  
    end 
    
    V_FP_stable = []
    V_FP_saddle = []
    V_FP_unstable = []
    I_saddle_ = []
    I_stable_ = []
    I_unstable_ = []
    Ca_FP_saddle = []
	Ca_FP_stable = []
	Ca_FP_unstable = []


    for i in range(I_l_lf,stop=I_h_lf,step=step_Illf_to_Ihlf)
        #println(range(I_l_lf,stop=I_h_lf,step=step_Illf_to_Ihlf))

        p_var_i = [i,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
        V_FP_est = V_FP_range[findfirst(x -> x>= i,I_FP_range)]
        Ca_FP_est = Ca_FP_range[findfirst(x -> x>= i,I_FP_range)]
		FP_l = nlsolve((F,x) ->nlsolve_membrane_grad_0_V_!(F,x,[p_var_i,p_fixed]), [V_FP_est,Ca_FP_est],ftol=1E-7,xtol=1E-7)
        res_h = nlsolve_membrane_grad_0_V_!([NaN,NaN],FP_l.zero,[p_var_i,p_fixed])
        #println(res_h)
        if res_h < 1e-6
            lambda = real.(eigvals(J_matrix(FP_l.zero[1],FP_l.zero[2],[p_var_i,p_fixed])) )
            if maximum(lambda) <=0
                push!(V_FP_stable,FP_l.zero[1])
                push!(Ca_FP_stable,FP_l.zero[2])
                push!(I_stable_,i)
            else
                if maximum(lambda[lambda.!=maximum(lambda)]) <=0
                    push!(V_FP_saddle,FP_l.zero[1])
                    push!(Ca_FP_saddle,FP_l.zero[2])
                    push!(I_saddle_,i)
                else
                    push!(V_FP_unstable,FP_l.zero[1])
                    push!(Ca_FP_unstable,FP_l.zero[2])
                    push!(I_unstable_,i)
                end
            end
        end

    end

    return V_FP_saddle,Ca_FP_saddle,I_saddle_,V_FP_stable,Ca_FP_stable,I_stable_,V_FP_unstable,Ca_FP_unstable,I_unstable_
end
function find_fixed_points(p,V_FP_range,IX,n_IX_to_ISN)
	## Computation of the SN from maximum I_FP observed when setting all gradients to 0 
	## in the conductance-based model for V_FP_range and computation of stability in the 
	## specific range of current "range(I_l_lf,stop=I_h_lf-1e-3,length=n_IX_to_ISN)"
	## Outputs may contain NaN for a better and easy result display in a plot

    p_var,p_fixed = p

	## Varying parameters
	I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak = p_var

    ## Fixed parameters
    Ca_i_0 = p_fixed[10]

	I_FP_range = []
    Ca_FP_range = []
	V_saddle = []
	Ca_saddle = []
    I_saddle = []
	V_stable = []
	Ca_stable = []
    I_stable = []
	V_unstable = []
	Ca_unstable = []
    I_unstable = []

	for v in V_FP_range
		solve_membrane_grad_0_I_ = nlsolve((F,x) ->nlsolve_membrane_grad_0_I_!(F,x,[p_var,p_fixed,v]), [-1.,Ca_i_0],ftol=1E-7,xtol=1E-7)
		push!(I_FP_range,solve_membrane_grad_0_I_.zero[1])
        push!(Ca_FP_range,solve_membrane_grad_0_I_.zero[2])  
		lambda = real.(eigvals(J_matrix(v,Ca_FP_range[end],p)) )
		if maximum(lambda) <=0
			push!(V_stable,v)
			push!(I_stable,solve_membrane_grad_0_I_.zero[1])
            push!(Ca_stable,solve_membrane_grad_0_I_.zero[2])
			push!(V_saddle,v)
			push!(I_saddle,NaN)
            push!(Ca_saddle,NaN)
			push!(V_unstable,v)
			push!(I_unstable,NaN)
            push!(Ca_unstable,NaN)
		else
			if maximum(lambda[lambda.!=maximum(lambda)]) <=0
				push!(V_saddle,v)
				push!(I_saddle,solve_membrane_grad_0_I_.zero[1])
                push!(Ca_saddle,solve_membrane_grad_0_I_.zero[2])
				push!(V_stable,v)
				push!(I_stable,NaN)
                push!(Ca_stable,NaN)
				push!(V_unstable,v)
				push!(I_unstable,NaN)
                push!(Ca_unstable,NaN)
			else
				push!(V_unstable,v)
				push!(I_unstable,solve_membrane_grad_0_I_.zero[1])
                push!(Ca_unstable,solve_membrane_grad_0_I_.zero[2])
				push!(V_saddle,v)
				push!(I_saddle,NaN)
                push!(Ca_saddle,NaN)
				push!(V_stable,v)
				push!(I_stable,NaN)
                push!(Ca_stable,NaN)
			end
		end
	end

	I_SN = maximum(I_FP_range[V_FP_range.<=-40])
	i_I_SN = findfirst(x -> x == I_SN, I_FP_range)
	V_SN = V_FP_range[i_I_SN]
    Ca_SN = Ca_FP_range[i_I_SN]

	V_FP_stable = []
	V_FP_saddle = []
	V_FP_unstable = []
	I_saddle_ = []
	I_stable_ = []
	I_unstable_ = []
    Ca_FP_saddle = []
	Ca_FP_stable = []
	Ca_FP_unstable = []

	I_l_lf = minimum([IX,maximum(I_stable)-0.5])
	I_h_lf = I_SN

	for i in range(I_l_lf,stop=I_h_lf-1e-3,length=n_IX_to_ISN)
		p_var_i = [i,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
		FP_l = nlsolve((F,x) ->nlsolve_membrane_grad_0_V_!(F,x,[p_var_i,p_fixed]), [V_SN-0.1,Ca_SN-1e-5],ftol=1E-7,xtol=1E-7)
		res_h = nlsolve_membrane_grad_0_V_!([NaN,NaN],FP_l.zero,[p_var_i,p_fixed])
		println(res_h)
		if maximum(res_h) < 1e-6
			lambda = real.(eigvals(J_matrix(FP_l.zero[1],FP_l.zero[2],[p_var_i,p_fixed])) )
			if maximum(lambda) <=0
				push!(V_FP_stable,FP_l.zero[1])
                push!(Ca_FP_stable,FP_l.zero[2])
				push!(I_stable_,i)
			else
				if maximum(lambda[lambda.!=maximum(lambda)]) <=0
					push!(V_FP_saddle,FP_l.zero[1])
                    push!(Ca_FP_saddle,FP_l.zero[2])
					push!(I_saddle_,i)
				else
					push!(V_FP_unstable,FP_l.zero[1])
                    push!(Ca_FP_unstable,FP_l.zero[2])
					push!(I_unstable_,i)
				end
			end
		end

	end

	#==
	# Display the results
    plotly()
	plt_V = Plots.plot(fontfamily="Computer Modern", xlabel=" ",ylabel="V",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(IX,round(I_SN+1)),size=(800,400))
    plot!(I_stable ,V_stable ,lc=:green,lw=2)
	plot!(I_saddle ,V_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(I_unstable,V_unstable , lc=:red, ls=:dash,lw=2)
    scatter!(I_SN .*ones(2), V_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)
	if n_IX_to_ISN>0
		scatter!(I_unstable_,V_FP_unstable,mc=:red,markerstrokewidth=0)
		scatter!(I_stable_,V_FP_stable,mc=:green,markerstrokewidth=0)
		scatter!(I_saddle_,V_FP_saddle,mc=:orange,markerstrokewidth=0)
	end

    plt_Ca = Plots.plot(fontfamily="Computer Modern", xlabel="I",ylabel="Ca",legend=:outertopright,legendfontsize=13,xlabelfontsize=15,ylabelfontsize=15,tickfontsize=12,xlims=(IX,round(I_SN+1)),ylims=(-0.005,(Ca_SN*10)),size=(800,400))
    plot!(I_stable ,Ca_stable ,lc=:green,lw=2)
	plot!(I_saddle ,Ca_saddle , lc=:orange, ls=:dot,lw=2)
	plot!(I_unstable,Ca_unstable , lc=:red, ls=:dash,lw=2)
    scatter!(I_SN .*ones(2), Ca_SN .*ones(2),mc=:black,markerstrokewidth=0,label=:none)
	plot!(legend=:none)
	if n_IX_to_ISN>0
		scatter!(I_unstable_,Ca_FP_unstable,mc=:red,markerstrokewidth=0)
		scatter!(I_stable_,Ca_FP_stable,mc=:green,markerstrokewidth=0)
		scatter!(I_saddle_,Ca_FP_saddle,mc=:orange,markerstrokewidth=0)
	end
	
    plt = Plots.plot(plt_V,plt_Ca,layout=(2,1),size=(600,600))
    display(plt)
    gr()
	==#

	if n_IX_to_ISN<=0
	    return I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable
    else
	    return I_SN,V_SN,Ca_SN,V_saddle,Ca_saddle,I_saddle,V_stable,Ca_stable,I_stable,V_unstable,Ca_unstable,I_unstable,V_FP_saddle,Ca_FP_saddle,I_saddle_,V_FP_stable,Ca_FP_stable,I_stable_,V_FP_unstable,Ca_FP_unstable,I_unstable_
    end
end


#find_fixed_points_ = find_fixed_points([p_var,p_fixed],range(-200.,stop=-20.,step=0.1),-4,0)