
include("../../membrane.jl")
include("../../fixedpoints.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

function load_jld(jld_name)
	data = []
	try 
		push!(data,load(jld_name))
	catch err
		push!(data ,NaN)
		println("Not found...")
	end
	return data[1]
end
function load_array(data,array_name)
	return data[array_name]
end

function load_dI_mat(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        dI = []
        I1 = []
        V_SN = []
        I2 = []
        push!(dI,load_array(jld,"dI_mat"))
        push!(I1,load_array(jld,"I1_mat"))
        push!(V_SN,load_array(jld,"VSN_mat"))
        push!(I2,load_array(jld,"I2_mat"))
    else
        dI = NaN
        I1 = NaN
        V_SN = NaN
        I2 = NaN
        println("Not found")
    end
    return dI[1],I1[1],I2[1],V_SN[1]
end

function σ_membrane!(du,u,p,t)
   
    ## Noise parameters
    sigma_V,sigma_mNa,sigma_hNa,sigma_mKDR,sigma_mCaLf,sigma_hCaLf,sigma_mCaLs,sigma_hCaLs,sigma_mKir,sigma_mKM,sigma_Ca = p[3]

    ## Gradients
	du[1] = sigma_V
	du[2] = sigma_mNa
	du[3] = sigma_hNa
	du[4] = sigma_mKDR
    du[5] = sigma_mCaLf
	du[6] = sigma_hCaLf
	du[7] = sigma_mCaLs
	du[8] = sigma_hCaLs
	du[9] = sigma_mKir
    du[10] = sigma_mKM
    du[11] = sigma_Ca

end


## Fixed parameters
	C = 1 #µF/cm²
	eNa = 50 # [mV]
	eK = -90 # [mV]
	eCa = 120 # [mV]
	eleak = -60.1 # [mV]
	Ca_o = 2 # Extracellular Ca concentration [mM]

## Intracellular calcium dynamics parameters
	k = 1.0e7 # [nm.cm^(-1)] 
	F = 96520 # Faraday default value [C/mol]
	d = 1 # [nm]
	Ca_i_0 = 5.0e-5 # Initial intracellular Ca concentration [mM]
	tau_Ca = 10 # [ms]

p_fixed = [C,eNa,eK,eCa,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca]

## Output of dI_grad tests computing the size of the bistability window 
pCaLs_range_KMCaLs = range(0.,stop=0.01725/1000,length=47)
gKM_range_KMCaLs = range(0.,stop=0.3,length=51)

jld_name = "ChannelUpdate/cluster/dI_noise/data/input/dI_KMCaLs"
dI_KMCaLs,I1_KMCaLs,I2_KMCaLs,VSN_KMCaLs = load_dI_mat(string(jld_name,"_1e-2_1e-5.jld"))
 

## Varying parameters
	I1bist = I1_KMCaLs[41,end]
	I2bist = I2_KMCaLs[41,end]
	I0 = I1_KMCaLs[41,end] +1e-3

	I(t) = I0  #+  2. *pulse(t, 15000, 15000+1000)  # [µA/cm^2]
	gNa = 30 #mS/cm²
	gKDR = 4 #mS/cm²
	gleak = 0.03268 #mS/cm²
	gKir = 0.           #mS/cm²
	gKM = gKM_range_KMCaLs[end]           #mS/cm²
	pCaLf = 0. /1000    # [cm/s]
	pCaLs = pCaLs_range_KMCaLs[41]   # [cm/s]

p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]


## Initial conditions for spiking state
	V_ic = -90. # [mV]
	Ca_ic = 5.0e-5 # [mM]
	ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

	## Noise-free Simulation
	prob = ODEProblem(membrane!,ic,(0.,80000.),[p_var,p_fixed])
	sol = solve(prob,dtmax=1,Rodas5P())
	sol_cb = solve(prob,callback=cb,save_everystep=false,save_start=false,save_end=false,Rodas5P())
	plot(sol.t,sol[1,:])
	display(Plots.scatter!(sol_cb.t,sol_cb[1,:]))
	#Plots.scatter(sol_cb.t[2:end],1000 ./(sol_cb.t[2:end]-sol_cb.t[1:end-1]),xlabel="t",ylabel="f",fontfamiliy="Computer Modern")
ic_cyc = sol[end]

f_eq = mean(1000 ./(sol_cb.t[sol_cb.t .>= sol_cb.t[end]*3/4][2:end]-sol_cb.t[sol_cb.t .>= sol_cb.t[end]*3/4][1:end-1]))


I(t) = I0 # [µA/cm^2]
p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

## Initial conditions for resting state
	V_ic = -70. # [mV]
	Ca_ic = 5.0e-5 # [mM]
	ic = [V_ic,mNa(V_ic),hNa(V_ic),mKDR(V_ic),mCaLf(V_ic),hCaLf(V_ic),mCaLs(V_ic),hCaLs(V_ic),mKir(V_ic),mKM(V_ic),Ca_ic]

	prob = ODEProblem(membrane!,ic,(0.,20000.),[p_var,p_fixed])
	sol_st = solve(prob,dtmax=1,Rodas5P())
	display(plot(sol_st.t,sol_st[1,:],xlabel="t",ylabel="V",fontfamiliy="Computer Modern"))
ic_st = sol_st[end]

