
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

 
jld_name = "tests/fI_curve_altgt/data/input/dI_KirKMCaLs"
dI_KirKMCaLs,I1_KirKMCaLs,I2_KirKMCaLs,VSN_KirKMCaLs = load_dI_mat(string(jld_name,".jld"))


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

## Varying parameters
I(t) = 0.  + 0. *pulse(t, 5000, 5000+500)  # [µA/cm^2]
gNa = 30 #mS/cm²
gKDR = 4 #mS/cm²
gleak = 0.03268 #mS/cm²
gKir = 0.           # mS/cm²
gKM = 0.           # mS/cm²
pCaLf = 0. /1000    # [cm/s]
pCaLs = 0. /1000   # [cm/s]

## Range of varying parameters
pCaLs_range_KirKMCaLs = range(0.,stop=0.01725/1000,length=47)[1:10:end]
gKir_range_KirKMCaLs = range(0.,stop=0.3,length=31)
gKM_range_KirKMCaLs = range(0.,stop=0.3,length=31)

pCaLs_range_KirKMCaLs_fI = pCaLs_range_KirKMCaLs[3:2:5]
gKir_range_KirKMCaLs_fI = gKir_range_KirKMCaLs[3:26:29]
gKM_range_KirKMCaLs_fI = gKM_range_KirKMCaLs[3:26:29]


I1_KirKMCaLs_fI = I1_KirKMCaLs[3:2:5,3:26:29,3:26:29]
I2_KirKMCaLs_fI = I2_KirKMCaLs[3:2:5,3:26:29,3:26:29]


p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
