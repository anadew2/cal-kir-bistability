
include("membrane_CaLf_replaced_by_CaLs.jl")
include("fixedpoints_CaLf_replaced_by_CaLs.jl")

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

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
I(t) = 37  + 2 *pulse(t, 2000, 2000+500)  # [µA/cm^2]
gNa = 30 #mS/cm²
gKDR = 4 #mS/cm²
gleak = 0.03268 #mS/cm²
gKir = 0.           # mS/cm²
gKM = 0.           # mS/cm²
pCaLf = 0. /1000    # [cm/s]
pCaLs = 0. /1000   # [cm/s]

## Range of varying parameters
pCaLs_range_KirCaLs = range(0.,stop=0.01725/1000 *3 ,length=24)
gKir_range_KirCaLs = range(0.,stop=1.,length=11)
pCaLs_range_KMCaLs = range(0.,stop=0.01725/1000 *3,length=24)
gKM_range_KMCaLs = range(0.,stop=1.5,length=16)

p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

## Estimated IX from simulation
IX_KirCaLs_1_1 = -0.15
IX_KirCaLs_1_end = -3
IX_KirCaLs_end_end = 9.
IX_KirCaLs_end_1 = 16

IX_KMCaLs_1_1 = -0.15
IX_KMCaLs_1_end = -3.
IX_KMCaLs_end_end = 36
IX_KMCaLs_end_1 = 37.5

## Estimated IX 
m_IX_Kir = (IX_KirCaLs_end_1-IX_KirCaLs_1_1)/gKir_range_KirCaLs[end]
m_IX_CaLs = (IX_KirCaLs_1_end-IX_KirCaLs_1_1)/pCaLs_range_KirCaLs[end]

IX_CaLs = IX_KirCaLs_1_1 .+ m_IX_CaLs .*pCaLs_range_KirCaLs
IX_KirCaLs = []
for ix in IX_CaLs
    IX = ix .+ m_IX_Kir .*0.75 .*gKir_range_KirCaLs
    push!(IX_KirCaLs,IX)
end

m_IX_KM = (IX_KMCaLs_end_1-IX_KirCaLs_1_1)/gKM_range_KMCaLs[end]
m_IX_CaLs = (IX_KMCaLs_1_end-IX_KMCaLs_1_1)/pCaLs_range_KMCaLs[end]

IX_CaLs = IX_KMCaLs_1_1 .+ m_IX_CaLs .*pCaLs_range_KMCaLs
IX_KMCaLs = []
for ix in IX_CaLs
    IX = ix .+ m_IX_KM .*0. .*gKM_range_KMCaLs .-20
    push!(IX_KMCaLs,IX)
end