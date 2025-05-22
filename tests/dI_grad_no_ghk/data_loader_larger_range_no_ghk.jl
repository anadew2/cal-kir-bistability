
include("membrane_no_ghk.jl")
include("fixedpoints_no_ghk.jl")

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
I(t) = 6  + 1 *pulse(t, 2000, 2000+500)  # [µA/cm^2]
gNa = 30 #mS/cm²
gKDR = 4 #mS/cm²
gleak = 0.03268 #mS/cm²
gKir = 0.           # [mS/cm^2]
gKM = 0.           # [mS/cm^2]
gCaLf = 0.          # [mS/cm^2]
gCaLs = 0.          # [mS/cm^2]

## Range of varying parameters
gCaLs_range_KirCaLs = range(0.,stop=0.3,length=16)
gKir_range_KirCaLs = range(0.,stop=0.5,length=6)
gCaLs_range_KMCaLs = range(0.,stop=0.4,length=21)
gKM_range_KMCaLs = range(0.,stop=1.5,length=16)

p_var = [I,gNa,gKDR,gKir,gKM,gCaLf,gCaLs,gleak]

## Estimated IX from simulation
IX_KirCaLs_1_1 = -0.15
IX_KirCaLs_1_end = -9.2
IX_KirCaLs_end_end = -4.2
IX_KirCaLs_end_1 = 5.9

IX_KMCaLs_1_1 = -0.15
IX_KMCaLs_1_end = -3.
IX_KMCaLs_end_end = 29.1
IX_KMCaLs_end_1 = 37.2

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
    IX = ix .+ m_IX_KM .*0.5 .*gKM_range_KMCaLs
    push!(IX_KMCaLs,IX)
end