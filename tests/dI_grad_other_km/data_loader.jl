
include("membrane._other_kmjl")
include("fixedpoints_other_km.jl")

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
k = 1.0e7 # [nm.cm^(-1)] #TYPO
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
pCaLs_range_KirCaLs = range(0.,stop=0.01725/1000,length=47)[1:2:end]
gKir_range_KirCaLs = range(0.,stop=0.3,length=51)[1:3:end]
pCaLs_range_KMCaLs = range(0.,stop=0.01725/1000,length=47)[1:2:end]
gKM_range_KMCaLs = range(0.,stop=0.3,length=51)[1:3:end]

p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

## Estimated IX 
IX_CaLs = -5 .* ones(length(pCaLs_range_KirCaLs))

IX_KirCaLs = []
for ix in IX_CaLs
    IX = ix .*ones(length(gKir_range_KirCaLs))
    push!(IX_KirCaLs,IX)
end

IX_KMCaLs = []
for ix in IX_CaLs
    IX = ix .*ones(length(gKM_range_KMCaLs))
    push!(IX_KMCaLs,IX)
end
