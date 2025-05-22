
include("membrane_outw.jl")
include("fixedpoints_outw.jl")

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
I(t) = 0.  # [µA/cm^2]
gNa = 30 #mS/cm²
gKDR = 4 #mS/cm²
gleak = 0.03268 #mS/cm²
gKir = 0.           # [mS/cm^2]
gKM = 0.           # [mS/cm^2]
pCaLf = 0. /1000    # [cm/s]
pCaLs = 0. /1000   # [cm/s]

## Range of varying parameters
pCaLs_KirCaLs_outw = range(0.,stop=0.01725/1000,length=24)[21]
gKir_KirCaLs_outw = range(0.,stop=0.5,length=51)
pCaLs_KMCaLs_outw = range(0.,stop=0.01725/1000,length=24)[21]
gKM_KMCaLs_outw = range(0.,stop=0.5,length=51)

p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]
