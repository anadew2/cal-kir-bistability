

heaviside(t)=0*(t<0)+1*(t>=0)
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

 
jld_name = "tests/LeFranc-multicomp-jl/fi_curve/data/input/dI_Kir_s_CaLs_d_compart_ext_"
dI_Kir_s_CaLs_d,I1_Kir_s_CaLs_d,I2_Kir_s_CaLs_d,f_min_Kir_s_CaLs_d = load_dI_mat(string(jld_name,".jld"))

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
    I(t) = 0/10^(6)   # [µA]

    p_fixed = [I,C,eNa,eK,eCa,eCaAN,eleak,Ca_o,k,F,d,Ca_i_0,tau_Ca,Ls,Ds,Ld,Dd,Ra]


## Varying parameters
    gNa_s = 30 #mS/cm²
    gKDR_s = 4 #mS/cm²
    gleak_s = 0.03268 #mS/cm²
    gKir_s = 0.02             # mS/cm² 
    gKM_s = 0.            # mS/cm²
    gKCa_s = 0.03           # mS/cm²
    gCaAN_s = 0.            # mS/cm²
    pCaLf_s = 0.15  /1000    # [cm/s]
    pCaLs_s = 0. /1000   # [cm/s]

    gNa_d = 0 #mS/cm²
    gKDR_d = 0 #mS/cm²
    gleak_d = 0.03 #mS/cm²
    gKir_d = 0.           # mS/cm²
    gKM_d = 0.            # mS/cm²
    gKCa_d = 0.0185         # mS/cm²
    gCaAN_d = 0.125         # mS/cm²
    pCaLf_d = 0. /1000    # [cm/s]
    pCaLs_d = 0.02  /1000   # [cm/s]

    p_var_d = [gNa_d,gKDR_d,gKir_d,gKM_d,gKCa_d,gCaAN_d,pCaLf_d,pCaLs_d,gleak_d]
    p_var_s = [gNa_s,gKDR_s,gKir_s,gKM_s,gKCa_s,gCaAN_s,pCaLf_s,pCaLs_s,gleak_s]


## Range of varying parameters
ind_CaLs_d = [3,16,29]
ind_Kir_s = [4,21,38]
ind_Kir_d = [4,21,38]
pCaLs_d_range_KirCaLs_fI = range(0.,stop=0.03/1000,length=31)[ind_CaLs_d[:]]
gKir_s_range_KirCaLs_fI = range(0.,stop=0.2,length=41)[ind_Kir_s[:]]
gKir_d_range_KirCaLs_fI = range(0.,stop=0.2,length=41)[ind_Kir_d[:]]


I1_Kir_s_CaLs_d_fI = I1_Kir_s_CaLs_d[ind_CaLs_d[:],ind_Kir_s[:]]
I2_Kir_s_CaLs_d_fI = I2_Kir_s_CaLs_d[ind_CaLs_d[:],ind_Kir_s[:]]


jld_name = "tests/LeFranc-multicomp-jl/fi_curve/data/input/dI_Kir_d_CaLs_d_compart_ext_"
dI_Kir_d_CaLs_d,I1_Kir_d_CaLs_d,I2_Kir_d_CaLs_d,f_min_Kir_d_CaLs_d = load_dI_mat(string(jld_name,".jld"))

I1_Kir_d_CaLs_d_fI = I1_Kir_d_CaLs_d[ind_CaLs_d[:],ind_Kir_d[:]]
I2_Kir_d_CaLs_d_fI = I2_Kir_d_CaLs_d[ind_CaLs_d[:],ind_Kir_d[:]]

println(gKir_d_range_KirCaLs_fI)
println(I1_Kir_d_CaLs_d_fI)
println(I2_Kir_d_CaLs_d_fI)

