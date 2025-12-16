using DifferentialEquations
using Distributions
using JLD 

function save_inputs_of_script_dI_mc(C_sampled,gleak_sampled,gNa_sampled,jld_name)
    save(jld_name, "C_sampled", C_sampled, "gleak_sampled", gleak_sampled, "gNa_sampled", gNa_sampled)
end
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
function load_inputs_of_script_dI_mc(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
		C_sampled = []
		gleak_sampled = []
		gNa_sampled = []
		push!(C_sampled,load_array(jld,"C_sampled"))
		push!(gleak_sampled,load_array(jld,"gleak_sampled"))
		push!(gNa_sampled,load_array(jld,"gNa_sampled"))
    else
		C_sampled = NaN
		gleak_sampled = NaN
		gNa_sampled = NaN
        println("Not found")
    end
    return C_sampled[1],gleak_sampled[1],gNa_sampled[1]
end


### GNA VARIABILITY ###

    var_level = 0.2
    #==
    # This block generates the random data, within the interval determined by the value of var_level
        C_var_mat = rand(Uniform(C-var_level*C,C+var_level*C),300)
        gleak_var_mat = rand(Uniform(gleak-var_level*gleak,gleak+var_level*gleak),300)
        gNa_var_mat = rand(Uniform(gnabar-var_level*gnabar,gnabar+var_level*gnabar),300)#rand(Uniform(gnabar-0.1*gnabar,gnabar+0.1*gnabar),100)

        inputs_jld_name = "tests/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
        save_inputs_of_script_dI_mc(C_var_mat,gleak_var_mat,gNa_var_mat, inputs_jld_name)
    ==#

    inputs_jld_name = "tests/dI_var_mc/data/input/inputs_of_script_dI_var_mc__var$(Int(var_level*100))_300.jld"
    C_var_mat,gleak_var_mat,gNa_var_mat = load_inputs_of_script_dI_mc(inputs_jld_name)


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
    I(t) = 0. # [µA/cm^2]
    gNa = 30 #mS/cm²
    gKDR = 4 #mS/cm²
    gleak = 0.03268 #mS/cm²
    gKir = 0.           # [mS/cm^2]
    gKM = 0.           # [mS/cm^2]
    pCaLf = 0. /1000    # [cm/s]
    pCaLs = 0. /1000   # [cm/s]

    p_var = [I,gNa,gKDR,gKir,gKM,pCaLf,pCaLs,gleak]

    ## Range of varying parameters
    pCaLs_range_KirCaLs = range(0.,stop=0.01725/1000,length=24)
    gKir_range_KirCaLs = range(0.,stop=1.5,length=16)
    pCaLs_range_KMCaLs = range(0.,stop=0.01725/1000,length=24)
    gKM_range_KMCaLs = range(0.,stop=1.5,length=16)

    gKir_KirCaLs = gKir_range_KirCaLs[4]
    gKM_KirCaLs = 0.
    pCaLs_KirCaLs = pCaLs_range_KirCaLs[21]

    gKir_KMCaLs = 0.
    gKM_KMCaLs = gKM_range_KMCaLs[4]
    pCaLs_KMCaLs = pCaLs_range_KMCaLs[21]

