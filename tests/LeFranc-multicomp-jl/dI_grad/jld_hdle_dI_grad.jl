#Functions used in script_dI_grad.jl with compartment.jl for jld files management

function save_dI(p_var,find_dI, jld_name)
    save(jld_name, "p_var", p_var, "dI", find_dI[1], "I1", find_dI[2], "I2", find_dI[3], "f_min", find_dI[4])
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
function load_dI(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        p_var = []
        dI = []
        I1 = []
        I2 = []
        f_min = []
        push!(p_var,load_array(jld,"p_var"))
        push!(dI,load_array(jld,"dI"))
        push!(I1,load_array(jld,"I1"))
        push!(I2,load_array(jld,"I2"))
        push!(f_min,load_array(jld,"f_min"))
    else
        p_var = NaN
        dI = NaN
        I1 = NaN
        I2 = NaN
        f_min = NaN
        println("Not found")
    end
    find_dI = []
    push!(find_dI,dI[1])
    push!(find_dI,I1[1])
    push!(find_dI,I2[1])
    push!(find_dI,f_min[1])
    return p_var,find_dI
end