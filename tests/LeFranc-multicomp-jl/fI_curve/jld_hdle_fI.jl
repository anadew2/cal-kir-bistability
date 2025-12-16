
function save_full_bif(p, find_full_bif, jld_name)
    save(jld_name, "p_fixed", p[1][2:end], "p_dend",p[2],"p_soma",p[3], "f_list", find_full_bif[1], "Vd_max_c", find_full_bif[2], "Vd_min_c", find_full_bif[3],  "Cad_max_c", find_full_bif[4], "Cad_min_c", find_full_bif[5],  "Vs_max_c", find_full_bif[6], "Vs_min_c", find_full_bif[7],  "Cas_max_c", find_full_bif[8], "Cas_min_c", find_full_bif[9])
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

function load_full_bif(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        p_fixed = []
        p_dend = []
        p_soma = []
        f_list = []
        Vd_max_c = []
        Vd_min_c = []
        Cad_max_c = []
        Cad_min_c = []
        Vs_max_c = []
        Vs_min_c = []
        Cas_max_c = []
        Cas_min_c = []
        push!(p_fixed,load_array(jld,"p_fixed"))
        push!(p_dend,load_array(jld,"p_dend"))
        push!(p_soma,load_array(jld,"p_soma"))
        push!(f_list,load_array(jld,"f_list"))
        push!(Vd_max_c,load_array(jld,"Vd_max_c"))
        push!(Vd_min_c,load_array(jld,"Vd_min_c"))
        push!(Cad_max_c,load_array(jld,"Cad_max_c"))
        push!(Cad_min_c,load_array(jld,"Cad_min_c"))
        push!(Vs_max_c,load_array(jld,"Vs_max_c"))
        push!(Vs_min_c,load_array(jld,"Vs_min_c"))
        push!(Cas_max_c,load_array(jld,"Cas_max_c"))
        push!(Cas_min_c,load_array(jld,"Cas_min_c"))
    else
        p_fixed = NaN
        p_dend = NaN
        p_soma = NaN
        f_list = NaN
        Vd_min_c = NaN
        Cad_max_c = NaN
        Cad_min_c = NaN
        Vs_max_c = NaN
        Vs_min_c = NaN
        Cas_max_c = NaN
        Cas_min_c = NaN
        println("Not found")
    end
    p = [vcat(0,p_fixed),p_dend,p_soma]
    find_full_bif = []
    push!(find_full_bif,f_list[1])
    push!(find_full_bif,Vd_max_c[1])
    push!(find_full_bif,Vd_min_c[1])
    push!(find_full_bif,Cad_max_c[1])
    push!(find_full_bif,Cad_min_c[1])
    push!(find_full_bif,Vs_max_c[1])
    push!(find_full_bif,Vs_min_c[1])
    push!(find_full_bif,Cas_max_c[1])
    push!(find_full_bif,Cas_min_c[1])
    return p,find_full_bif
end
function load_dI_mat(jld_name)
    jld = load_jld(jld_name)
    if jld != NaN
        dI = []
        I1 = []
        f_min = []
        I2 = []
        push!(dI,load_array(jld,"dI_mat"))
        push!(I1,load_array(jld,"I1_mat"))
        push!(f_min,load_array(jld,"f_min_mat"))
        push!(I2,load_array(jld,"I2_mat"))
    else
        dI = NaN
        I1 = NaN
        f_min = NaN
        I2 = NaN
        println("Not found")
    end
    return dI[1],I1[1],I2[1],f_min[1]
end
