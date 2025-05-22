using JLD

include("../../membrane.jl")
include("../../fixedpoints.jl")

function save_full_bif(p_var, find_full_bif, jld_name)
    save(jld_name, "p_var", p_var, "f_list", find_full_bif[1], "V_max_c", find_full_bif[2], "V_min_c", find_full_bif[3])
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
        p_var = []
        f_list = []
        V_max_c = []
        V_min_c = []
        push!(p_var,load_array(jld,"p_var"))
        push!(f_list,load_array(jld,"f_list"))
        push!(V_max_c,load_array(jld,"V_max_c"))
        push!(V_min_c,load_array(jld,"V_min_c"))
    else
        p_var = NaN
        f_list = NaN
        V_max_c = NaN
        V_min_c = NaN
        println("Not found")
    end
    find_full_bif = []
    push!(find_full_bif,f_list[1])
    push!(find_full_bif,V_max_c[1])
    push!(find_full_bif,V_min_c[1])
    return p_var,find_full_bif
end

include("data_loader.jl")


#

jld_name = "ChannelUpdate/cluster/fI_curve/data/jld_1/fI_KirCaLs__"

f_list_KirCaLs_fI = []
V_max_c_KirCaLs_fI = []
V_min_c_KirCaLs_fI = []
I_list_KirCaLs_fI = []
I_list_KirCaLs_st = []

n_batch = 1
stp_I_list_batch_1=0.01
stp_I_list_batch_2=0.1

for i_pCaLs_range_KirCaLs in 1:3
    f_list_ = []
    V_max_c_ = []
    V_min_c_ = []
    I_list_ = []
    I_list_st_ = []
    for i_gKir_range_KirCaLs in 1:3
        f_list__ = []
        V_max_c__ = []
        V_min_c__ = []
        I_list__ = []
        I_list_st__ = []
        for i_batch in 1:n_batch
            jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KirCaLs)_i_gKir_$(i_gKir_range_KirCaLs)_i_batch_$(i_batch).jld")
            #println(jld_name_i)
            try
                p_var_,find_full_bif_= load_full_bif(jld_name_i)
                push!(f_list__,reverse(find_full_bif_[1]))
                push!(V_max_c__,reverse(find_full_bif_[2]))
                push!(V_min_c__,reverse(find_full_bif_[3]))
                
                I_list_batch_1 = collect((I1_KirCaLs_fI[i_pCaLs_range_KirCaLs,i_gKir_range_KirCaLs]-0.5):stp_I_list_batch_1:(I1_KirCaLs_fI[i_pCaLs_range_KirCaLs,i_gKir_range_KirCaLs]+0.5))
                I_list_batch_2 = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I2_KirCaLs_fI[i_pCaLs_range_KirCaLs,i_gKir_range_KirCaLs]+0.5))
                I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2)
            
                ind_batch = Int(ceil(length(I_list_batch_)/n_batch))
                I_list_batch = I_list_batch_[( (Int(1+ ind_batch*(i_batch-1))) : Int(minimum([ind_batch + ind_batch*(i_batch-1),length(I_list_batch_)]))) ]
                push!(I_list__,I_list_batch)
                
                I_list_st = I_list_batch[I_list_batch.<=I2_KirCaLs_fI[i_pCaLs_range_KirCaLs,i_gKir_range_KirCaLs]]
                push!(I_list_st__,I_list_st)
            catch err
                jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KirCaLs)_i_gKir_$(i_gKir_range_KirCaLs)_i_batch_$(i_batch).jld")
                println(jld_name_i)
                println("ouh")
            end       
        end
        
        #
        push!(f_list_,reduce(vcat,f_list__))
        push!(V_max_c_,reduce(vcat,V_max_c__))
        push!(V_min_c_ ,reduce(vcat,V_min_c__))
        push!(I_list_,reduce(vcat,I_list__))
        push!(I_list_st_,reduce(vcat,I_list_st__))
        ==#
        println(f_list_)
    end
    push!(f_list_KirCaLs_fI,f_list_)
    push!(V_max_c_KirCaLs_fI,V_max_c_)
    push!(V_min_c_KirCaLs_fI,V_min_c_)
    push!(I_list_KirCaLs_fI,I_list_)
    push!(I_list_KirCaLs_st,I_list_st_)
end

==#

## KM ##
#
jld_name = "ChannelUpdate/cluster/fI_curve/data/jld_1/fI_KMCaLs__"
f_list_KMCaLs_fI = []
V_max_c_KMCaLs_fI = []
V_min_c_KMCaLs_fI = []
I_list_KMCaLs_fI = []
I_list_KMCaLs_st = []

n_batch = 1
stp_I_list_batch_1=0.01
stp_I_list_batch_2=0.1

for i_pCaLs_range_KMCaLs in 1:3
    f_list_ = []
    V_max_c_ = []
    V_min_c_ = []
    I_list_ = []
    I_list_st_ = []
    for i_gKM_range_KMCaLs in 1:3
        f_list__ = []
        V_max_c__ = []
        V_min_c__ = []
        I_list__ = []
        I_list_st__ = []
        for i_batch in 1:n_batch
            jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KMCaLs)_i_gKM_$(i_gKM_range_KMCaLs)_i_batch_$(i_batch).jld")
            #println(jld_name_i)
            try
                p_var_,find_full_bif_= load_full_bif(jld_name_i)
                push!(f_list__,reverse(find_full_bif_[1]))
                push!(V_max_c__,reverse(find_full_bif_[2]))
                push!(V_min_c__,reverse(find_full_bif_[3]))
                
                I_list_batch_1 = collect((I1_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]-0.5):stp_I_list_batch_1:(I1_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]+0.5))
                I_list_batch_2 = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I2_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]+0.5))
                I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2)
            
                ind_batch = Int(ceil(length(I_list_batch_)/n_batch))
                I_list_batch = I_list_batch_[( (Int(1+ ind_batch*(i_batch-1))) : Int(minimum([ind_batch + ind_batch*(i_batch-1),length(I_list_batch_)]))) ]
                push!(I_list__,I_list_batch)
                
                I_list_st = I_list_batch[I_list_batch.<=I2_KMCaLs_fI[i_pCaLs_range_KMCaLs,i_gKM_range_KMCaLs]]
                push!(I_list_st__,I_list_st)
            catch err
                jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KMCaLs)_i_gKM_$(i_gKM_range_KMCaLs)_i_batch_$(i_batch).jld")
                println(jld_name_i)
            end       
        end
        push!(f_list_,reduce(vcat,f_list__))
        push!(V_max_c_,reduce(vcat,V_max_c__))
        push!(V_min_c_ ,reduce(vcat,V_min_c__))
        push!(I_list_,reduce(vcat,I_list__))
        push!(I_list_st_,reduce(vcat,I_list_st__))
        ==#
        println(f_list_)
    end
    push!(f_list_KMCaLs_fI,f_list_)
    push!(V_max_c_KMCaLs_fI,V_max_c_)
    push!(V_min_c_KMCaLs_fI,V_min_c_)
    push!(I_list_KMCaLs_fI,I_list_)
    push!(I_list_KMCaLs_st,I_list_st_)
end

Plots.plot(reverse.(I_list_KMCaLs_fI),reverse.(f_list_KMCaLs_fI),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)
#scatter!(reverse(I_list_KMCaLs_st),reverse(0 .*I_list_KMCaLs_st),markershape=:rect,markerstrokewidth=0.,ms=2,legend=:outertopright)

Plots.plot(reverse.(I_list_KirCaLs_fI),reverse.(f_list_KirCaLs_fI),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)
#scatter!(reverse(I_list_KirCaLs_st),reverse(0 .*I_list_KirCaLs_st),markershape=:rect,markerstrokewidth=0.,ms=2,legend=:outertopright)

plt = Plots.plot()
for i in 1:3
    for j in 1:3
        try 
        plot!(plt,I_list_KirCaLs_fI[i][j],(f_list_KirCaLs_fI[i][j]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)
        display(plt)
        catch err 
            println("i=$i , j=$j")
        end
    end
end
