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




jld_name = "tests/fI_curve_altgt/data/jld_3/fI_KirKMCaLs__"

f_list_KirKMCaLs_fI = []
V_max_c_KirKMCaLs_fI = []
V_min_c_KirKMCaLs_fI = []
I_list_KirKMCaLs_fI = []
I_list_KirKMCaLs_st = []

n_batch = 1
stp_I_list_batch_1=0.01
stp_I_list_batch_2=0.1
stp_I_list_batch_3=0.25
stp_I = 0.5
stp_I_2= 3

for i_pCaLs_range_KirKMCaLs in 1:2
    f_list_ = []
    V_max_c_ = []
    V_min_c_ = []
    I_list_ = []
    I_list_st_ = []
    for i_gKir_range_KirKMCaLs in 1:2
        f_list__ = []
        V_max_c__ = []
        V_min_c__ = []
        I_list__ = []
        I_list_st__ = []
        for i_gKM_range_KirKMCaLs in 1:2
            #if i_gKir_range_KirKMCaLs+i_gKM_range_KirKMCaLs>2
                f_list___ = []
                V_max_c___ = []
                V_min_c___ = []
                I_list___ = []
                I_list_st___ = []
                for i_batch in 1:n_batch
                    jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KirKMCaLs)_i_Kir_$(i_gKir_range_KirKMCaLs)_i_gKM_$(i_gKM_range_KirKMCaLs)_i_batch_$(i_batch).jld")
                    #println(jld_name_i)
                    try
                        p_var_,find_full_bif_= load_full_bif(jld_name_i)
                        push!(f_list___,reverse(find_full_bif_[1]))
                        push!(V_max_c___,reverse(find_full_bif_[2]))
                        push!(V_min_c___,reverse(find_full_bif_[3]))
                        
                            
                            if I2_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]-I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs] >stp_I
                                I_list_batch_1 = collect((I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]-stp_I):stp_I_list_batch_1:(I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]+stp_I))
                                I_list_batch_2a = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I2_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]))
                                I_list_batch_2b = collect((I2_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]):stp_I_list_batch_2:(I2_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]+stp_I))
                                I_list_batch_3 = collect((I_list_batch_2b[end]+stp_I_list_batch_3):stp_I_list_batch_3:(maximum([10,I_list_batch_2b[end]+stp_I_2])))
                                
                                I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2a,I_list_batch_2b,I_list_batch_3)
                            else
                                I_list_batch_1 = collect((I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]-stp_I):stp_I_list_batch_1:(I1_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]+stp_I))
                                I_list_batch_2 = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I_list_batch_1[end]+10*stp_I_list_batch_2))
                                I_list_batch_3 = collect((I_list_batch_2[end]+stp_I_list_batch_3):stp_I_list_batch_3:(maximum([10,I_list_batch_1[end]+stp_I_2])))
                                
                                I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2,I_list_batch_3)
                            end
                    
                        ind_batch = Int(ceil(length(I_list_batch_)/n_batch))
                        I_list_batch = I_list_batch_[( (Int(1+ ind_batch*(i_batch-1))) : Int(minimum([ind_batch + ind_batch*(i_batch-1),length(I_list_batch_)]))) ]
                        push!(I_list___,I_list_batch)
                        
                        I_list_st = I_list_batch[I_list_batch.<=I2_KirKMCaLs_fI[i_pCaLs_range_KirKMCaLs,i_gKir_range_KirKMCaLs,i_gKM_range_KirKMCaLs]]
                        push!(I_list_st___,I_list_st)
                    catch err
                        jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_range_KirKMCaLs)_i_Kir_$(i_gKir_range_KirKMCaLs)_i_gKM_$(i_gKM_range_KirKMCaLs)_i_batch_$(i_batch).jld")
                        println(jld_name_i)
                        println("ouh")
                    end       
                end
                push!(f_list__,reduce(vcat,f_list___))
                push!(V_max_c__,reduce(vcat,V_max_c___))
                push!(V_min_c__ ,reduce(vcat,V_min_c___))
                push!(I_list__,reduce(vcat,I_list___))
                push!(I_list_st__,reduce(vcat,I_list_st___))
            #end
        end
        println("in the kir loop : l(f_list__) = $(length(f_list__))")
        #
        push!(f_list_,f_list__)
        push!(V_max_c_,V_max_c__)
        push!(V_min_c_ ,V_min_c__)
        push!(I_list_,I_list__)
        push!(I_list_st_,I_list_st__)
        ==#
    end
    println("in the cal loop : l(f_list_) =$(length(f_list_))")

    push!(f_list_KirKMCaLs_fI,f_list_)
    push!(V_max_c_KirKMCaLs_fI,V_max_c_)
    push!(V_min_c_KirKMCaLs_fI,V_min_c_)
    push!(I_list_KirKMCaLs_fI,I_list_)
    push!(I_list_KirKMCaLs_st,I_list_st_)
end
size(f_list_KirKMCaLs_fI[1][1][1])

==#

Plots.plot(reverse.(I_list_KirKMCaLs_fI),reverse.(f_list_KirKMCaLs_fI),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)
Plots.plot(reverse.(I_list_KirKMCaLs_fI[1]),reverse.(f_list_KirKMCaLs_fI[1]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)
Plots.plot(reverse.(I_list_KirKMCaLs_fI[2]),reverse.(f_list_KirKMCaLs_fI[2]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)

#scatter!(reverse(I_list_KirKMCaLs_st),reverse(0 .*I_list_KirKMCaLs_st),markershape=:rect,markerstrokewidth=0.,ms=2,legend=:outertopright)

plt = Plots.plot()
for i in 1:2
    for j in 1:2
        for k in 1:2
            
            try 
            plot!(plt,I_list_KirKMCaLs_fI[i][j][k],(f_list_KirKMCaLs_fI[i][j][k]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Paired_8,label="CaL = $(i); Kir = $(j); KM = $(k)")
            display(plt)
            catch err 
                println("i=$i , j=$j , k=$k")
            end
        end
    end
end
plot!(legend=:outertopright)