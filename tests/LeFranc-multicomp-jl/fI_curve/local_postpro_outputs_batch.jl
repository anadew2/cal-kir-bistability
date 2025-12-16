using JLD

include("../gates_compart.jl")
include("../currents.jl")
include("../compartment.jl")
include("../fixedpoints.jl")

include("jld_hdle_fI.jl")
include("data_loader.jl")
include("fct_fI_curve.jl")

n_batch = 1
stp_I_list_batch_1= 0.1/10^6 #1/10^6 #0.1/10^6 #[µA]
stp_I_list_batch_2= 5/10^6 #[µA]
stp_I = 50/10^6 #[µA]
stp_I_2 = 75/10^6 #[µA]


jld_name = "tests/LeFranc-multicomp-jl/fI_curve/data/jld_3/fI_Kir_s_CaLs_d_compart_ext_"

f_list_Kir_s_CaLs_d_fI = []
Vd_max_c_Kir_s_CaLs_d_fI = []
Vd_min_c_Kir_s_CaLs_d_fI = []
Cad_max_c_Kir_s_CaLs_d_fI = []
Cad_min_c_Kir_s_CaLs_d_fI = []
Vs_max_c_Kir_s_CaLs_d_fI = []
Vs_min_c_Kir_s_CaLs_d_fI = []
Cas_max_c_Kir_s_CaLs_d_fI = []
Cas_min_c_Kir_s_CaLs_d_fI = []
I_list_Kir_s_CaLs_d_fI = []
I_list_Kir_s_CaLs_st = []

for i_pCaLs_d in 1:3
    f_list_ = []
    Vd_max_c_ = []
    Vd_min_c_ = []
    Cad_max_c_ = []
    Cad_min_c_ = []
    Vs_max_c_ = []
    Vs_min_c_ = []
    Cas_max_c_ = []
    Cas_min_c_ = []
    I_list_ = []
    I_list_st_ = []

    for i_gKir_s in 1:3
        f_list__ = []
        Vd_max_c__ = []
        Vd_min_c__ = []
        Cad_max_c__ = []
        Cad_min_c__ = []
        Vs_max_c__ = []
        Vs_min_c__ = []
        Cas_max_c__ = []
        Cas_min_c__ = []
        I_list__ = []
        I_list_st__ = []
        
        for i_batch in 1:n_batch
            jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_d)_i_Kir_$(i_gKir_s)_i_batch_$(i_batch).jld")
            #println(jld_name_i)
            try
                p_var_,find_full_bif_= load_full_bif(jld_name_i)
                push!(f_list__,reverse(find_full_bif_[1]))
                push!(Vd_max_c__,reverse(find_full_bif_[2]))
                push!(Vd_min_c__,reverse(find_full_bif_[3]))
                push!(Cad_max_c__,reverse(find_full_bif_[4]))
                push!(Cad_min_c__,reverse(find_full_bif_[5]))
                push!(Vs_max_c__,reverse(find_full_bif_[6]))
                push!(Vs_min_c__,reverse(find_full_bif_[7]))
                push!(Cas_max_c__,reverse(find_full_bif_[8]))
                push!(Cas_min_c__,reverse(find_full_bif_[9]))
                
                    if I2_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]-I1_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s] >stp_I
                        I_list_batch_1 = collect((I1_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]-stp_I):stp_I_list_batch_1:(I1_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]+stp_I))
                        I_list_batch_2a = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(I2_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]))
                        I_list_batch_2b = collect((I2_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]):stp_I_list_batch_2:(minimum([275/10^6,I2_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]+stp_I+stp_I_2])))
                        I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2a,I_list_batch_2b)
                    else
                        I_list_batch_1 = collect((I1_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]-stp_I):stp_I_list_batch_1:(I1_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]+stp_I))
                        I_list_batch_2 = collect((I_list_batch_1[end]+stp_I_list_batch_2):stp_I_list_batch_2:(minimum([275/10^6,I_list_batch_1[end]+stp_I_list_batch_2+4*stp_I_2])))
                        I_list_batch_ = vcat(I_list_batch_1,I_list_batch_2)
                    end
                    

                    ind_batch = Int(ceil(length(I_list_batch_)/n_batch))
                    I_list_batch = I_list_batch_[( (Int(1+ ind_batch*(i_batch-1))) : Int(minimum([ind_batch + ind_batch*(i_batch-1),length(I_list_batch_)]))) ]

                push!(I_list__,I_list_batch)
                
                I_list_st = I_list_batch[I_list_batch.<=I2_Kir_s_CaLs_d_fI[i_pCaLs_d,i_gKir_s]]
                push!(I_list_st__,I_list_st)
            catch err
                jld_name_i = string(jld_name,"i_pCaLs_$(i_pCaLs_d)_i_Kir_$(i_gKir_s)_i_batch_$(i_batch).jld")
                println(jld_name_i)
                println("ouh")
            end       
        end
        println("---------------- i=$i_gKir_s   ;   j=$i_pCaLs_d ---------------- ")
        println("in the kir loop : l(f_list__) = $(length(f_list__))")
        println("in the kir loop : l(f_list_) = $(length(reduce(vcat,f_list__)))")
        println("in the kir loop : l(I_list_) = $(length(reduce(vcat,I_list__)))")
        
        #
        push!(f_list_,reduce(vcat,f_list__))
        push!(Vd_max_c_,reduce(vcat,Vd_max_c__))
        push!(Vd_min_c_ ,reduce(vcat,Vd_min_c__))
        push!(Cad_max_c_,reduce(vcat,Cad_max_c__))
        push!(Cad_min_c_ ,reduce(vcat,Cad_min_c__))
        push!(Vs_max_c_,reduce(vcat,Vs_max_c__))
        push!(Vs_min_c_ ,reduce(vcat,Vs_min_c__))
        push!(Cas_max_c_,reduce(vcat,Cas_max_c__))
        push!(Cas_min_c_ ,reduce(vcat,Cas_min_c__))
        push!(I_list_,reduce(vcat,I_list__))
        push!(I_list_st_,reduce(vcat,I_list_st__))
        ==#
    end
    println("in the cal loop : l(f_list_) =$(length(f_list_))")

    push!(f_list_Kir_s_CaLs_d_fI,f_list_)
    push!(Vd_max_c_Kir_s_CaLs_d_fI,Vd_max_c_)
    push!(Vd_min_c_Kir_s_CaLs_d_fI,Vd_min_c_)
    push!(Cad_max_c_Kir_s_CaLs_d_fI,Cad_max_c_)
    push!(Cad_min_c_Kir_s_CaLs_d_fI,Cad_min_c_)
    push!(Vs_max_c_Kir_s_CaLs_d_fI,Vs_max_c_)
    push!(Vs_min_c_Kir_s_CaLs_d_fI,Vs_min_c_)
    push!(Cas_max_c_Kir_s_CaLs_d_fI,Cas_max_c_)
    push!(Cas_min_c_Kir_s_CaLs_d_fI,Cas_min_c_)
    push!(I_list_Kir_s_CaLs_d_fI,I_list_)
    push!(I_list_Kir_s_CaLs_st,I_list_st_)
end
size(f_list_Kir_s_CaLs_d_fI[1][1][1])

==#

Plots.plot(reverse.(I_list_Kir_s_CaLs_d_fI).*10^6,reverse.(f_list_Kir_s_CaLs_d_fI),markershape=:circle,markerstrokewidth=0.,markerstrokecolor=:white,ms=3,palette=:tab10,legend=:none,ylabel="f",xlabel="I",fontfamily="Computer Modern")
#Plots.plot(reverse.(I_list_Kir_s_CaLs_d_fI[1]),reverse.(f_list_Kir_s_CaLs_d_fI[1]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)
#Plots.plot(reverse.(I_list_Kir_s_CaLs_d_fI[2]),reverse.(f_list_Kir_s_CaLs_d_fI[2]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:Set3_9)

plt = Plots.plot()
for i in 1:3
    for j in 1:3
        println("i=$i ; j=$j")
            
            try 
            plot!(plt,I_list_Kir_s_CaLs_d_fI[i][j].*10^6,(f_list_Kir_s_CaLs_d_fI[i][j]),markershape=:circle,markerstrokewidth=0.,ms=3,palette=:tab10,label="CaLs_d = $(i); Kir_s = $(j)",xlabel="I [pA]",fontfamily="Computer Modern")
            catch err 
                println("Err : i=$i , j=$j")
            end
        
    end
end
plot!(legend=:topleft)

