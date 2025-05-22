V_ss = -120:0.1:0 
dV_ss = 0.01

include("../../currents.jl")

ICaLs_ss = ICaLs.(1,mCaLs.(V_ss),hCaLs.(V_ss),V_ss,Ca_o,Ca_o)
Plots.plot(V_ss,ICaLs_ss)

IKir_ss = IKir.(1,mKir.(V_ss),V_ss,eK)
Plots.plot(V_ss,IKir_ss)

IKM_ss = IKM.(1,mKM.(V_ss),V_ss,eK)
Plots.plot!(V_ss,IKM_ss)


dICaLs_ss_dV = (1/(2*dV_ss)) .*(ICaLs.(1,mCaLs.(V_ss.+dV_ss),hCaLs.(V_ss.+dV_ss),V_ss.+dV_ss,Ca_o,Ca_o) -ICaLs.(1,mCaLs.(V_ss.-dV_ss),hCaLs.(V_ss.-dV_ss),V_ss.-dV_ss,Ca_o,Ca_o))
Plots.plot(V_ss,dICaLs_ss_dV)

dIKir_ss_dV = (1/(2*dV_ss)) .*(IKir.(1,mKir.(V_ss.+dV_ss),V_ss.+dV_ss,eK) -IKir.(1,mKir.(V_ss.-dV_ss),V_ss.-dV_ss,eK))
Plots.plot(V_ss,dIKir_ss_dV)

dIKM_ss_dV = (1/(2*dV_ss)) .*(IKM.(1,mKM.(V_ss.+dV_ss),V_ss.+dV_ss,eK) -IKM.(1,mKM.(V_ss.-dV_ss),V_ss.-dV_ss,eK))
Plots.plot(V_ss,dIKM_ss_dV)


include("../dI_outw/currents_outw.jl")

IKir_outw_ss = IKir.(1,mKir.(V_ss),V_ss,eK)
Plots.plot(V_ss,IKir_outw_ss)

IKM_outw_ss = IKM.(1,mKM.(V_ss),V_ss,eK)
Plots.plot!(V_ss,IKM_outw_ss)


include("../dI_inw/currents_inw.jl")

IKir_inw_ss = IKir.(1,mKir.(V_ss),V_ss,eK)
Plots.plot(V_ss,IKir_inw_ss)

IKM_inw_ss = IKM.(1,mKM.(V_ss),V_ss,eK)
Plots.plot!(V_ss,IKM_inw_ss)


include("../../currents.jl")