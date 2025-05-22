include("gates_CaLf_replaced_by_CaLs.jl")

function INa(gNa,mNa,hNa,V,eNa)
    INa = gNa*mNa^3*hNa*(V-eNa)     # [µA/cm^2]
    return INa
end
function IKDR(gKDR,mKDR,V,eK)
    IKDR = gKDR*mKDR^4*(V-eK)       # [µA/cm^2]
    return IKDR
end
function ICaLf(gCaLf,mCaLf,hCaLf,V,Ca,Ca_o)
    ICaLf = gCaLf*mCaLf^2*hCaLf*(V-eCa) 
    return ICaLf
end
function ICaLs(gCaLs,mCaLs,hCaLs,V,Ca,Ca_o)
    ICaLs = gCaLs*mCaLs*hCaLs*(V-eCa)
    return ICaLs
end
function IKir(gKir,mKir,V,eK)
    IKir = gKir*mKir*(V-eK)         # [µA/cm^2]
    return IKir
end
function IKM(gKM,mKM,V,eK)
    IKM = gKM*mKM*(V-eK)            # [µA/cm^2]
    return IKM
end
function Ileak(gleak,V,eleak)
    Ileak = gleak*(V-eleak)         # [µA/cm^2]
    return Ileak 
end
