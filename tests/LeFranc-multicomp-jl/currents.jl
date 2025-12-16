
# ----------------------------------------------- GHK -----------------------------------------------
function ghk_LeFranc( V, ci, co) #v(mV), ci(mM), co(mM), z) 
    FARADAY = 96520 #default value ( https://www.neuron.yale.edu/neuron/static/docs/units/units.html)
    R       = 8.3134
    celsius = 36

    z = (1e-3)*2*FARADAY*V/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	ghk = (.001)*2*FARADAY*(eci - eco)  

    return ghk
end
function geq_ghk_LeFranc( ci, co) #v(mV), ci(mM), co(mM), z) 
    FARADAY = 96520 #default value ( https://www.neuron.yale.edu/neuron/static/docs/units/units.html)
    R       = 8.3134
    celsius = 36

    z = (1e-3)*2*FARADAY*V/(R*(celsius+273.15))
	z/(exp(z) - 1)
	eco = co*efun(z)
	eci = ci*efun(-z)
	ghk = (.001)*2*FARADAY*(eci - eco)  

    return ghk
end

# ----------------------------------------------- Currents -----------------------------------------------
function INa(gNa,mNa,hNa,V,eNa)
    INa = gNa*(mNa^3)*hNa*(V-eNa)     # [µA/cm^2]
    return INa
end
function IKDR(gKDR,mKDR,V,eK)
    IKDR = gKDR*mKDR^4*(V-eK)       # [µA/cm^2]
    return IKDR
end
function ICaLf(pCaLf,mCaLf,hCaLf,V,Ca,Ca_o)
    ICaLf = pCaLf*(mCaLf^2)*hCaLf*ghk_LeFranc(V, Ca, Ca_o) *1000  # [µA/cm^2] instead of [mA/cm^2]
    return ICaLf
end
function ICaLs(pCaLs,mCaLs,hCaLs,V,Ca,Ca_o)
    ICaLs = pCaLs*mCaLs*hCaLs*ghk_LeFranc(V, Ca, Ca_o) *1000    # [µA/cm^2] instead of [mA/cm^2]
    return ICaLs
end
function IKir(gKir,mKir,V,eK)
    IKir = gKir*mKir*(V-eK)         # [µA/cm^2]
    return IKir
end
function ICaAN(gCaAN,mCaAN,V,eCaAN)
    ICaAN = gCaAN*(mCaAN^2)*(V-eCaAN)     # [µA/cm^2]
    return ICaAN
end
function IKCa(gKCa,mKCa,V,eK)
    IKCa = gKCa*mKCa*(V-eK)     # [µA/cm^2]
    return IKCa
end
function IKM(gKM,mKM,V,eK)
    IKM = gKM*mKM*(V-eK)            # [µA/cm^2]
    return IKM
end
function Ileak(gleak,V,eleak)
    Ileak = gleak*(V-eleak)         # [µA/cm^2]
    return Ileak 
end


function IKA_Medlock(gKA,mKA,hKA,V,eK)
    IKA = gKA*mKA*hKA*(V-eK)       # [µA/cm^2]
    return IKA
end
function ICaT_Feng2019(pCaT,mCaT,hCaT,V,Ca,Ca_o)
    ICaT = pCaT*mCaT^2*hCaT*ghk_LeFranc(V, Ca, Ca_o) *1000  # [µA/cm^2] instead of [mA/cm^2]
    return ICaT
end