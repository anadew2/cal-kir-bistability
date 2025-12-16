#-------------------------------------------------------------------------------
#                                 CONSTANTS
#
#-----------------------------------------------------------------------------

function mNa(v)
	vtraub = -50 
	v2 = v - vtraub
	a = 0.32 * vtrap_lf(13-v2, 4)
	b = 0.28 * vtrap_lf(v2-40, 5)
	m_inf = a / (a + b)
	return m_inf
end
function tau_mNa(v)
	celsius =36
	tadj = 3.0 ^ ((celsius-36)/ 10 )
	vtraub =  -50 
	v2 = v - vtraub
	a = 0.32 * vtrap_lf(13-v2, 4)
	b = 0.28 * vtrap_lf(v2-40, 5)
	tau_m = 1 / (a + b) / tadj

	return tau_m
end
function hNa(v)
	vtraub =  -50 
	v2 = v - vtraub
	a = 0.128 * exp((17-v2)/18)
	b = 4 / ( 1 + exp((40-v2)/5) )
	h_inf = a / (a + b)

	return h_inf
end
function tau_hNa(v)
	celsius =36
	tadj = 3.0 ^ ((celsius-36)/ 10 )
	vtraub =  -50 
	v2 = v - vtraub
	a = 0.128 * exp((17-v2)/18)
	b = 4 / ( 1 + exp((40-v2)/5) )
	tau_h = 1 / (a + b) / tadj
	
	return tau_h
end
function mKDR(v)
	vtraub =  -50 
	v2 = v - vtraub

	a = 0.032 * vtrap_lf(15-v2, 5)
	b = 0.5 * exp((10-v2)/40)
	n_inf = a / (a + b)

	return n_inf
end
function tau_mKDR(v)
	celsius =36
	tadj = 3.0 ^ ((celsius-36)/ 10 )
	vtraub =  -50 
	v2 = v - vtraub

	a = 0.032 * vtrap_lf(15-v2, 5)
	b = 0.5 * exp((10-v2)/40)
	tau_n = 1 / (a + b) / tadj

	return tau_n
end
function vtrap_lf(x,y)
	if (abs(x/y) < 1e-6)
		vtrap = y*(1 - x/y/2)
	else
		vtrap = x/(exp(x/y)-1)
	end
end

function mCaLf(v)
	## Based on LeFrancLeMasson2010 --> m_inf = -0.0012 + (1.0029/(1+exp((-(v+14.3907))/3.1029  #paper values without the typo (V1/2~-14, not the opposite)
	m_inf = 1/(1+exp(-(v+17.5)/4.3))
	return m_inf
end
function tau_mCaLf(v)
    taufactor=0.5
	v = v+65
	a = 1*efun(.1*(25-v))
	b = 4*exp(-v/18)
	tau_m = taufactor/(a + b)
    return tau_m
end
function hCaLf(v)
    zshift=0 
	zpente=0
    h_inf =   1 / (1+exp((v+(14+zshift))/(zpente+4.03)))
    return h_inf
end
function tau_hCaLf(v)
    tau_h =  1500
    return tau_h
end

function mCaLs(v)
	## Based on LeFrancLeMasson2010 --> m_inf =  -0.0048+(1.0257/(1+exp(-(v-(-20.4565))/4))^0.4731) #paper values
	m_inf = 1/(1+exp(-(v+25.4)/6))
	return m_inf
end
function tau_mCaLs(v)
    taufactor = 160
	v = v+65
	a = 1*efun(.1*(25-v))
	b = 4*exp(-v/18)
	tau_m = (taufactor/(a + b))
    return tau_m
end
function hCaLs(v)
    zshift= 0 #0	
	zpente=0
    h_inf =   1 / (1+exp((v+(14+zshift))/(zpente+4.03)))
    return h_inf
end
function tau_hCaLs(v)
    tau_h =  10000
    return tau_h
end

function efun(z)
	if (abs(z) < 1e-4)
		efun = 1 - z/2
	else
		efun = z/(exp(z) - 1)
	end
end

#-----------------------------------------------------------------------------------------------------
function mKir(v)
	vhalf=-65 #(mV)
	zslope=10 #(mV)
	m_inf = 1/(1+ exp((v-vhalf)/zslope))
	return m_inf
end
function tau_mKir(v)
	tau_kir=1 #(ms)
	tau_m = tau_kir
    return tau_m
end



#-----------------------------------------------------------------------------------------------------
# ican2.mod
	function mCAN(ca)
		cac	= 0.9
		beta = 0.002
		alpha2 = beta * (ca/cac)^2

		m_inf = alpha2 / (alpha2 + beta)
		return m_inf
	end
	function tau_mCAN(ca)
		tau_m = 1000
		return tau_m
	end

#-----------------------------------------------------------------------------------------------------
# sk.mod
	function mKCa(ca)
		k= 0.5 #[mM]

		m_inf = ca/(ca+k)	

		return m_inf
	end
	function tau_mKCa(ca)
		tau_m = 10
		return tau_m
	end

#-------------------------------------------------------------------------------
#                                  KM GATE
#
#-------------------------------------------------------------------------------
function mKM(v)
	## Inspired from Kv7.2 Micelli #mK72inf_Micelli
    m_inf = 1 ./ (1 .+ exp.((-76.1-v)./25.7))
    return m_inf
end
function tau_mKM(v)
	## Inspired from Kv7.2 Micelli #taumK72_Micelli
    tau = 103 #wiegthed average
    return tau
end