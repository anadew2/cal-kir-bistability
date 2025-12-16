#-------------------------------------------------------------------------------
#                            SOLVER FUNCTIONS
#
#-------------------------------------------------------------------------------
function condition(u,t,integrator) # Event when event_f(u,t) == 0
  u[14]  
end

function affect!(integrator)
end

cb = ContinuousCallback(condition,affect!,nothing,save_positions=(true,false));