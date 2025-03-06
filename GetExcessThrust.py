def GetExcessThrust(UEFC, opt_vars, AR, S):
    
    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

    # Get excess thrust: maximum thrust - required thrust
    N = opt_vars[0]
    taper = opt_vars[1]
    tau = opt_vars[2]
    UEFC.taper = taper
    UEFC.tau = tau
    
    V    = UEFC.flight_velocity(opt_vars, AR, S)
    Tmax = UEFC.maximum_thrust(V)
    Treq = UEFC.required_thrust(opt_vars, AR, S)
    
    return Tmax - Treq
