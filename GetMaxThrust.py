import numpy as np

def GetMaxThrust(UEFC, V):

    # You need to finish this file
    rho = UEFC.rho

    ct0 = 0.2093
    ct1 = -0.2484
    ct2 = -0.1386 #

    Tmax_static = 2               # DO NOT CHANGE, Maximum thrust desired at static conditions
    Rprop       = 0.1016          # Propeller radius (m)
    Aprop       = np.pi*Rprop**2          # Propeller disk area

    Omega  = np.sqrt(Tmax_static/(0.5*rho*Rprop**2*Aprop*ct0))

    # calculate Lambda, CT, and then Tmax from given variables
    Lambda = V/(Omega*Rprop) # Advance ratio

    CT = ct0+ct1*Lambda+ct2*Lambda**2 # Thrust coefficient

    Tmax = CT*0.5*rho*((Omega*Rprop)**2)*Aprop

    return Tmax

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    CLOSE_TOL = 1E-10
    aircraft = UEFC()
    aircraft.mpay_g = 300. # set payload mass (g)
    aircraft.tau = 0.12
    aircraft.taper = 0.7
    aircraft.dbmax = 0.10
    V = 10.0
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 0.6425093736549011, CLOSE_TOL)

    V = 12.0
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 0.30493918745386056, CLOSE_TOL)

    V = 0.0
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 2.0, CLOSE_TOL)

    V = 7.3
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 1.0632935228084999, CLOSE_TOL)
    print(f"==> All GetMaxThrust tests have passed!")


if __name__ == "__main__":
    tests()