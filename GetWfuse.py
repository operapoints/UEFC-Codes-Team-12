import numpy as np
def GetWfuse(UEFC, AR, S):

    # You need to finish this file

    # YOU MAY NEED TO ADJUST THE CONSTANTS TO BETTER FIT YOUR ESTIMATED
    # AIRPLANE.  PROBABLY, THIS WOULD INVOLVE SETTING NEW VALUES OF
    # mfusel, mfuseS, b0, and S0

    # Calculate fuselage weight from UEFC parameters and S and AR
    mfuse0 = np.nan  # fixed mass (kg)
    mfusel = np.nan  # span (length) dependent mass (kg)
    mfuseS = np.nan  # wing area dependent mass (kg)

    SPV = np.nan   # Wing area for which mfusel and mfuseS were calculated (m^2)
    bPV = np.nan    # Wingspan for which mfusel and mfuseS were calculated (m)

    b = UEFC.wing_dimensions(AR, S)["Span"]
    g = UEFC.g

    # Calculate Wfuse from the given variables
    Wfuse = np.nan

    return Wfuse

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
    AR = 10
    S = 0.225
    Wfuse = GetWfuse(aircraft, AR, S)
    assert check_close(Wfuse, 2.65851, CLOSE_TOL)

    AR = 12
    S = 0.225
    Wfuse = GetWfuse(aircraft, AR, S)
    assert check_close(Wfuse, 2.714688994695082, CLOSE_TOL)

    AR = 12
    S = 0.7
    Wfuse = GetWfuse(aircraft, AR, S)
    assert check_close(Wfuse, 3.745653247040947, CLOSE_TOL)
    print(f"==> All GetWfuse tests have passed!")


if __name__ == "__main__":
    tests()