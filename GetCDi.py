import numpy as np

def GetCDi(UEFC, opt_vars, AR, S):

    # You need to finish this file

    # Calculate induced drag coefficient
    CL  = UEFC.lift_coefficient(opt_vars, AR, S)
    e   = UEFC.span_efficiency(opt_vars, AR, S)

    # calculate CDi from the given variables
    CDi = (CL**2)/(np.pi*e*AR)

    return CDi

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    aircraft = UEFC()
    aircraft.mpay_g = 300. # set payload mass (g)
    aircraft.tau = 0.12
    aircraft.taper = 0.7
    aircraft.dbmax = 0.10
    AR = 10
    S = 0.225
    opt_vars = np.array([1.1]) # load factor
    CDi = GetCDi(aircraft, opt_vars, AR, S)
    CLOSE_TOL = 1E-10
    assert check_close(CDi, 0.02697067562224389, CLOSE_TOL)

    aircraft = UEFC()
    aircraft.mpay_g = 280. # set payload mass (g)
    aircraft.tau = 0.10
    aircraft.taper = 0.4
    aircraft.dbmax = 0.05
    AR = 11
    S = 0.7
    opt_vars = np.array([1.06]) # load factor
    CDi = GetCDi(aircraft, opt_vars, AR, S)
    assert check_close(CDi, 0.010385481202078137, CLOSE_TOL)
    print(f"==> All GetCDi tests have passed!")

if __name__ == "__main__":
    tests()
