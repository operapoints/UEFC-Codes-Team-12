import numpy as np

def GetCD(UEFC, opt_vars, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

    # Profile drag
    CDp = UEFC.profile_drag_coefficient(opt_vars, AR, S)

    # Induced drag coefficient
    CDi = UEFC.induced_drag_coefficient(opt_vars, AR, S)

    # Fuselage drag model
    CDfuse = UEFC.fuse_drag_coefficient(opt_vars, AR, S)

    # Payload drag coefficient increment
    CDpay = UEFC.payload_drag_coefficient(opt_vars, AR, S)

    # Total drag coefficient
    CD = CDfuse + CDp + CDi + CDpay

    return {
            "Total": CD,
            "Breakdown": {
                    "Fuselage": CDfuse,
                    "Wing":     CDp,
                    "Payload":  CDpay,
                    "Induced":  CDi,
                    }
            }

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
    cd_breakdown = GetCD(aircraft, opt_vars, AR, S)
    CLOSE_TOL = 1E-10
    assert check_close(cd_breakdown["Total"], 0.06970350362524477, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Fuselage'], 0.017777777777777778, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Wing'], 0.024955050225223097, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Payload'], 0.0, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Induced'], 0.02697067562224389, CLOSE_TOL)

    aircraft = UEFC()
    aircraft.mpay_g = 280. # set payload mass (g)
    aircraft.tau = 0.10
    aircraft.taper = 0.4
    aircraft.dbmax = 0.05
    AR = 11
    S = 0.7
    opt_vars = np.array([1.06]) # load factor
    cd_breakdown = GetCD(aircraft, opt_vars, AR, S)
#     print(cd_breakdown)
    assert check_close(cd_breakdown["Total"], 0.044494068599582604, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Fuselage'], 0.011746031746031746, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Wing'], 0.022362555651472726, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Payload'], 0.0, CLOSE_TOL)
    assert check_close(cd_breakdown["Breakdown"]['Induced'], 0.010385481202078137, CLOSE_TOL)
    print(f"==> All GetCD tests have passed!")

if __name__ == "__main__":
    tests()