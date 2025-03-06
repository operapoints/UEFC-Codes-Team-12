import numpy as np
from scipy.optimize import minimize
from GetUEFC        import UEFC
from calc_obj       import calc_obj

#opt_vars: 0:N 1:taper 2:tau
def obj_fcn(opt_vars):
    N = opt_vars[0]
    taper = opt_vars[1]
    tau = opt_vars[2]
    


def opt_obj(UEFC, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FUNCTION FOR THIS PROBLEM

    # Determine the maximum objective function (velocity)
    # achievable for an airplane with the inputted values of AR, S.

    # Optimization variables: N

    # Maximize objective: minimize negative of objective (ORIGINAL)
    #obj_fcn = lambda opt_vars: -calc_obj(UEFC, opt_vars, AR, S)

    # Modification (otherwise solution not always found)
    obj_fcn = lambda opt_vars: -calc_obj(UEFC, opt_vars, AR, S) * (1./3.)
    # optimization variable is the load factor of the turn.
    N_initialGuess = 1.25
    N_lowerBound   = 1.0001
    N_upperBound   = 5.0

    taper_initialGuess = 0.7
    taper_lowerBound = 0.4
    taper_upperBound = 1

    tau_initialGuess = 0.1
    tau_lowerBound = 0.08
    tau_upperBound = 0.12

    # R_initialGuess = 6.0
    # R_lowerBound   = 0.1
    # R_upperBound   = 12.5

    # mpay_initialGuess = 10.
    # mpay_lowerBound   = 0.01
    # mpay_upperBound   = 1000.

    initialGuess = [N_initialGuess,taper_initialGuess,tau_initialGuess]
    bounds       = ((N_lowerBound,    N_upperBound),
                    (taper_lowerBound,taper_upperBound),
                    (tau_lowerBound,tau_upperBound))

    # Constraint format is different, depending on algorithm.
    method = "SLSQP"
    if method == "SLSQP":

        # Excess thrust (max - required thrust) must be positive.
        T_constraint = {
                "type": "ineq",
                "fun":  lambda opt_vars: UEFC.excess_thrust(opt_vars, AR, S)
                }

        # Wingtip deflection must be less than the maximum allowed value.
        db_constraint = {
                "type": "ineq",
                "fun":  lambda opt_vars: UEFC.dbmax \
                - UEFC.wing_tip_deflection(opt_vars, AR, S)
                }

        # Lift coefficient must be less than the maximum allowed cruise value.
        CL_constraint = {
                "type": "ineq",
                "fun": lambda opt_vars: UEFC.CLdes \
                - UEFC.lift_coefficient(opt_vars, AR, S)
            }

    else:
        raise AttributeError("Optimization method " + method \
                             + " not recognized.")

    constraints = [T_constraint, db_constraint, CL_constraint]

    # try:
    result = minimize(fun=obj_fcn, x0=initialGuess, bounds=bounds,
                      constraints=constraints, method=method,
                      options={"maxiter": 20000})

    success = result.success
    # print(success)
    # except:  # Optimizer failed
    #     result  = None
    #     success = False

    if success:
        opt_vars_maxObj = result.x  # Variables that maximize objective
        obj_max         = calc_obj(UEFC, opt_vars_maxObj, AR, S)

    else:  # If optimizer fails
        opt_vars_maxObj = np.zeros(np.size(initialGuess))
        obj_max         = 0

    return opt_vars_maxObj, obj_max, success


if __name__ == "__main__":

    # Simple test case. Feel free to modify this part of the file.
    aircraft = UEFC()

    S  = (0.1 + 0.2)/2 * 1.5  # m^2
    AR = 1.5**2 / S
    aircraft.CLdes = 0.75
    aircraft.mpay_g = 0.0
    print(aircraft.weight(1.1, AR, S)["Total"])
    opt_vars_maxobj, obj_max, success = opt_obj(aircraft, AR, S)

    print()
    print("Aspect ratio: %0.4f"     % AR)
    print("Wing area:    %0.4f m^2" % S)

    print()
    print("Load factor:  %0.3f"   % opt_vars_maxobj[0])
    # print("Turn radius:  %0.2f m" % opt_vars_maxobj[1])
    # print("Payload mass: %0.0f g" % opt_vars_maxobj[2])

    print("Objective: %0.0f m/s" % obj_max)

