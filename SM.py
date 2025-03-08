import numpy as np


## Relevant Attributes for Plane Vanilla ##

c = 0.15 # avg cord
c_h = 0.06 # avg cord of the horizontal tail
b = 1.5 # units: m
b_h = 0.4 # units: m
f_e = 0.6
l_h = 0.65 # units: m
CLw_nom = 0.65
CMw_nom = -0.15

# Planform Area:
# S = c*b
# S_h = c_h*b_h

# Aspect Ratio:
# AR = b**2/S
# AR_h = b_h**2/S_h

# AoA slope:
# a_h = (2*np.pi)/(1 + (2/AR_h))
# a_w = (2*np.pi)/(1 + (2/AR))

# Tail Volume Coefficient:
# V_h = (S_h * l_h)/(S*c)

## Calculate Xcg/c ##
def CalcXcg_c( CMw_nom, CLw_nom):
    CM = CMw_nom
    CL = CLw_nom

    Xcg_c = 1/4 - (CM/CL)
    return Xcg_c

## Calculate X_np/c ##
def CalcXnp_c(c , b , c_h, b_h , l_h):

    S = c * b
    S_h = c_h * b_h

    AR = (b**2)/S
    AR_h = (b_h**2)/S_h

    a_h = (2*np.pi)/(1 + (2/AR_h))
    a_w = (2*np.pi)/(1 + (2/AR))

    V_h = (S_h * l_h)/(S*c)
    Xnp_c = ((1/4)*(a_w/a_h) + V_h*(1 + c/(4*l_h)))/((a_w/a_h) + V_h*(c/l_h))

    return Xnp_c

Xcg = CalcXcg_c(CMw_nom, CLw_nom)
print(f"Xcg/c = {Xcg}")
print("##################")

Xnp = CalcXnp_c(c, b, c_h, b_h, l_h)
print(f"Xnp/c = {Xnp}")
