## Relevant Attributes for Plane Vanilla ##

c = 0.15 # cord
c_h = 0.06 # cord of the horizontal tail
b = 1.5 # units: m
b_h = 0,4 # units: m
f_e = 0.6
l_h = 0.65 # units: m
CLw_nom = 0.65
CMw_nom = -0.15

S = c*b
S_h = c_h*b_h
## Calculate Xcg/c ##

Xcg_c = 1/4 - (CMw_nom/CLw_nom)
print(f"Xcg / avg c = {Xcg_c}")

## Calculate X_np/c ##
Xnp_c = NotImplemented



