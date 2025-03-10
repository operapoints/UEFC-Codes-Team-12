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
S = c*b
S_h = c_h*b_h

# Aspect Ratio:
AR = b**2/S
AR_h = b_h**2/S_h

# AoA slope:
a_h = (2*np.pi)/(1 + (2/AR_h))
a_w = (2*np.pi)/(1 + (2/AR))

# Tail Volume Coefficient:
V_h = (S_h * l_h)/(S*c)

## Calculate Xcg/c: SM.1/.2 ##
def CalcXcg_c_nom( CMw_nom, CLw_nom):
    CM = CMw_nom
    CL = CLw_nom

    Xcg_c = 1/4 - (CM/CL)
    return Xcg_c

## Calculate X_np/c SM.3/4##
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

## Return Values for Alpha with respect to the Alpha trim (SM.5) ##

def Alpha(alpha_e_trim,a_w,a_h,V_h,c,l_h):
    ''' Calculating alpha in terms of alpha_e_trim'''


## Calc Xcg From Weight + Locations of Components SM.6 ##
def CalcXcg(Comp_Weight_list):
    ## Units are in cm and grams
    ## Elements in Comp_Weight_list are given as location, weight

    weight_sum = 0.
    COM_track = 0.
    for e in Comp_Weight_list:
        COM_track += e[0]*e[1]
        weight_sum += e[1]

    COM = COM_track/weight_sum
    return COM

def Find_nom_payload_x(Comp_Weight_list, Xcg_nom):
    weight_sum = 0.
    COM_track = 0.
    for e in Comp_Weight_list:
        COM_track += e[0]*e[1]
        weight_sum += e[1]

    m_pay = 300
    x_pay = (Xcg_nom * (weight_sum + m_pay) - COM_track)/m_pay
    return x_pay

Xcg_nom = CalcXcg_c_nom(CMw_nom, CLw_nom)
print(f"Xcg/c = {Xcg_nom}")
print("##################")

Xnp = CalcXnp_c(c, b, c_h, b_h, l_h)
print(f"Xnp/c = {Xnp}")
print("##################")

Motor_Prop_Housing = (-16., 75.)
battery = (5.,75.)
servos = (20., 15.)
radio = (0., 8.)
landing_gear = (-13. , 12.)
fuselage_stick = (28., 36.)
pushrods_housing_wiring = (45., 24.)
wing = (3., 90.)
wing_mount = (7., 25.)
tail_surfaces = (75., 20.)

PV_Comp_Weight_list = [Motor_Prop_Housing, battery, servos, radio, landing_gear, fuselage_stick, pushrods_housing_wiring, wing, wing_mount, tail_surfaces]
print(f"Xcg from components = {CalcXcg(PV_Comp_Weight_list)} cm")
print(f"Xcg/c = {CalcXcg(PV_Comp_Weight_list)/c}")
print("##################")
x_pay_nom = Find_nom_payload_x(PV_Comp_Weight_list, Xcg_nom)
print(f"Payload location for nominal Xcg = {x_pay_nom} cm")
