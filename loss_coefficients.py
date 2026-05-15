from math import exp, log

def kc_lam(L_over_D, Re, sigma):
    # 4L/D / Re > 0.05
    if 4*(L_over_D)/Re < 0.05:
        print("Warning, 4L/D / Re below 0.05")
    x = 4*(L_over_D)/Re
    offset = 1.19*exp(-10.65*x**0.597) # some function that scales it from 4(L/d)/Re = inf
    kc = (1.08 - 0.41*sigma) - offset
    return kc

def kc_turb(Re, sigma):
    # Re > 3000
    offset = 0.14*(1-exp(-0.00136*(Re-3000)**0.622))
    kc = (0.54-0.39*sigma) - offset
    return kc

def ke_turb(Re, sigma):
    offset = 0.12/(1+0.42*(log(Re/3000))**2.05)
    ke = (1+0.85*sigma)*(1-sigma)**2.35 - sigma*offset
    return ke