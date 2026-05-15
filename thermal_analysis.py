from parameters import *
import numpy as np

# shell side
def handoutThermalCoefficient(ReSh,ReTu,shape='triangle'):

    if shape == 'triangle':
        c = 0.2
    elif shape == 'square':
        c= 0.15
    else:
        c= 0.2

    Nui = 0.023 * ReTu**0.8 * Pr**0.3
    Nuo = c * ReSh**0.6 * Pr**0.3

    hi = Nui * k_w / di
    ho = Nuo * k_w / do

    Hinv = 1/hi + (di * np.log(do/di))/(2*k_tube) + 1/ho * di/do
    H = 1/Hinv

    return H

def tempIteratorLMTD(H,A,passes,mdot1,mdot2,T1in=Tcold_in,T2in=Thot_in):
    # guess initial
    T1out = 25
    T2out = T2in - mdot1/mdot2 * (T1out - T1in)

    def LMTD(T1in=T1in,T1out=T1out,T2in=T2in,T2out=T2out):
        dT1 = T2in - T1out
        dT2 = T2out - T1in
        lmtd = (dT1 - dT2) / np.log(dT1/dT2)
        return lmtd

    def passCorrection(T1in=T1in,T1out=T1out,T2in=T2in,T2out=T2out,passes=1):
        P1 = (T1out-T1in)/(T2in-T1in)
        R = (T2in-T2out)/(T1out-T1in)

        Pdash = ((1-P1*R)/(1-P1))**passes
        P = (Pdash-1)/(Pdash-R)

        if R == 1:
            numerator = np.sqrt(2) * P
            denominator = (1 - P) * np.log((2 - P * (2 - np.sqrt(2))) / (2 - P * (2 + np.sqrt(2))))
        else:
            numerator = np.sqrt(R**2 + 1) * np.log((1 - P) / (1 - P * R))
            denominator = (R - 1) * np.log((2 - P * (R + 1 - np.sqrt(R**2 + 1))) /(2 - P * (R + 1 + np.sqrt(R**2 + 1))))

        F = numerator / denominator

        return F

    Qe = mdot1 * c_p * (T1out - T1in)
    Qlmtd = H * A * passCorrection(passes=passes) * LMTD()

    while abs(Qe-Qlmtd)/Qe > 10e-2:
        T1out = T1in + Qlmtd/(mdot1*c_p)
        T2out = T2in - mdot1/mdot2 * (T1out - T1in)

        Qe = mdot1 * c_p * (T1out - T1in)
        Qlmtd = H * A * passCorrection() * LMTD()

    return T1out,T2out,Qe

H = 3920
N = 13
passes = 1
A = N * 0.35 * np.pi * di
mdot1 = 0.48
mdot2 = 0.33

print(tempIteratorLMTD(H,A,passes,mdot1,mdot2))

ReTu = 11000
ReSh = 15000

print(handoutThermalCoefficient(ReSh,ReTu))