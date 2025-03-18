import numpy as np
import matplotlib.pyplot as plt
import math
#-----------------------------

def OpenShuntStub(Zo, Z_L, f, eps_eff):

    pm = np.array([1,-1])
    R_L = np.real(Z_L)
    X_L = np.imag(Z_L)

    num = X_L + pm * np.sqrt(R_L * (((Zo - R_L) ** 2 + X_L ** 2) / Zo))
    den = R_L - Zo

    t = num/den

    k = (2*np.pi*f*np.sqrt(eps_eff))/3e8

    dSeries = (1/k)*np.arctan(t)

    lShunt = 0


    return dSeries, lShunt







if __name__ == "__main__":

    eps_eff = 3.27
    dSeries, lShunt = OpenShuntStub(50, (75-66.314j), 2e9, eps_eff)

    print("Series length:",dSeries*1e3, "Shunt length:",lShunt )