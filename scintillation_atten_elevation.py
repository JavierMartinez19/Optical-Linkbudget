import numpy as np
from math import exp, log, sin, pi, degrees
from scipy.stats import norm
from scipy.special import erfinv
from scipy.integrate import quad
import matplotlib.pyplot as plt

class Scint_atten():
    def __init__(self) -> None:
        pass
    
    def Cn2(self, z, elevation, hgs):
        h = hgs + z * np.sin(elevation)
        return pow(z,5/6)*0.00594*pow(10/27,2)*pow((pow(10,-5)*(h)),10)*exp(-h/1000)+2.7e-16*exp(-h/1500)+1.7e-14*exp(-h/100)
    
    def scintillation_index(self, wavelength, hgs, distance, elevation):
        
        Cn2, _ = quad(self.Cn2, 0, distance, args=(elevation,hgs))
        scintillation_index = (19.2/pow(wavelength,7/6))*Cn2
        
        return scintillation_index
    
    def scintillation_loss(self, scintillation_index):
        var = log(1+scintillation_index)
        mean = -var/2
        sample = exp(np.random.normal(mean, pow(var,0.5)))
        
        pthr = norm.cdf(sample)
        
        loss = -4.343*(erfinv(2*pthr-1)*pow(2*log(sample+1),0.5)-0.5*log(sample+1))
        
        if loss < -50:  loss = -50 # loss = -inf => problems
        
        return loss
        
        
if __name__ == "__main__":
    
    sa = Scint_atten()
    
    wavelength = 1550e-9 #m
    scintillation_index = 0.1
    elevation = np.linspace(0, pi/2, 1000)
    hgs = 0 #m
    h_atm = 20000 #m
    
    scintillation_loss = []
    scintillation_index = []
    elevation_degrees = []
    for alt in elevation:
        if np.sin(alt) != 0:  # Evitar división por cero cuando theta es 0 grados
            elevation_degrees.append(degrees(alt))
            distance = (h_atm-hgs)/sin(alt)
            scintillation_index.append(sa.scintillation_index(wavelength, hgs, distance, alt))
            #scintillation_loss.append(sa.scintillation_loss(scintillation_index))
        else: pass
    
    plt.figure()
    plt.plot(elevation_degrees, scintillation_index, label="scintillation index")
    plt.yscale('log')
    plt.xlim(0, 90)
    plt.title("Scintillation Index")
    plt.xlabel("Elevation [º]")
    plt.ylabel("Scintillation index")
    
    plt.show()
    
    