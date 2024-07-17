from math import pi, exp, log
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
from scipy.integrate import simps
from scipy.stats import lognorm

#C²n Tenerife: 9e-16
#C²n Vigo: 1e-15
#C²n Sierra Nevada : 4e-16
#C²n Atacama: 10e-17
#C²n Hawai 9e-16

class Lognormal():
    def __init__(self) -> None:
        pass
        
    def lognormal(self, i, var_i):
        p = []
        var = log(1+var_i)
        mean = -0.5*var
        for value in i:
            if value != 0:
                y = (1/(value*pow(2*pi*var, 0.5)))*exp(-(pow(log(value)-mean,2)/(2*var)))
                p.append(y)
                
            else: 
                y = 0
                p.append(y)
                          
        return p
    
    def losses(self, pthr, sci_index):
        if type(sci_index) == np.float64:
            losses = -4.343*(erfinv(2*pthr-1)*pow(2*log(sci_index+1),0.5)-0.5*log(sci_index+1))
        else:
            losses = []
            for value in sci_index:
                losses.append(-4.343*(erfinv(2*pthr-1)*pow(2*log(value+1),0.5)-0.5*log(value+1)))
        return losses

if __name__ == "__main__":
    
    ln = Lognormal() 
    i = np.linspace(0, 3, 1000)
    vars_i_orig = np.array([0.0012056286162777739, 0.1, 0.2, 0.4])

    for j in range(1):
        vars_i = pow(vars_i_orig,j+1)
        
        plt.figure()
        for var_i in vars_i:
            
            p = ln.lognormal(i, var_i)

            i2 = np.linspace(0,1, 1000)
            p2 = ln.lognormal(i2, var_i)
            area = simps(p2, i2)
            losses = ln.losses(area, var_i)
            print("Area: " + str(area))
            print("σ²: " + str(var_i))
            print("p_max: " + str(max(p)))
            print("irradiance: " + str(i[p.index(max(p))]))
            print()
            plt.plot(i, p, label="σ²={:.4f}".format(var_i))
            #plt.plot(i2, p2, label="σ²={:.2f}".format(var_i))
           
        plt.axvline(x=1, color = "red", linestyle='--', label='Mean')    
        plt.title("Irradiance: lognormal distribution")
        plt.xlabel("Irradiance [W/m²]")
        plt.ylabel("Probability")
        plt.legend()
        print("Losses: "+str(losses)+" dB\n")
        
    # plt.figure()
    # for var_i in vars_i:
    #     areas = np.linspace(0, 1, 1000)
    #     perdidas = ln.losses(areas, var_i)
    #     plt.plot(areas, perdidas, label="σ²={:.4f}".format(var_i))
    # plt.title("Losses")
    # plt.xlabel("pthr [%]")
    # plt.ylabel("Losses [dB]")
    # plt.legend()
    # print("Losses: "+str(losses)+" dB\n")
    
    # var_i = np.linspace(0, 4, 1000)
    # pthr = 1e-6
    # pthr = 0.5
    
    # y = ln.losses(pthr, var_i)
    # plt.figure()
    # plt.plot(var_i, y, label="pthr={}".format(pthr))
    # plt.title("Losses")
    # plt.xlabel("Scintillation index")
    # plt.ylabel("Losses [dB]")
    # plt.legend()
    
    
    # print("\n\n\n")
    # sample = np.random.lognormal(1, 0.31) #create a samples with lognormal distribution
    # print("sample: " + str(sample))
    # size = 100000
    # indices = np.arange(size)

    # # Añadir etiquetas y título
    # plt.xlabel('Índice')
    # plt.ylabel('Valor Lognormal')
    # plt.title('Puntos Generados a partir de una Distribución Lognormal')
    # plt.legend()
    
    # shape, loc, scale = lognorm.fit(sample, floc=0)
    
    # p = lognorm.cdf(sample, shape, loc, scale)
    
    # print("p: " + str(p))
    
    plt.show()