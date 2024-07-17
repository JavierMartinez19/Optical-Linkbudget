from mycalendar import MyCalendar
from math import pow,pi,log10,log,tan,exp,log2,radians,sin,cos,degrees,asin,atan
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button
from scipy.stats import norm
from scipy.special import erfinv
from scipy.integrate import quad

# Utilities

#Natural to dB
def nat2db(nat):
    return 10*log10(nat)

#dB to natural
def db2nat(db):
    return pow(10, db/10)

#Arcsec to rad
def arcsec2rad(arcsec):
    return arcsec*(2*pi)/1296000

#Rad to arcsec
def rad2arcsec(rad):
    return (rad*1296000)/(2*pi)

#Geodetic coordinates to Cartesian coordinates
def geodetic2cartesian(geodetic):
    lat, long, h = radians(geodetic[0]), radians(geodetic[1]), geodetic[2]
    #slat, long, h = geodetic[0], geodetic[1], geodetic[2]
    
    semimajor_axis = Re/1e3
    minor_semiaxis = Re/1e3
    
    ecc = pow((pow(semimajor_axis,2)-pow(minor_semiaxis,2))/pow(semimajor_axis,2),0.5)
    
    n = semimajor_axis/pow(1-(pow(ecc,2)*pow(sin(lat),2)),0.5)
    
    x = (n+h)*cos(lat)*cos(long)
    y = (n+h)*cos(lat)*sin(long)
    z = (n*(1-pow(ecc,2))+h)*sin(lat)
    
    cartesian =(x,y,z)
    
    return cartesian
    

def interpol(value, y, x):
    position = 0
    if len(x)==len(y):
        aux = abs(y[0]-value)
        for i in range(len(x)):
            if abs(y[i]-value) < aux:
                position = x[i]
            aux = abs(y[i]-value)      
    return position

#CONSTANTS

#Earth radius
Re = 6378137 #m

#Boltzmann's constant
k = 1.38e-23

#Reference temperature
To = 290 #K

#Plank's constant
h = 6.623e-34
        
#Electronic charge
q = -1.602e-19

#Speed of light
c = 299792458 #m/s

#Minimum SNR required
snr_min = 3.00

#Snow parameters
b_dry_snow = 1.38
b_wet_snow = 0.72
a_dry_snow = 5.42e-5 + 5.4958776
a_wet_snow = 1.023e-4 + 3.7855466
snow_rate = None #mm/h

#Rain parameters
rain_rate = None #mm/h

#Coordinates
coord_teleco = (42.16996139529612, -8.68781081995837) #latitude, longitude
coord_tenerife = (28.299637836987383, -16.511022084597542)

class LinkBudget(): 
    def __init__(
        
        self
        , elevation_angle = 10 #grades
        , wavelength=1550e-9 #m
        , height=10e3 #m
        , jitter=10 #aiming divergence arcseconds
        , baud=0.31e9 #bps
        , tx_power=30 #dBm
        , tx_diameter=2e-2 #m
        , tx_optic_eff=0.9 #%
        , rx_diameter=1 #m
        , rx_obscuration_rat=0.2 #%
        , rx_optic_eff=0.9 #%
        , rx_filter_loss= 0.8 
        , rx_sensitivity=-24 #dBm
        ): 
        
        self.figures = []
        self.info = []
        self.mycalendar = MyCalendar()
        # Channel
        self.wavelength = wavelength
        self.height = height
        self.elevation_angle_rad = radians(elevation_angle)
        self.baud_rate = baud
        self.jitter_rad = arcsec2rad(jitter)
        #self.free_space_losses
        #self.atmosphere_absorption_losses

        # Transmitter
        self.tx_power_dBm = tx_power
        self.tx_diameter = tx_diameter
        self.tx_area = pi*pow(self.tx_diameter/2,2)
        self.tx_optic_eff = tx_optic_eff
        #self.tx_gain_dB
        #self.tx_divergence_rad
        self.tx_spot_m = self.tx_diameter/1e2
        
        # Receiver
        self.rx_diameter = rx_diameter
        self.rx_area = pi*pow(self.rx_diameter/2,2)
        self.rx_obscuration_rat = rx_obscuration_rat
        self.rx_optic_eff = rx_optic_eff
        #self.rx_gain_dB
        #self.rx_spot_m
        self.rx_filter_loss = rx_filter_loss
        self.rx_sensitivity_dBm = rx_sensitivity
        
        #parameters
          
        #Bandwidth
        self.B = 1e9 #Hz
      
        #Detector quantum efficiency (0.5 - 0.9)
        self.eta = 0.9
        
        self.eta_min = 0.5
        
        self.eta_max = 0.9
        
        #Multiplication Factor M (10 - 30)
        self.M = 10
        
        #excess noise factor arising due to random nature of multiplication factor: Fano factor F * η ().
        self.F = 0.984*self.eta
        
        #Background Noise Power
        self.Pb = 0
        
        #Dark current
        self.Id = 100e-7 #A
        
        #Bulk dark current
        self.Idb = 0.8*self.Id
        
        #Surface dark current
        self.Ids = self.Id - self.Idb 
        
        #Load resistance
        self.Rl = 50 #ohms
        
        #Transmittance
        self.transmittance = 0.972
        
        #Visibility range (km)
        self.V = 25
        
        self.fog_thickness = 1 #km
        
        self.clouds = (0, "Clear", "Link OK")
        
        #Height of ground station
        self.hgs=2267 #m
        
        self.scintillation_losses = 0 #dB
        
    def jitter_gauss(self, num_data):
        mean = 12
        #mean = 10
        sigma = 0.1 #standard deviation
        #sigma = 0
        gaussian_arcsec = np.abs(np.random.normal(mean, sigma, num_data)) #positive values
        gaussian_rad = []
        for value in gaussian_arcsec: 
            gaussian_rad.append(arcsec2rad(value))
           
        return gaussian_rad, mean, sigma
        
    def distance(self):
        return (Re+self.hgs)*(pow((pow(Re+self.height,2)/pow((Re+self.hgs),2)) - pow(cos(self.elevation_angle_rad),2),0.5) - sin(self.elevation_angle_rad))
        
    def field(self):
        if self.distance() > ((0.8*self.wavelength)/(2*pi)) and self.distance() > ((2*pow(self.rx_diameter,2))/self.wavelength):
            field = "far-field"
        else:
            field = "near-field"    
        return field  
        
    def transmittance_db(self):
        return nat2db(self.transmittance)
            
    def free_space_loss(self):
        return nat2db(pow(((self.wavelength) / (4*pi*self.distance())), 2))
    
    def free_space_loss_sat(self, dist):
        return nat2db(pow(((self.wavelength) / (4*pi*dist)), 2))
     
    # Requires fortran
    def atmosphere_transmitance_losses(self):
        pass
        # return nat2db(atmo.transmitance(wavelength=self.wavelength, elevation_deg=10))
       
    # #wavenumber   
    def wn(self):
        return (2*pi)/self.wavelength       

    def tx_divergence_rad(self):
        return (4*self.wavelength) / (pi*self.tx_diameter)
    
    def tx_gain_dB(self): 
        return nat2db(pow(pi*self.tx_diameter/self.wavelength, 2))
    
    def rx_gain_dB(self):
        return nat2db(pow(pi*self.rx_diameter/self.wavelength, 2)\
                      *(1-pow(self.rx_obscuration_rat, 2)))

    def rx_spot_m(self):
        return tan(self.tx_divergence_rad()/2) * (self.distance()) + self.tx_diameter
    
    def rx_power_dBm(self):
        print("Gtx: ", self.tx_gain_dB())
        print("Perdidas espacio: ", self.free_space_loss())
        print("perdidas pointing", self.pointing_loss_dB())
        print("Grx: ", self.rx_gain_dB())
        return self.tx_power_dBm + nat2db(self.tx_optic_eff) + self.tx_gain_dB() \
        + self.free_space_loss() + self.pointing_loss_dB() \
        + self.rx_gain_dB() + nat2db(self.eta) + self.transmittance_db() - self.fog_attenuation() + self.scintillation_losses
        
    def rx_power_sat_dBm(self, dist):
        return self.tx_power_dBm + nat2db(self.tx_optic_eff) + self.tx_gain_dB() \
        + self.free_space_loss_sat(dist) + self.pointing_loss_dB() \
        + self.rx_gain_dB() + nat2db(self.eta) + self.transmittance_db() - self.fog_attenuation() + self.scintillation_losses  
        
    def rx_power_w(self):
        return db2nat(self.rx_power_dBm() - 30)
    
    def rx_power_sat_w(self, dist):
        return db2nat(self.rx_power_sat_dBm(dist) - 30)

    def link_margin(self):
        return self.rx_power_dBm() - self.rx_sensitivity_dBm
    
    #Noise spectral density
    def No(self):
        return k*To
    
    #Channel capacity
    def capacity(self):
        return self.baud_rate*log2(1+(self.rx_power_w()/(self.No()*self.baud_rate))) #bps
    
    def capacity_sat(self, dist,snr):
        #c = self.baud_rate*log2(1+(self.rx_power_sat_w(dist)/(self.No()*self.baud_rate))) #bps
        #print("Capacidad: ", c)
        #print("Prx: ", self.rx_power_sat_w(dist))
        
        c = self.baud_rate*log2(1+pow(10,(snr/10)))
        return c
    
    def operation_freq(self):    
        return c/self.wavelength
    
    def pointing_loss_dB(self):
        try:
            return nat2db(exp(-8*pow(self.jitter_rad, 2)/pow(self.tx_divergence_rad(), 2)))
        except ValueError:
            return -500
        
    #Turbulence attenuation
    def dist_atm(self, hgs, elev): #atm: 20 km
        return ((Re/1000)+hgs)*(pow((pow((Re/1000)+20,2)/pow(((Re/1000)+hgs),2)) - pow(cos(elev),2),0.5) - sin(elev)) 
    
    def scintillation_index(self, hgs, distance, elevation):
        hgs = hgs*1e3
        distance = distance*1e3
        Cn2, _ = quad(self.Cn2, 0, distance, args=(elevation,hgs))
        scintillation_index = (19.2/pow(self.wavelength,7/6))*Cn2
        
        return scintillation_index
    
    def Cn2(self, z, elevation, hgs):
        h = hgs + z * np.sin(elevation)
        return pow(z,5/6)*0.00594*pow(10/27,2)*pow((pow(10,-5)*(h)),10)*exp(-h/1000)+2.7e-16*exp(-h/1500)+1.7e-14*exp(-h/100)
    
    def scintillation_loss(self, scintillation_index):
        var = log(1+scintillation_index)
        mean = -var/2
        sample = exp(np.random.normal(mean, pow(var,0.5)))
        
        pthr = norm.cdf(sample)
        
        loss = -4.343*(erfinv(2*pthr-1)*pow(2*log(sample+1),0.5)-0.5*log(sample+1))
        
        if loss < -50:  loss = -50 # loss = -inf => problems
        
        return loss

    def Ro(self):
        return (self.eta*q)/(h*self.operation_freq())
    
    def snr_pin(self):
        return pow((self.Ro()*self.rx_power_w()),2)/  \
        (2*q*self.baud_rate*(self.Ro()*self.rx_power_w()+self.Ro()*self.Pb+self.Id)+(4*k*To*self.baud_rate/self.Rl))
        
    def snr_apd(self):
        return pow(self.M*self.Ro()*self.rx_power_w(),2)/ \
        ((2*q*self.baud_rate*(self.Ro()*self.rx_power_w()+self.Ro()*self.Pb+self.Idb)*pow(self.M,2)*self.F+self.Ids)+(4*k*To*self.baud_rate/self.Rl))
        
    def snr_pin_sat(self, dist):
        return pow((self.Ro()*self.rx_power_sat_w(dist)),2)/  \
        (2*q*self.baud_rate*(self.Ro()*self.rx_power_sat_w(dist)+self.Ro()*self.Pb+self.Id)+(4*k*To*self.baud_rate/self.Rl))   
        
    def snr_pin_sat_2023(self, dist, visibility):
        self.V = visibility
        return pow((self.Ro()*self.rx_power_sat_w(dist)),2)/  \
        (2*q*self.baud_rate*(self.Ro()*self.rx_power_sat_w(dist)+self.Ro()*self.Pb+self.Id)+(4*k*To*self.baud_rate/self.Rl))    
        
    def dist_gs_sat(self, coord_gs, coord_sat): 
        return (Re/1e3+coord_gs[2])*(pow((pow(Re/1e3+coord_sat[2],2)/pow((Re/1e3+coord_gs[2]),2))\
                - pow(cos(radians(self.elevation_angle_sat(coord_gs, coord_sat))),2),0.5) - sin(radians(self.elevation_angle_sat(coord_gs, coord_sat))))
    
    def elevation_angle_sat(self, coord_gs, coord_sat):
        cartesian_gs = geodetic2cartesian(coord_gs)
        cartesian_sat = geodetic2cartesian(coord_sat)
        
        #vector from gs to sat
        vector = tuple(a - b for a, b in zip(cartesian_sat, cartesian_gs))
        vector_norm = np.linalg.norm(vector)
        elevation_angle_rad = asin((cartesian_sat[2]-cartesian_gs[2])/vector_norm)
        return degrees(elevation_angle_rad)
    
    def azimuth(self, coord_gs, coord_sat):
        coord_gs = (coord_gs[0], coord_gs[1], self.hgs/1e3)
        cartesian_gs = geodetic2cartesian(coord_gs)
        cartesian_sat = geodetic2cartesian(coord_sat)
        
        azimuth = atan((cartesian_sat[1]-cartesian_gs[1])/(cartesian_sat[0]-cartesian_gs[0]))
        
        return azimuth
        
    def los(self, coord_gs, coord_sat): #LINE OF SIGHT
        elevation_angle = self.elevation_angle_sat(coord_gs, coord_sat)
        if elevation_angle >=5:
            los = True
        else:
            los = False
        return los             
                 
    # ------------------------------------------------------------------ WEATHER EFFECTS ------------------------------------------------------------  
        
    def specific_attenuation(self): #attenuation per unit length expressed in dB/km 
        return (1/(self.distance()*1e-3))*10*log10(db2nat(self.tx_power_dBm)/self.rx_power_w())
    
    # size distribution coefficient of scattering (Kim model)
    def p(self): #v: visibility range [km]
        if self.V > 50:
            value = 1.6
        elif self.V<=50 and self.V>6:
            value = 1.3
        elif self.V<=6 and self.V>1:
            value = 0.16*self.V + 0.34
        elif self.V<=1 and self.V>0.5:
            value = self.V-0.5
        else: 
            value = 0      
            
        return value
    
    def specific_attenuation_fog(self):
        return (3.91/self.V)*pow(self.wavelength*1e9/550,-self.p())
    
    def fog_attenuation(self):
        return self.specific_attenuation_fog()*self.fog_thickness
    
    def specific_attenuation_snow_visibility(self):
        return 58/self.V

    def specific_attenuation_snow_dry(self):
        return a_dry_snow*pow(snow_rate,b_dry_snow)
    
    def specific_attenuation_snow_wet(self):
        return a_wet_snow*pow(snow_rate,b_wet_snow)
    
    def specific_attenuation_rain_rate(self):
        return 1.067*pow(rain_rate,0.67)
    
    def specific_attenuation_rain_visibility(self):
        return 2.8/self.V

    def print_all(self):
        print("--- System parameters ---")
        print(f"Type of field: {self.field()}")
        print("Wavelength: {} nm".format(self.wavelength*1e9))
        print("Link distance: {:.3f} km".format(self.distance()*1e-3))
        print("Link margin @ {:.0f} Mbps: {:.2f} dB".format(self.baud_rate/1e6, self.link_margin())) 
        print()
        print("--- Transmitter parameters ---")
        print("Transmit power: {:.2f} dBm".format(self.tx_power_dBm))
        print("Optical efficiency losses: {:.2f} dB".format(nat2db(self.tx_optic_eff)))
        print("Aperture diameter: {} cm".format(self.tx_diameter*1e2))
        print("Telescope Gain: {:.2f} dB".format(self.tx_gain_dB()))
        print("Beam size on TX: {} m".format(self.tx_diameter))
        print("Divergence angle: {:.4f}º / {:.6f} rad".format(self.tx_divergence_rad()*180/pi, self.tx_divergence_rad()))
        print("Pointing jitter:  {:.4f}º / {:.6f} rad".format(self.jitter_rad*180/pi, self.jitter_rad))
        print()
        print("--- Channel parameters ---")
        print(f"Fog Thickness: {self.fog_thickness} km")
        print("Free space losses: {:.3f} dB".format(self.free_space_loss()))
       # print("Atmosphere transmitance losses: {:.3f} dB".format(self.atmosphere_transmitance_losses()))
        print("Pointing losses: {:.2f} dB".format(self.pointing_loss_dB()))
        print("Spot size on receiver: {:.2f} m".format(self.rx_spot_m()))
        print()
        print("--- Receiver parameters ---")
        print("Aperture diameter: {} cm".format(self.rx_diameter*1e2))
        print("Obscuration ratio: {:.3f} ".format(self.rx_obscuration_rat))
        print("Telescope Gain: {:.2f} dB".format(self.rx_gain_dB()))
        print("Receiver sensitivity: {} dBm".format(self.rx_sensitivity_dBm))
        print("Optical efficiency losses: {:.3f} dB".format(nat2db(self.rx_optic_eff)))
        print("Filter losses: {:.3f} dB".format(nat2db(self.rx_filter_loss)))
        print("Received power: {:.2f} dBm".format(self.rx_power_dBm()))
        print()
        print("Specific attenuation: {:.3f} dB/km".format(self.specific_attenuation()))
        print()
        #print("Eb/No: {:.2f} dB".format(nat2db(self.eb_no())))
        self.eta = self.eta_min
        print("SNR for:")
        print("η={:.2f}".format(self.eta))
        self.M = 10
        print("SNR pin: {:.2f} dB".format(nat2db(self.snr_pin())))
        print("SNR apd (M = 10): {:.2f} dB".format(nat2db(self.snr_apd())))
        self.M = 30
        print("SNR apd (M = 30): {:.2f} dB".format(nat2db(self.snr_apd())))
        self.eta = self.eta_max   
        print("η={:.2f}".format(self.eta)) 
        self.M = 10
        print("SNR pin: {:.2f} dB".format(nat2db(self.snr_pin())))
        print("SNR apd (M = 10): {:.2f} dB".format(nat2db(self.snr_apd())))
        self.M = 30
        print("SNR apd (M = 30): {:.2f} dB".format(nat2db(self.snr_apd())))
        print()
        
    # ------------------------------------------------------------------ GRAPHICS ------------------------------------------------------------ 
        
    def graph_snr_ptx(self):
        aux = self.tx_power_dBm
        x = np.linspace(self.tx_power_dBm-10,self.tx_power_dBm+20,400)
        y1 = []
        y2 = []
        for value in x:
            self.tx_power_dBm = value
            self.eta = self.eta_min
            y1.append(nat2db(self.snr_pin()))
            self.eta = self.eta_max
            y2.append(nat2db(self.snr_pin()))

        x_eta_min = np.interp(snr_min, y1, x)
        x_eta_max = np.interp(snr_min, y2, x)
        plt.plot(x,y1, label="η={:.2f}".format(self.eta_min), color="blue")
        plt.plot(x,y2, label="η={:.2f}".format(self.eta_max), color="orange")
        if x_eta_min > x[0] and x_eta_min < x[len(x)-1]:
            plt.scatter(x_eta_min, snr_min,  color="blue", marker="x", label="{:.2f}".format(x_eta_min))
        if x_eta_max > x[0] and x_eta_max < x[len(x)-1]:    
            plt.scatter(x_eta_max, snr_min,  color="orange", marker="x", label="{:.2f}".format(x_eta_max))
        plt.axhline(y=snr_min, color="red", linestyle="--")
        if self.fog_thickness == 0:
            plt.title("SNR as function of transmitted power")
        else:
            plt.title("SNR as function of transmitted power")     
        plt.xlabel("Ptx [dBm]")
        plt.ylabel("SNR [dB]")
        info = str("Distance: {:.2f} km\nJitter {} arcsecs\nBaudrate: {} Mbps".format(self.distance()/1e3, rad2arcsec(self.jitter_rad), self.baud_rate/1e6))
        plt.legend()
        ax = plt.gca()
        extra_text = plt.text(-0.2, 1.05,info, transform=ax.transAxes, fontsize=10)
        ax.add_artist(extra_text)
        plt.grid()
        self.tx_power_dBm = aux
        return 
    
    def graph_snr_jitter(self):
        aux = self.jitter_rad
        
        x = np.linspace(0,30,400)
        y1 = []
        y2 = []
        for value in x:
            self.jitter_rad = arcsec2rad(value)
            self.eta = self.eta_min
            y1.append(nat2db(self.snr_pin()))
            self.eta = self.eta_max
            y2.append(nat2db(self.snr_pin()))

        x_eta_min = interpol(snr_min, y1, x)
        x_eta_max = interpol(snr_min, y2, x)
        plt.plot(x,y1, label="η={:.2f}".format(self.eta_min), color="blue")
        plt.plot(x,y2, label="η={:.2f}".format(self.eta_max), color="orange")
        if x_eta_min > x[0] and x_eta_min < x[len(x)-1]:
            plt.scatter(x_eta_min, snr_min,  color="blue", marker="x", label="{:.2f}".format(x_eta_min))
        if x_eta_max > x[0] and x_eta_max < x[len(x)-1]:    
            plt.scatter(x_eta_max, snr_min,  color="orange", marker="x", label="{:.2f}".format(x_eta_max))
        plt.axhline(y=snr_min, color="red", linestyle="--")
        if self.fog_thickness == 0:
            plt.title("SNR as function of jitter")
        else:
            plt.title("SNR as function of jitter")     
        plt.xlabel("Jitter [arcsecs]")
        plt.ylabel("SNR [dB]")
        info = str("Distance: {:.2f} km\nTransmitted power {} dBm\nBaudrate: {} Mbps".format(self.distance()/1e3, self.tx_power_dBm, self.baud_rate/1e6))
        plt.legend()
        ax = plt.gca()
        extra_text = plt.text(-0.2, 1.05,info, transform=ax.transAxes, fontsize=10)
        ax.add_artist(extra_text)
        plt.grid()
        self.jitter_rad = aux
        return
    
    def graph_snr_elevation_angle(self):
        aux = self.elevation_angle_rad
        x = np.linspace(10,90,180)
        iteration = 0
        y1 = []
        y2 = []
        for value in x:
            
            self.elevation_angle_rad = radians(value)
            if iteration == 0:
                max_dist = self.distance()
            self.eta = self.eta_min
            y1.append(nat2db(self.snr_pin()))
            self.eta = self.eta_max
            y2.append(nat2db(self.snr_pin()))
            iteration = iteration + 1

        x_eta_min = np.interp(snr_min, y1, x)
        x_eta_max = np.interp(snr_min, y2, x)
        plt.plot(x,y1, label="η={:.2f}".format(self.eta_min), color="blue")
        plt.plot(x,y2, label="η={:.2f}".format(self.eta_max), color="orange")
        if x_eta_min > x[0] and x_eta_min < x[len(x)-1]:
            plt.scatter(x_eta_min, snr_min,  color="blue", marker="x", label="{:.2f}".format(x_eta_min))
        if x_eta_max > x[0] and x_eta_max < x[len(x)-1]:    
            plt.scatter(x_eta_max, snr_min,  color="orange", marker="x", label="{:.2f}".format(x_eta_max))
        plt.axhline(y=snr_min, color="red", linestyle="--")
        if self.fog_thickness == 0:
            plt.title("SNR as function of elevation angle")
        else:
            plt.title("SNR as function of elevation angle")  
        plt.xlabel("Elevation angle [degrees]")
        plt.ylabel("SNR [dB]")
        info = str("Distance: {:.1f} km (90º) / {:.1f} km (10º)\nTransmitted power {} dBm\nJitter {} arcsecs\nBaudrate: {} Mbps" \
            .format(self.distance()/1e3, max_dist/1e3, self.tx_power_dBm, rad2arcsec(self.jitter_rad), self.baud_rate/1e6))
        plt.legend()
        ax = plt.gca()
        extra_text = plt.text(-0.2, 1.2,info, transform=ax.transAxes, fontsize=10)
        ax.add_artist(extra_text)
        plt.grid()
        self.elevation_angle_rad = aux 
        return
    
    def graph_snr_tx_diameter(self):
        value_bot_min = False
        value_top_min = False
        value_bot_max = False
        value_top_max = False
        x = np.linspace(1e-2,3e-2,400)
        aux = self.tx_diameter
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        x_div = []
        iteration = 0
        for value in x:
            try:
                self.tx_diameter = value
                self.eta = self.eta_min
                y1.append(nat2db(self.snr_pin()))
                x1.append(value*100)
                if value_bot_min == False:
                    if y1[iteration] <= snr_min: 
                        value_bot_min = True                   
                if value_top_min == False:
                    if y1[iteration] >= snr_min:
                        value_top_min = True
                self.eta = self.eta_max
                y2.append(nat2db(self.snr_pin()))
                x2.append(value*100)
                if value_bot_max == False:
                    if y2[iteration] <= snr_min: 
                        value_bot_max = True                   
                if value_top_max == False:
                    if y2[iteration] >= snr_min:
                        value_top_max = True
                x_div.append((180*self.tx_divergence_rad())/pi)
                iteration = iteration +1
            except ValueError:
                break
        
        x_eta_min = interpol(snr_min, y1, x)
        x_eta_max = interpol(snr_min, y2, x)
        plt.plot(x1,y1, label="η={:.2f}".format(self.eta_min), color="blue")
        plt.plot(x2,y2, label="η={:.2f}".format(self.eta_max), color="orange")
        if x_eta_min > x[0] and x_eta_min < x[len(x)-1] and value_bot_min == True and value_top_min == True:
            plt.scatter(x_eta_min*100, snr_min,  color="blue", marker="x", label="{:.2f}".format(x_eta_min*100))
        if x_eta_max > x[0] and x_eta_max < x[len(x)-1] and value_bot_max == True and value_top_max == True:    
            plt.scatter(x_eta_max*100, snr_min,  color="orange", marker="x", label="{:.2f}".format(x_eta_max*100))
        plt.axhline(y=snr_min, color="red", linestyle="--")
        if self.fog_thickness == 0:
            plt.title("SNR as function of diameter in tx")
        else:
            plt.title("SNR as function of diameter in tx")  
        plt.xlabel("Tx diameter [cm]")
        plt.ylabel("SNR [dB]")
        info = str("Transmitted power {} dBm\nDistance: {:.2f} km\nJitter {} arcsecs\nBaudrate: {} Mbps" \
            .format(self.tx_power_dBm, self.distance()/1e3, rad2arcsec(self.jitter_rad), self.baud_rate/1e6))
        plt.legend()
        ax = plt.gca()
        extra_text = plt.text(-0.2, 1.3,info, transform=ax.transAxes, fontsize=10)
        ax.add_artist(extra_text)
        plt.grid()
        self.tx_diameter = aux
        return
    
    def graph_snr_sat(self, data, coord_gs, max_passes, vis, cl, thr_elev): #data-> 0: height, 1: alt, 2: azimuth, 3: distance, 4: number_pass, 5: time, 6: num_passes, 7: alt_min    
        aux = self.jitter_rad
        k = 0
        for j in range(len(coord_gs)):
            if len(data[j]) > 0:
                number_pass = data[j][0][4]
                max_passes = min(data[j][0][6], max_passes)
                print("Number of passes with elevation > {}º in {}: {}\n".format((data[j][0][7]), coord_gs[j][3], data[j][0][6]))
                i = 0
                for _ in range(max_passes):
                    distance = []
                    visibility = []
                    clouds = []
                    t = []
                    seconds = 0
                    snr_sat_1 = []
                    snr_sat_2 = []
                    date = []
                    azimuth = []
                    alt = []
                    hour = []
                    capacity = []
                    agreggate_capacity = 0
                    turbulence_atten = []
                    try:
                        number_pass = data[j][i][4]
                        while number_pass == data[j][i][4]:
                            distance.append(data[j][i][3]*1e3)
                            visibility.append(vis[j][i])
                            clouds.append(cl[j][i])
                            t.append(self.mycalendar.float2mins(seconds/60))
                            date.append(data[j][i][5])
                            hour.append(data[j][i][5][11:19])
                            azimuth.append(data[j][i][2])
                            alt.append(data[j][i][1])
                            i = i + 1
                            seconds = seconds + 10
                    except Exception as e:
                        pass  
                    
                    jitter_rad_rand, mean, sigma = self.jitter_gauss(len(distance))
                    
                    for iter in range(len(distance)):
                        self.jitter_rad = jitter_rad_rand[iter]
                        self.V = visibility[iter]
                        self.clouds = (clouds[iter], self.cloud_cover(clouds[iter])[0], self.cloud_cover(clouds[iter])[1])
                        dist_atm = self.dist_atm(coord_gs[j][2], np.radians(alt[iter]))
                        scintillation_index = self.scintillation_index(coord_gs[j][2], dist_atm, radians(alt[iter]))
                        self.scintillation_losses = self.scintillation_loss(scintillation_index)
                        turbulence_atten.append(self.scintillation_losses)
                        self.eta = self.eta_min
                        snr_sat_1.append(nat2db(self.snr_pin_sat(distance[iter])))
                        #capacity.append(self.capacity_sat(distance[iter], snr_sat_1))
                        #agreggate_capacity = agreggate_capacity + self.capacity_sat(distance[iter])
                        self.eta = self.eta_max
                        snr_sat_2.append(nat2db(self.snr_pin_sat(distance[iter])))
                        capacity.append(self.capacity_sat(distance[iter], snr_sat_2[len(snr_sat_2)-1]))
                        agreggate_capacity = agreggate_capacity + capacity[len(capacity)-1]
                    max_alt = max(alt)  
                    index_max_alt = alt.index(max_alt)
                    fig = plt.figure(figsize=(8, 6))
                    ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.8])  # [left, bottom, width, height]
                    # First plot
                    ax1.plot(t,snr_sat_1, label="η={:.2f}".format(self.eta_min), color="blue")
                    ax1.plot(t,snr_sat_2, label="η={:.2f}".format(self.eta_max), color="orange")
                    ax1.axhline(y=snr_min, color="red", linestyle="--")
                    ax1.set_title(f"Date: {date[0][0:10]}. GS: {str(coord_gs[j][0])[0:10]}º, {str(coord_gs[j][1])[0:10]}º ({coord_gs[j][3]})")
                    x_location = [0, (len(t)-1)/5, 2*(len(t)-1)/5, 3*(len(t)-1)/5, 4*(len(t)-1)/5, len(t)-1]
                    ax1.set_xticks(x_location)
                    legend1= ax1.legend(loc="upper left")
                    plt.gca().add_artist(legend1)
                    ax1.grid()
                    ax1.set_xlabel("Time", fontsize = 13)
                    ax1.set_ylabel("SNR [dB]", fontsize = 13) 
                    ax2 = ax1.twiny() 
                    x_location = [0, (len(hour)-1)/5, 2*(len(hour)-1)/5, 3*(len(hour)-1)/5, 4*(len(hour)-1)/5, len(hour)-1]
                    ax2.set_xticks(x_location)
                    line = ax2.plot(hour,snr_sat_1)[0]
                    ax2.tick_params(axis='x', labelsize=10)
                    ax2.set_xlabel("Hour (UTC)", fontsize = 13)
                    line.set_visible(False)
                    
                    # Second plot
                    ax2 = plt.axes([0.6, 0.6, 0.32, 0.32])  # [left, bottom, width, height]
                    x_location = [0, index_max_alt, len(date)-1]
                    y1_location = np.linspace(0, 360, 7)
                    ax2.set_ylim(0,360)
                    ax2.set_xticks(x_location)
                    ax2.set_yticks(y1_location)
                    ax2.plot(hour, azimuth, label="Azimuth", color="green")   
                    ax2.set_xlabel("Hour (UTC)", fontsize = 13)
                    ax2.grid()
                    ax2.tick_params(axis='y', colors='green', labelsize=12, which='both', direction='in', length=6, width=2)
                    ax2.tick_params(axis='x', labelsize=10)
                    ax2.set_ylabel("Azimuth [degrees]", color = "green", fontsize = 13) 
                    ax2 = ax2.twinx() 
                    y2_location = np.linspace(0, 90, 7)
                    ax2.set_ylim(0,90)
                    ax2.set_yticks(y2_location)
                    ax2.plot(hour, alt, label="Elevation", color="red")
                    ax2.tick_params(axis='y', colors='red', labelsize=12, which='both', direction='in', length=6, width=2)
                    ax2.set_ylabel("Elevation [degrees]", color = "red", fontsize = 13) 
                    
                    # Third plot
                    #Config polar plot
                    ax3 = plt.axes([0.6, 0.085, 0.35, 0.35], polar = True)
                    directions = ['N', 'E', 'S', 'W']
                    ang_rads = np.linspace(0, 2*np.pi, len(directions), endpoint=False)
                    ax3.set_theta_direction(-1)  # Config angle direction
                    ax3.set_theta_offset(np.pi/2)  # Config start angle
                    ax3.set_xticks(ang_rads)
                    ax3.set_xticklabels(directions)
                    rad_tags = []
                    ax3.set_yticklabels(rad_tags)
                    
                    for iter in range(len(alt)):
                        alt[iter] = radians(90 - alt[iter])
                        azimuth[iter] = radians(azimuth[iter])
                        
                    ax3.plot(azimuth, alt, color = "red")
                    ax3.text(azimuth[0], alt[0], "AOS: "+hour[0], weight='bold', color = "green", fontsize=10)
                    ax3.text(azimuth[int(index_max_alt)], alt[int(index_max_alt)], "TCA: "+ hour[int(index_max_alt)], weight='bold', color = "green", fontsize=10)
                    ax3.text(azimuth[len(azimuth)-1], alt[len(alt)-1], "LOS: " + hour[len(hour)-1], weight='bold', color = "green", fontsize=10)
                    v = str(self.V)
                    cc = str(self.clouds[0])
                    message1 = str(self.clouds[1])
                    message2 = str(self.clouds[2])
                    string = f"Visibility: "+ v +" km\n\nCloud Cover: "+ cc + "% (" + message1 + ", " + message2 + ")\n\n"
                    text_legend = string + "AOS: Acquisition Of Signal\n\nLOS: Loss Of Signal\n\nTCA: Time Of Closest Approach\n\nMean capacity: {:.4f} Gbps\n\nMax capacity: {:.4f} Gbps\n\nAgreggate capacity: {:.4f} Gb\n\nFog attenuation: {:.4f} dB\n\nJitter ~ N({}, {}) arcsecs\n\nMax scintillation attenuation: {:.4f} dB\n\nMin scintillation attenuation: {:.4f} dB".format((sum(capacity)/len(capacity))/1e9, np.max(capacity)/1e9,  agreggate_capacity/1e9,self.fog_attenuation(),mean, sigma, -1*np.min(turbulence_atten), -1*np.max(turbulence_atten))
                    # number_pass = number_pass + 1

                    # Conect button
                    self.info.append(text_legend)
                    button_ax = fig.add_axes([0.55, 0.1, 0.1, 0.075])
                    button = Button(button_ax, 'Info')
                    button.on_clicked(lambda _, text=self.info[k]: self.show_info(text)) 
                    button.label.set_fontsize(12)            
                    self.figures.append((fig, button))
                    k = k + 1
            else:
                print("There are no passes above {}º in {}\n".format(thr_elev, coord_gs[j][3]))   
        self.jitter_rad = aux 
        return
    
    def cloud_cover(self, value):
        data = None
        message = None
        if value <= 12.5: 
            data = "Clear" 
            message = "link OK"
        elif value <= 37.5: 
            data = "Mostly Clear" 
            message = "probably link OK"
        elif value <= 62.5: 
            data = "Partly Cloudy"
            message = "link may fail" 
        elif value <= 87.5: 
            data = "Mostly Cloudy"
            message = "probably no link"
        else: 
            data = "Cloudy" 
            message = "no link"
        return data, message
    
    #Functions for button
    def show_info(self, text):
        plt.figure()
        plt.text(0.5, 0.5, text, ha='center', va='center', fontsize=12)
        plt.axis('off')
        plt.show()
        return
    
    def graphics_ptx_jitter(self):
        plt.figure()
        plt.subplot(1,2,1)
        self.graph_snr_ptx()
        plt.subplot(1,2,2)
        self.graph_snr_jitter()
        return
    
    def graphics_leo(self):
        plt.figure()
        plt.subplot(2,2,1)
        self.graph_snr_ptx()
        plt.subplot(2,2,2)
        self.graph_snr_jitter()
        plt.subplot(2,2,3)
        self.graph_snr_elevation_angle()
        plt.subplot(2,2,4)
        self.graph_snr_tx_diameter()
        plt.tight_layout()
        return
        
    def graph_visibility(self):
        fig, ax = plt.subplots()
        aux = self.V
        x = np.linspace(0,0.6,100) #km
        y_fog = []
        y_snow = []
        y_rain = []
        for value in x:
            if value != 0:
                self.V = value
                y_fog.append(self.specific_attenuation_fog())
                y_snow.append(self.specific_attenuation_snow_visibility())
                y_rain.append(self.specific_attenuation_rain_visibility())
            else: 
                y_fog.append(None)
                y_snow.append(None)
                y_rain.append(None)
                
        ax.plot(x,y_fog, label="Fog", color="green")
        ax.plot(x,y_snow, label="Snow", color="orange")
        ax.plot(x,y_rain, label="Rain", color="blue")   
        ax.set_ylim(0,400)     
        plt.title(f"Specific attenuation due to weather conditions (λ = {1e9*self.wavelength} nm)")
        plt.xlabel("Visibility range [km]")
        plt.ylabel("Specific attenuation [dB/km]")
        ax.legend()
        ax.grid()
        self.V = aux
        return    
    
    def table_fog(self):
        fig, ax = plt.subplots()
        aux_thickness = self.fog_thickness
        aux_visibility = self.V
        rows = 10
        cols = 10
        x = np.linspace(0.5,9.5,cols) #visibility
        y = np.linspace(0.1,1,rows) #fog thickness
        matrix_np = np.zeros((rows,cols), dtype = float)
        for i in range(rows):
            self.V = x[i]
            for j in range(cols):
                self.fog_thickness = y[j]                 
                matrix_np[rows-1-j][i] = "{:.2f}".format(self.fog_attenuation()) 
        tabla = ax.table(cellText=matrix_np, loc="center", cellLoc="center", colWidths=[0.1]*cols)
        tabla.auto_set_font_size(False)
        tabla.set_fontsize(11)    
        tabla.scale(1, 1.5)
        for i in range(rows):
            for j in range(cols):
                tabla[(i,j)].set_height(1/rows)
                
        ax.set_xticks(x)   
        ax.set_yticks(y-0.05) 
        ax.set_xticklabels(x)  
        ax.set_yticklabels(np.round(y,1))     
        ax.set_title("Losses [dB]")       
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")
        ax.set_xlabel("Visibility [km]")
        ax.set_ylabel("Fog thickness [km]")      
        self.fog_thickness = aux_thickness
        self.V = aux_visibility
        return

if __name__ == "__main__":
    # print("------ GROUND LINK: ------\n")
    # lb_ground = LinkBudget( wavelength=1550e-9 #m
    #     , elevation_angle = 90 #grades                   
    #     , height=10e3 #m
    #     , jitter=10 #aiming divergence arcseconds
    #     , baud=0.31e9 #bps
    #     , tx_power = 5 #dBm
    #     , tx_diameter=2e-2 #m
    #     , tx_optic_eff=0.9 #%
    #     , rx_diameter=20e-2 #m
    #     , rx_obscuration_rat=0.2 #%
    #     , rx_optic_eff=0.9 #%
    #     , rx_filter_loss= 0.8 
    #     , rx_sensitivity=-24)

    # lb_ground.print_all()
    # lb_ground.graphics_ptx_jitter()

    print("------ LEO ORBIT: ------\n")
    lb_leo = LinkBudget( wavelength=1550e-9 #m
        , elevation_angle = 30 #grades                
        , height=700*1e3 #m
        , jitter=10 #aiming divergence arcseconds
        , baud=0.31e9 #bps
        , tx_power = 20 #dBm
        , tx_diameter=1.5e-2 #m
        , tx_optic_eff=0.9 #%
        , rx_diameter=1 #m
        , rx_obscuration_rat=0.2 #%
        , rx_optic_eff=0.9 #%
        , rx_filter_loss= 0.8 
        , rx_sensitivity=-24)
    lb_leo.print_all()
    lb_leo.graphics_leo()
    lb_leo.graph_visibility()
    
    #CASE 1 KM FOG IN EARTH SURFACE
    lb_leo.fog_thickness = 1 #km
    lb_leo.print_all()
    lb_leo.graphics_leo()
    lb_leo.table_fog()
   
    coord_gs = (coord_tenerife[0], coord_tenerife[1], lb_leo.hgs/1e3)
    coord_gs_cartesian = geodetic2cartesian(coord_gs)

    lb_leo.fog_thickness = 0 #km
     
    print("Geodetic coordinates gs:")
    print("Latitude: {}º".format(coord_gs[0]))
    print("Longitude: {}º".format(coord_gs[1]))
    print("Height: {} km".format(coord_gs[2]))
    print("Cartesian coordinates gs:")
    print("x: {}".format(coord_gs_cartesian[0]))
    print("y: {}".format(coord_gs_cartesian[1]))
    print("z: {}".format(coord_gs_cartesian[2]))
    
    
    plt.show()      
    
    