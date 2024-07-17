from mycalendar import MyCalendar
from plotmap import PlotMap 
from weather import Weather
from skyfield.api import load, wgs84
from datetime import datetime, timedelta
from sgp4.api import SGP4_ERRORS
from math import pi
from linkbudget import LinkBudget 
import matplotlib.pyplot as plt
import numpy as np

def eci2geodesic(eci):
    #Calculate the magnitude of the position vector (distance from the center of the Earth)
    magnitude = pow(pow(eci[0],2) + pow(eci[1],2) + pow(eci[2],2),0.5)
    lat = 180/pi * (eci[1] / magnitude)
    long = 180/pi * (eci[0] / magnitude)
    height = magnitude - Re
    
    return lat, long, height

Re = 6378.137 #km

class Prop():
    
    def __init__(self):
        self.mycalendar = MyCalendar()
        self.plotmap = PlotMap()
        self.linkbudget = LinkBudget()
        self.weather = Weather()
     
    def parse_tle(self, url, coord_gs, thr_elev):
        #Loading weather data from 2023
        self.weather.read_json()
        
        days_forecast = 2        
        coord_pos = []
        coord_neg = []
        data_pos = [] #height, alt, azimuth, distance, number_pass, time
        visibility = []
        clouds = []
        
        #TLE AQUISITION
        
        elev_min = 5
        
        satellites = load.tle_file(url)
        #satellites = load.tle_file(url, reload = True)
        by_number = {sat.model.satnum: sat for sat in satellites}
        satellite =  by_number[int(url[50:55])]
        print("Satellite Data\nEpoch: ", satellite.epoch.utc_strftime())
        
        time_start = datetime.utcnow()
        time_end = time_start + timedelta(days=days_forecast)
        ts = load.timescale()
        iter = 0
        
        #Location
        time = ts.utc(time_start.year, time_start.month, time_start.day, time_start.hour, time_start.minute, time_start.second)
        geocentric = satellite.at(time)
        lat, lon = wgs84.latlon_of(geocentric)
        height = wgs84.height_of(geocentric).km
        latitude = lat.degrees
        longitude = lon.degrees
        coord_init = (latitude, longitude)
        print("Date: {} UTC".format(time_start))
        print("Coordinates: {:.4f}º, {:.4f}º".format(latitude, longitude))
        print("Height: {:.4f} km".format(height))
        print("\nGround Stations Data")
        
        #Next passes file
        with open("Next_passes.txt", "w") as file:
            file.write(str("Epoch: "+satellite.epoch.utc_strftime()+"\n"))
        
        #Loop for gs list
        for gs in coord_gs:
            #Weather 2023
            if gs[3] == "Vigo": self.weather.graph_weather(self.weather.df_vigo_2023,(gs[0], gs[1]), "Vigo")
            elif gs[3] == "Tenerife": self.weather.graph_weather(self.weather.df_tenerife_2023,(gs[0], gs[1]), "Tenerife")
            elif gs[3] == "Sierra Nevada": self.weather.graph_weather(self.weather.df_sn_2023,(gs[0], gs[1]), "Sierra Nevada")
            elif gs[3] == "Atacama": self.weather.graph_weather(self.weather.df_atacama_2023,(gs[0], gs[1]), "Atacama")
            elif gs[3] == "Hawai": self.weather.graph_weather(self.weather.df_hawai_2023,(gs[0], gs[1]), "Hawai")
            else: pass
            
            time = ts.utc(time_start.year, time_start.month, time_start.day, time_start.hour, time_start.minute, time_start.second)
            coord_pos.append([])
            coord_neg.append([])
            data_pos.append([])
            visibility.append([])
            clouds.append([])
            bluffton = wgs84.latlon(gs[0], gs[1])
            t1 = ts.utc(time_end.year, time_end.month, time_end.day, time_end.hour, time_end.minute, time_end.second)
            
            #Visibility
            self.weather.download_json(f"https://api.open-meteo.com/v1/forecast?latitude={gs[0]}&longitude={gs[1]}&hourly=cloud_cover,visibility&forecast_days=3", f"weather_{gs[3]}.json")

            #Look Up Table Weather
            lut = {}
            list_times = self.weather.json_data["hourly"]["time"]

            for i in range(len(list_times)):
                list_times[i] = list_times[i][0:13] 
                
            list_values = []
            for i in range(len(self.weather.json_data["hourly"]["cloud_cover"])):
                list_values.append((self.weather.json_data["hourly"]["cloud_cover"][i], self.weather.json_data["hourly"]["visibility"][i]))    
    
            for times, values in zip(list_times, list_values):
                lut[times] = values
            
            # Search start and end of passes over a latlon location 
            t, events = satellite.find_events(bluffton, time, t1, altitude_degrees=elev_min)
            event_names = 'rise above {}°'.format(elev_min), 'TCA', 'set below {}°'.format(elev_min), "culminate"
            difference = satellite - bluffton
            topocentric = difference.at(time)
            alt, azimuth, distance = topocentric.altaz()
            prev_alt = alt.degrees
            print("Coordinates: {:.4f}º, {:.4f}º ({})".format(gs[0], gs[1], gs[3]))
            print("Height: {:.3f} km".format(gs[2]))
            print('Elevation: {:.4f}º'.format(alt.degrees))
            print('Azimuth: {:.4f}º'.format(azimuth.degrees))
            print('Distance: {:.4f} km\n'.format(distance.km))
            number_pass = 0
            while time.tt < ts.utc(time_start.year, time_start.month, time_start.day+days_forecast, time_start.hour, time_start.minute, time_start.second).tt:
                geocentric = satellite.at(time)
                lat, lon = wgs84.latlon_of(geocentric)
                height = wgs84.height_of(geocentric).km
                difference = satellite - bluffton
                topocentric = difference.at(time)
                alt, azimuth, distance = topocentric.altaz()
                alt = alt.degrees
                azimuth = azimuth.degrees
                distance = distance.km
                
                if alt > elev_min and prev_alt < elev_min:
                    number_pass = number_pass + 1 
                
                prev_alt = alt
                
                if alt >= elev_min:

                    coord_pos[iter].append([lat.degrees, lon.degrees, number_pass])
                    data_pos[iter].append([height, alt, azimuth, distance, number_pass, time.utc_strftime(), 0, thr_elev])
                    if str(time.utc_strftime()[0:13]) in lut:
                        visibility[iter].append(lut[str(time.utc_strftime()[0:13])][1])
                        clouds[iter].append(lut[str(time.utc_strftime()[0:13])][0])
                    else: 
                        visibility[iter].append(24)
                        clouds[iter].append(0)  
                        
                else:
                    coord_neg[iter].append([lat.degrees, lon.degrees])
    
                time = time + timedelta(seconds=10)   
            
            #check elev min
            rows = len(data_pos)-1
            data_elev = []
            delete_passes = []
            num_passes = 0
            num_pass = 1
            
            for j in range(len(data_pos[rows])):
                if data_pos[rows][j][4] <= num_pass:
                    data_elev.append(data_pos[rows][j][1])
                else: 
                    if max(data_elev) < thr_elev:
                        delete_passes.append(num_pass)
                    else:
                        num_passes = num_passes+1   
                    num_pass = num_pass+1 
                    data_elev = []   
            
            #last iteration
            if max(data_elev) < thr_elev:
                delete_passes.append(num_pass)
            else:
                num_passes = num_passes+1
            num_pass = num_pass+1 
            data_elev = [] 
            
            iter_check = 0
            for j in range(len(data_pos[rows])):
                if data_pos[rows][iter_check][4] in delete_passes:
                    del data_pos[rows][iter_check]
                    iter_check = iter_check-1
                
                else:
                    data_pos[rows][iter_check][6] = num_passes
                iter_check = iter_check+1     
                
            number_pass = 1    
            #Save next passes     
            with open("Next_passes.txt", "a") as file:
                file.write("\nCoordinates Ground Station: {:.4f}º, {:.4f}º ({})\n".format(gs[0], gs[1], gs[3]))
                for ti, event in zip(t, events):
                    name = event_names[event]
                    if name == "rise above {}°".format(elev_min):
                        file.write(f"Pass: {number_pass}\n"+str(ti.utc_strftime('%Y %b %d %H:%M:%S')+" "+name+"\n"))
                        number_pass = number_pass + 1
                    else:
                        file.write(str(ti.utc_strftime('%Y %b %d %H:%M:%S')+" "+name+"\n"))
        
            iter = iter + 1   
        return coord_init, coord_pos, coord_neg, data_pos, visibility, clouds  
          
if __name__ == "__main__":
    weather = Weather()
    mp = Prop()          
    pm = PlotMap()
    lb = LinkBudget()
    coord_tenerife = (28.299637836987383, -16.511022084597542, 2.267, "Tenerife")
    coord_vigo = (42.16996139529612, -8.68781081995837, 0.4, "Vigo")
    coord_sn = (37.06285060758391, -3.3859696887312576, 2.896, "Sierra Nevada")
    coord_atacama = (-23.07277599786815, -67.9804469733406, 5.058, "Atacama")  
    coord_hawai = (19.82605380110838, -155.4719893582595, 4.205, "Hawai")
    
    max_passes = 3
    thr_elev = 5 #º
    
    coord_gs = [coord_vigo, coord_tenerife, coord_sn, coord_atacama, coord_hawai]
    #coord_gs = [coord_tenerife]
    
    coord_init, coord_pos, coord_neg, data_pos, visibility, clouds = mp.parse_tle('https://celestrak.org/NORAD/elements/gp.php?CATNR=25994&FORMAT=tle', coord_gs, thr_elev)
    
    lb.graph_snr_sat(data_pos, coord_gs, max_passes, visibility, clouds, thr_elev)
    pm.plotmap("Terra", coord_gs, coord_init, coord_pos, coord_neg)
    
    plt.show()