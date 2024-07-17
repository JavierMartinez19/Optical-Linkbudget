from skyfield.api import load, wgs84, EarthSatellite
from datetime import datetime, timedelta
import pandas as pd
from mycalendar import MyCalendar
from weather import Weather
from linkbudget import LinkBudget
from math import log10

def nat2db(nat):
    return 10*log10(nat)

class TLE_2023():
    
    def __init__(self) -> None:
        self.mycalendar = MyCalendar()
        
        self.df_vigo = pd.DataFrame(columns=["Rise above 10º", "TCA", "Set below 10º", "Link duration", "Cloud cover", "Visibility", "SNR rise (dB)", "SNR TCA (dB)", "SNR set (dB)", "Cloud info"])
        self.df_tenerife = pd.DataFrame(columns=["Rise above 10º", "TCA", "Set below 10º", "Link duration","Cloud cover", "Visibility", "SNR rise (dB)", "SNR TCA (dB)", "SNR set (dB)", "Cloud info"])
        self.df_sn = pd.DataFrame(columns=["Rise above 10º", "TCA", "Set below 10º", "Link duration","Cloud cover", "Visibility", "SNR rise (dB)", "SNR TCA (dB)", "SNR set (dB)", "Cloud info"])
        self.df_atacama = pd.DataFrame(columns=["Rise above 10º", "TCA", "Set below 10º", "Link duration","Cloud cover", "Visibility", "SNR rise (dB)", "SNR TCA (dB)", "SNR set (dB)", "Cloud info"])
        self.df_hawai = pd.DataFrame(columns=["Rise above 10º", "TCA", "Set below 10º", "Link duration", "Cloud cover", "Visibility", "SNR rise (dB)", "SNR TCA (dB)", "SNR set (dB)", "Cloud info"])
        self.weather = Weather()
        self.lb = LinkBudget()
        
    def epoch_tle(self, tle_line_1):
        epoch = float(tle_line_1[18:32])
        year, month, day, hours, minutes, seconds = self.mycalendar.frac2date(epoch)
        epoch_str = str(year).zfill(2) + "-" + str(month).zfill(2) + "-" + str(day).zfill(2) + " " + \
            str(hours).zfill(2) + ":" + str(minutes).zfill(2) + ":" + str(seconds).zfill(2) +" UTC"
        return year, month, day, hours, minutes, seconds, epoch_str
    
    def save_dataframe(self, tles):
        #Clear dataframes
        self.df_tles = pd.DataFrame(columns=["Epoch", "TLE"])
        self.df_vigo.drop(self.df_vigo.index, inplace=True)
        self.df_tenerife.drop(self.df_tenerife.index, inplace=True)
        self.df_sn.drop(self.df_sn.index, inplace=True)
        self.df_atacama.drop(self.df_atacama.index, inplace=True)
        self.df_hawai.drop(self.df_hawai.index, inplace=True)
        
        #Read TLEs from file
        with open(tles, 'r') as f:
            tle_lines = f.readlines()
            
        #Name of the saved file
        saved_file = tles[0:len(tles)-4]+"_epoch.txt"
        
        #Process each pair of TLE lines
        with open(saved_file, "w") as f:
            for i in range(0, len(tle_lines), 2):
                line1 = tle_lines[i].strip()
                line2 = tle_lines[i+1].strip()
                epoch_str = self.epoch_tle(line1)[6]
                prev_epoch = None
                
                if i > 0:
                    prev_epoch = self.df_tles["Epoch"][len(self.df_tles)-1]
                    
                if i == 0 or prev_epoch != epoch_str:    
                    
                    f.write(epoch_str+"\n"+line1+"\n"+line2+"\n")
                
                    #Add new TLEs to the dataframe
                    new_data = pd.DataFrame({'Epoch': [epoch_str], 'TLE': [line1+"\n"+line2]})
            
                    self.df_tles = pd.concat([self.df_tles, new_data], ignore_index=True) 
                
        self.df_tles.to_csv(tles[0:len(tles)-4]+"_dataframe.csv", index=False)        
          
        return  
    
    def passes_2023(self, name_sat, coord_gs):
        #Read weather json
        self.weather.read_json()
        
        #Convert date cols in datetime format
        self.weather.df_vigo_2023["datetime"] = pd.to_datetime(self.weather.df_vigo_2023["datetime"])
        self.weather.df_tenerife_2023["datetime"] = pd.to_datetime(self.weather.df_tenerife_2023["datetime"])
        self.weather.df_sn_2023["datetime"] = pd.to_datetime(self.weather.df_sn_2023["datetime"])
        self.weather.df_atacama_2023["datetime"] = pd.to_datetime(self.weather.df_atacama_2023["datetime"])
        self.weather.df_hawai_2023["datetime"] = pd.to_datetime(self.weather.df_hawai_2023["datetime"])
        
        print("Passes "+name_sat)
        event_names = 'rise above 10º', 'TCA', 'set below 10º'
        #One iteration for each gs
        for gs in coord_gs:
            print("Iteration: "+gs[3])
            aos = []
            tca = []
            los = []
            link_duration = []
            cloud_cover = []
            visibility = []
            snr_rise = []
            snr_tca = []
            snr_set = []
            agg_capacity = []
            cloud_info = []
            counter = 3 
            for i in range(0,len(self.df_tles)):
                #Read lines in DataFrame
                line1 = self.df_tles["TLE"][i][0:69]
                line2 = self.df_tles["TLE"][i][70:len(self.df_tles["TLE"][i])]
            
                #Read this and next epoch
                this_epoch =  self.df_tles["Epoch"][i]
                this_epoch = datetime(int(this_epoch[0:4]), int(this_epoch[5:7]), int(this_epoch[8:10]), int(this_epoch[11:13]), int(this_epoch[14:16]),int(this_epoch[17:19]))
                next_epoch = None
                if i < len(self.df_tles)-1: #enters if it is not the last iteration
                    next_epoch = self.df_tles["Epoch"][i+1]
                    next_epoch = datetime(int(next_epoch[0:4]), int(next_epoch[5:7]), int(next_epoch[8:10]), int(next_epoch[11:13]), int(next_epoch[14:16]),int(next_epoch[17:19]))
                    
                else: #if this is the last iteration, next_epoch = 31/12/2023 23:59:59
                    next_epoch = datetime(2023, 12, 31, 23, 59, 59) 
        
                #TLE acquisition
                ts = load.timescale()
                satellite = EarthSatellite(line1, line2, name_sat, ts)
                time_start = ts.utc(this_epoch.year, this_epoch.month, this_epoch.day, this_epoch.hour, this_epoch.minute, this_epoch.second)
                time_end = ts.utc(next_epoch.year, next_epoch.month, next_epoch.day, next_epoch.hour, next_epoch.minute, next_epoch.second)            

                #Search next passes
                bluffton = wgs84.latlon(gs[0], gs[1])
                t, events = satellite.find_events(bluffton, time_start, time_end, altitude_degrees=10.0)

                for ti, event in zip(t, events):
                    #save date in datetime format
                    name = event_names[event]
                    date_utc_str = ti.utc_strftime('%Y %b %d %H:%M:%S')
                    date_datetime = datetime.strptime(date_utc_str, "%Y %b %d %H:%M:%S")
                    date_datetime_partial = datetime(date_datetime.year, date_datetime.month, date_datetime.day)
                    
                    #Save events in arrays
                    if name == "rise above 10º":
                        aos.append(date_datetime)
                        #Find row in the weather dataframe
                        if gs[3] == "Vigo":
                            diff = abs(self.weather.df_vigo_2023['datetime'] - date_datetime_partial)
                            #print(diff)
                            index = diff.idxmin()
                            #print("Indice: "+str(index))
                            try:
                                clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                if vis == 0: vis = 0.1
                                #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                            except:
                                if date_datetime.hour < 23:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                                else:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
        
                            cloud_cover.append(clouds)
                            visibility.append(vis)

                        elif gs[3] == "Tenerife":
                            diff = abs(self.weather.df_tenerife_2023['datetime'] - date_datetime_partial)
                            index = diff.idxmin()
                            try:
                                clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                if vis == 0: vis = 0.1
                                #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                            except:
                                if date_datetime.hour < 23:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                                else:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                    
                            cloud_cover.append(clouds)
                            visibility.append(vis)
                            
                        elif gs[3] == "Sierra Nevada":
                            diff = abs(self.weather.df_sn_2023['datetime'] - date_datetime_partial)
                            index = diff.idxmin()
                            try:
                                clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                if vis == 0: vis = 0.1
                                #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                            except:
                                if date_datetime.hour < 23:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                                else:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                    
                            cloud_cover.append(clouds)
                            visibility.append(vis)
                            
                        elif gs[3] == "Atacama":
                            diff = abs(self.weather.df_atacama_2023['datetime'] - date_datetime_partial)
                            index = diff.idxmin()
                            try:
                                clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                if vis == 0: vis = 0.1
                                #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                            except:
                                if date_datetime.hour < 23:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                                else:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                    
                            cloud_cover.append(clouds)
                            visibility.append(vis)
                        else:
                            diff = abs(self.weather.df_hawai_2023['datetime'] - date_datetime_partial)
                            index = diff.idxmin()
                            try:
                                clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                if vis == 0: vis = 0.1
                                #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                            except:
                                if date_datetime.hour < 23:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                            
                                else:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                    
                            cloud_cover.append(clouds)
                            visibility.append(vis)
                        snr = self.snr(date_datetime, satellite, bluffton, cloud_cover[len(cloud_cover)-1], visibility[len(visibility)-1])[0]
                        snr_rise.append(snr)    
                        
                    elif name == "TCA":
                        tca.append(date_datetime)
                        snr, data_cc, message_cc = self.snr(date_datetime, satellite, bluffton, cloud_cover[len(cloud_cover)-1], visibility[len(visibility)-1])
                        snr_tca.append(snr)  
                        cloud_info.append(data_cc + ", " + message_cc)
                        
                    elif name == "set below 10º":
                        los.append(date_datetime)
                        snr = self.snr(date_datetime, satellite, bluffton, cloud_cover[len(cloud_cover)-1], visibility[len(visibility)-1])[0]
                        snr_set.append(snr)  
                        link_duration.append(str((date_datetime - aos[len(aos)-1]).seconds//60).zfill(2) + ":" + str((date_datetime - aos[len(aos)-1]).seconds%60).zfill(2)) 
                        
                    else: pass   
                    
                    #control in case of fail in aos
                    if counter != 0:
                        counter = counter -1
                    else:
                        counter = 3
                        if len(aos) < len(tca):
                            delta = timedelta(minutes=int(link_duration[len(link_duration)-1][0:2]), seconds=int(link_duration[len(link_duration)-1][3:5]))
                            aos.append(los[len(los)-1]-delta)
                            #Find row in the weather dataframe
                            if gs[3] == "Vigo":
                                diff = abs(self.weather.df_vigo_2023['datetime'] - date_datetime_partial)
                                index = diff.idxmin()
                                try:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                except:
                                    if date_datetime.hour < 23:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                    else:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                        
                                cloud_cover.append(clouds)
                                visibility.append(vis)
                            elif gs[3] == "Tenerife":
                                diff = abs(self.weather.df_tenerife_2023['datetime'] - date_datetime_partial)
                                index = diff.idxmin()
                                try:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                    if vis == 0: vis = 0.1
                                
                                except:
                                    if date_datetime.hour < 23:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                        if vis == 0: vis = 0.1
                                
                                    else:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        
                                cloud_cover.append(clouds)
                                visibility.append(vis)
                            elif gs[3] == "Sierra Nevada":
                                diff = abs(self.weather.df_sn_2023['datetime'] - date_datetime_partial)
                                index = diff.idxmin()
                                try:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                except:
                                    if date_datetime.hour < 23:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                    else:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                        
                                cloud_cover.append(clouds)
                                visibility.append(vis)
                            elif gs[3] == "Atacama":
                                diff = abs(self.weather.df_atacama_2023['datetime'] - date_datetime_partial)
                                index = diff.idxmin()
                                try:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                except:
                                    if date_datetime.hour < 23:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                    else:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                        
                                cloud_cover.append(clouds)
                                visibility.append(vis)
                            else:
                                diff = abs(self.weather.df_hawai_2023['datetime'] - date_datetime_partial)
                                index = diff.idxmin()
                                try:
                                    clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["cloudcover"]
                                    vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour]["visibility"]
                                    if vis == 0: vis = 0.1
                                    #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                except:
                                    if date_datetime.hour < 23:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour+1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                
                                    else:
                                        clouds = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["cloudcover"]
                                        vis = self.weather.df_vigo_2023.iloc[index]["hours"][date_datetime.hour-1]["visibility"]
                                        if vis == 0: vis = 0.1
                                        #print("Nubes "+ str(date_datetime) + ": " + str(clouds)+ "\n\n")
                                        
                                cloud_cover.append(clouds)
                                visibility.append(vis)
                            snr = self.snr(date_datetime, satellite, bluffton, cloud_cover[len(cloud_cover)-1], visibility[len(visibility)-1])[0]
                            snr_rise.append(snr)
                        else:
                            pass       
                    
            #Check lists
            if len(tca) > len(los):
                los.append("Lost Of Signal in 2024")
                link_duration.append("Lost Of Signal in 2024")
                snr_set.append("Lost Of Signal in 2024")  
        
            print("len aos: "+str(len(aos)))
            print("len tca: "+str(len(tca)))
            print("len los: "+str(len(los)))
            print("link duration: "+str(len(link_duration)))
            print("cloud cover: "+str(len(cloud_cover)))
            print("visibility: "+str(len(visibility)))
            print("snr rise: "+str(len(snr_rise)))
            print("snr tca: "+str(len(snr_tca)))
            print("snr_set: "+str(len(snr_set))) 
            print("cloud info: "+str(len(cloud_info)))  
            print()  
                
            #Save data in the DataFrames
            for j in range(len(aos)):
                new_data = pd.DataFrame({"Rise above 10º": [aos[j]], "TCA": [tca[j]], "Set below 10º": [los[j]], "Link duration": [link_duration[j]], "Cloud cover": [cloud_cover[j]],\
                    "Visibility": [visibility[j]], "SNR rise (dB)": [snr_rise[j]], "SNR TCA (dB)": [snr_tca[j]], "SNR set (dB)": [snr_set[j]], "Cloud info": [cloud_info[j]]})
                if gs[3] == "Vigo":
                    self.df_vigo = pd.concat([self.df_vigo, new_data], ignore_index=True)
                elif gs[3] == "Tenerife":
                    self.df_tenerife = pd.concat([self.df_tenerife, new_data], ignore_index=True)
                elif gs[3] == "Sierra Nevada":
                    self.df_sn = pd.concat([self.df_sn, new_data], ignore_index=True)
                elif gs[3] == "Atacama":
                    self.df_atacama = pd.concat([self.df_atacama, new_data], ignore_index=True)
                elif gs[3] == "Hawai":
                    self.df_hawai = pd.concat([self.df_hawai, new_data], ignore_index=True)            

            print("Done")    
            #Save data in csv
            self.df_vigo.to_csv("passes_"+name_sat+"_vigo_2023_dataframe.csv", index=False)
            self.df_tenerife.to_csv("passes_"+name_sat+"_tenerife_2023_dataframe.csv", index=False)
            self.df_sn.to_csv("passes_"+name_sat+"_sn_2023_dataframe.csv", index=False)
            self.df_atacama.to_csv("passes_"+name_sat+"_atacama_2023_dataframe.csv", index=False)
            self.df_hawai.to_csv("passes_"+name_sat+"_hawai_2023_dataframe.csv", index=False)
       
        return
    
    def snr(self, time, sat, bluffton, cloud_cover, visibility):
        difference = sat - bluffton
        ts = load.timescale()
        time = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second)
        topocentric = difference.at(time)
        alt, azimuth, distance = topocentric.altaz()
        alt = alt.degrees
        azimuth = azimuth.degrees
        distance = distance.km
        snr = str(nat2db(self.lb.snr_pin_sat_2023(distance, visibility)))[0:7]
        data_cc, message_cc = self.lb.cloud_cover(cloud_cover)
        
        if snr[0] == "-": 
            snr = "None"
        
        return snr, data_cc, message_cc
    
if __name__ == "__main__":    
    coord_tenerife = (28.299637836987383, -16.511022084597542, 2.267, "Tenerife")
    coord_vigo = (42.16996139529612, -8.68781081995837, 0.4, "Vigo")
    coord_sn = (37.06285060758391, -3.3859696887312576, 2.896, "Sierra Nevada")
    coord_atacama = (-23.07277599786815, -67.9804469733406, 5.058, "Atacama")  
    coord_hawai = (19.82605380110838, -155.4719893582595, 4.205, "Hawai")
    coord_gs = [coord_vigo, coord_tenerife, coord_sn, coord_atacama, coord_hawai]
    t23 = TLE_2023()
    
    #Lume-1
    t23.save_dataframe("lume-1_2023.txt")
    name_sat = "Lume-1"
    t23.passes_2023(name_sat, coord_gs)
    
    #Terra
    t23.save_dataframe("terra_2023.txt")
    name_sat = "Terra"
    t23.passes_2023(name_sat, coord_gs)
    
    
    