import json, requests
import pandas as pd
import matplotlib.pyplot as plt

class Weather():
    
    def __init__(self):
        self.json_data = None
        self.json_vigo_2023 = None
        self.json_tenerife_2023 = None
        self.json_sn_2023 = None
        self.json_atacama_2023 = None
        self.json_hawai_2023 = None
        self.df_vigo_2023 = None
        self.df_tenerife_2023 = None
        self.df_sn_2023 = None
        self.df_atacama_2023 = None
        self.df_hawai_2023 = None
    
    #Load current data
    def download_json(self, url, saved_file):
            response = requests.get(url)

            if response.status_code == 200:
                self.json_data = response.json()

                with open(saved_file, 'w') as json_file:
                    for i in range(len(self.json_data["hourly"]["visibility"])):
                        self.json_data["hourly"]["visibility"][i] = self.json_data["hourly"]["visibility"][i] / 1000 #km
                        self.json_data["hourly"]["time"][i] = self.json_data["hourly"]["time"][i][0:10] + " " + self.json_data["hourly"]["time"][i][11:] + ":00 UTC"
                    json.dump(self.json_data, json_file, indent=2)
            else:
                print(f"Error downloading weather json. Status code: {response.status_code}") 
            return 
        
    #Load data 2023    
    def read_json(self):
        with open("weather_vigo_2023.json") as f:
            self.json_vigo_2023 = json.load(f)
        self.df_vigo_2023 = pd.json_normalize(self.json_vigo_2023["days"]) 
        with open("weather_tenerife_2023.json") as f:
            self.json_tenerife_2023 = json.load(f)  
        self.df_tenerife_2023 = pd.json_normalize(self.json_tenerife_2023["days"])      
        with open("weather_sn_2023.json") as f:
            self.json_sn_2023 = json.load(f)
        self.df_sn_2023 = pd.json_normalize(self.json_sn_2023["days"])    
        with open("weather_atacama_2023.json") as f:
            self.json_atacama_2023 = json.load(f)
        self.df_atacama_2023 = pd.json_normalize(self.json_atacama_2023["days"])    
        with open("weather_hawai_2023.json") as f:
            self.json_hawai_2023 = json.load(f)
        self.df_hawai_2023 = pd.json_normalize(self.json_hawai_2023["days"])    
        return
    
    def graph_weather(self, df, coord_gs, name):
        df['datetime'] = pd.to_datetime(df['datetime'])
        init = pd.Timestamp('2023-03-01')
        end = pd.Timestamp('2023-03-31')
        df_partial = df[(df['datetime'] >= init) & (df['datetime'] <= end)]
        
        fig, ax = plt.subplots()
        ax.set_title("Cloud cover and visibility: {:.4f}º, {:.4f}º ({})".format(coord_gs[0], coord_gs[1], name))
        ax.plot(df_partial["datetime"], df_partial["cloudcover"], label = "Cloudcover", color = "blue")
        ax.set_xlabel("Date", fontsize = 13)
        ax.tick_params(axis='y', colors='blue', labelsize=12, which='both', direction='in', length=6, width=2)
        ax.tick_params(axis='x', labelsize=10)
        ax.set_ylabel("Cloudcover [%]", color = "blue", fontsize = 13) 
        ax = ax.twinx()      
        ax.plot(df_partial["datetime"], df_partial["visibility"], label = "Visibility", color = "red")
        ax.tick_params(axis='y', colors='red', labelsize=12, which='both', direction='in', length=6, width=2)
        ax.set_ylabel("Visibility [km]", color = "red", fontsize = 13) 
        
        return
                    
if __name__ == "__main__":
    
    weather = Weather()
    weather.download_json("https://api.open-meteo.com/v1/forecast?latitude=28.299637836987383&longitude=-16.511022084597542&hourly=visibility&forecast_days=1", "weather.json")
    elem = weather.json_data["hourly"]["visibility"][3]
    weather.read_json()
    print(weather.df_vigo_2023)
    print()
    print(weather.df_vigo_2023["cloudcover"][0])
    