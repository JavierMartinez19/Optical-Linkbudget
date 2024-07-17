import json, csv, requests
from datetime import datetime,timedelta
from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS
from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import ITRS

class SatLoc():
    def __init__(self):
        self.coordinates = None
        self.sat_name = "LUME 1"
        self.timestamp = None
        self.json_name = 'ejemplo_local.json'
        
    def main(self):
        # Name of the file to add satellite position
        #csvfile = "P_AOCS_MAG_0_Y2022.csv"
        tlefile = "lume-1_2022_TLE.json"
        #json_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=json'
        saved_file = self.json_name
        #self.download_json(json_url, saved_file)
        #tlefile=saved_file
        sat_name = self.sat_name
        
        self.coordinates = self.add_position(tlefile, datetime.utcnow(), sat_name)
        
    # def download_json(self, url, save_path):
    #     response = requests.get(url)

    #     if response.status_code == 200:
    #         json_data = response.json()

    #         with open(save_path, 'w') as json_file:
    #             json.dump(json_data, json_file, indent=None)
            
    #         print(f"JSON downloaded and saved in '{save_path}'.")
    #     else:
    #         print(f"ERROR. Status code: {response.status_code}") 
        
    #     with open(save_path, 'r') as json_file:
    #         datos_satelites = json.load(json_file) 
        
    #     if self.sat_name in datos_satelites:       
    #         # Obtener los TLE del satélite deseado
    #         tle_satelite_deseado = datos_satelites[self.sat_name]
            
    #         with open(save_path, 'w') as json_file:
    #             json.dump(json_data, json_file, indent=None)
            

    def fetch_tle_json(self, timestamp, tle_file='ejemplo_local.json'):
        #print("Reading TLE file") 
        with open('{}'.format(tle_file), 'r') as tle:
            data = json.load(tle)

        # Initialize variables to track the closest epoch
        closest_diff = timedelta.max

        for item in data:
            epoch_str = item.get('EPOCH')
            if epoch_str:
                epoch = datetime.fromisoformat(epoch_str)
                diff = abs(timestamp - epoch)

                # Check if the current epoch is closer than the previous closest epoch
                if diff < closest_diff:
                    closest_diff = diff
                    closest_item = item

        return closest_item.get('TLE_LINE1'), closest_item.get('TLE_LINE2'), closest_item.get('EPOCH')


    def get_satellite_coordinates(self, timestamp, tle1,tle2):
        satellite = Satrec.twoline2rv(tle1, tle2)

        # Get the current timestamp
        t = Time(timestamp)

        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s

        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])
        
        # Calculate the position vector of the satellite at the current time
        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
        itrs_geo = teme.transform_to(ITRS(obstime=t))

        location = itrs_geo.earth_location

        return location.geodetic

    def add_position(self, tle_file, timestamp, sat_name):
        self.timestamp = datetime.strptime(str(timestamp)+"Z", "%Y-%m-%d %H:%M:%S.%fZ")

        tle1, tle2, epoch = self.fetch_tle_json(timestamp, tle_file=tle_file)
        position = self.get_satellite_coordinates(timestamp, tle1, tle2)

        latitude = position.lat.deg
        longitude = position.lon.deg
        height = position.height.value
        #position = "Latitude: "+ str(latitude) + "\nLongitude: " + str(longitude) + "\nHeight: " + str(height)
        # print(sat_name)
        # print(timestamp)
        # print(position)
        
        return latitude, longitude, height*1e3

    def print_all(self):
        print(self.sat_name)
        print(self.timestamp)
        print(f"Latitude: {self.coordinates[0]}")
        print(f"Longitude: {self.coordinates[1]}")
        print("Height: {:.2f} km".format(self.coordinates[2]/1e3))
        print()
        
if __name__ == "__main__":
    sc = SatLoc()
    sc.main()