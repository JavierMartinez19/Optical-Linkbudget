import matplotlib.pyplot as plt
import os 
import geopandas as gpd
import seaborn as sns
import earthpy as et 

class PlotMap():
    
    def __init__(self):
        self.lat = None
        self.long = None
       
    
    def plotmap(self, sat_name, coord_gs, coord_sat_init, coord_pos, coord_neg): 
        
        lat_pos = []
        long_pos = []
        lat_neg = []
        long_neg = []

        for i in range(len(coord_pos)):
            for j in range(len(coord_pos[i])):
                lat_pos.append(coord_pos[i][j][0])
                long_pos.append(coord_pos[i][j][1])
        
        for i in range(len(coord_neg)):
            for j in range(len(coord_neg[i])):
                lat_neg.append(coord_neg[i][j][0])
                long_neg.append(coord_neg[i][j][1])
            
        os.chdir(os.path.join(et.io.HOME, 'earth-analytics'))
        
        sns.set_theme(font_scale=1.5)
        sns.set_style("whitegrid")
        
        worldBound_path = os.path.join("data", "spatial-vector-lidar", "global", "ne_110m_land", "ne_110m_land.shp")
        worldBound = gpd.read_file(worldBound_path)
        
        # add_gs = []
        # for i in range(len(coord_gs)):
        #     add_gs.append([])
        #     add_gs[i].append([coord_gs[i][1], coord_gs[i][0]])
        
        # add_gs = np.array(add_gs)
        # gs_locations = [Point(xy) for xy in add_gs]
        # gs_locations
        
        # gs_locations = gpd.GeoDataFrame(gs_locations, columns=['geometry'], crs=worldBound.crs)
        # gs_locations.head(2)
        
        # Set working dir & get data
        
        figure1, ax = plt.subplots(figsize=(12,8))
        
        worldBound.plot(figsize=(10,5),color='black',ax=ax)
        ax.scatter(long_neg, lat_neg, s=1, color = "red")
        ax.scatter(long_pos, lat_pos, s=1, color = "green")
        ax.scatter(coord_sat_init[1], coord_sat_init[0], s=120, color = "black", marker='*')
        ax.scatter(coord_sat_init[1], coord_sat_init[0], s=60, color = "yellow", marker='*', label = "Satellite: "+sat_name)
        
        for gs in coord_gs:
            ax.scatter(gs[1], gs[0], s=120, color = "black", marker='*')
            ax.scatter(gs[1], gs[0], s=60, color = "cyan", marker='*')
        # gs_locations.plot(ax=ax, color='black', marker='*', markersize=90),
        # gs_locations.plot(ax=ax, color='cyan', marker='*', markersize=45, label = "gs")
        ax.set(xlabel="Longitude (degrees)",ylabel="Latitude (degrees)", title="Global Map - Geographic Coordinate System")
        ax.set_axisbelow(True)
        
        #Leyend
        ax.scatter([], [], color='cyan', marker='*', label='Ground Stations')
        ax.scatter([], [], color='green', linestyle='-', label='LOS: True')
        ax.scatter([], [], color='red', linestyle='-', label='LOS: False')
        ax.legend(fontsize = 12)
        
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        
        return