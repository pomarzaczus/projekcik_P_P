# import the necessary packages
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import copy
import geopandas as gpd
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Polygon, Wedge
import matplotlib.patches as mpatches
import time
from datetime import datetime
from sgp4.api import Satrec, WGS72
from math import radians, cos
from skyfield.api import Topos, load
from sgp4.api import Satrec
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.api import EarthSatellite
# Function to generate TLE lines
#from orbital import  KeplerianElements,earth

import math
import ephem
from datetime import timedelta
#https://www.satcatalog.com/component/x-band-transmitter-payload-txc-200/
#https://www.satcatalog.com/component/swift-ktx-transmitter/
#około 12W zakładając 550cm2 powierchni i 100% osiwtlenia
#kamera
#dane

# wymiary obszaru na ziemi który będą pokrywały kamery (km)
height = 115 
width = 30
angle = -5
# Constants for the orbit
distance_between_sat = 10 #dystans w stopniach pomiędzy satelitami
travel_time = 20 #11 czas podruży w sekundach do natępnego miesca zrobienia obrazu
testtime = 4320 #8000 czas trwania badania rzeczywisty testtime*travel_time w sek
dist =  22#22.5 #min dysytans w sekundach jak później od siembie zostały wystartowane 90/4 czyli oddalenie od siebie na mapie nie wiem czy to ma znaczen9e powinno
dist_sek = 30 #sek dysytans w sekundach jak później od siembie zostały wystartowane 
year =2024  #czas rozpoczecia badania
month = 5
day = 18
a = 12  #hour
b= 0    #minutes
c=0     #seckonds
przez = 1         # przeźrosczystość
altitude = 620  # km
mean_motion = 15.079  # revs per day (based on altitude)
eccentricity = 0.0  # Circular orbit
inclination = 90.0  # Polar orbit
arg_perigee = 0.0  # Argument of perigee
mean_anomaly = 0.0  # Mean anomaly
twilight = 0 #-12 * ephem.degree #stopień w którym uznajemy że jest zmierzch/ ciemno w danym miejscu
# Generate TLEs
satellite_number = 10001


def km_to_degrees_latitude(distance_km):
    return distance_km / 111

def km_to_degrees_longitude(distance_km, latitude_deg):
    #print(math.cos(math.radians(latitude_deg)))
    #if latitude_deg > 0:
    return abs(distance_km / (111 * math.cos(math.radians(latitude_deg))))


epoch = datetime(year, month, day,a,b,c)

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

fig, ax = plt.subplots(figsize=(30, 20))

patches = []

plt.ion() 



base =world.plot(
    ax=ax,
    color="lightgray",
    edgecolor="black",
    alpha=0.5
)

colors = 100 * np.random.rand(len(patches))
p = PatchCollection(patches, alpha=0.4)


ax.set_xticks([])
ax.set_yticks([])



plt.title("Basic Map of World with GeoPandas")
plt.show()



ts = load.timescale()
def calculate_time(interval):
    interval*=travel_time 
    #print(year,month,day,a,b,c)
    temp1 = datetime(year,month,day,a,b,c)
    temp2 = temp1 + timedelta(seconds =interval)
    #print(temp2.time())
    return temp2.day,temp2.hour,temp2.minute,temp2.second

def rotate_point(point, center, angle):
    # Convert angle to radians
    angle_rad = np.deg2rad(angle)
    
    # Translate point back to origin
    translated_point = point - center
    
    # Calculate the new coordinates
    rotated_x = translated_point[0] * np.cos(angle_rad) - translated_point[1] * np.sin(angle_rad)
    rotated_y = translated_point[0] * np.sin(angle_rad) + translated_point[1] * np.cos(angle_rad)
    
    # Translate point back to its original location
    final_point = np.array([rotated_x, rotated_y]) + center
    
    return final_point

def get_iss_location(interval,satellite):
    #print(interval)
    #dzien,godz,min,sek =  day,a,b,c
    dzien,godz,min,sek =  calculate_time(interval)
    t = ts.utc(2024, 1, dzien, godz, min, sek)
    #t = ts.utc(2021, 12, dzien, godz, min)
    #print(t)
    geocentric = satellite.at(t)



    position = wgs84.geographic_position_of(geocentric)
    latitude = float(position.latitude.degrees)
    longitude = float(position.longitude.degrees)

    return latitude, longitude


def date_to_tle_format(date):
    # Extract year and day of year
    year = date.year % 100  # Last two digits of the year
    day_of_year = date.timetuple().tm_yday
    # Calculate the fractional part of the day
    fraction_of_day = (date.hour + date.minute / 60 + date.second / 3600) / 24
    # Combine day of year and fractional day
    tle_day = day_of_year + fraction_of_day
    # Format as TLE epoch string
    tle_epoch = f"{year:02}{tle_day:012.8f}"
    return tle_epoch




def generate_tle_elements(satellite_number, epoch, inclination, raan, eccentricity, arg_perigee, mean_anomaly, mean_motion):



    raanek =  ['0','0','0','.','0','0','0','0']
    temp = str(raan)
    xxxxx = 0
    if raan < 10:
        xxxxx =2
        raanek[xxxxx] = temp[0]
        xxxxx=4
        for a in range(2,len(temp)):
            raanek[xxxxx] = temp[len]
            xxxxx+=1
    if raan < 100 and raan >= 10:
        xxxxx =1
        raanek[xxxxx] = temp[0]
        raanek[xxxxx+1] = temp[1]
        xxxxx=4
        for a in range(3,len(temp)):
            raanek[xxxxx] = temp[len]
            xxxxx+=1

    if raan > 100:
        xxxxx =0
        raanek[xxxxx] = temp[0]
        raanek[xxxxx+1] = temp[1]
        raanek[xxxxx+2] = temp[2]
        xxxxx=4
        for a in range(4,len(temp)):
            raanek[xxxxx] = temp[len]
            xxxxx+=1
    
    #print()
     #{epoch.strftime('%y%j.%f')[:14]}00                           
    tle1 = f"1 {satellite_number:05d}U 21001A   {date_to_tle_format(epoch)}  .00000000  00000-0  00000-0 0  9991"
    tle2 = f"2 {satellite_number:05d} {inclination:8.4f} {"".join(raanek)} 0000000 000.0000 000.0000 {mean_motion:11.8f}    01"
    
    return tle1, tle2

# Parameters of time 

colors = [
    "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
    "#800000", "#008000", "#000080", "#808000", "#800080", "#008080",
    "#C0C0C0", "#808080", "#9999FF", "#993366", "#FFFFCC", "#CCFFFF",
    "#660066", "#FF8080", "#0066CC", "#CCCCFF", "#000080", "#FF00FF",
    "#FFFF00", "#00FFFF", "#800080", "#800000", "#008080", "#0000FF",
    "#00CCFF", "#CCFFFF", "#CCFFCC", "#FFFF99", "#99CCFF", "#FF99CC",
    "#CC99FF", "#FFCC99", "#3366FF", "#33CCCC", "#99CC00", "#FFCC00",
    "#FF9900", "#FF6600", "#666699", "#969696", "#003366", "#339966",
    "#003300", "#333300", "#993300", "#993366", "#333399", "#333333",
    "#000000", "#330000", "#660000", "#990000", "#CC0000", "#FF0000",
    "#003300", "#006600", "#009900", "#00CC00"
]


tle_list = []
def generateTle(hour,min,sek,satellite_number):
    a1 = hour
    b1= min
    c1= sek
    longitude = 0
    for s in range(0,16):  # 16 different longitudes
        longitude+=distance_between_sat
        
        raan = longitude 
        for latitude in range(0,4):  # 4 different latitudes
            b1+=dist
            c1+=dist_sek
            #print(a+ int(b/60),int(b)%60,c%60)
            epoch = datetime(year, month, day,a1+ int(b1/60),int(b1)%60,c1%60)
            #print(epoch)
            tle1, tle2 = generate_tle_elements(satellite_number, epoch, inclination, raan, eccentricity, arg_perigee, mean_anomaly, mean_motion)
                #print(tle1, tle2)
            tle_list.append((tle1, tle2))
            satellite_number += 1
        a1 = hour
        b1= min
        c1= sek
generateTle(a,b,c,satellite_number)
# Print TLEs
satelites = []

def gelenerateSsateltietes():
    i=1
    for line1, line2 in tle_list:
        
        #print(line1, line2)
        satelites.append(EarthSatellite(line1, line2, 'ISS kutas '+str(i), ts))
        i+=1


gelenerateSsateltietes()


# month = 5
# day = 18
timee = datetime.now().time()
s = ephem.Sun()
gatech = ephem.Observer()
x = 0
#dzien, godz,minu,seku= copy.copy(day),copy.copy(a),copy.copy(b),copy.copy(c)
#print(year,month,dzien,godz,minu,seku)
day_time_picture = []
prev = []
for _ in range(64):
    prev.append(0)
    day_time_picture.append(0)
lat, len = 0,0

while True:
    #time.sleep(1)
    x+=1
    for xdddd in range(0,64):
       
        dzien, godz,minu,seku= calculate_time(x)
   
        lat, len = get_iss_location(x,satelites[xdddd])
                #print(lat,len)
            #print(a,lat, len)
        gatech.lon, gatech.lat = str(len), str(lat)
        gatech.date = str(year)+'/'+str(month)+'/'+str(dzien)+" "+ str(godz)+":"+str(minu)+':'+str(seku)
        # print(gatech.date)
        
        # print(year,month,day,a,b,c)
        s.compute(gatech)
        if s.alt < twilight:
            continue
        day_time_picture[xdddd]+=1
        longitude_deg2 = km_to_degrees_longitude(height,lat)
        longitude_deg = km_to_degrees_longitude(width,lat)
        lat= lat + longitude_deg2/2
        len = len + longitude_deg/2
        if longitude_deg2 > 10:
            continue
            #print("----------")
        
        
        if prev[xdddd] > lat:
            rectangle = np.array([(len ,lat-longitude_deg2),(len,lat),(len-longitude_deg,lat),(len-longitude_deg,lat-longitude_deg2)])
            center = np.mean(rectangle, axis=0)
            rotated_points = np.array([rotate_point(point, center, angle) for point in rectangle])
            ax.add_patch(mpatches.Polygon(rotated_points, color = "red",alpha = przez))
        else:
            rectangle = np.array([(len ,lat-longitude_deg2),(len,lat),(len-longitude_deg,lat),(len-longitude_deg,lat-longitude_deg2)])
            center = np.mean(rectangle, axis=0)
            rotated_points = np.array([rotate_point(point, center, -1*angle) for point in rectangle])
            ax.add_patch(mpatches.Polygon(rotated_points, color = "red",alpha = przez))
        prev[xdddd] = lat
        # rectangle = np.array([(len ,lat-longitude_deg2),(len,lat),(len-longitude_deg,lat),(len-longitude_deg,lat-longitude_deg2)])
        # ax.add_patch(mpatches.Polygon(rectangle, color = "red",alpha = 1))
    # fig.canvas.draw()
    # fig.canvas.flush_events()  
    if x% 100 == 0:
        print(str(x/testtime*100)+" %")
    if x % testtime == 0:
        #fig.canvas.draw()
        #fig.canvas.flush_events()   
        print((np.array(day_time_picture)/(x) ) *100)
        plt.savefig('iamge_'+str(height)+'_'+str(width)+'_'+str(angle)+'_'+str(distance_between_sat)+'_'+str(travel_time)+'_'+str(testtime)+'_'+str(dist)+'_'+str(dist_sek)+'_'+str(year)+'_'+str(month)+'_'+str(day)+'_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(inclination)+'.png',dpi = 420)
        #input()
        break
#/(x*travel_time) *100
#procent zdjęc w świetle po tygodniu
# [68.56545455 68.66181818 68.86727273 68.75818182 66.39454545 66.59636364
#  66.72545455 66.45454545 60.27090909 60.50181818 60.41454545 60.19454545
#  57.67818182 57.89454545 57.76727273 57.55090909 56.74       56.96181818
#  56.82545455 56.6        56.78       57.00181818 56.86727273 56.64
#  57.82       58.04363636 57.92181818 57.69818182 60.65818182 60.89090909
#  60.80727273 60.57636364 67.19090909 67.42545455 67.47454545 67.28727273
#  67.82       67.94       68.12545455 68.06545455 61.01272727 61.
#  61.21818182 61.28545455 57.93454545 57.86       58.08909091 58.20909091
#  56.77272727 56.69636364 56.91454545 57.05454545 56.65272727 56.56727273
#  56.79090909 56.93090909 57.48545455 57.41090909 57.63090909 57.76181818
#  59.87090909 59.84181818 60.06181818 60.15272727]