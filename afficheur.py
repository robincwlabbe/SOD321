import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

instances = '40'
with open("C:/Users/User/Desktop/3A/M2/SOD321/SOD321/instance_"+instances+"_1.txt") as f:
    lines = f.readlines()
x = pd.read_csv("C:/Users/User/Desktop/3A/M2/SOD321/SOD321/results_"+instances+".csv", header = None)

ids = []
for i in range(len(x)):
    for j in range(len(x)):
        if x.iloc[i,j] == 1:
            ids.append([i+1,j+1])

ids = pd.DataFrame(ids)


    
longitude = []
latitude = []

for i in range(10, len(lines)):
    xy = lines[i].split()
    longitude.append(int(xy[0]))
    latitude.append(int(xy[1]))


i_0 = int(lines[1])
i_f = int(lines[2])
path = [i_0]
while len(path)<=len(ids):
    path.append(ids[ids.iloc[:,0]==path[-1]].iloc[0,1])  

regions = []
for r in lines[6].split():
    regions.append(int(r))

path_latitude = [latitude[i-1] for i in path]
path_longitude = [longitude[i-1] for i in path]


scatter = plt.scatter(longitude,latitude,c=regions)
labels_regions = [regions[i] for i in sorted(np.unique(regions, return_index=True)[1])]
plt.legend(handles=scatter.legend_elements()[0],
           labels=labels_regions)
plt.plot(path_longitude,path_latitude,label = "path")

plt.grid()


for i in range(len(x)):
    plt.text(x = longitude[i]+0.5 ,y = latitude[i]+0.5, s = i+1)
plt.xlabel("longitude")
plt.ylabel("latitude")
plt.show()
