# Scraping csv of stadiums / arenas from pdf

import numpy as np
import pandas as pd
import re
from PyPDF2 import PdfFileReader

file = open("stadiums.txt","r")
content = file.read()

lineList = content.split("\n")[:-1]
townList = np.array([])
stateList = np.array([])
nameList = np.array([])
capacityList = np.array([])

for line in lineList:

    town = line.split(" (")[0]
    state = line.split(")-")[0][-2:]
    name = line.split("-")[-1].split(" (")[0]
    capacity = line.split(" (")[-1][:-1]
    
    townList = np.append(townList, town)
    stateList = np.append(stateList, state)
    nameList = np.append(nameList, name)
    try:
        capacityList = np.append(capacityList, int(capacity.replace(",","")))
    except:
        capacityList = np.append(capacityList, 8000)
        print("DEBUG: Missing Capacity")

data = {}
data["town"] = townList
data["state"] = stateList
data["name"] = nameList
data["capacity"] = capacityList

df = pd.DataFrame(data)
df.to_csv("stadiums.csv")
