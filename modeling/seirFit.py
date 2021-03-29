# Class to fir seirsplus model checkpoints to coronavirus projection data

import numpy as np
import pandas as pd
from seirsplus.models import *
import sys
import matplotlib.pyplot as plt
import copy
from datetime import datetime

class Model:

    # Initializing params and ihme data
    def __init__(self):

        # CSV of projected case data
        # SOURCE: http://www.healthdata.org/covid/data-downloads
        self.ihmePath = sys.path[0] + "/ihme-covid19/2021_03_25/reference_hospitalization_all_locs.csv"
        self.ihmeData = pd.read_csv(self.ihmePath) # 4 second read time

        self.usData = self.ihmeData[self.ihmeData.location_name == "United States of America"]
        self.usPop = 328000000

        # self.data = self.usData # Initializing first runthrough as usData
        
        # CSV of state names and populations
        # SOURCE: https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-total.html
        self.stateData = pd.read_csv("/home/justinmiller/devel/vaccine-distribution/modeling/stateData.csv")
        self.stateData = self.stateData.drop(index = [0,1,2,3,4,13,56])
        self.stateList = self.stateData.NAME
        self.statePop = self.stateData.POPESTIMATE2019

        # CSV of total vaccinated by state
        # SOURCE: https://ourworldindata.org/covid-vaccinations
        self.vaccineData = pd.read_csv("/home/justinmiller/devel/vaccine-distribution/modeling/owid-covid-data.csv")
        self.usVaccineData = self.vaccineData[self.vaccineData["location"] == "United States"]
        self.vaccinationStart = 358 # Start of reported fully vaccinated (01/09/21) 
        self.vaccineDates = self.usVaccineData.date.values[self.vaccinationStart:]
        self.fullyVaccinated = self.usVaccineData.people_fully_vaccinated.values[self.vaccinationStart:]
        
        # left setting nans
        for i in np.where(np.isnan(self.fullyVaccinated))[0]:
            self.fullyVaccinated[i] = self.fullyVaccinated[i-1]
        
        
        self.infections = self.usData["est_infections_mean"].values.astype(int)
        self.dates = self.usData["date"].values # Stored as str
        self.dateList = [datetime.strptime(date, '%Y-%m-%d').date() for date in self.dates]
        
        self.initBeta = 1 # Not too important in the long run
        self.descentDay = 1 # initBeta handles first checkpoint (day 2)
        self.checkPoints = {"t":[1],"beta":[self.initBeta]}

        self.betaInterval = 10 # How often to change beta value


    def runSeir(self, runLength, checkPoints):

       
        self.seir = SEIRSModel(beta = self.initBeta, sigma = 1/5.2, gamma = 1/12.39, initN = 328000000, initI=33)
        self.seir.run(runLength, checkpoints = checkPoints)
   
        
        # Normalizing seir output
        self.seirI = self.seir.numI[-10 * runLength:][::10]
        self.histI = self.infections[:runLength]
        
    def descent(self):

        # Params
        targetError = 0.0005
        a = .04 # this is the change in beta rate over day
        aDecreaseFactor = 0.5
        
        # Inits
        lastSign = 0
        lastD = 0
        error = 1


        # Init seir arrays
        testCheckpoint = self.checkPoints
        testCheckpoint["t"].append(self.descentDay)
        testCheckpoint["beta"].append(testCheckpoint["beta"][-1])
            

        # Rebounding descent until acceptable error
        while abs(error) > targetError:

           

            self.runSeir(self.descentDay + self.betaInterval, testCheckpoint)
            
            d = self.seirI[-1] - self.histI[-1] # getting distance
            # print(d)
            signD = abs(d) / d
            #print(signD)

            overshot = signD * -1 == lastSign
            

            lastSign = signD
            
            
            # Break condition
            error = d / self.histI[-1]

            
            if abs(error) < targetError:

                self.checkPoints = {}
                self.checkPoints["t"] = testCheckpoint["t"]
                self.checkPoints["beta"] = testCheckpoint["beta"]

                print(str(self.descentDay) + ": " + str(testCheckpoint["beta"][-1])[:5])

            elif not overshot:

                # If stuck in ascent swap direction
                if (d > lastD and lastD != 0):
                    testCheckpoint["beta"][-1] = testCheckpoint["beta"][-1] - signD * a
                    
                else:
                    testCheckpoint["beta"][-1] = testCheckpoint["beta"][-1] - signD * a

                tempCheckpoint = {}
                tempCheckpoint["t"] = testCheckpoint["t"]
                tempCheckpoint["beta"] = testCheckpoint["beta"]
                
                
            elif overshot:

                a = a * aDecreaseFactor
                testCheckpoint["beta"][-1] = testCheckpoint["beta"][-1] - signD * a
                tempCheckpoint = {}
                tempCheckpoint["t"] = testCheckpoint["t"]
                tempCheckpoint["beta"] = testCheckpoint["beta"]
                
            lastD = d
            
    def runModelMatcher(self):

        while self.descentDay < 417: # March 27th 2021 from start of pandemic data
            
            self.descent()
            self.descentDay = self.descentDay + self.betaInterval

            
    def createUnvaccinatedCheckpoints(self):

        self.unvaxCheckpoints = {}
        self.unvaxCheckpoints["t"] = self.checkPoints["t"]
        self.unvaxCheckpoints["beta"] = self.checkPoints["beta"]

        for timeIdx in range(len(self.checkPoints["t"])):

            if self.checkPoints["t"][timeIdx] >= self.vaccinationStart:

                adjTime = self.checkPoints["t"][timeIdx] - self.vaccinationStart
                vaccinations = self.fullyVaccinated[adjTime]

                vaxProp = (vaccinations / self.usPop)
                self.unvaxCheckpoints["beta"][timeIdx] = self.checkPoints["beta"][timeIdx] / (1 / 1-vaxProp)
        
    def createDiffVaccinatedCheckpoints(self, vaxEfficiency):

        self.diffVaxCheckpoints = {}
        self.diffVaxCheckpoints["t"] = self.checkPoints["t"]
        self.diffVaxCheckpoints["beta"] = self.checkPoints["beta"]

        for timeIdx in range(len(self.checkPoints["t"])):

            if self.checkPoints["t"][timeIdx] >= self.vaccinationStart:

                adjTime = self.checkPoints["t"][timeIdx] - self.vaccinationStart
                vaccinations = self.fullyVaccinated[adjTime]

                vaccinations = vaccinations * vaxEfficiency / 0.78 # Rate of shots in arms SOURCE: https://www.npr.org/sections/health-shots/2021/01/28/960901166/how-is-the-covid-19-vaccination-campaign-going-in-your-state
                
                vaxProp = (vaccinations / self.usPop)
                self.diffVaxCheckpoints["beta"][timeIdx] = self.checkPoints["beta"][timeIdx] * (1 / 1-vaxProp)
        

    def deathWithVaccine(self, seir, vaxDeathRate, vaxEfficiency = 0.78):

        defaultDeathRate = 0.01
        
        infections = seir.numI[::10]
        dateList = self.dateList[:len(infections)]

        vaxDataStart = self.vaccinationStart - 5 # diff between vax and case dates
        
        deathList = infections[:vaxDataStart] * 0.01 

        for vaxIdx in range(len(self.fullyVaccinated)):

            deathRate = (self.usPop * defaultDeathRate - self.fullyVaccinated[vaxIdx] * (vaxEfficiency / 0.78) * vaxDeathRate) / (self.usPop - self.fullyVaccinated[vaxIdx]) # Ratio calculation
            #print(deathRate)
            deathList = np.append(deathList, infections[vaxDataStart + vaxIdx] * deathRate)

        return deathList
        
# Instantiating
model = Model()
model.runModelMatcher()
vaxCheckpoints = copy.deepcopy(model.checkPoints)
vaxSeir = copy.deepcopy(model.seir)
#print(vaxCheckpoints)



model.createUnvaccinatedCheckpoints()
unVaxCheckpoints = copy.deepcopy(model.checkPoints)
model.runSeir(model.descentDay, model.unvaxCheckpoints)
unVaxSeir = copy.deepcopy(model.seir)
#print(unVaxCheckpoints)

model.createDiffVaccinatedCheckpoints(1)
bestVaxCheckpoints = copy.deepcopy(model.checkPoints)
model.runSeir(model.descentDay, model.diffVaxCheckpoints)
bestVaxSeir = copy.deepcopy(model.seir)
#print(bestVaxCheckpoints)

input()

# Plotting historical to model
plt.plot(model.dateList[:len(model.histI)], model.histI, label = "Historical Data")
plt.plot(model.dateList[:len(model.histI)], model.seirI, label = "SEIRS Model fit to Transmission Rate")
plt.legend()
plt.xlabel("Date")
plt.ylabel("Cases (Retrospective Projection)")
plt.title("COVID-19 Daily Case Modeling")
plt.show()

"""
# Plot of transmission
plt.plot(model.dateList[:len(vaxSeir.numI[::10])][::10][1:],model.checkPoints["beta"][1:])
plt.xlabel("Date")
plt.ylabel("Rate of Transmission")
plt.title("Dynamics of the COVID-19 Rate of Transmission")
plt.show()
"""

# Plotting model with vaccine administrations
plt.plot(model.dateList[:len(vaxSeir.numI[::10])], vaxSeir.numI[::10], label = "Actual Vaccinations (78% Administered)")
plt.plot(model.dateList[:len(vaxSeir.numI[::10])], bestVaxSeir.numI[::10], label = "Best Case Vaccinations (100% Administered)")
plt.plot(model.dateList[:len(vaxSeir.numI[::10])], unVaxSeir.numI[::10], label = "Without Vaccinations")
plt.legend()
plt.xlabel("Date")
plt.ylabel("Daily Cases (Modelled)")
plt.title("COVID-19 Case Model w/ Vaccine Administration")
plt.show()


# Death with vaccine
def plotDeaths(title, vaxDeathRate = 0.04):

    caseToDeathDelay = 21
    deathListActual = model.deathWithVaccine(vaxSeir, vaxDeathRate)
    plt.plot(model.dateList[:len(vaxSeir.numI[::10]) + caseToDeathDelay][1 + caseToDeathDelay:], deathListActual, label = "Actual Vaccinations (78% Administered)")
    deathListBestCase = model.deathWithVaccine(vaxSeir, vaxDeathRate, 1)
    plt.plot(model.dateList[:len(bestVaxSeir.numI[::10]) + caseToDeathDelay][1 + caseToDeathDelay:], deathListBestCase, label = "Best Case Vaccinations (100% Administered)")
    plt.plot(model.dateList[:len(unVaxSeir.numI[::10]) + caseToDeathDelay][caseToDeathDelay:], unVaxSeir.numI[::10] * 0.01, label = "Without Vaccinations") # 0.01 is pretty locked in as defaultDeathRate
    plt.legend()
    plt.xlabel("Date")
    plt.ylabel("Daily Deaths (Modelled)")
    plt.title(title)
    plt.show()

# plotDeaths("COVID-19 Death Model w/ RANDOM Vaccine Administration", 0.01)
plotDeaths("COVID-19 Death Model w/ Targeted Vaccine Administration")
