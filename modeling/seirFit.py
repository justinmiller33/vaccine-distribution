# Class to fir seirsplus model checkpoints to coronavirus projection data

import numpy as np
import pandas as pd
from seirsplus.models import *
import sys
import matplotlib.pyplot as plt

class Model:

    # Initializing params and ihme data
    def __init__(self):

        self.ihmePath = sys.path[0].replace("\\modeling","\\IHME_DATA_2021_03_11/reference_hospitalization_all_locs.csv")
        self.ihmeData = pd.read_csv(self.ihmePath) # 4 second read time

        self.usData = self.ihmeData[self.ihmeData.location_name == "United States of America"]
        self.usPop = 328000000

        # self.data = self.usData # Initializing first runthrough as usData
        
        # CSV of state names and populations
        self.stateData = pd.read_csv("C:/devel/vaccine-distribution/stateData.csv")
        self.stateData = self.stateData.drop(index = [0,1,2,3,4,13,56])
        self.stateList = self.stateData.NAME
        self.statePop = self.stateData.POPESTIMATE2019

        self.vaccineData = pd.read_csv("C:/devel/vaccine-distribution/coronavirus-data-explorer.csv")
        self.usVaccineData = self.vaccineData[self.vaccineData["Country name"] == "United States"]
        self.vaccinationStart = 358
        self.vaccineDates = self.usVaccineData.Day.values[self.vaccinationStart:]
        self.fullyVaccinated = self.usVaccineData.people_fully_vaccinated.values[self.vaccinationStart:]
        
        # left setting nans
        for i in np.where(np.isnan(self.fullyVaccinated))[0]:
            self.fullyVaccinated[i] = self.fullyVaccinated[i-1]
        
        
        self.infections = self.usData["est_infections_lower"].values.astype(int)
        self.dates = self.usData["date"].values # Stored as str

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

        while self.descentDay < 400:
            
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

                vaxProp = vaccinations / self.usPop
                self.unvaxCheckpoints["beta"][timeIdx] = self.checkPoints["beta"][timeIdx] / (1 / 1-vaxProp)
        

        
# Instantiating
model = Model()
model.runModelMatcher()
vaxCheckpoints = model.checkPoints
vaxSeir = model.seir

model.createUnvaccinatedCheckpoints()
unvaxCheckpoints = model.checkPoints
unvaxSeir = model.seir

# model.runSeir(model.descentDay, model.unvaxCheckpoints)

