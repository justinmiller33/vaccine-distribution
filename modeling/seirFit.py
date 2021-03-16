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


        self.infections = self.usData["est_infections_lower"].values.astype(int)
        self.dates = self.usData["date"].values # Stored as str

        self.initBeta = 1 # Not too important in the long run
        self.descentDay = 1 # initBeta handles first checkpoint (day 2)
        self.checkPoints = {"t":[1],"beta":[self.initBeta]}

        self.betaInterval = 2 # How often to change beta value


    def runSeir(self, runLength, checkPoints):

        print("here")
        self.seir = SEIRSModel(beta = self.initBeta, sigma = 1/5.2, gamma = 1/12.39, initN = 300000, initI=33)
        print("model")
        self.seir.run(runLength, checkpoints = checkPoints)
        print("ran")
        
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

            print(testCheckpoint["beta"])

            self.runSeir(self.descentDay + self.betaInterval, testCheckpoint)
            
            d = self.seirI[-1] - self.histI[-1] # getting distance
            #print(d)
            signD = abs(d) / d
            #print(signD)

            overshot = signD * -1 == lastSign
            # print(overshot)

            lastSign = signD
            lastD = d
            
            # Break condition
            error = d / self.histI[-1]

            
            if abs(error) < targetError:

                self.checkPoints = {}
                self.checkPoints["t"] = testCheckpoint["t"]
                self.checkPoints["beta"] = testCheckpoint["beta"]

            elif not overshot:

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
                
            print(testCheckpoint["beta"])
            
    def runModelMatcher(self):

        while self.descentDay < 100:
            
            self.descent()
            self.descentDay = self.descentDay + self.betaInterval

            
            
        

        
# Instantiating
model = Model()

    

