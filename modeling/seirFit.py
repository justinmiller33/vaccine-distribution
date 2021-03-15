# Class to fir seirsplus model checkpoints to coronavirus projection data

import numpy as np
import pandas as pd
from seirsplus.models import *
import sys
import matplotlib.pyplot as plt

class Model:

    # Initializing params and ihme data
    def __init__(self):

        self.ihmePath = sys.path[0].replace("/modeling","/IHME_DATA_2021_03_11/reference_hospitalization_all_locs.csv")
        self.ihmeData = pd.read_csv(self.ihmePath) # 4 second read time

        self.usData = self.ihmeData[self.ihmeData.location_name == "United States of America"]


        self.infections = self.usData["est_infections_lower"].values.astype(int)
        self.dates = self.usData["date"].values # Stored as str

        self.initBeta = 1 # Not too important in the long run
        self.descentDay = 1
        self.checkPoints = {"t":[1],"beta":[self.initBeta]}

        self.betaInterval = 1 # How often to change beta value


    def runSeir(self, runLength, checkPoints):

        
        self.seir = SEIRSModel(beta = self.initBeta, sigma = 1/5.2, gamma = 1/12.39, initN = 300000000, initI=33)
        self.seir.run(runLength, checkpoints = checkPoints)

        # Normalizing seir output
        self.seirI = self.seir.numI[1::10]
        self.histI = self.infections[:runLength]

    
    def descent(self):

        # Params
        targetError = 0.05
        a = .04 * self.betaInterval # this is the change in beta rate over day
        
        # Inits
        lastD = 0
        error = 1

        # Rebounding descent until acceptable error
        while abs(error) < targetError:
            
            testCheckpoint = self.checkpoints
            testCheckpoint["t"] = testCheckpoint["t"].append(self.descentDay)
            testCheckpoint["beta"] = testCheckpoint["beta"].append(testCheckpoint["beta"][-1])
            

            self.runSeir(self.descentDay, testCheckPoint)

            d = self.seirI[-1] - self.histI[-1]
            delta = d - lastD 

            lastD = d

            # Error condition to decide rerun
            error = d / self.histI

    def runModelMatcher(self):

        self.descentDay = self.descentDay + self.betaInterval
        self.descent()
            
            
        

        
# Instantiating
model = Model()

    

