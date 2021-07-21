# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 15:35:02 2021

@author: elozon
"""
import pandas as pd
import numpy as np
import os

from InputSetup import FFInputCreation as inputCreation
from ParameterManipulation import createTSParamFiles, createFFParamFiles, copyNominalFiles
from stochasticTurbulenceTools_mod import stochasticTurbulence

#Taken from the CaseCreation.ipynb

nSeeds = 1
HubHt = 148.84
D=246.0
nTurbs = 4
Ranges = False
workdir=os.getcwd()

indexes = ['TI', 'Vhub', 'Shear', 'T1x-locs', 'T1y-locs', 'T2x-locs', 'T2y-locs', 'T3x-locs', 'T3y-locs', 'T4x-locs', 'T4y-locs']
TI = ['A']
Vhub = [11.4]
Shear = [0.1]
locs = [-800, -1600,-800, 0 , 800, 0, 800, 1600]
data = TI + Vhub +Shear +locs
if Ranges == False:
    ParamVals = pd.DataFrame(data, index = indexes)
    nCases = int(ParamVals.shape[1])
else:
    #writeParamVals()
    OrigVals = pd.read_csv('../SampleFiles/ParamRanges.csv', names=['Var','MinVar','MaxVar','NomVar'])
    OrigVals = OrigVals['NomVar']

caseNames=['None']*nCases
for case in range(nCases):
    if case < 10:
        case = '0'+str(case)
    caseNames[int(case)] = 'Case{0}'.format(case)
    
ParamVals.columns=caseNames



os.chdir(workdir)

createTSParamFiles(caseNames,nSeeds,ParamVals,D,HubHt,'TSParams.txt')

def writeTimeSeriesFile(fileOut,yloc,zloc,u,v,w,time):
    import math
    """ Write a TurbSim primary input file, 

    """

    print('Writing {0}'.format(fileOut))
    # --- Writing TurbSim user-defined time series file
    with open(fileOut, 'w') as f:
        f.write('--------------TurbSim v2.00.* User Time Series Input File-----------------------\n')
        f.write('     Time series input from Experimental Data\n')
        f.write('--------------------------------------------------------------------------------\n')
        f.write('          3 nComp - Number of velocity components in the file\n')
        f.write('          1 nPoints - Number of time series points contained in this file (-)\n')
        f.write('          1 RefPtID - Index of the reference point (1-nPoints)\n')
        f.write('     Pointyi Pointzi ! nPoints listed in order of increasing height\n')
        f.write('       (m)     (m)\n')
        f.write('       {0}   {1}\n'.format(yloc,zloc))
        f.write('--------Time Series-------------------------------------------------------------\n')
        f.write('Elapsed Time            Point01u                Point01v                Point01w\n')
        f.write('       (s)             (m/s)                   (m/s)                   (m/s)\n')
        for i in range(time.shape[0]):
            f.write('{:.2f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(time[i],u[i],v[i],w[i]))

nSeeds=1
yloc = 0.0
zloc = HubHt
for case in ['Case00']:#caseNames:
    xlocs=['None']*nTurbs
    ylocs=['None']*nTurbs
    for wt in range(nTurbs):
        xlocs[wt]=float(ParamVals[case][3+wt*2])
        ylocs[wt]=float(ParamVals[case][4+wt*2])
    for seed in range(nSeeds):
        abspath = workdir.format(case,seed)
        
        TSpathLow  = abspath+'\Low'
        print(TSpathLow)
        IFdata = stochasticTurbulence(D,prefix=TSpathLow)
        IFdata.readBTS('.',HubHt)
        meanu = IFdata.u[:,IFdata.jHub,IFdata.kHub].mean()
        Width = IFdata.dY*(IFdata.nY-1)
        lowTime = np.arange(0, IFdata.nSeconds, IFdata.dT)
        
        for wt in range(nTurbs):
            tstart=int(xlocs[wt]/meanu/IFdata.dT+3.0*D) ## This 3D is based on the default in this scripts to start the FFarm domain 3D upstream of the most upstream turbine. Modify if needed

            tmp = lowTime.shape[0]-tstart

            TurbLoc_rel=ylocs[wt]+Width/2.

            fileOut = abspath+'\\USRTimeSeries_T{0}.txt'.format(wt)

            uvel = np.zeros(lowTime.shape[0])
            vvel = np.zeros(lowTime.shape[0])
            wvel = np.zeros(lowTime.shape[0])

            uvel[:tmp] = IFdata.u[tstart:,IFdata.y2j(TurbLoc_rel),IFdata.kHub]
            vvel[:tmp] = IFdata.v[tstart:,IFdata.y2j(TurbLoc_rel),IFdata.kHub]
            wvel[:tmp] = IFdata.w[tstart:,IFdata.y2j(TurbLoc_rel),IFdata.kHub]

            uvel[tmp:] = IFdata.u[:tstart,IFdata.y2j(TurbLoc_rel),IFdata.kHub]
            vvel[tmp:] = IFdata.v[:tstart,IFdata.y2j(TurbLoc_rel),IFdata.kHub]
            wvel[tmp:] = IFdata.w[:tstart,IFdata.y2j(TurbLoc_rel),IFdata.kHub]


            writeTimeSeriesFile(fileOut,yloc,zloc,uvel,vvel,wvel,lowTime)
            
#createFFParamFiles(caseNames,nSeeds,ParamVals,HubHt,D,'FFarmParams.txt')