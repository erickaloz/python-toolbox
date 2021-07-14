def copyTSNominalFiles(caseNames,nSeeds,nomDataPath,windName):
    import os, subprocess
    nCases = len(caseNames)
    for b in range(nCases):
        name2=caseNames[b]
        if not os.path.exists(name2):
            os.makedirs(name2)
        os.chdir(name2)
        for seed in range(0,nSeeds):
            name3='Seed_'+str(seed)
            if not os.path.exists(name3):
                os.makedirs(name3)
            os.chdir(name3)
            print(os.getcwd())

            NomWindFiles1=str(nomDataPath)+os.sep+str(windName)+'_'+str(name2)+'.txt'#'Turbsim_'+str(name2)+'.txt'
            print('WindFile: ',NomWindFiles1)
            subprocess.call(["cp", NomWindFiles1, str(windName)+'.txt'])
            os.chdir('../')
            
        os.chdir('../')

def copyFFNominalFiles(caseNames,nSeeds,nTurbs,nomDataPath,rootname):
    import os
    from shutil import copyfile
    nCases = len(caseNames)
    for b in range(nCases):
        name2=caseNames[b]
        if not os.path.exists(name2):
            os.makedirs(name2)
        os.chdir(name2)
        #src=str(nomDataPath)+os.sep+str(rootname)+'ElastoDyn.dat'
        #dst='.'+os.sep+str(rootname)+'ElastoDyn.dat'
        #copyfile(src,dst)

        for wt in range(nTurbs):
            #src=str(nomDataPath)+os.sep+str(rootname)+'ServoDyn.'+str(wt)+'.dat'
            #dst='.'+os.sep+str(rootname)+'ServoDyn.'+str(wt)+'.dat'
            #copyfile(src, dst)
            
            src=str(nomDataPath)+os.sep+str(rootname)+'WT'+str(wt+1)+'.fst'
            dst='.'+os.sep+str(rootname)+'WT'+str(wt+1)+'.fst'
            copyfile(src, dst)

        for seed in range(nSeeds):
            name3='Seed_'+str(seed)
            if not os.path.exists(name3):
                os.makedirs(name3)
            os.chdir(name3)
            src='../../TurbSim_mod.bts'
            copyfile(src,'./TurbSim_mod.bts')
            print(os.getcwd())
            os.chdir('../')
            
        os.chdir('../')
        
def createTSParamFiles(caseNames,nSeeds,data,D,HubHt,fileOut):
    import numpy as np
    import os, subprocess
    import pandas as pd
    import random
    from InputSetup import TSInputCreation as inputCreation
    
    nCases = len(caseNames)
    nTurbs=int(data[3:].shape[0]/2)
    seedNum = ['None']*nSeeds
    
    for seed in range(nSeeds):
        seedNum[seed] = random.randint(-2147483648,2147483647)
        
    for b in range(nCases):
        name2=caseNames[b]
        os.chdir(name2)
        
        xlocs=['None']*nTurbs
        ylocs=['None']*nTurbs
        for turb in range(nTurbs):
            xlocs[turb]=data[caseNames[b]][3+turb*2]
            ylocs[turb]=data[caseNames[b]][4+turb*2] 
        xmin=min(xlocs); xmax=max(xlocs)
        ymin=min(ylocs); ymax=max(ylocs)
        for seed in range(nSeeds):
            name3='Seed_'+str(seed)
            os.chdir(name3)
            print(os.getcwd())
            print('Writing TurbSim input file: ', fileOut)
            with open(fileOut, 'w') as f:
                f.write('Seed\t{0}\n'.format(seedNum[seed]))
                f.write('D\t{0}\n'.format(D))
                f.write('HubHt\t{0}\n'.format(HubHt))
                f.write('Vhub\t{0}\n'.format(data[name2][1]))
                f.write('TI\t{0}\n'.format(data[name2][0]))
                f.write('Shear\t{0}\n'.format(data[name2][2]))
                f.write('xlocs\t{0}\t{1}\n'.format(xmin,xmax))
                f.write('ylocs\t{0}\t{1}\n'.format(ymin,ymax))
            
            os.chdir('TurbSim')
            inputCreation('..',nTurbs)
            
            os.chdir('../..')
            
        os.chdir('../')


def createFFParamFiles(caseNames,nSeeds,data,HubHt,D,fileOut):
    import os, subprocess
    import pandas as pd
    
    import sys
    from InputSetup import FFInputCreation as inputCreation
    sys.path.append('/home/kshaler/PostProcessing/General')
    from stochasticTurbulenceTools_mod import stochasticTurbulence
    
    nCases = len(caseNames)
    nTurbs = int(data[3:].shape[0]/3)
    for b in range(nCases):
        name2=caseNames[b]
        os.chdir(name2)
        
        xlocs=['None']*nTurbs
        ylocs=['None']*nTurbs
        zlocs=['None']*nTurbs
        for turb in range(nTurbs):
            xlocs[turb]=data[caseNames[b]][turb*3]
            ylocs[turb]=data[caseNames[b]][1+turb*3]
            zlocs[turb]=data[caseNames[b]][2+turb*3]

        ReadTS = True

        ReadINP = False

        #if caseNames[b] == 'Case000':
        if ReadTS == True:
            if ReadINP == True:
                TSFile = os.path.join('.'+os.sep+Case.prefix+'.inp')

            else:
                print(os.getcwd())
                HubHt_HighTS = 148.5
                
                IFdata = stochasticTurbulence(D,prefix='Low')
                IFdata.readBTS('Seed_0/TurbSim', HubHt)

                IFdata.kHub = IFdata.z2k(HubHt)

                Vhub = IFdata.u[:,IFdata.jHub,IFdata.kHub].mean()
                
                IFdata = stochasticTurbulence(D,prefix='Low')
                HubHt_LowTS = IFdata.RefHt 
                IFdata.readBTS('Seed_0/TurbSim', HubHt_LowTS)

                IFdata.kHub = IFdata.z2k(HubHt_LowTS)

                Vhub_Low = IFdata.u[:,IFdata.jHub,IFdata.kHub].mean()
                
                IFdata = stochasticTurbulence(D,prefix='HighT1')
                IFdata.readBTS('Seed_0/TurbSim', HubHt_HighTS)

                IFdata.kHub = IFdata.z2k(HubHt_HighTS)

                Vhub_High = IFdata.u[:,IFdata.jHub,IFdata.kHub].mean()
        else:
            print('TurbSim parameters must be entered directly.')

        for seed in range(nSeeds):
            name3='Seed_'+str(seed)
            os.chdir(name3)
            print(os.getcwd())
            print('Writing FFarm input file: ', fileOut)
            with open(fileOut, 'w') as f:
                f.write('D\t{0}\n'.format(D))
                f.write('HubHt\t{0}\n'.format(HubHt))
                f.write('high_extent_X\t1.2\n')
                f.write('high_extent_Y\t1.2\n')
                f.write('xlocs')
                for turb in range(nTurbs):
                    f.write('\t{:.3f}'.format(xlocs[turb]))
                f.write('\n')
                f.write('ylocs')
                for turb in range(nTurbs):
                    f.write('\t{:.3f}'.format(ylocs[turb]))
                f.write('\n')
                f.write('zlocs')
                for turb in range(nTurbs):
                    f.write('\t{:.3f}'.format(zlocs[turb]))
                    
            inputCreation('.','../NREL5MW.T','TurbSim/',Vhub,Vhub_Low,Vhub_High)
            
            os.chdir('../')
            
        os.chdir('../')
        
def copyNominalFiles(NumSeeds,NumTurbs):
    import os
    import subprocess
    
    NomDataPath='/scratch/kshaler/FFarmSetup/Nominal'
    TurbineDataPath=os.getcwd()+'/.'
    
    for wt in range(NumTurbs):
        fstFile=NomDataPath+'/NREL5MW_ElastoDyn.{0}.dat'.format(wt)
        subprocess.call(["cp", fstFile, TurbineDataPath])

        fstFile=NomDataPath+'/NREL5MW_ServoDyn.{0}.dat'.format(wt)
        subprocess.call(["cp", fstFile, TurbineDataPath])

    for seed in range(NumSeeds):
        name='Seed_'+str(seed)
        if not os.path.exists(name):
            os.makedirs(name)
        os.chdir(name)
        print(os.getcwd())

        TurbineDataPath=os.getcwd()+'/.'

        os.makedirs('TurbSim')

        for wt in range(NumTurbs):
            fstFile=NomDataPath+'/NREL5MW_WT{0}.fst'.format(wt+1)
            subprocess.call(["cp", fstFile, TurbineDataPath])

        fstFile=NomDataPath+'/FFarm.fstf'
        subprocess.call(["cp", fstFile, TurbineDataPath])

        fstFile=NomDataPath+'/InflowWind.dat'
        subprocess.call(["cp", fstFile, TurbineDataPath])

        os.chdir('../')
        
