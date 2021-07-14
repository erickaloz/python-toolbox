def TSInputCreation(Path,totTurbs):

    import os,sys
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import sys
    sys.path.append('/home/kshaler/PostProcessing/General')
    from TurbSimCaseCreation import TSCaseCreation, WriteTSFile

    plt.rc("font",family="serif")
    plt.rc("font",size=14)

    zbot=1  ####This is the default value; can be changed if you wish
    ParamsFile = Path+'/TSParams.txt'

    Pars=pd.read_csv(ParamsFile, delimiter='\t',header=None,index_col=0,nrows=8,names=['Var','Val'])
    xlocs_tmp=pd.read_csv(ParamsFile, delimiter='\t',header=None,skiprows=6,index_col=0,nrows=1)
    ylocs_tmp=pd.read_csv(ParamsFile, delimiter='\t',header=None,skiprows=7,index_col=0,nrows=1)

    nTurbs = int(xlocs_tmp.shape[1])

    xlocs=['None']*nTurbs
    ylocs=['None']*nTurbs
    for turb in range(nTurbs):
        xlocs[turb] = xlocs_tmp[turb+1][0]
        ylocs[turb] = ylocs_tmp[turb+1][0]

    D=np.float(Pars['Val']['D'])
    RefHt=np.float(Pars['Val']['HubHt'])
    Vhub=np.float(Pars['Val']['Vhub'])
    TI=Pars['Val']['TI']
    Shear=Pars['Val']['Shear']

    CaseLow = TSCaseCreation('TurbSim_Low')
    CaseLow.Turb(D,RefHt,fmax=0.2) ### if fmax or cmax differ from NREL-5MW value, specify that here, e.g., CaseLow.Turb(D,RefHt,fmax=0.642,cmax=6.0)
    CaseLow.turbLocs(xlocs,ylocs)
    CaseLow.Inflow(Pars)

    CaseWT1 = TSCaseCreation('TurbSim_WT')
    CaseWT1.Turb(D,RefHt,fmax=0.2) ### if fmax or cmax differ from NREL-5MW value, specify that here, e.g., CaseLow.Turb(D,RefHt,fmax=0.642,cmax=6.0)
    CaseWT1.turbLocs(xlocs,ylocs)
    CaseWT1.Inflow(Pars)

    CaseLow.discretization_low(Vhub,Shear)
    CaseLow.domainSize_low(zbot)

    CaseWT1.discretization_high(Vhub,Shear)
    CaseWT1.domainSize_high(zbot)

    
    OrigFile = os.path.join(CaseLow.prefix+'.inp')
    NewFile = os.path.join('Low.inp')

    WriteTSFile(OrigFile, NewFile, CaseLow, 0,'low')
    
    for wt in range(totTurbs):  
        OrigFile = os.path.join(CaseWT1.prefix+'.inp')
        NewFile = os.path.join('HighT{0}.inp'.format(wt+1))
        WriteTSFile(OrigFile, NewFile, CaseWT1, wt,'high')

    fig = plt.figure(figsize=(6,5))
    ax  = fig.add_subplot(111,aspect="equal")

    xminLow = min(CaseLow.x)
    xmaxLow = max(CaseLow.x)+7*CaseLow.D

    yminLow = -CaseLow.Width/2
    ymaxLow = CaseLow.Width/2

    # TurbLocs
    for wt in range(nTurbs):
        ax.plot(CaseLow.x[wt],CaseLow.y[wt],'x',ms=8,mew=2,label="WT{0}".format(wt+1))

    # low-res box
    ax.plot([xminLow,xmaxLow,xmaxLow,xminLow,xminLow],
            [yminLow,yminLow,ymaxLow,ymaxLow,yminLow],'--k',lw=2,label='Low')

    # high-res boxes
    for wt in range(nTurbs):
        xmin = CaseWT1.x[wt]-0.6*CaseWT1.D
        xmax = CaseWT1.x[wt]+0.6*CaseWT1.D

        ymin = CaseWT1.y[wt]-CaseWT1.Width/2
        ymax = CaseWT1.y[wt]+CaseWT1.Width/2

        ax.plot([xmin,xmax,xmax,xmin,xmin],
                [ymin,ymin,ymax,ymax,ymin],'-',lw=2,label="WT{0}".format(wt+1))


    plt.legend(bbox_to_anchor=(1.05,1.015),frameon=False)

    ax.set_xlabel("x-location [m]")
    ax.set_ylabel("y-location [m]")

    fig.tight_layout
    fig.savefig('TSLayout.pdf',bbox_to_inches='tight',dpi=500)
    
def FFInputCreation(Path,tpath,TSpath,Vhub,Vhub_Low=0.0,Vhub_High=0.0):

    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    from stochasticTurbulenceTools_mod import stochasticTurbulence
    from FFarmCaseCreation import FFarmCaseCreation, WriteFFarmFile
    import matplotlib.pylab as pl
    import numpy as np
    
    plt.rc("font",family="serif")
    plt.rc("font",size=14)

    ParamsFile = Path+'/FFarmParams.txt'
    Pars=pd.read_csv(ParamsFile, delimiter='\t',header=None,index_col=0,nrows=4,names=['Var','Val'])
    xlocs_tmp=pd.read_csv(ParamsFile, delimiter='\t',header=None,skiprows=4,index_col=0,nrows=1)
    ylocs_tmp=pd.read_csv(ParamsFile, delimiter='\t',header=None,skiprows=5,index_col=0,nrows=1)
    zlocs_tmp=pd.read_csv(ParamsFile, delimiter='\t',header=None,skiprows=6,index_col=0,nrows=1)

    nTurbs = xlocs_tmp.shape[1]

    xlocs=['None']*nTurbs
    ylocs=['None']*nTurbs
    zlocs=['None']*nTurbs
    for turb in range(nTurbs):
        xlocs[turb] = xlocs_tmp[turb+1][0]
        ylocs[turb] = ylocs_tmp[turb+1][0]
        zlocs[turb] = zlocs_tmp[turb+1][0]

    D=Pars['Val']['D']
    HubHt=Pars['Val']['HubHt']
    high_extent_X=Pars['Val']['high_extent_X']
    high_extent_Y=Pars['Val']['high_extent_Y']

    Case = FFarmCaseCreation(nTurbs,'FFarm')
    Case.Turb(D,HubHt,tpath)

    ## If ReadTS=True, TurbSim parameters will be read directly from the TurbSim.inp or TurbSim.bts file
    ReadTS = True

    ## if ReadINP = True, TurbSim parameters will be read directly from the TurbSim.inp file. Otherwise, TurbSim.bts file will be used.
    ReadINP = False
    
    if ReadTS == True:
        if ReadINP == True:
            TSFile = os.path.join(TSpath+os.sep+Case.prefix+'.inp')

        else:
            IFdata_Low = stochasticTurbulence(Case.D,prefix='Low')
            IFdata_Low.readBTS(TSpath, HubHt)
            IFdata_High = stochasticTurbulence(Case.D,prefix='HighT1')
            IFdata_High.readBTS(TSpath, 148.5)
    else:
        print('TurbSim parameters must be entered directly.')

    Case.turbLocs(xlocs,ylocs,zlocs)
    Case.discretization(Vhub)
    Case.highResDomain(IFdata_High,Vhub_High,Case.dX_High_desired,high_extent_X,high_extent_Y)
    Case.lowResDomain(IFdata_Low,Case.dX_Low_desired,Vhub,Vhub_Low)

    OrigFile = os.path.join(Path+os.sep+Case.prefix+'.fstf')
    NewFile = os.path.join(Path+os.sep+Case.prefix+'_mod.fstf')
    WriteFFarmFile(OrigFile, NewFile, Case, NewFile=False)


    fig = plt.figure(figsize=(13.5,10))
    ax  = fig.add_subplot(111,aspect="equal")

    xmax_low = Case.X0_Low+Case.dX_Low*Case.nX_Low
    ymax_low = Case.Y0_Low+Case.dY_Low*Case.nY_Low

    colors=pl.cm.tab20b(np.linspace(0,1,12))
    
    # low-res box
    ax.plot([Case.X0_Low,xmax_low,xmax_low,Case.X0_Low,Case.X0_Low],
            [Case.Y0_Low,Case.Y0_Low,ymax_low,ymax_low,Case.Y0_Low],'--k',lw=2,label='Low')

    # high-res boxes
    for wt in range(Case.nTurbs):
        xmax_high = Case.X0_High[wt]+Case.dX_High*Case.nX_High
        ymax_high = Case.Y0_High[wt]+Case.dY_High*Case.nY_High
        ax.plot([Case.X0_High[wt],xmax_high,xmax_high,Case.X0_High[wt],Case.X0_High[wt]],
                [Case.Y0_High[wt],Case.Y0_High[wt],ymax_high,ymax_high,Case.Y0_High[wt]],
                '-',
                label="HighT{0}".format(wt+1),color=colors[wt])
        ax.plot(Case.x[wt],Case.y[wt],'x',ms=8,mew=2,label="WT{0}".format(wt+1),color=colors[wt])

    plt.legend(bbox_to_anchor=(1.05,1.015),frameon=False)

    ax.set_xlabel("x-location [m]")
    ax.set_ylabel("y-location [m]")

    fig.tight_layout
    fig.savefig(Path+'/FFarmLayout.pdf',bbox_to_inches='tight',dpi=500)
