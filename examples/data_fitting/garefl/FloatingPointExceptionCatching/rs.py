#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Reflscript for ga_refl 9-Dec-2010 F.H.
# Modified 4-Sep-2012 F.H.

#issues: ISIS data files cannot contain headerline without specifier as they normally do, this will
#        lead to MC not starting up

import math, os, random, re, sys, string, subprocess, time, operator
import scipy, scipy.stats, scipy.special, shutil, numpy

class CFitInteractor():
    def __init__(self):
        pass

    #The sparse parameter loads only a fraction of the data, if sparse is larger or equal than 1 this equates to the
    #number of lines loaded. If sparse is smaller than one this equates to a probability that any line will be stored
    def fnLoadSingleColumns(self, sFileName, data=None, exceptions=None, header=True, headerline=None, LoadList=None,
                            sparse=0):

        file  = open(sFileName,"r")
        content=file.readlines()
        file.close()

        if data is None:
            data={}

        if header:                                                                  #if headerline in file, read it from there
            splitheaderline=content[0].split()
            content=content[1:]
        else:                                                                       #else use the one given as an attribute
            splitheaderline=headerline

        for i,columnname in enumerate(splitheaderline):
            if LoadList is None or (columnname in LoadList):
                data[columnname]=[]

        j=0
        random.seed()
        for line in content:
            if sparse==0 or (sparse>=1 and j<sparse) or (sparse<1 and random.random()<sparse):
                splitline=line.split()
                if splitline:
                    for i,column in enumerate(splitline):
                        if LoadList is None or (splitheaderline[i] in LoadList):
                            try:
                                data[splitheaderline[i]].append(float(splitline[i]))
                            except:
                                data[splitheaderline[i]].append(splitline[i])
                j+=1

        return data


    #Loads single row data into a dictionary ["rowname",'Values',[data]]
    #it appends the data found in the file to the provided list 'data'
    #it will skip rows with names which are either in the exception list
    #or appends data with name extension "_2nd" that are already present in the data list
    #The sparse parameter loads only a fraction of the data, if sparse is larger or equal than 1 this equates to the
    #number of lines loaded. If sparse is smaller than one this equates to a probability that any line will be stored
    def fnLoadSingleColumnsIntoStatDict(self, sFileName, data=None, exceptions=None, header=True, headerline=None,
                                        LoadList=None, sparse=0):

        if not data: data = {}
        file  = open(sFileName,"r")
        content=file.readlines()
        file.close()

        if data is None:
            data={}

        if header:                                                                  #if headerline in file, read it from there
            splitheaderline=content[0].split()
            content=content[1:]
        else:                                                                       #else use the one given as an attribute
            splitheaderline=headerline

        for i,columnname in enumerate(splitheaderline):
            if LoadList is None or (columnname in LoadList):
                data[columnname]={}
                data[columnname]['Values']=[]

        j=0
        random.seed()
        for line in content:
            if sparse==0 or (sparse>=1 and j<sparse) or (sparse<1 and random.random()<sparse):
                splitline=line.split()
                if splitline:
                    for i,column in enumerate(splitline):
                        if LoadList is None or (splitheaderline[i] in LoadList):
                            data[splitheaderline[i]]['Values'].append(float(splitline[i]))
                j+=1


        return data

    def fnLoadSingleRows(self, sFileName, data=None, exceptions=None):

        file  = open(sFileName,"r")
        content=file.readlines()
        file.close()

        if data is None:
            data={}

        for line in content:
            splitline=line.split()
            if splitline[0] not in exceptions:
                data[splitline[0]]=[]
                for entry in range(1, len(splitline)):
                    data[splitline[0]].append(splitline[entry])

        return data


    #saves all data out to a file
    def fnSaveSingleColumns(self, sFilename, data):

        file=open(sFilename,"w")

        for element in data:
            file.write(element+" ")
        file.write("\n")

        for i in range(len(data[data.keys()[0]])):
            for element in data:
                file.write(str(data[element][i])+" ")
            file.write("\n")

        file.close()
        #saves all data out to a file

    def fnSaveSingleColumnsFromStatDict(self, sFilename, data):

        file=open(sFilename,"w")

        for element in data:
            file.write(element+" ")
        file.write("\n")

        for i in range(len(data[data.keys()[0]]['Values'])):
            for element in data:
                file.write(str(data[element]['Values'][i])+" ")
            file.write("\n")

        file.close()


#Garefl methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
#error analysis cannot be found.
#The MCMC directory is called 'MCMC'

class CGaReflInteractor(CFitInteractor):
    def __init__(self):
        CFitInteractor.__init__(self)
        pass

    def fnLoadStatData(self,dSparse=0):                                          #Load a file like iSErr or SErr into
                                                                                #the self.liStatResult object
        sStatResultHeader=''
        liStatResult=[]
        sFileName='isErr.dat'

        try:
            if os.path.isfile('isErr.dat') and os.path.isfile('sErr.dat'):
                print '-------------------------------'
                print 'Found isErr.dat and sErr.dat ?!'
                print 'Load isErr.dat as default.'
                print '-------------------------------'
            elif os.path.isfile('isErr.dat'):                                     #check which type of MC output present
                print 'Found isErr.dat\n'
            elif os.path.isfile('sErr.dat'):
                print 'Found sErr.dat\n'
                sFileName='sErr.dat'

            diStatRawData= {'Parameters': self.fnLoadSingleColumnsIntoStatDict(sFileName,sparse=dSparse)}
            return diStatRawData

        except IOError:
            print 'Could not load '+sFileName+'. \n'
            sys.exit(1)





#Refl1D methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
#error analysis are found.
#The MCMC directory is called 'MCMC'
#The refl1d script name has to be run.py.

class CRefl1DInteractor(CFitInteractor):
    def __init__(self):
        CFitInteractor.__init__(self)
        pass

    #LoadStatResults returns a list of variable names, a logP array, and a numpy.ndarray
    #[values,var_numbers].

    def fnLoadStatData(self,dSparse=0):
        import dream.state
        state = dream.state.load_state("MCMC/run")
        state.mark_outliers() # ignore outlier chains

        #load Parameter
        data=self.fnLoadSingleColumns("MCMC/run.par",header=False,headerline=['parname','bestfitvalue'])
        lParName=[]; vars=[]
        i=0
        for parname in data['parname']:
            lParName.append(parname)
            vars.append(i)
            i+=1

        points, logp = state.sample(None,vars,None)

        diStatRawData={'Parameters': {}}
        diStatRawData['Parameters']['Chisq']={}                                #TODO: Work on better chisq handling
        diStatRawData['Parameters']['Chisq']['Values']=[]
        for parname in lParName:
            diStatRawData['Parameters'][parname]={}
            diStatRawData['Parameters'][parname]['Values']=[]

        random.seed()
        for j in range(len(points[:,0])):
            if dSparse==0 or (dSparse>1 and j<dSparse) or (dSparse<1 and random.random()<dSparse):
                diStatRawData['Parameters']['Chisq']['Values'].append(logp[j])
                for i,parname in enumerate(lParName):
                    if 'rho_' in parname or 'background' in parname:         #TODO: this is a hack because Paul does not scale down after scaling up
                        points[j,i]*=1E-6
                    diStatRawData['Parameters'][parname]['Values'].append(points[j,i])

        self.fnSaveSingleColumnsFromStatDict('sErr.dat',diStatRawData['Parameters'])


        return diStatRawData


class CReflPar:

    def __init__(self, cGaReflInteractor, cRefl1DInteractor):
        """

        """
        self.diParameters={}                                                    #Dictionary with all the parameters
    # self.diParameters data structure:
    # dictionary: sparameter : dictionary
    #                          'number'    : int                                # order of initialization in ga_refl
    #                          'lowerlimit'  : float                            # constraint range lower limit
    #                          'upperlimit'  : float                            # constraint range upper limit
    #                          'value'  : float                                 # absolute value
    #                          'error'  : float                                 # error derived from covariance matrix
    #                          'relval' : float                                 # relative data value between 0 and 1 in terms of the constraints
    #                          'variable: string                                # associated ga_refl variable
        self.diMolgroups={}
        self.diStatResults={}                                                   # dictionary of statistical results for various routines
                                                                                # from fnAnalyzeStatFile
                                                                                #   'Parameters' key contains dictionary
                                                                                #       parameter names are keys contain dictionary
                                                                                #           'Values' key contains ordered list of all MC parameters
                                                                                #           'LowPerc' contains lower percentile
                                                                                #           'Median' contains median
                                                                                #           'HighPerc' contains higher percentile
                                                                                #   'MaxParameterLength' contains length of longest parameter name
                                                                                # from fnRecreateStatFile
                                                                                #   'NumberOfStatValues' number of MC points
                                                                                #   'nSLDProfiles' key contains list of lists of lists
                                                                                #      of all profiles 'nSLDProfiles'[MCiteration][model][z,rho]
                                                                                # from fnLoadMolgroups
                                                                                #   'Molgroups' contains information about molecular groups
                                                                                #     of all iterations 'Molgroups'[MCiteration] {molgroupname {'zaxis' [z]
                                                                                #                                                               'areaxis' [area]
                                                                                #                                                               'nslaxis' [nsl]
                                                                                #                                                               'headerdata' {'property' value}
                                                                                #                                                               'property1' value
                                                                                #                                                               'property2' value}}

        self.liStatResult=[]                                                    #a list that contains the isErr.dat or sErr.dat file line by line
        self.sStatResultHeader=''                                               #headerline from isErr.dat or sErr.dat
        self.sStatFileName=''                                                   #Name of the statistical File
        self.chisq=0.
        self.fMolgroupsStepsize=0.
        self.iMolgroupsDimension=0
        self.fMolgroupsNormArea=0
                                                                                #check for system and type of setup file
        if os.path.exists("setup.c"):
            self.setupfilename="setup.c"
        else:
            self.setupfilename="setup.cc"

        self.cGaReflInteractor=cGaReflInteractor
        self.cRefl1DInteractor=cRefl1DInteractor



#-------------------------------------------------------------------------------

    def fnAnalyzeStatFile(self, fConfidence,sparse=0):                          #summarizes stat file

        self.fnLoadStatData(sparse)                                             #Load data from file into list
        self.fnLoadParameters()                                                 #Load Parameters for limits



        if fConfidence>1:
            fConfidence=1
        if fConfidence<0:
            fConfidence=scipy.special.erf(-1*fConfidence/scipy.sqrt(2))

        fLowerPercentileMark=100.0*(1-fConfidence)/2
        fHigherPercentileMark=(100-fLowerPercentileMark)

        iNumberOfMCIterations=self.diStatResults['NumberOfStatValues']                                                  #how many iterations already done
        iMaxParameterNameLength=self.diStatResults['MaxParameterLength']
        print 'Analysis of current MC simulation ...'
        print 'Number of iterations: %(ni)d' % {'ni':iNumberOfMCIterations}

        fMaxConvergence=0
        iHistoryLength=5

        for element in sorted(self.diParameters.keys(),key= lambda sParameter: self.diParameters[sParameter]['number']):

            fLowPerc=scipy.stats.scoreatpercentile\
            (self.diStatResults['Parameters'][element]['Values'],fLowerPercentileMark)           #Calculate Percentiles
            fMedian =scipy.stats.scoreatpercentile\
            (self.diStatResults['Parameters'][element]['Values'],50.)
            fHighPerc=scipy.stats.scoreatpercentile\
            (self.diStatResults['Parameters'][element]['Values'],fHigherPercentileMark)

            self.diStatResults['Parameters'][element]['LowPerc']=fLowPerc
            self.diStatResults['Parameters'][element]['Median']=fMedian
            self.diStatResults['Parameters'][element]['HighPerc']=fHighPerc

            flowerlimit=self.diParameters[element]['lowerlimit']
            fupperlimit=self.diParameters[element]['upperlimit']
            temp=abs(fupperlimit-flowerlimit)

            sGraphOutput='['
            itemp1=int((fLowPerc-flowerlimit)/temp*10+0.5)
            itemp2=int((fMedian-flowerlimit)/temp*10+0.5)
            itemp3=int((fHighPerc-flowerlimit)/temp*10+0.5)

            for i in range(11):
                s1=' '
                if itemp1==i or itemp3==i: s1='|'
                if itemp2==i:
                    if s1=='|': s1='+'
                    else: s1='-'
                sGraphOutput += s1
            sGraphOutput += ']'

            if (fLowPerc-flowerlimit)<temp*0.01:
                self.diStatResults['Parameters'][element]['LowerLimitCollision']=True
                sGraphOutput='#'+sGraphOutput[1:]
            else:
                self.diStatResults['Parameters'][element]['LowerLimitCollision']=False

            if (fupperlimit-fHighPerc)<temp*0.01:
                self.diStatResults['Parameters'][element]['UpperLimitCollision']=True
                sGraphOutput=sGraphOutput[:-1]+'#'
            else:
                self.diStatResults['Parameters'][element]['UpperLimitCollision']=False


            if iNumberOfMCIterations > iHistoryLength:                          #at least five iterations

                fLowPercHistory=[]
                fMedianHistory=[]
                fHighPercHistory=[]

                for i in range((iHistoryLength-1)*(-1),0):
                    fLowPercHistory.append(scipy.stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:i],fLowerPercentileMark))
                    fMedianHistory.append(scipy.stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:i],50))
                    fHighPercHistory.append(scipy.stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:i],fHigherPercentileMark))

                fLowPercHistory.append(fLowPerc)
                fMedianHistory.append(fMedian)
                fHighPercHistory.append(fHighPerc)

                fLowPercAverage=numpy.average(fLowPercHistory)
                fMedianAverage=numpy.average(fMedianHistory)
                fHighPercAverage=numpy.average(fHighPercHistory)

                fSigma=[]
                for i in range(iHistoryLength):
                    fSigma.append(abs(fLowPercHistory[i]-fLowPercAverage)/temp)
                    fSigma.append(abs(fMedianHistory[i]-fMedianAverage)/temp)
                    fSigma.append(abs(fHighPercHistory[i]-fHighPercAverage)/temp)

                fMaxConvergence=max(fMaxConvergence,max(fSigma))

                fLowPercConv=abs(fLowPerc-scipy.stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:-1],fLowerPercentileMark))/temp
                fMedianConv=abs(fMedian-scipy.stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:-1],50))/temp
                fHighPercConv=abs(fHighPerc-scipy.stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:-1],fHigherPercentileMark))/temp

                self.diStatResults['Parameters'][element]['LowPercConv']=fLowPercConv
                self.diStatResults['Parameters'][element]['MedianConv']=fMedianConv
                self.diStatResults['Parameters'][element]['HighPercConv']=fHighPercConv

            else:
                fMaxConvergence=1e6
                fLowPercConv=fMedianConv=fHighPercConv=0

            sPrintString='%(el)'
            sPrintString += str(iMaxParameterNameLength)
            sPrintString += 's  %(sg)s  [%(ll)10.4g,%(ul)10.4g]  [%(lp)10.4g(%(lpc).3f), %(m)10.4g(%(mc).3f), ' \
                            '%(hp)10.4g(%(hpc).3f)] (-%(ld)10.4g, +%(hd)10.4g)'

            print  sPrintString % \
            {'el':element,'ll': flowerlimit, 'ul':fupperlimit,'lp': fLowPerc,
             'lpc': fLowPercConv,'ld': (fMedian-fLowPerc), 'm': fMedian,
             'mc': fMedianConv, 'hd': (fHighPerc-fMedian), 'hp':fHighPerc,
             'hpc': fHighPercConv, 'sg':sGraphOutput}

        self.diStatResults['Convergence']=fMaxConvergence
        print 'Maximum deviation from average over last %(iHL)d iterations: %(maxc).4f' % \
        {'iHL':iHistoryLength, 'maxc':fMaxConvergence}
        print 'Confidence level: %(fCL).4f' % {'fCL':fConfidence}





#-------------------------------------------------------------------------------

    def fnBackup(self):
        if not os.path.isdir('rsbackup'):                                       #create backup dir
            pr=subprocess.Popen(["mkdir",     "rsbackup"])
            pr.wait()
        pr=subprocess.Popen(["cp","pop.dat",  "rsbackup/"])                       #save all data that will
        pr.wait()
        pr=subprocess.Popen(["cp","par.dat",  "rsbackup/"])                       #change during fit
        pr.wait()
        pr=subprocess.Popen(["cp","covar.dat","rsbackup/"])
        pr.wait()
        subprocess.call("cp fit*.dat rsbackup/", shell=True)
        pr=subprocess.Popen(["cp", "fit", "rsbackup/"])
        pr.wait()
        subprocess.call("cp model*.dat rsbackup/", shell=True)
        pr=subprocess.Popen(["cp", "pop_bak.dat", "rsbackup/"])
        pr.wait()
        subprocess.call("cp profile*.dat rsbackup/", shell=True)
        subprocess.call("cp "+self.setupfilename+" rsbackup/", shell=True)
        pr=subprocess.Popen(["cp", "setup.o", "rsbackup/"])
        pr.wait()

#-------------------------------------------------------------------------------

    def fnCalculateMolgroupProperty(self,fConfidence):

        import scipy.stats

        try:
            self.diStatResults=self.fnLoadObject('StatDataPython.dat')

        except IOError:
            print 'Failure to calculate values from StatDataPython.dat.'
            print 'Recreate statistical data from sErr.dat.'
            self.fnRecreateStatistical()


        diResults= {}
        sMolgroups=self.diStatResults['Molgroups'][0].keys()
        for sMolgroup in sMolgroups:                                                                                    #create results for individual molgroups
            diResults[sMolgroup+'_COM']=[]
            diResults[sMolgroup+'_INT']=[]
            diResults[sMolgroup+'_AVG']=[]

        for i in range(len(self.diStatResults['Molgroups'])):                                                           #cycle over all MC iterations
            mgdict=self.diStatResults['Molgroups'][i]
            for sMolgroup in mgdict.keys():                                                                             #cycle over all molgroups

                if sum(mgdict[sMolgroup]['areaaxis']):
                    fCOM,fInt=numpy.average(mgdict[sMolgroup]['zaxis'],weights=mgdict[sMolgroup]['areaaxis'],returned=True) #Calculate Center of Mass and Integral
                    fInt=numpy.sum(mgdict[sMolgroup]['areaaxis'])*(mgdict[sMolgroup]['zaxis'][1]-mgdict[sMolgroup]['zaxis'][0])
                else:
                    fCOM=1E5; fInt=0                                                                                        #total integral in volume per surface area (unit: Å)
                                                                                                                        #taking into account grid spacing of z-axis
                fInt=fInt/mgdict['normarea']['areaaxis'][0]

                for j in range(len(mgdict[sMolgroup]['areaaxis'])):
                    if mgdict[sMolgroup]['areaaxis'][j]:
                        fCOM=fCOM-mgdict[sMolgroup]['zaxis'][j]                                                         #normalize COM to start of molecular group
                        break

                fAvg=numpy.average(mgdict[sMolgroup]['areaaxis'])

                diResults[sMolgroup+'_COM'].append(fCOM)
                diResults[sMolgroup+'_INT'].append(fInt)
                diResults[sMolgroup+'_AVG'].append(fAvg)

                                                                                                                        #percentage water in sub-membrane space for tBLM
            if 'tether' in sMolgroups and 'tetherg' in sMolgroups \
                and 'normarea' in sMolgroups and 'bME' in sMolgroups \
                and 'vf_bilayer' in sMolgroups:

                fVolSubmembrane=mgdict['normarea']['areaaxis'][0]*(float(mgdict['tether']['headerdata']['l'])+
                 float(mgdict['tetherg']['headerdata']['l']))*self.diStatResults['Parameters']['vf_bilayer']['Values'][i]
                fVolComponents=float(mgdict['bME']['headerdata']['vol'])*float(mgdict['bME']['headerdata']['nf'])+\
                               float(mgdict['tether']['headerdata']['vol'])*float(mgdict['tether']['headerdata']['nf'])+\
                               float(mgdict['tetherg']['headerdata']['vol'])*float(mgdict['tetherg']['headerdata']['nf'])+\
                               float(mgdict['headgroup1']['headerdata']['vol'])*float(mgdict['headgroup1']['headerdata']['nf'])
                j=2                                                                                                     #look for other headgroups
                while True:
                    if 'headgroup1_'+str(j) in sMolgroups:
                        fVolComponents+=float(mgdict['headgroup1_'+str(j)]['headerdata']['vol'])*\
                                        float(mgdict['headgroup1_'+str(j)]['headerdata']['nf'])
                    else:
                        break
                    j+=1
                fWaterFraction=1-fVolComponents/fVolSubmembrane
                if not 'WaterFracSubmembrane' in diResults.keys():
                    diResults['WaterFracSubmembrane']=[]
                diResults['WaterFracSubmembrane'].append(fWaterFraction)

                fTotalTetherLength=float(mgdict['tether']['headerdata']['l'])+float(mgdict['tetherg']['headerdata']['l'])
                if not 'fTotalTetherLength' in diResults.keys():
                    diResults['fTotalTetherLength']=[]
                diResults['fTotalTetherLength'].append(fTotalTetherLength)

                fTotalLipid1Length=float(mgdict['lipid1']['headerdata']['l'])+float(mgdict['methyl1']['headerdata']['l'])
                if not 'fTotalLipid1Length' in diResults.keys():
                    diResults['fTotalLipid1Length']=[]
                diResults['fTotalLipid1Length'].append(fTotalLipid1Length)

                fTotalLipid2Length=float(mgdict['lipid2']['headerdata']['l'])+float(mgdict['methyl2']['headerdata']['l'])
                if not 'fTotalLipid2Length' in diResults.keys():
                    diResults['fTotalLipid2Length']=[]
                diResults['fTotalLipid2Length'].append(fTotalLipid2Length)

                fAreaPerLipid2=float(mgdict['lipid2']['headerdata']['vol'])/float(mgdict['lipid2']['headerdata']['l'])
                if not 'fAreaPerLipid2' in diResults.keys():
                    diResults['fAreaPerLipid2']=[]
                diResults['fAreaPerLipid2'].append(fAreaPerLipid2)

                fTetherDensity=float(mgdict['tether']['headerdata']['nf'])/mgdict['normarea']['areaaxis'][0]
                if not 'fTetherDensity' in diResults.keys():
                    diResults['fTetherDensity']=[]
                diResults['fTetherDensity'].append(fTetherDensity)


                #output
        if fConfidence>1:
            fConfidence=1
        if fConfidence<0:
            fConfidence=scipy.special.erf(-1*fConfidence/scipy.sqrt(2))

        fLowerPercentileMark=100.0*(1-fConfidence)/2
        fHigherPercentileMark=(100-fLowerPercentileMark)

        File=open("CalculationResults.dat","w")
        for element in diResults:

            fLowPerc=scipy.stats.scoreatpercentile(diResults[element],fLowerPercentileMark)                             #Calculate Percentiles
            fMedian =scipy.stats.scoreatpercentile(diResults[element],50.)
            fHighPerc=scipy.stats.scoreatpercentile(diResults[element],fHigherPercentileMark)

            sPrintString='%(el)s'
            sPrintString += '  [%(lp)10.4g, %(m)10.4g, %(hp)10.4g] (-%(ld)10.4g, +%(hd)10.4g)'

            File.write(sPrintString % \
            {'el':element,'lp': fLowPerc,'ld': (fMedian-fLowPerc), 'm': fMedian,'hd': (fHighPerc-fMedian), 'hp':fHighPerc})
            File.write('\n')

        File.close()


#-------------------------------------------------------------------------------

    def fnCheckFit(self):                                                       #checks fit directory for certain conditions

        if os.path.isfile('isErr.dat'):                                         #check which type of MC output present
            self.sStatFileName='isErr.dat'
        elif os.path.isfile('isErr.dat'):
            self.sStatFileName='sErr.dat'
        else:
            self.sStatFileName=''

        try:
            self.fnLoadStatData()                                                #Load data from file into list
            self.fnLoadParameters()                                              #Load Parameters for limits

            if sorted(self.diParameters.keys())!=sorted(self.diStatResults['Parameters'].keys()):
                print '----------------------------------------------------'
                print '----------------------------------------------------'
                print 'setup.c and Stat File do not agree -----------------'
                print 'backing up Stat File -------------------------------'
                print '----------------------------------------------------'
                print '----------------------------------------------------'

                i=2
                while True:
                    if os.path.isfile(self.sStatFileName[:-4]+str(i)[1:]+".bak"):
                        pass
                    else:
                        pr=subprocess.Popen(["mv ",self.sStatFileName,self.sStatFileName[:-4]+str(i)[1:]+".bak"])
                        pr.wait()
                        self.sStatFileName=''
                        self.liStatResult=[]
                        self.sStatResultHeader=''
                        break

                    i+=1

        except IOError:
            pass




#-------------------------------------------------------------------------------

    def fnContourData(self,sname,dZGrid,dRhoGrid, dAreaGrid=0.1):

        def fnMap2Array(liArray,liDimensions,dx1,dy1,dx2,dy2):                  #map nSLD profile onto array
                                                                                #by drawing lines between (dz1,drho1)
                                                                                #and (dz2,drho2)

            if not liArray:                                                     #if array is empty create with first point
                liArray.append([0])
                liDimensions[0]=dx1
                liDimensions[1]=dx1
                liDimensions[3]=dy1
                liDimensions[4]=dy1


            while dx1>liDimensions[1]:                                          #check if array has to be extended
                for i in range(len(liArray)):                                   #and extend it if necessary
                    liArray[i].append(0)
                liDimensions[1]=liDimensions[1]+liDimensions[2]
            while dx1<liDimensions[0]:
                for i in range(len(liArray)):
                    liArray[i].insert(0,0)
                liDimensions[0]=liDimensions[0]-liDimensions[2]
            while dy1>liDimensions[4]:
                li=[0]*len(liArray[0])
                liArray.append(li)
                liDimensions[4]=liDimensions[4]+liDimensions[5]
            while dy1<liDimensions[3]:
                li=[0]*len(liArray[0])
                liArray.insert(0,li)
                liDimensions[3]=liDimensions[3]-liDimensions[5]
            while dx2>liDimensions[1]:
                for i in range(len(liArray)):
                    liArray[i].append(0)
                liDimensions[1]=liDimensions[1]+liDimensions[2]
            while dx2<liDimensions[0]:
                for i in range(len(liArray)):
                    liArray[i].insert(0,0)
                liDimensions[0]=liDimensions[0]-liDimensions[2]
            while dy2>liDimensions[4]:
                li=[0]*len(liArray[0])
                liArray.append(li)
                liDimensions[4]=liDimensions[4]+liDimensions[5]
            while dy2<liDimensions[3]:
                li=[0]*len(liArray[0])
                liArray.insert(0,li)
                liDimensions[3]=liDimensions[3]-liDimensions[5]


            iXStart=int((dx1-liDimensions[0])/liDimensions[2])                  #calculate coordinates from (z,rho)
            iXEnd=int((dx2-liDimensions[0])/liDimensions[2])                    #data points
            iYStart=int((dy1-liDimensions[3])/liDimensions[5])
            iYEnd=int((dy2-liDimensions[3])/liDimensions[5])

            diffX=iXEnd-iXStart                                                 #how many indizes do go?
            diffY=iYEnd-iYStart

            if abs(diffX)>abs(diffY):                                           #which direction more steps?
                diff=abs(diffX)
            else:
                diff=abs(diffY)

            if diff==0:                                                         #treat single point differently because
                liArray[iYStart][iXStart]=liArray[iYStart][iXStart]             #of division by zero error -> do nothing!
            else:
                fStepX=float(diffX)/float(diff)                                 #calculate stepsize for each direction
                fStepY=float(diffY)/float(diff)

                i=0                                                             #draw by increasing field occupancy
                while i<diff:                                                   #and thus create surface/contour plot
                    iX=iXStart+int(round(fStepX*float(i),0))
                    iY=iYStart+int(round(fStepY*float(i),0))
                    liArray[iY][iX] += 1
                    i += 1



        self.fnRecreateStatistical()                                            #Load Parameters and profiles from stat data
        iNumberOfModels=self.fnGetNumberOfModelsFromSetupC()                    #how many profiles to analyze?


        liContourArray=[]                                                       #initiate array
        liContourArrayDimensions=[]                                             #array dimensions
        i=0                                                                     #result for contour plot of profiles
        while i<iNumberOfModels:                                                #list of lists for each profile
            liContourArray=liContourArray+[[]]
            liContourArrayDimensions=liContourArrayDimensions+[[0.,0.,dRhoGrid,0.,0.,dZGrid]]
                                                                                #dimension format ymin, ymax, ystep, xmin, xmax, xstep
            i=i+1

        for iteration in self.diStatResults['nSLDProfiles']:                    #cycle through all individual stat. results
                i=0
                for model in iteration:                                         #cycle through all models
                    for l in range(len(model[0])):                              #extract nSLD profile data point by point
                        dz  =round(model[0][l]/dZGrid)*dZGrid                   #round to grid precision
                        drho=round(model[1][l]/dRhoGrid)*dRhoGrid               #round nSLD to grid precision
                        if l!=0:
                            fnMap2Array(liContourArray[i],
                                liContourArrayDimensions[i],
                                drhoold,dzold,drho, dz)

                        dzold=dz
                        drhoold=drho
                    i += 1


        i=0
        print 'Processing data for output ...'
        while i<iNumberOfModels:                                                #loop through all models

            print 'Model %i: %i x %i' % (i, len(liContourArray[i][0]),
                len(liContourArray[i]))
            sFileName='Cont_nSLD_Array'+str(i)+'.dat'                                 #write out array
            file = open(sFileName,"w")
            for line in liContourArray[i]:
                sLine=''
                for element in line:
                    sLine = sLine+str(element)+' '
                sLine = sLine + '\n'
                file.write(sLine)
            file.close()

            dZMin=liContourArrayDimensions[i][3]
            dZMax=liContourArrayDimensions[i][4]
            dZStep=liContourArrayDimensions[i][5]
            dRhoMin=liContourArrayDimensions[i][0]
            dRhoMax=liContourArrayDimensions[i][1]
            dRhoStep=liContourArrayDimensions[i][2]

            dZ=dZMin                                                            #write out x-dimension wave
            sFileName='Cont_nSLD_DimZ'+str(i)+'.dat'                                  #dimension wave has one point extra for Igor
            file = open(sFileName,"w")
            while (dZ<=dZMax+dZStep):
                sLine=str(round(dZ/dZGrid)*dZGrid)+'\n'
                file.write(sLine)
                dZ=dZ+dZStep
            file.close()

            dRho=dRhoMin                                                        #write out y-dimension wave
            sFileName='Cont_nSLD_DimRho'+str(i)+'.dat'                                #dimension wave has one point extra for Igor
            file = open(sFileName,"w")
            while (dRho<=dRhoMax+dRhoStep):
                sLine=str(round(dRho/dRhoGrid)*dRhoGrid)+'\n'
                file.write(sLine)
                dRho=dRho+dRhoStep
            file.close()


            i=i+1

        if 'Molgroups' in self.diStatResults:

            liContourArray=[]														      #initiate array
            liContourArrayDimensions=[]												      #array dimensions
            for _ in self.diStatResults['Molgroups'][0]:					              #iterate through molecular group names
                liContourArray=liContourArray+[[]]
                liContourArrayDimensions=liContourArrayDimensions+[[0.,0.,dAreaGrid,0.,0.,dZGrid]]
                                                                                          #dimension format ymin, ymax, ystep, xmin, xmax, xstep

            for iteration in self.diStatResults['Molgroups']:					          #cycle through all individual stat. results
                    i=0;dareaold=0;dzold=0
                    for molgroup in iteration:											  #cycle through all molecular groups
                        for l in range(len(iteration[molgroup]['zaxis'])):							  #extract area profile data point by point
                            dz	=round(iteration[molgroup]['zaxis'][l]/dZGrid)*dZGrid					#round to grid precision
                            darea=round(iteration[molgroup]['areaaxis'][l]/dAreaGrid)*dAreaGrid				#round nSLD to grid precision
                            if l!=0:
                                fnMap2Array(liContourArray[i],
                                    liContourArrayDimensions[i],
                                    dareaold,dzold,darea,dz)

                            dzold=dz
                            dareaold=darea
                        i+=1

            i=0
            for molgroup in self.diStatResults['Molgroups'][0]:							  #loop through all models

                print '%s %i: %i x %i' % (molgroup,i,len(liContourArray[i][0]),
                    len(liContourArray[i]))
                sFileName='Cont_'+molgroup+'_Array'+'.dat'						  #write out array
                file = open(sFileName,"w")
                for line in liContourArray[i]:
                    sLine=''
                    for element in line:
                        sLine = sLine+str(element)+' '
                    sLine = sLine + '\n'
                    file.write(sLine)
                file.close()

                dZMin=liContourArrayDimensions[i][3]
                dZMax=liContourArrayDimensions[i][4]
                dZStep=liContourArrayDimensions[i][5]
                dAreaMin=liContourArrayDimensions[i][0]
                dAreaMax=liContourArrayDimensions[i][1]
                dAreaStep=liContourArrayDimensions[i][2]

                dZ=dZMin															      #write out x-dimension wave
                sFileName='Cont_'+molgroup+'_DimZ'+'.dat'							  #dimension wave has one point extra for Igor
                file = open(sFileName,"w")
                while (dZ<=dZMax+dZStep):
                    sLine=str(round(dZ/dZGrid)*dZGrid)+'\n'
                    file.write(sLine)
                    dZ=dZ+dZStep
                file.close()

                dArea=dAreaMin														      #write out y-dimension wave
                sFileName='Cont_'+molgroup+'_DimArea'+'.dat'						  #dimension wave has one point extra for Igor
                file = open(sFileName,"w")
                while (dArea<=dAreaMax+dAreaStep):
                    sLine=str(round(dArea/dAreaGrid)*dAreaGrid)+'\n'
                    file.write(sLine)
                    dArea=dArea+dAreaStep
                file.close()


                i=i+1

#-------------------------------------------------------------------------------


    def fnGetChiSq(self):                                                       #export chi squared
        return self.chisq
#-------------------------------------------------------------------------------
    def fnGetNumberOfModelsFromSetupC(self):

        file = open(self.setupfilename,"r")                                     #open setup.c
        data = file.readlines()
        file.close()
        smatch=re.compile(r'define\s+MODELS\s+(.+?)\n',re.IGNORECASE | re.VERBOSE)
        for line in data:                                                       #search through setup.c
            if smatch.search(line):                                             #searching for MODELS constant
                i=smatch.search(line).group(1)
                return(int(i))
        return(0)

#-------------------------------------------------------------------------------

    def fnGetParameterValue(self,sname):                                        #export absolute parameter value
        return self.diParameters[sname]['value']                                #for given name
#-------------------------------------------------------------------------------

    def fnGetSortedParNames(self):                                              #return a list of sorted parameter
        litest=self.diParameters.keys()                                         #names, the sort criterium is defined
        litest.sort(self.fndiParametersNumberSort)                              #in fndiParametersNumberSort and it is
        return litest                                                           #the parameter number
#-------------------------------------------------------------------------------

    def fnGetSortedParValues(self):                                             #the same as above but it returns
        litest=self.diParameters.keys()                                         #a list of sorted parameter values
        litest.sort(self.fndiParametersNumberSort)
        lvalue=[]
        for parameter in litest:
            lvalue.append(str(self.diParameters[parameter]['value']))
        return lvalue
#-------------------------------------------------------------------------------

    def fnGetTaggedParameters(self):                                            #returns a list of the name and the
        file = open(self.setupfilename)                                              #range + stepsize information of parameters
        data = file.readlines()                                                 #which are tagged for displacement error
        file.close()                                                            #analysis
        output=[]
        for line in data:
            if '!rstag' in line:
                smatch=re.compile(r'pars_add\(pars.*?\"(.+?)\"\s*,.+?,(.+?),(.+?)\).+!rstag\s+(.+?)\s+(.+?)\s+(.+?)\s+!',
                                  re.IGNORECASE | re.VERBOSE)
                output.append(smatch.search(line).groups())
        return output
#-------------------------------------------------------------------------------


    def fnLoadAndPrintPar(self):
        self.fnLoadParameters()
        self.fnLoadCovar()
        self.fnPrintPar()
#-------------------------------------------------------------------------------

    def fnLoadCovar(self):                                                      #loads a covariance matrix,
                                                                                #calculates the errors and stores
                                                                                #it into self.diParameters
                                                                                #the function fnLoadParameters
                                                                                #must have been already carried out

        try:
            file = open('covar.dat')
            data = file.readlines()
            file.close()
            for i in range(len(string.split(data[1]))):                         #number of columns in second row
                for parameter in self.diParameters.keys():                      #search for parameter number and
                    if self.diParameters[parameter]['number']==i:               #retrieve value
                        fvalue=self.diParameters[parameter]['value']
                        break
                ferror=float(string.split(data[i+1])[i])                        #the first row contains a comment
                if ferror<0:
                    ferror=0
                ferror=math.sqrt(ferror)*fvalue                                 #calculate error
                self.diParameters[parameter]['error']=ferror
        except IOError:
            for parameter in self.diParameters.keys():
                self.diParameters[parameter]['error']=float(0)                  #fill all errors with zero

#-------------------------------------------------------------------------------

    def fnLoadFileListAndChangeToLocal(self):                                   #scans the setup.c file and creates a
        file = open(self.setupfilename,"r")                                          #list of the filenames
        data = file.readlines()                                                 #of the reflectivity data files, it also
        file.close()                                                            #modifies setup.c in this way that it
                                                                                #loads copies of reflectivity dat files
                                                                                #located in the working directory with the
                                                                                #file ending .mce, those files are the
                                                                                #modified files by the statistical error
                                                                                #analysis
        newdata = []
        filelist = []
        smatch1 = re.compile(r'fit_data.+?\"(.+?)\"',
                             re.IGNORECASE | re.VERBOSE)
        smatch2 = re.compile(r'(fit_data.+?\").+?(\")',
                             re.IGNORECASE | re.VERBOSE)
        smatch3=re.compile(r'\s*//',re.IGNORECASE | re.VERBOSE)
        for line in data:
            if ('fit_data' in line) and (not smatch3.match(line)):              #scan for file loading but not comment lines
                filelist.append(smatch1.search(line).group(1))                  #append path+filename to filelist
                newdata.append(smatch2.sub(r'\1'+
                              os.path.split(filelist[-1])[-1]+r'.mce\2',line))  #modifiy setup.c line with new filename
            else:
                newdata.append(line)                                            #no file loading -> do not modify setup.c line
        file=open(self.setupfilename,"w")                                            #write setup.c
        file.writelines(newdata)
        file.close()
        return filelist
#-------------------------------------------------------------------------------

    def fnLoadObject(self,sFileName):

        import pickle

        File=open(sFileName,"r")
        object=pickle.load(File)
        File.close()

        return object


#-------------------------------------------------------------------------------

    def fnLoadMolgroups(self):

        self.diMolgroups={}
        li=[]
        #load mol.dat

        file=open('mol.dat')
        data=file.readlines()
        file.close()

        i=0
        while 1:
            tdata=string.split(data[i])                                         #read header that contains molgroup data
            self.diMolgroups[tdata[1]]={}
            self.diMolgroups[tdata[1]].update({'headerdata':{}})
            self.diMolgroups[tdata[1]]['headerdata'].update({'Type':tdata[0]})
            self.diMolgroups[tdata[1]]['headerdata'].update({'ID':tdata[1]})
            j=2
            while 1:
                self.diMolgroups[tdata[1]]['headerdata'].update({tdata[j]:tdata[j+1]})
                j += 2
                if j==len(tdata):
                    break

            i+=2                                                                #skip header line for data columns
            zax=li[:]; areaax=li[:]; nslax=li[:]
            self.diMolgroups[tdata[1]].update({'zaxis':zax,'areaaxis':areaax,
                                              'nslaxis':nslax})

            while 1:
                if i>=len(data):
                    break
                tline=string.split(data[i])
                if tline:
                    self.diMolgroups[tdata[1]]['zaxis'].append(float(tline[0]))
                    self.diMolgroups[tdata[1]]['areaaxis'].append(float(tline[1]))
                    self.diMolgroups[tdata[1]]['nslaxis'].append(float(tline[2]))
                    i += 1
                else:
                    break

            i += 1
            if i>=len(data):
                break



#-------------------------------------------------------------------------------

    def fnLoadParameters(self):                                                 #loads the last row's parameters, and
                                                                                #ranges from par.dat and stores them into
                                                                                #self.diParameters; parameter names are
                                                                                #read from setup.c, since par.dat truncates
                                                                                #parameter names over 17 characters

                                                                                #after reading in the parameters check for
                                                                                #definitions of GAUSSIAN,
                                                                                #and check vs. active fit parameters
                                                                                #load them into self.molgroups dictionary
        self.diParameters={}

        iBurn=4000                                                              #check wether an existing MCMC exists
        bMCMCexists=False
        for i in range(1,9):
            iBurn *= 2
            if os.path.isdir('MCMC_'+str(iBurn)+'_500'):
                print 'Found '+'MCMC_'+str(iBurn)+'_500 \n'
                bMCMCexists=True
                break

        if bMCMCexists:
            file = open('MCMC_'+str(iBurn)+'_500/run.par')
            data=file.readlines()
            tParValues=[]
            for line in data:
                tParValues.append(float(string.split(line)[-1]))
        else:
            try:
                file = open('par.dat')
                data = file.readlines()
            except:
                print '--------------------------------------'
                print 'No par.dat found - Initializing fit.'
                print '--------------------------------------'
                self.fnMake()
                pr=subprocess.Popen(["nice","./fit","-S","-n","1"])                  # initial genetic run
                pr.wait()
                if os.path.isfile('par.dat'):
                    file = open('par.dat')
                    data = file.readlines()
                else:
                    print '--------------------------------------'
                    print 'Could not start up fit.'
                    print '--------------------------------------'
                    raise StandardError, 'Could not start up fit.'

            file.close()
            while data[-1][0]=='#':                                                 # delete comments in last lines
                data=data[:-1]
            tParValues=string.split(data[-1])[2:]                                   # best fit values from last row

        file = open(self.setupfilename)
        setupdata = file.readlines()
        file.close()
                                                                                # reminder .+? match all character
        smatch=re.compile(r'pars_add\(pars.*?\"(.+?)\"\s*,.+?\((.+?)\)\s*,(.+?),(.+?)\)',
                                  re.IGNORECASE | re.VERBOSE)
        smatch2=re.compile(r'\s*//',re.IGNORECASE | re.VERBOSE)
        for i in range(len(tParValues)):                                        #iterate over columns of par.dat
            for j in range(len(setupdata)):                                     #scan setup.c for the first parameter
                setupline=setupdata[j]                                          #  initialization and fetch that name
                if smatch.search(setupline) and (not smatch2.match(setupline)): #  plus its fit ranges
                    sParName   =      smatch.search(setupline).group(1)
                    sVarName   =      smatch.search(setupline).group(2)
                    flowerlimit=float(smatch.search(setupline).group(3))
                    fupperlimit=float(smatch.search(setupline).group(4))
                    del setupdata[j]                                            #delete the line where parameter found
                    break
            else:
                print '--------------------------------------'
                print 'Parameters do not match in setup file and par.dat'       #no parameter in setup.c left!
                print '--------------------------------------'
                raise StandardError,'Mismatch between setup file and par.dat.'
            if sParName in self.diParameters:
                print '--------------------------------------'
                print 'The parameter %s is defined twice in the garefl setup file.' % sParName
                print 'This is not supported by rs.py.'
                print '--------------------------------------'
                raise StandardError,'Doubled parameter in setup file.'
            ipar=int(i)                                                         #parameter number
            frelval=float(tParValues[i])                                        #relative par value is stored in file
            fvalue=(fupperlimit-flowerlimit)*frelval+flowerlimit                #calculate absolute par value
            self.diParameters[sParName]={'number'    :ipar,
                                         'variable'  :sVarName,
                                         'lowerlimit':flowerlimit,
                                         'upperlimit':fupperlimit,
                                         'relval'    :frelval,
                                         'value'     :fvalue}                   #save entry in diParameters[]


        self.chisq=float(string.split(data[-1])[1])                             #chi squared from second column last row

        bAbsoluteParameters=False
        for sParName in self.diParameters.keys():
            if self.diParameters[sParName]['relval']>1:
                bAbsoluteParameters=True
        if bAbsoluteParameters:
            for sParName in self.diParameters.keys():
                self.diParameters[sParName]['value']=self.diParameters[sParName]['relval']
                self.diParameters[sParName]['relval']=(self.diParameters[sParName]['value']-self.diParameters[sParName]['lowerlimit'])/(self.diParameters[sParName]['upperlimit']-self.diParameters[sParName]['lowerlimit'])



#-------------------------------------------------------------------------------

    def fnLoadStatData(self,sparse=0):

        cInteractor=self.cGaReflInteractor
        if os.path.isfile('sErr.dat') or os.path.isfile('isErr.dat'):
            pass
        elif os.path.isdir('MCMC'):
            cInteractor=self.cRefl1DInteractor

        self.diStatResults=cInteractor.fnLoadStatData(sparse)

        iMaxParameterNameLength=0
        for parname in self.diStatResults['Parameters'].keys():                 #cycle through all parameters
            if len(parname)>iMaxParameterNameLength:                            #determine length of longest parameter name for displaying
                iMaxParameterNameLength=len(parname)
        self.diStatResults['MaxParameterLength']=iMaxParameterNameLength

        self.diStatResults['NumberOfStatValues']=\
        len(self.diStatResults['Parameters'][self.diStatResults['Parameters'].keys()[0]]['Values'])



#-------------------------------------------------------------------------------

    def fnMake(self):                                                           #make setup.c and print sys output
        pr=subprocess.Popen(["rm","-f","setup.o"])
        pr.wait()
        pr=subprocess.Popen(["rm","-f","fit"])
        pr.wait()
        pr=subprocess.Popen(["make"])
        pr.wait()

#-------------------------------------------------------------------------------

    def fnModifyAndCopyFiles(self,filelist):                                    #reads original reflectivity files
                                                                                #and modifies them with normal deviates

        random.seed()                                                           #initialize random number generator
        for reflfile in filelist:                                               #iterate over all refl files
            file=open(reflfile)
            data=file.readlines()
            file.close()
            newdata=[]
            iFileType=0
            for line in data:
                if '#'in line:                                                  #do not modify comment lines
                    newdata.append(line)
                else:
                    columns=line.split()                                        #access columns
                    if not iFileType:                                         #Stuart's Mod for 3 Column i.e. NIST
                        if len(columns)==3:
                            iFileType=1
                        elif len(columns)==4:
                            iFileType=2                                         #Stuart's Mod for 4 Column i.e. ISIS
                        else:
                            iFileType=99                                        #unrecognized format
                    if not len(columns):
                        print 'Found empty line in data file.'
                    if len(columns)==3 and iFileType==1:
                        try:
                            fvalue=float(columns[1])                            #reflectivity
                            ferror=float(columns[2])                            #error
                            columns[1]=str(fvalue+random.normalvariate(0,1)*ferror) #modify reflectivity
                            newline=''
                            for column in columns:                              #glue columns together to one line
                                newline=newline+column+' '
                            newline=newline[:-1]+'\n'                           #add newline to end of line
                            newdata.append(newline)                             #add line to new data file
                        except:
                            print '-----------------------------------'
                            print 'Data file %s corrupt.' % (reflfile)
                            print 'File was identified being NIST type'
                            print '-----------------------------------'
                            raise StandardError, 'Corrupt Data file.'
                    elif len(columns)==4  and iFileType==2:
                        try:
                            fvalue=float(columns[2])                            #reflectivity
                            ferror=float(columns[3])                            #error
                            columns[2]=str(fvalue+random.normalvariate(0,1)*ferror)     #modify reflectivity
                            newline=''
                            for column in columns:                              #glue columns together to one line
                                newline=newline+column+' '
                            newline=newline[:-1]+'\n'                           #add newline to end of line
                            newdata.append(newline)                             #add line to new data file
                        except:
                            print '-----------------------------------'
                            print 'Data file %s corrupt.' % (reflfile)
                            print 'File was identified being ISIS type.'
                            print '-----------------------------------'
                            raise StandardError, 'Corrupt Data file.'

                    else:
                        print '-----------------------------------------------------'
                        print 'Filetype not recognized or contains errors: %s' % (reflfile)
                        print '-----------------------------------------------------'
                        raise StandardError, 'Corrupt Data file.'

            file=open(os.path.split(reflfile)[-1]+'.mce',"w")                   #write modified file into .mce file in
            file.writelines(newdata)                                            #working directory
            file.close()


#-------------------------------------------------------------------------------
    def fnnSLDEnvelopes(self,fGrid,fSigma,sname,shortflag=False, iContrast=-1):

        def fnInterpolateData(xdata,ydata,fMin,fMax,fGrid):

            f=fMin                                                              #target start
            ix=0                                                                #data start
            liInterpolated=[[],[]]                                              #interpolation result

            while f<=fMax:                                                      #main loop
                #print f,fMin,fMax,fGrid,ix,xdata[ix],len(xdata)
                if f<xdata[0]:                                                  #fill area where no data available with first value
                    liInterpolated[0].append(f)
                    liInterpolated[1].append(ydata[0])
                    f=f+fGrid
                elif f>xdata[-1]:                                               #fill remaining cells with last value
                    liInterpolated[0].append(f)
                    liInterpolated[1].append(ydata[-1])
                    f=f+fGrid
                else:                                                           #at least one data point surpassed by f
                    while (f>xdata[ix]) and (ix<(len(xdata)-1)):                #searching first data point past f
                        ix=ix+1
                    if f<xdata[ix]:                                             #was there a data point past f?
                        LowerX=ix-1                                             #calculate data surrounding f
                        UpperX=ix
                        fDataGrid=xdata[UpperX]-xdata[LowerX]
                        fInterpolate=ydata[LowerX]*((f-xdata[LowerX])/fDataGrid #do the interpolation
                            )+ydata[UpperX]*((xdata[UpperX]-f)/fDataGrid)
                        liInterpolated[0].append(f)
                        liInterpolated[1].append(fInterpolate)
                        f=f+fGrid
                    elif f==xdata[ix]:                                          #no interpolation needed
                        liInterpolated[0].append(xdata[ix])
                        liInterpolated[1].append(ydata[ix])
                        f=f+fGrid

            return liInterpolated


        def fnStoreEnvelope(liStoredEnvelopes,liStoredEnvelopeHeaders,envelope,percentile):
            iposition=len(liStoredEnvelopes)/2
            liStoredEnvelopes.insert(iposition,envelope[1])
            liStoredEnvelopes.insert(iposition,envelope[0])
            liStoredEnvelopeHeaders.insert(iposition,str(1-percentile/2))
            liStoredEnvelopeHeaders.insert(iposition,str(percentile/2))

        def fnSaveEnvelopes(liStoredEnvelopes,liStoredEnvelopeHeaders,iModel,fMin,fMax,fGrid,fSigma):
            file=open('Envelopes'+str(iModel)+'.dat',"w")

            if fSigma==0:                                                       #save all computed envelopes

                liSaveEnvelopeHeaders=liStoredEnvelopeHeaders
                liSaveEnvelopes=liStoredEnvelopes


            else:                                                               #save only multiples of fSigma in Sigma

                liSaveEnvelopeHeaders=[]
                liSaveEnvelopes=[]

                fmult=0.
                while 1:

                    fConfidence=scipy.special.erf(fmult/scipy.sqrt(2))          #upper and lower Percentiles for fmult*sigma
                    fLowerPerc=(1-fConfidence)/2
                    fUpperPerc=1-(1-fConfidence)/2

                    fStepSize=1/float(len(liStoredEnvelopeHeaders))             #positions in envelopes list
                    iLowerPerc=int(math.floor(fLowerPerc/fStepSize))
                    iUpperPerc=int(math.ceil(fUpperPerc/fStepSize))

                    if (iLowerPerc==0) or (iUpperPerc>=len(liStoredEnvelopeHeaders)-1):
                        break

                    liSaveEnvelopeHeaders.insert(0,'minus'+str(fmult)+'sigma')
                    liSaveEnvelopes.insert(0,liStoredEnvelopes[iLowerPerc])

                    if iUpperPerc <> iLowerPerc:
                        liSaveEnvelopeHeaders.append('plus'+str(fmult)+'sigma')
                        liSaveEnvelopes.append(liStoredEnvelopes[iUpperPerc])

                    fmult=fmult+fSigma

                liSaveEnvelopeHeaders.insert(0,'LowerEnvelope')
                liSaveEnvelopes.insert(0,liStoredEnvelopes[0])
                liSaveEnvelopeHeaders.append('UpperEnvelope')
                liSaveEnvelopes.append(liStoredEnvelopes[len(liStoredEnvelopeHeaders)-1])




            file.write("z ")
            for element in liSaveEnvelopeHeaders:
                if fSigma==0:
                    file.write("p"+element+" ")
                else:
                    file.write(element+" ")
            file.write("\n")

            f=fMin
            for i in range(len(liSaveEnvelopes[0])):
                file.write(str(f)+" ")
                for element in liSaveEnvelopes:
                    file.write(str(element[i])+" ")
                file.write("\n")
                f=f+fGrid


            file.close()

        def fnCalculateEnvelope(profilelist):

            LowerEnvelope=numpy.array(profilelist[0][1])
            UpperEnvelope=numpy.array(profilelist[0][1])

            for iteration in profilelist:
                LowerEnvelope=numpy.minimum(LowerEnvelope,numpy.array(iteration[1]))
                UpperEnvelope=numpy.maximum(UpperEnvelope,numpy.array(iteration[1]))

            fArea=numpy.sum(numpy.subtract(UpperEnvelope,LowerEnvelope))

            envelope=[LowerEnvelope.tolist(),UpperEnvelope.tolist()]

            return fArea,envelope



        self.fnRecreateStatistical()                                            #Load Parameters and profiles from stat data
        iNumberOfModels=self.fnGetNumberOfModelsFromSetupC()                    #how many profiles to analyze?

        if iContrast == -1:
            iModelStart=0
            iModelEnd=len(self.diStatResults['nSLDProfiles'][0])
        else:
            iModelStart=iContrast
            iModelEnd=iContrast+1

        for iModel in range(iModelStart,iModelEnd):                             #do the analysis separately for each model
            print 'Analyzing model %i ...' % (iModel)

            liStoredEnvelopes=[]                                                #List of stored envelopes
            liStoredEnvelopeHeaders=[]                                          #and their headers

            fMin=0;fMax=0                                                       #initializing boundaries
            profilelist=[]                                                      #extracting all profiles related to the actual
            for iteration in self.diStatResults['nSLDProfiles']:                #model
                profilelist.append(iteration[iModel][:])
                if fMin>profilelist[-1][0][0]:
                    fMin=profilelist[-1][0][0]
                if fMax<profilelist[-1][0][-1]:
                    fMax=profilelist[-1][0][-1]
            fMax=math.floor((fMax-fMin)/fGrid)*fGrid+fMin                       #make fMax compatible with fGrid and fMin

            print 'Rebinning data...'
            for i in range(len(profilelist)):
                profilelist[i]=fnInterpolateData(profilelist[i][0],profilelist[i][1],fMin,fMax,fGrid)



            iNumberOfProfiles=len(profilelist)

            if not shortflag:
                for iPercentile in range(iNumberOfProfiles):

                    print 'Calculating %f percentile...' % (1-float(iPercentile)/float(iNumberOfProfiles))

                    fArea,liEnvelope=fnCalculateEnvelope(profilelist)                   #calculate envelope
                    fnStoreEnvelope(liStoredEnvelopes,liStoredEnvelopeHeaders,
                                    liEnvelope,float(iPercentile)/float(iNumberOfProfiles))     #store envlope in list

                    if iPercentile<>iNumberOfProfiles-1:
                        iMaxReduction=0
                        for i in range(len(profilelist)):                               #eliminate the profile with the largest reduction
                                                                                        #in area
                            profiletestlist=profilelist[:]
                            profiletestlist.pop(i)
                            fTestArea,fTestEnvelope=fnCalculateEnvelope(profiletestlist)
                            if fArea>fTestArea:
                                iMaxReduction = i
                                fArea=fTestArea

                        profilelist.pop(iMaxReduction)                                  #eliminate profile
            else:
                fSigmaStep=0.1
                iPercentile=0
                fArea,liEnvelope=fnCalculateEnvelope(profilelist)                           #calculate envelope
                fnStoreEnvelope(liStoredEnvelopes,liStoredEnvelopeHeaders,
                                liEnvelope,float(iPercentile)/float(iNumberOfProfiles))     #store envlope in list
                iPercentile+=1

                while iPercentile<iNumberOfProfiles:

                    print 'Short: Calculating %f percentile...' % (1-float(iPercentile)/float(iNumberOfProfiles))

                    lScoring=[]                                                         #list of profile scores
                    for i in range(len(profilelist)):                                   #eliminate the profile with the largest reduction
                                                                                        #in area
                        profiletestlist=profilelist[:]
                        profiletestlist.pop(i)
                        fTestArea,fTestEnvelope=fnCalculateEnvelope(profiletestlist)
                        lScoring.append([i,fTestArea])
                        #print lScoring

                    lScoring=sorted(lScoring, key=operator.itemgetter(1))               #sort by lowest area

                    fConf=(1-float(iPercentile)/float(iNumberOfProfiles))               #momentary confidence level
                    fSig=scipy.special.erfinv(fConf)*scipy.sqrt(2)                      #related sigma value
                    fSigT=math.floor(fSig/fSigmaStep)*fSigmaStep                              #get to the lower sigma step
                    fConfT=scipy.special.erf(fSigT/scipy.sqrt(2))                        #target confidence level
                                                                                        #number of profiles to eliminate
                    iElimination=iNumberOfProfiles-iPercentile-int(fConfT*iNumberOfProfiles)

                    print "iPercentile %i iNumberOfProfiles %i iElimination %i fConf %e fSig %e fSigT %e fConfT %e" % (iPercentile, iNumberOfProfiles, iElimination, fConf, fSig, fSigT, fConfT)

                    lScoring=sorted(lScoring[0:iElimination], key=operator.itemgetter(0), reverse=True)

                    for i in range(iElimination):
                        profilelist.pop(lScoring[i][0])                                  #delete profiles starting with highest indices
                        fArea,liEnvelope=fnCalculateEnvelope(profilelist)                #calculate envelope
                        fnStoreEnvelope(liStoredEnvelopes,liStoredEnvelopeHeaders,
                                liEnvelope,float(iPercentile)/float(iNumberOfProfiles))  #store envlope in list
                        iPercentile+=1

                    #raw_input("Press Enter to continue...")

            fnSaveEnvelopes(liStoredEnvelopes,liStoredEnvelopeHeaders,iModel,
                fMin,fMax,fGrid,fSigma)

#-------------------------------------------------------------------------------


    def fndiParametersNumberSort(self,a,b):                                     #sort criterium for dParameters is the
                                                                                #parameter number
        return cmp(self.diParameters[a]['number'],self.diParameters[b]['number'])

#-------------------------------------------------------------------------------

    def fnPlotMolgroups(self, plotname):

        from matplotlib.font_manager import fontManager, FontProperties

        font = FontProperties(size='x-small')

        self.fnLoadParameters()                                                 #Load Parameters and modify setup.cc
        self.fnBackup()                                                         #Backup setup.c, and other files
        try:
            liParameters=self.diParameters.keys()                               #get list of parameters from setup.c/par.dat
            liParameters.sort(self.fndiParametersNumberSort)                    #sort by number of appereance in setup.c

            liAddition=[]
            for parameter in liParameters:                                      #cycle through all parameters
                liAddition.append(('%s = %s;\n'%                                #change setup.c to quasi fix all parameters
                    (self.diParameters[parameter]['variable'],self.diParameters[parameter]['value'])))
            self.fnWriteConstraint2SetupC(liAddition)                               #write out
            subprocess.call(["rm","-f","setup.o"])
            subprocess.call(["rm","-f","mol.dat"])
            subprocess.call(["sync"])                                           #synchronize file system
            time.sleep(1)                                                       #wait for system to clean up
            self.fnMake()                                                       #compile changed setup.c
            subprocess.call(["./fit", "-o"])                                    #write out profile.dat and fit.dat
            subprocess.call(["sync"])                                           #synchronize file system
            time.sleep(1)                                                       #wait for system to clean up

        finally:
            self.fnRemoveBackup()

        self.fnLoadMolgroups()                                                  #Load Molgroups into self.diMolgroups

        file=open('areatab.dat','w')                                            #save all molgroupdata in table for loading
                                                                                #into Igor

        file.write('z ')                                                        #write header line
        for element in self.diMolgroups:
            file.write(self.diMolgroups[element]['headerdata']['ID']+' ')
        file.write("summol water waterperc")
        file.write('\n')

        element=self.diMolgroups.keys()[0]
        datalength=len(self.diMolgroups[element]['zaxis'])

        for i in range(datalength):
            file.write(str(self.diMolgroups[element]['zaxis'][i])+' ')
            summe=0; normarea=0
            for element in self.diMolgroups:
                if element!='normarea':
                    summe=summe+self.diMolgroups[element]['areaaxis'][i]
                else:
                    normarea=self.diMolgroups[element]['areaaxis'][i]
                file.write(str(self.diMolgroups[element]['areaaxis'][i])+' ')
            file.write(str(summe)+' '+str(normarea-summe)+' '+str((normarea-summe)/normarea))
            file.write('\n')

        file.close()

        plt.figure(1, figsize=(14,10))

        plt.subplot(221)

        areasum=[]                                                              #add up area
        nSLsum=[]                                                               #calculate nSL
        zax=[]
        for element in self.diMolgroups:
            area=[]
            nSL=[]
            zax=self.diMolgroups[element]['zaxis']
            stepsize=float(zax[1]-zax[0])
            length=len(zax)
            for i in range(length):
                    area.append(self.diMolgroups[element]['areaaxis'][i])
                    nSL.append(self.diMolgroups[element]['nslaxis'][i]*1e4)
            plt.subplot(221)
            plt.plot(zax,area,label=element)
            plt.subplot(222)
            plt.plot(zax,nSL,label=element)

            if element!='normarea':                                             #do not sum up normarea indicator
                if not areasum:
                    for i in range(length):
                        areasum.append(0.)
                        nSLsum.append(0.)
                for i in range(length):
                    areasum[i]=areasum[i]+area[i]
                    nSLsum[i]=nSLsum[i]+nSL[i]

            if element=='normarea':                                             #Get normarea for nSLD calculations
                normarea=self.diMolgroups[element]['areaaxis'][0]

        plt.subplot(221)
        plt.plot(zax,areasum,label='Sum')
        plt.subplot(222)
        plt.plot(zax,nSLsum,label='Sum')


        nSLDSum=[]                                                             #calculate nSLD sum
        nSLDtVolFracSum=[]                                                     #calculate nSLD times vol frac sum
        for i in range(len(areasum)):
            nSLDSum.append(0)
            nSLDtVolFracSum.append(0)


        for element in self.diMolgroups:
            nSLDtVolFrac=[]                                                     #calculate nSLD times volfrac, at the moment times volfrac is not used
            if element!='normarea':                                             #Get normarea for nSLD calculations
                for i in range(len(self.diMolgroups[element]['areaaxis'])):
                    area=self.diMolgroups[element]['areaaxis'][i]
                    if area:
                        nSLD=self.diMolgroups[element]['nslaxis'][i]/(stepsize*area)*1e6
                    else:
                        nSLD=0
                    nSLDtVolFrac.append(nSLD)#*(area/normarea))
                    nSLDtVolFracSum[i]=nSLDtVolFracSum[i]+nSLDtVolFrac[i]
                plt.subplot(223)
                plt.plot(zax,nSLDtVolFrac,label=element)
        plt.subplot(223)
        #plt.plot(zax,nSLDtVolFracSum,label='nSLDtVolFracSum')

        for j in range(4):                                                      #loop over contrast mixtures
            if j==0:
                fBulknSLD=6.34
            if j==1:
                fBulknSLD=4.00
            if j==2:
                fBulknSLD=0.00
            if j==3:
                fBulknSLD=-0.56

            for i in range(length):                                             #calculate nSLD for several cases
                nSLDSum[i]=nSLsum[i]*1E2/(stepsize*normarea)+fBulknSLD*(1-(areasum[i]/normarea))
            plt.subplot(224)
            plt.plot(zax,nSLDSum,label='nSLDsum CM'+str(fBulknSLD))




        plt.subplot(221)
        plt.ylabel('area / Ang+2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0 ,prop=font)

        plt.subplot(222)
        plt.ylabel('nSL / 1e-4 Ang')
        plt.xlabel('z / Ang')
        plt.legend(loc=0 ,prop=font)

        plt.subplot(223)
        plt.ylabel('nSLD')# * volfrac / 1e-6 Ang-2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0 ,prop=font)


        plt.subplot(224)
        plt.ylabel('nSLD / 1e-6 Ang-2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0 ,prop=font)

        plt.suptitle('%s \n\n Area and nSLD Profile' % (plotname))
        plt.savefig('%s_mol.png' % (plotname), format='png')
        plt.show()
        #plt.close()
#-------------------------------------------------------------------------------

    def fnPlotFit(self, plotname):

        from matplotlib.font_manager import fontManager, FontProperties

        font = FontProperties(size='x-small')

        plt.figure(1, figsize=(14,10))

        iCounter=0
        while 1:
            sfilename='fit'+str(iCounter)+'.dat'
            if os.path.isfile(sfilename):
                file = open(sfilename,'r')
                data=file.readlines()
                file.close()
                data=data[1:]

                k=0; l=0
                qlist=[]; dqlist=[]
                Rlist=[]; dRlist=[]; fitlist=[]; fitRFlist=[]
                RFlist=[]; dRFlist=[]; reslist=[]
                resplus=[]; resminus=[]
                for line in data:
                    splitline=line.split()
                    qlist.append(float(splitline[0]))
                    dqlist.append(float(splitline[1]))
                    Rlist.append(float(splitline[2]))
                    dRlist.append(float(splitline[3]))
                    fitlist.append(float(splitline[4]))
                    RFlist.append(float(splitline[2])*pow(float(splitline[0]),4))
                    dRFlist.append(float(splitline[3])*pow(float(splitline[0]),4))
                    fitRFlist.append(float(splitline[4])*pow(float(splitline[0]),4))
                    reslist.append((float(splitline[2])-float(splitline[4]))*pow(float(splitline[3]),-1))
                    resplus.append(1)
                    resminus.append(-1)

                plt.subplot(221)
                plt.errorbar(qlist,Rlist,yerr=dRlist,xerr=dqlist,fmt='.')
                plt.semilogy(qlist,fitlist,label='fit'+str(iCounter))
                plt.xlim(xmin=-0.01)

                plt.subplot(222)
                plt.errorbar(qlist,RFlist,yerr=dRFlist,xerr=dqlist,fmt='.')
                plt.semilogy(qlist,fitRFlist,label='fit'+str(iCounter))
                plt.xlim(xmin=-0.01)

                plt.subplot(223)
                plt.plot(qlist,reslist,label='fit'+str(iCounter))
                plt.plot(qlist,resplus,'r')
                plt.plot(qlist,resminus,'r')

                iCounter=iCounter+1

            else:
                break

        iCounter=0
        while 1:
            sfilename='profile'+str(iCounter)+'.dat'
            if os.path.isfile(sfilename):
                file = open(sfilename,'r')
                data=file.readlines()
                file.close()
                data=data[1:]

                k=0; l=0
                zlist=[]; rholist=[]
                for line in data:
                    splitline=line.split()
                    zlist.append(float(splitline[0]))
                    rholist.append(float(splitline[1])*1e6)

                plt.subplot(224)
                plt.plot(zlist,rholist,label='profile'+str(iCounter))

                iCounter=iCounter+1

            else:
                break

        plt.subplot(221)
        plt.ylabel('Reflectivity / R')
        plt.xlabel('q / Ang-1')
        plt.legend(loc=0 ,prop=font)

        plt.subplot(222)
        plt.ylabel('R*Q^4')
        plt.xlabel('q / Ang-1')
        plt.legend(loc=0 ,prop=font)

        plt.subplot(223)
        plt.ylabel('Normalized Residuals/ (R -F(q))/dR')
        plt.xlabel('q / Ang-1')
        plt.legend(loc=0 ,prop=font)

        plt.subplot(224)
        plt.ylabel('nSLD / 1e-6 Ang-2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0 ,prop=font)

        plt.suptitle('%s \n\n Reflectivity data and fit, residuals and profile' % (plotname))
        plt.savefig('%s_fit.png' %(plotname), format='png')
        plt.show()
        #plt.close()
#-------------------------------------------------------------------------------

    def fnPrintPar(self):                                                       #prints parameters and their errors
                                                                                #from the covariance matrix on the screen

        litest=self.diParameters.keys()
        litest.sort(self.fndiParametersNumberSort)
        for parameter in litest:
            fRange=(self.diParameters[parameter]['upperlimit']
                            -self.diParameters[parameter]['lowerlimit'])
            fLowLim=self.diParameters[parameter]['lowerlimit']
            fValue=self.diParameters[parameter]['value']
            sRangeIndicator=''
            for i in range(10):                                                 #creates the par range ascii overview
                if ((fValue >= float(i)/10*fRange+fLowLim) and
                    (fValue < float(i+1)/10*fRange+fLowLim)):
                    sRangeIndicator += '|'
                else:
                    if (fValue == float(i+1)/10*fRange+fLowLim) and (i==9):
                        sRangeIndicator += '|'
                    else:
                        sRangeIndicator += '.'
            print '%2i %25s  %s %15g +/- %g in [%g,%g]' % (self.diParameters[parameter]['number'],
                                                 parameter,sRangeIndicator,fValue,
                                                 self.diParameters[parameter]['error'],
                                                 fLowLim,self.diParameters[parameter]['upperlimit'])
        print 'Chi squared: %g' % self.chisq

#-------------------------------------------------------------------------------

    def fnPullMolgroup(self, liMolgroupNames,dSparse=0):

        """
        Function recreates statistical data and extracts only area profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        """

        self.fnRecreateStatistical(sparse=dSparse)

        di={'zaxis' : self.diStatResults['Molgroups'][0][self.diStatResults['Molgroups'][0].keys()[0]]['zaxis']}

        for i,iteration in enumerate(self.diStatResults['Molgroups']):

            sumprofile=[]
            for molgroup in liMolgroupNames:
                try:
                    sumprofile = [i+j for i,j in zip(sumprofile,iteration[molgroup]['areaaxis'])] if \
                    sumprofile else iteration[molgroup]['areaaxis']
                except:
                    print'Molecular group %s does not exist.' % molgroup
                    sys.exit()

            di['iter%i' % i]=sumprofile

        diStat={'zaxis' : [], 'm2sigma' : [], 'msigma' : [], 'mean' : [], 'psigma' : [], 'p2sigma' : []}
        for i in range(len(di[di.keys()[0]])):
            liOnePosition=[iteration[i] for key,iteration in di.iteritems() if key != 'zaxis']
            stat=[scipy.stats.scoreatpercentile(liOnePosition,percentile) for percentile in [2.3,15.9,50,84.1,97.7]]
            diStat['zaxis'].append(str(di['zaxis'][i]))
            diStat['m2sigma'].append(stat[0])
            diStat['msigma'].append(stat[1])
            diStat['mean'].append(stat[2])
            diStat['psigma'].append(stat[3])
            diStat['p2sigma'].append(stat[4])

        cInteractor=self.cGaReflInteractor
        cInteractor.fnSaveSingleColumns('pulledmolgroups.dat',di)
        cInteractor.fnSaveSingleColumns('pulledmolgroupsstat.dat',diStat)

#-------------------------------------------------------------------------------
    def fnRecreateStatistical(self,bRecreateMolgroups=True,sparse=0):           #Recreates profile and fit data
                                                                                #associated with stat file

        self.fnLoadParameters()                                                 #Load Parameters into self.diParameters
        self.fnLoadStatData(sparse)                                             #Load Results from statistical analysis
        iNumberOfModels=self.fnGetNumberOfModelsFromSetupC()                    #how many profiles to analyze?
        j=0
        self.diStatResults['nSLDProfiles']=[]                                   #delete list of all nSLD profiles
        self.diStatResults['Molgroups']=[]                                      #delete list of all molecular groups

        for iteration in range(self.diStatResults['NumberOfStatValues']):       #cycle through all individual stat. results
            try:
                self.fnBackup()                                                 #Backup setup.c, and other files
                self.diStatResults['nSLDProfiles'].append([])                   #appends a new list for profiles for the
                                                                                #current MC iteration
                liParameters=self.diParameters.keys()                           #get list of parameters from setup.c/par.dat
                liParameters.sort(self.fndiParametersNumberSort)                #sort by number of appereance in setup.c
                bConsistency=True
                for element in liParameters:
                    if element not in self.diStatResults['Parameters'].keys():
                        bConsistency=False
                if bConsistency:                                                #check for consistency
                    liAddition=[]
                    for parameter in liParameters:                              #cycle through all parameters
                        liAddition.append(('%s = %s;\n'%                        #change setup.c to quasi fix all parameters
                            (self.diParameters[parameter]['variable'],          #to the result of the stat. analysis
                            self.diStatResults['Parameters'][parameter]['Values'][iteration])))
                    self.fnWriteConstraint2SetupC(liAddition)                   #write out
                    #raw_input("Please enter ...")
                    self.fnMake()                                               #compile changed setup.c
                    print 'Processing parameter set %i.\n \n' % (j)
                    self.fnWriteOutGareflModel()                                #write out profile.dat and fit.dat

                    i=0                                                         #store profile data in liContourProfile
                    while i<iNumberOfModels:
                        self.diStatResults['nSLDProfiles'][-1].append(([],[]))  #adding another profile for the model
                        sFileName='profile'+str(i)+'.dat'                       #Load profile#.dat
                        file=open(sFileName)
                        data=file.readlines()
                        file.close()
                        data = data[1:]                                         #delete header

                        for line in data:                                       #extract nSLD profile data line by line
                            splitline=line.split()
                            dz=float(splitline[0])
                            drho=float(splitline[1])
                            self.diStatResults['nSLDProfiles'][-1][-1][0].append(dz)
                            self.diStatResults['nSLDProfiles'][-1][-1][1].append(drho)
                        i += 1

                    if bRecreateMolgroups:
                        subprocess.call(["./fit", "-o"])                        #write out mol.dat
                        self.fnLoadMolgroups()                                  #populate self.diMolgroups
                        self.diStatResults['Molgroups'].append(self.diMolgroups)#append molgroup information to self.diStatResults


                else:
                    print 'Statistical error data and setup file do not match'
                    self.fnRemoveBackup()
                    raise ''


            finally:
                self.fnRemoveBackup()
                j += 1

        self.fnSaveObject(self.diStatResults,'StatDataPython.dat')              #save stat data to disk


#-------------------------------------------------------------------------------


    def fnReplaceParameterLimitsInSetup(self, sname, flowerlimit, fupperlimit): #scans setup.c file for parameter with
        file = open(self.setupfilename,'r+')                                    #name sname and replaces the lower and
        data = file.readlines()                                                 #upper fit limits by the given values
        file.close()
        smatch=re.compile(r'(pars_add\(pars.*?\"'+sname+
                          '.+?,.+?,).+?(,).+?(\))',
                          re.IGNORECASE | re.VERBOSE)
        newdata=[]
        for line in data:
            newdata.append(smatch.sub(r'\1 '+str(flowerlimit)+r'\2 '
                           +str(fupperlimit)+r'\3',line))

        file = open(self.setupfilename,'w')
        file.writelines(newdata)
        file.close()
#-------------------------------------------------------------------------------

    def fnRemoveBackup(self):                                                   #deletes the backup directory
        self.fnRestoreBackup()
        subprocess.call(['rm','-rf','rsbackup'])
#-------------------------------------------------------------------------------

    def fnRestoreBackup(self):                                                  #copies all files from the backup directory
        if os.path.isfile('rsbackup/pop.dat'):
            subprocess.call(['cp','rsbackup/pop.dat','.'])                       #back to the working directory
        if os.path.isfile('rsbackup/par.dat'):
            subprocess.call(['cp','rsbackup/par.dat', '.'])
        if os.path.isfile('rsbackup/covar.dat'):
            subprocess.call(['cp','rsbackup/covar.dat','.'])
        if os.path.isfile('rsbackup/fit*.dat'):
            subprocess.call(['cp','rsbackup/fit*.dat','.'])
        if os.path.isfile('rsbackup/fit'):
            subprocess.call(['cp','rsbackup/fit','.'])
        if os.path.isfile('rsbackup/fit*.dat'):
            subprocess.call(['cp','rsbackup/model*.dat','.'])
        if os.path.isfile('rsbackup/pop_bak.dat'):
            subprocess.call(['cp','rsbackup/pop_bak.dat','.'])
        if os.path.isfile('rsbackup/profile*.dat'):
            subprocess.call(['cp','rsbackup/profile*.dat','.'])
        if os.path.isfile('rsbackup/'+self.setupfilename):
            subprocess.call(['cp','rsbackup/'+self.setupfilename,'.'])
        if os.path.isfile('rsbackup/setup.o'):
            subprocess.call(['cp','rsbackup/setup.o','.'])
#-------------------------------------------------------------------------------


    def fnRestoreFileList(self,filelist):                                       #not used
        file = open(self.setupfilename,"r")
        data = file.readlines()
        file.close()
        newdata = []
        smatch = re.compile(r'(fit_data.+?\").+?(\")',re.IGNORECASE | re.VERBOSE)
        for line in data:
            if 'void constr_models' in line:
                newdata.append(smatch.sub(r'\1'+filelist[0]+r'\2',line))
                pr=subprocess.Popen(['rm','-f',os.path.split(filelist[0])[-1]+'.mce'])
                pr.wait()
                del filelist[0]
            else:
                newdata.append(line)
        file=open(self.setupfilename,"w")
        file.writelines(newdata)
        file.close()
#-------------------------------------------------------------------------------

    def fnSaveObject(self,object,sFileName):

        import pickle

        File=open(sFileName,"w")
        pickle.dump(object,File)
        File.close()


#-------------------------------------------------------------------------------

    def fnStatTable(self,sTableName, fConfidence):

        def fnTexFormatf(fLow,fMed,fHigh):

            def fnDetPrec(fF):
                if math.fabs(fF)<1e-4:                                          #applies to nSLDs
                    fF *= 1E6
                if fF>0:                                                        #determine precision
                    fPrec=math.ceil(math.log10(fF)*(-1))
                else: fPrec=0
                if fPrec>0:                                                     #takes care of numbers like fF=0.0095
                    iPrec=int(fPrec)
                    if round(fF,iPrec)==round(fF,iPrec-1):                      #which should be rounded to 0.01 and
                        fPrec -= 1                                              #not to 0.010
                return fF,fPrec

            fLowDiff,fLowPrec=fnDetPrec(fMed-fLow)
            fHighDiff,fHighPrec=fnDetPrec(fHigh-fMed)
            fMed,fMedPrec=fnDetPrec(fMed)
            fPrec=(max(fLowPrec,fHighPrec))+1.0
            iPrec=int(fPrec)

            fLowDiff=round(fLowDiff+0.5*math.pow(10,(-1)*(fPrec+1.0)),iPrec)           #conservative rounding
            fHighDiff=round(fHighDiff+0.5*math.pow(10,(-1)*(fPrec+1.0)),iPrec)
            fMed=round(fMed,iPrec)
            return fLowDiff,fMed,fHighDiff,iPrec


        self.fnAnalyzeStatFile(fConfidence)                                     #analyze stat data and
                                                                                #populate self.diStatresults

        file=open(sTableName,'r')                                             #load in template
        template=file.readlines()
        file.close()

        table=[]                                                                #table to be created

        for line in template:
            splitline=line.split()
            for i,phrasetex in enumerate(splitline):                            #look at each string in template
                phrase=phrasetex.replace('\\','')                               #remove all '\'
                if phrase in self.diStatResults['Parameters']:                                #if it resembles a paramter name -> replace
                    fMedian=self.diStatResults['Parameters'][phrase]['Median']
                    fLowPerc=self.diStatResults['Parameters'][phrase]['LowPerc']
                    fHighPerc=self.diStatResults['Parameters'][phrase]['HighPerc']
                    fLowPercDiff,fMedianDiff,fHighPercDiff,iPrec=fnTexFormatf(fLowPerc,fMedian,fHighPerc)
                    sPrec='%(#).'+str(iPrec)+'f'
                    temp='$'                                                    #Latex format
                    temp += (sPrec % {'#':fMedianDiff})
                    temp += '_{'
                    temp=temp+'-'+(sPrec % {'#':fLowPercDiff})
                    temp += '}^{'
                    temp=temp+'+'+(sPrec % {'#':fHighPercDiff})
                    temp += '}$'
                    splitline[i]=temp
            table.append(' '.join(splitline)+'\n')


        file=open("StatTable.tex",'w')                                          #save table to file
        file.writelines(table)
        file.close()

#-------------------------------------------------------------------------------
    def fnTruncateParDat(self):
        if os.path.exists("par.dat"):                                           #check if something is there to truncate
            file=open("par.dat",'r')
            data=file.readlines()
            file.close()

            iFirstDataLine=-1
            iLastDataLine=-1

            for i in range(int(len(data)/2)):                                   #find position of first and last line
                if (data[i][0]!='#') and (iFirstDataLine==-1):                  #with data
                    iFirstDataLine=i
                if (data[len(data)-1-i][0]!='#') and (iLastDataLine==-1):
                    iLastDataLine=len(data)-1-i
                if (iFirstDataLine!=-1) and (iLastDataLine!=-1):
                    break

            if (iFirstDataLine!=-1) and (iLastDataLine!=-1):                    #cut out everything between those two
                data1=data[:iFirstDataLine+1]                                   #lines
                data2=data[iLastDataLine:]
                file=open("par.dat",'w')
                file.writelines(data1)
                file.writelines(data2)
                file.close()


#-------------------------------------------------------------------------------

    def fnWriteConstraint2SetupC(self,liExpression):                            #Writes a list of expressions at the beginning of
                                                                                #the constraint section in setup.c
        file = open(self.setupfilename,"r")                                     #open setup.c
        data = file.readlines()
        file.close()
        for i in range(len(data)):                                              #need integer index for slicing
            if 'void constr_models(' in data[i]:                                #searching for function declaration
                newdata=data[:i+2]+liExpression+data[i+2:]                      #inserting Expression two lines later,
                                                                                #after the initial '{' of the function

                break                                                           #break loop, only one insertion
        file = open(self.setupfilename,"w")                                     #write changed setup.c
        file.writelines(newdata)
        file.close()

#-------------------------------------------------------------------------------

    def fnWriteOutGareflModel(self):
        pr=subprocess.Popen(["./fit", "-g"])                                     #write out profile.dat and fit.dat
        pr.wait()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


def Auto(convergence=0.001):                                                    # automatic fit
    ReflPar.fnMake()
    subprocess.call(['rm','-f','par.dat'],stdout=open(os.devnull,"w"))
    subprocess.call(['rm','-f','pop.dat'],stdout=open(os.devnull,"w"))
    subprocess.call(['rm','-f','covar.dat'],stdout=open(os.devnull,"w"))
    subprocess.call(['nice','./fit','-eS','-n','21'],stdout=open(os.devnull,"w"))
    subprocess.call(['cp','pop_bak.dat','pop.dat'])                             # copy newest population into pop.dat
    print 'Genetic run, approximate roughness'
    ReflPar.fnLoadAndPrintPar()
    fOldChiSq=ReflPar.fnGetChiSq()

    while 1:                                                                    # genetic runs until chisq<20 or not improving
        subprocess.call(['nice','./fit','-peS','-n','51'],stdout=open(os.devnull,"w"))
        subprocess.call(['cp','pop_bak.dat','pop.dat'])
        print 'Genetic run, approximate roughness'
        ReflPar.fnLoadAndPrintPar()
        if (ReflPar.fnGetChiSq()<10 or (fOldChiSq-ReflPar.fnGetChiSq())<convergence):
            subprocess.call(['cp','pop_bak.dat','pop.dat'])
            break
        fOldChiSq=ReflPar.fnGetChiSq()
    AutoFinish(convergence)                                                     # AutoFinish takes over
    ReflPar.fnTruncateParDat()

#-------------------------------------------------------------------------------

def AutoFinish(convergence=0.001):
                                                                                # automatic fit starting with Amoeba, approximate roughness
    ReflPar.fnLoadParameters()
    fOldChiSq=ReflPar.fnGetChiSq()

    while 1:                                                                    #Amoeba, approximate roughness
        subprocess.call(['nice','./fit','-peaS'],stdout=open(os.devnull,"w"))   #until chisq<10 or no improvement
        print 'Amoeba, approximate roughness'
        ReflPar.fnLoadAndPrintPar()
        subprocess.call(['cp','pop_bak.dat','pop.dat'])                         # copy newest population into pop.dat
        if ReflPar.chisq<5 or (fOldChiSq-ReflPar.fnGetChiSq()<convergence):
            break
        fOldChiSq=ReflPar.fnGetChiSq()
    AutoFinish2(convergence)                                                    #AutoFinish2 takes over

#-------------------------------------------------------------------------------


def AutoFinish2(convergence=0.001):                                             # automatic fit, only local minimum refinement
    ReflPar.fnLoadParameters()

    fOldChiSq=1E24
    fBlockChiSq=5E23
    fTempChiSq=1E23
    iGeneticIsUseless=False
    iGeneticSkipCount=0
    iFinishCounter=0

    while 1:

        if iGeneticIsUseless:                                                  # skip genetic algorithm after it proved to be
            iGeneticSkipCount += 1                                            # useless and give it a chance every 5 iterations
            print ' '
            print 'Genetic algorithm skipped'

        subprocess.call(['cp','pop.dat','pop_rspy.dat'])
        while (not iGeneticIsUseless) and iGeneticSkipCount<6:
            iGeneticIsUseless=False
            iGeneticSkipCount=0
            subprocess.call(['nice','./fit','-pe','-n','21'],stdout=open(os.devnull,"w"))
            print 'Genetic run, correct roughness'
            ReflPar.fnLoadAndPrintPar()
            fTempChiSq2=ReflPar.fnGetChiSq()
            if fTempChiSq-fTempChiSq2 < convergence:
                if fTempChiSq2 < fTempChiSq:
                    fTempChiSq=fTempChiSq2
                    subprocess.call(['cp','pop_bak.dat','pop.dat'])
                else:
                    pass
                    #iGeneticIsUseless=True
                break
            else:
                fTempChiSq=fTempChiSq2
                subprocess.call(['cp','pop_bak.dat','pop.dat'])
        if fTempChiSq>fBlockChiSq:
                    subprocess.call(['cp','pop_rspy.dat','pop.dat'])
                    fTempChiSq=fBlockChiSq


        subprocess.call(['cp','pop.dat','pop_rspy.dat'])
        while 1:
            subprocess.call(['nice','./fit','-pea'],stdout=open(os.devnull,"w"))
            print 'Amoeba, correct roughness'
            ReflPar.fnLoadAndPrintPar()
            fTempChiSq2=ReflPar.fnGetChiSq()
            if fTempChiSq-fTempChiSq2 < convergence:
                if fTempChiSq2 < fTempChiSq:
                    fTempChiSq=fTempChiSq2
                    subprocess.call(['cp','pop_bak.dat','pop.dat'])
                break
            else:
                subprocess.call(['cp','pop_bak.dat','pop.dat'])
                fTempChiSq=fTempChiSq2
        if fTempChiSq>fBlockChiSq:
            subprocess.call(['cp','pop_rspy.dat','pop.dat'])
            fTempChiSq=fBlockChiSq

        subprocess.call(['cp','pop.dat','pop_rspy.dat'])
        while 1:
            print 'Calling LM fit ..'
            subprocess.call(['nice','./fit','-pel'],stdout=open(os.devnull,"w"))
            print 'Levenberg-Marquardt, correct roughness'
            ReflPar.fnLoadAndPrintPar()
            fTempChiSq2=ReflPar.fnGetChiSq()
            if fTempChiSq-fTempChiSq2 < convergence:
                if fTempChiSq2 < fTempChiSq:
                    fTempChiSq=fTempChiSq2
                    subprocess.call(['cp','pop_bak.dat','pop.dat'])
                break
            else:
                subprocess.call(['cp','pop_bak.dat','pop.dat'])
                fTempChiSq=fTempChiSq2
        if fTempChiSq>fBlockChiSq:
            subprocess.call(['cp','pop_rspy.dat','pop.dat'])
            fTempChiSq=fBlockChiSq


        print 'old ChiSq: %g new ChiSq: %g' % (fOldChiSq, fTempChiSq)
        if (fOldChiSq-fTempChiSq)<convergence:
            iFinishCounter+=1
        else:
            iFinishCounter=0

        if iFinishCounter==2:
            break

        fOldChiSq=fTempChiSq

#-------------------------------------------------------------------------------
def AvgProfile():

    iCounter=0
    while 1:
        sfilename='ContDimRho'+str(iCounter)+'.dat'
        if os.path.isfile(sfilename):
            file = open(sfilename,'r')
            data=file.readlines()
            file.close()
            rho=[]
            for i in range(len(data)):
                rho.append(float(data[i]))

            sfilename='ContDimZ'+str(iCounter)+'.dat'
            file = open(sfilename,'r')
            data=file.readlines()
            file.close()
            z=[]
            for i in range(len(data)):
                z.append(float(data[i]))

            ztemp=[]                                                                                #correct for Igor-style array axis notation
            for i in range(len(z)-1):
                ztemp.append((z[i]+z[i+1])/2)
            z=ztemp

            sfilename='ContArray'+str(iCounter)+'.dat'
            file = open(sfilename,'r')
            data=file.readlines()
            file.close()
            arr=[]
            for i in range(len(data)):
                rhoarr=[]
                tdata=string.split(data[i])
                for j in range(len(tdata)):
                    rhoarr.append(float(tdata[j]))
                arr.append(rhoarr)
            #print len(tdata), len(rho), len(z), len(data)

            median=[]; maxlikely=[]; lowperc=[]; highperc=[]
            for i in range(len(z)):
                cumulative=numpy.cumsum(arr[i])
                mediancum=cumulative[-1]*0.5
                lowperccum=cumulative[-1]*0.17
                highperccum=cumulative[-1]*0.83
                for j in range(len(arr[i])-1):
                    if cumulative[j] <= mediancum <= cumulative[j+1]:
                        frac=(mediancum-cumulative[j])/(cumulative[j+1]-cumulative[j])
                        median.append(rho[j]*frac+rho[j-1]*(1-frac))
                    if cumulative[j] <= lowperccum <= cumulative[j+1]:
                        frac=(lowperccum-cumulative[j])/(cumulative[j+1]-cumulative[j])
                        lowperc.append(rho[j]*frac+rho[j-1]*(1-frac))
                    if cumulative[j] <= highperccum <= cumulative[j+1]:
                        frac=(highperccum-cumulative[j])/(cumulative[j+1]-cumulative[j])
                        highperc.append(rho[j]*frac+rho[j-1]*(1-frac))
                    if max(arr[i])==arr[i][j]:
                        maxlikely.append(rho[j])

            sfilename='AvgProfile'+str(iCounter)+'.dat'
            file=open(sfilename,"w")
            file.write('z    maxlikely    median    lowperc     highperc \n')
            for i in range(len(z)):
                file.write(str(z[i])+' '+str(maxlikely[i])+' '+str(median[i])+' '+str(lowperc[i])+' '+str(highperc[i])+'\n')
            file.close()
            iCounter += 1
        else:
            break


#-------------------------------------------------------------------------------

def fContour(dZGrid=0.5, dRhoGrid=1e-8, sStatMode='is'):                                         #Calculate contour plot data from SErr.dat
    ReflPar.fnContourData(sStatMode+'Err.dat',dZGrid, dRhoGrid)

#-------------------------------------------------------------------------------
def fCalculateMolgroups(fConfidence):
    ReflPar.fnCalculateMolgroupProperty(fConfidence)
#-------------------------------------------------------------------------------
def fnDetermineFitSoftware():

    if os.path.isfile('run.py'):
        print 'Refl1D setup identified.'
        return 'refl1d'
    else:
        print 'Garefl setup identified.'
        return 'garefl'
#-------------------------------------------------------------------------------
def fnDetermineStatFile():

    if os.path.isfile('isErr.dat'):
        return 'isErr.dat'
    elif os.path.isfile('sErr.dat'):
        return 'sErr.dat'
    else:
        return ''
#-------------------------------------------------------------------------------

def DErr(convergence=0.001):                                                    #function has not yet been thoroughly tested
                                                                                #Error analysis by stepwise parameter displacement,
                                                                                #fixing, and fitting within a defined range
    def fnWriteToFile(sParameterName, fTestValue, fChiSq):                      #append new data set to file
        try:
            file = open('Error_'+sParameterName+'.dat','r')                     #if file exists, open and read data
            data = file.readlines()
            file.close()
        except IOError:
            data = []                                                           #otherwise start with an empty data file
        newdata=data[:]                                                         #copy dat into new object
        newdata.append(str(fTestValue)+"  "+str(fChiSq)+'\n')                   #append line with the data parameter to store
        file = open('Error_'+sParameterName+'.dat','w')                         #create new file and write out the new data
        file.writelines(newdata)
        file.close()

    ReflPar.fnLoadAndPrintPar()                                                 #print latest fit parameters
    tTaggedParameters=ReflPar.fnGetTaggedParameters()                           #get a list of all tagged parameters and test ranges
    ReflPar.fnBackup()                                                          #backup working directory
    try:
        for tParameter in tTaggedParameters:                                    #iterate through all parameters
            fValue=ReflPar.fnGetParameterValue(tParameter[0])                   #get fitted par value, fit constraints,
            sParameterName=tParameter[0]
            fLowerLimit=float(tParameter[1])
            fUpperLimit=float(tParameter[2])
            fTestLLimit=float(tParameter[3])                                    #test ranges
            fTestULimit=float(tParameter[4])
            fTestStep=float(tParameter[5])
            for iTestPoint in range(int(fTestLLimit/fTestStep),                 #iterate through test points, which are centered
                                    int(fTestULimit/fTestStep)):                #around the fit value
                fTestValue=fValue+float(iTestPoint)*fTestStep
                if fTestValue<fLowerLimit or fTestValue>fUpperLimit:
                    continue
                print sParameterName, fTestValue
                ReflPar.fnReplaceParameterLimitsInSetup(sParameterName,         #replace fit constraint in setup.c with a range, which
                                    0.9999*fTestValue, 1.0001*fTestValue)       #is
                pr=subprocess.Popen(["make"])
                pr.wait()
                AutoFinish(convergence)                                         #fit using the stored minimum
                pr=subprocess.Popen(["cp","pop_bak.dat","pop.dat"])
                pr.wait()
                AutoFinish2(convergence)                                        #refit, because sometimes the fit gets trapped
                pr=subprocess.Popen(["cp","pop_bak.dat","pop.dat"])
                pr.wait()
                ReflPar.fnLoadParameters()
                fnWriteToFile(sParameterName, fTestValue, ReflPar.fnGetChiSq())
            ReflPar.fnRestoreBackup()
    finally:
        ReflPar.fnRemoveBackup()                                                #in case of crash restore working dir

#-------------------------------------------------------------------------------

def fEnvelope(fZGrid=0.1,fSigma=0,sStatMode='is',bShortFlag=False, iContrast=-1):   #Calculate nSLD envelopes plot data from SErr.dat
    ReflPar.fnnSLDEnvelopes(fZGrid,fSigma,sStatMode+'Err.dat',bShortFlag, iContrast)

#-------------------------------------------------------------------------------


def fMCMultiCore(iIterations=1000, fMCConvergence=0.01, iConvergence=0.01,
    fConfidence=0.9546, sMode='is',iNodes=1, iIterationsPerCall='10',sJobID=''):#Monte Carlo for Multicore systems

    def fCleanUpDirectory(sDir):
            subprocess.call(['rm', '-rf',sDir])                                 #remove rs directory

    def fCreateRsDir(sMode,lSubProcesses):
        while 1:                                                                #create random name directory
            sRandomName=str(int(random.random()*1000000000000))                 #and check if the name already
            sDirName='../rspy'+sRandomName                                      #exists as a Multicore DirName
            iNameExists=0
            for element in lSubProcesses:
                if sDirName==element[1]:
                    iNameExists=1
            if not iNameExists==1:
                break

        subprocess.call('mkdir '+sDirName,shell=True)

        subprocess.call('cp * '+sDirName+'/',shell=True)                        #copy whole main directory over
        subprocess.call('rm -f '+sDirName+'/'+sMode+'Err.dat',shell=True)       #not including the stat data calculated
                                                                                 #so far

        if os.path.isfile(sDirName+'/'+'StatFinished'):
            subprocess.call(["rm","-f",sDirName+'/'+"StatFinished"])
        if os.path.isfile(sDirName+'/'+'StatAborted'):
            subprocess.call(["rm","-f",sDirName+'/'+"StatAborted"])
        return sDirName


    def fKillAllProcesses():                                                    #kills all running processes
        for item in lSubProcesses:
            try:
                os.kill(int(item[0]),15)
            except:
                pass
            print 'Delete directories ...'
            fCleanUpDirectory(item[1])


    lSubProcesses=[]
    iDoneIterations=0
    iChange=0
    fTimeAverage=0.
    iExitFlag=0
    #ReflPar.fnCheckFit()



    try:

        while (iDoneIterations < iIterations) and (iExitFlag==0):
                                                                                #start new subprocesses
            while (iDoneIterations<(iIterations-len(lSubProcesses)*iIterationsPerCall)) and (len(lSubProcesses)<iNodes):
                sDirName = fCreateRsDir(sMode,lSubProcesses)                    #sDirName is name of directory for Multicore architecture
                                                                                #and a Filepointer for the PBS architecture

                pid=str(subprocess.Popen([sDirName+'/rs.py',
                    '-'+sMode,str(iIterationsPerCall),str(iConvergence)],
                    cwd=sDirName, stdout=open(os.devnull,"w"), stderr=open(os.devnull,"w")).pid)
                lSubProcesses.append((pid,sDirName, time.time()))
                iChange=1
                time.sleep(2)                                                   #wait for ga_refl random number generator


            if iChange==1:                                                      #is there any change to report?
                iChange=0
                fActualTime=time.time()
                fProjectedTime=0
                fInProcessTime=0
                for element in lSubProcesses:
                    fInProcessTime=fInProcessTime+(fActualTime-element[2])      #add up time already spent in not finished proc.

                fProjectedTime=(iIterations-iDoneIterations)*fTimeAverage       #total cpu time for remaining iterations
                fProjectedTime -= fInProcessTime                                #credit for total cpu time alrady done
                fProjectedTime /= len(lSubProcesses)                            #get real time

                if fProjectedTime<0:
                    fProjectedTime=0
                lTD=time.gmtime(fProjectedTime)
                lTA=time.gmtime(fTimeAverage)
                print ''
                print time.ctime(fActualTime)
                print '-%s Monte Carlo Error analysis using %i processes.' % (sMode, iNodes)
                print 'Computing %i iterations per call' % iIterationsPerCall
                print '%i of %i iterations done.' % (iDoneIterations, iIterations)
                print ''
                print 'Process list:'
                for i in range(len(lSubProcesses)):
                    lTS=time.gmtime(fActualTime-lSubProcesses[i][2])
                    print '%i: PID %s in directory %s running for %id %ih %imin %is' % (i+1, lSubProcesses[i][0],
                        lSubProcesses[i][1], lTS[2]-1, lTS[3], lTS[4], lTS[5])
                if fTimeAverage>0:
                    print ''
                    print 'Average time per iteration: %id %ih %imin %is' % (lTA[2]-1, lTA[3], lTA[4], lTA[5])
                    print 'Projected finish in %id %ih %imin %is' % (lTD[2]-1, lTD[3], lTD[4], lTD[5])
                    print 'on %s' % (time.ctime(fActualTime+fProjectedTime))
                print ''
                print ''

            time.sleep(30)                                                      #wait before checking for finished subprocesses

            while 1:
                iFinishedProcessFound=0

                for i in range(len(lSubProcesses)):                             #check for finished sub processes
                    sDirName=lSubProcesses[i][1]
                    pid=lSubProcesses[i][0]

                    if os.path.isfile(sDirName+'/'+'StatFinished'):             #look up if process is finished
                        if os.path.isfile(sMode+'Err.dat'):
                            file=open(sDirName+'/'+sMode+'Err.dat','r')         #get statistical data from rs subdir
                            data=file.readlines()
                            file.close()
                            data=data[1:]                                       #delete headerline
                            iFailureCounter=0
                            while 1:
                                try:
                                    file=open(sMode+'Err.dat','a')              #append new data to file in root dir
                                    file.writelines(data)
                                    file.close()
                                    break
                                except:
                                    print '(i)sErr.dat is in use. Wait for 2s'
                                    iFailureCounter += 1
                                    time.sleep(2)
                                    if iFailureCounter==15:
                                        print 'Cannot append to (i)sErr.dat -> Abort.'
                                        break
                        else:
                            subprocess.call(['cp',sDirName+'/'+sMode+'Err.dat','.'])


                        fCleanUpDirectory(sDirName)

                        iDoneIterations=iDoneIterations+iIterationsPerCall
                        iChange=1
                        iFinishedProcessFound=1
                        fDeltaTime=time.time()-lSubProcesses[i][2]
                        fTimeAverage=fTimeAverage*(float(iDoneIterations-
                            iIterationsPerCall))/float(iDoneIterations)+fDeltaTime/float(iDoneIterations)
                        del lSubProcesses[i]                                    #remove entry from list

                        try:
                            ReflPar.fnAnalyzeStatFile(fConfidence)              #see if convergence criterium for whole MC had been reached
                            if ReflPar.diStatResults['Convergence']<=fMCConvergence:
                                print 'MC has converged ..'
                                iExitFlag=1
                        except:
                            print 'Analysis failed...'


                        break                                                   #because we changed len(lSubProcesses)

                    if os.path.isfile(sDirName+'/'+'StatAborted'):              #look up if process is finished
                        fCleanUpDirectory(sDirName)
                        if (time.time()-lSubProcesses[i][2])<180:
                            iExitFlag=1
                            print '=========Multicore Error========='
                            print 'Process termination within 3 min.'
                            print '================================='
                        del lSubProcesses[i]                                    #remove entry from list
                        iChange=1
                        iFinishedProcessFound=1
                        break                                                   #because we changed len(lSubProcesses)



                if (iFinishedProcessFound==0) or (len(lSubProcesses)==0) or (iExitFlag==1):
                    break

        if not sJobID=='':
            file=open(sJobID,'w')                                               #indicate Multicore is finished if required
            file.write('Multicore finished \n')                                 #requirement comes from presence of a
            file.close()                                                        #job id as given by a PBS caller


    finally:
        print 'Exiting MC'
        print iExitFlag, iDoneIterations, iIterations
        fKillAllProcesses()


#-------------------------------------------------------------------------------
def fMCMC(iMaxIterations=1024000, liMolgroups=['protein'], fSparse=0):

    while True:

        if not os.path.isfile('run.py'):                                              #make sure there is a run.py
            file=open('run.py','w')
            file.write('from bumps.fitproblem import FitProblem\n')
            file.write('from refl1d import garefl\n')
            file.write('from refl1d.names import *\n')
            file.write('\n')
            file.write("problem = garefl.load('model.so')\n")
            file.write('\n')
            file.write('problem.penalty_limit = 50\n')
            file.close()

        iBurn=4000                                                                  #check wether an existing MCMC exists
        bMCMCexists=False
        for i in range(1,9):
            iBurn=iBurn*2
            if os.path.isdir('MCMC_'+str(iBurn)+'_500'):
                print 'Found '+'MCMC_'+str(iBurn)+'_500 \n'
                bMCMCexists=True
                break

        lCommand=['refl1d','run.py','--fit=dream','--parallel']
        if bMCMCexists:
            if iBurn>=iMaxIterations:
                print 'Maximum number of MCMC iterations reached\n'
                break                                                               #end
            lCommand.append('--resume=MCMC_'+str(iBurn)+'_500')
            lCommand.append('--store=MCMC_'+str(iBurn*2)+'_500')
            lCommand.append('--burn='+str(iBurn))
        else:
            lCommand.append('--init=lhs')
            lCommand.append('--store=MCMC_8000_500')
            lCommand.append('--burn=8000')
            iBurn=4000

        lCommand.append('--steps=500')

        subprocess.call(lCommand)                                                   #run MCMC

        os.rename('MCMC_'+str(iBurn*2)+'_500','MCMC')                               #create sErr.dat
        if os.path.isfile('isErr.dat'):
            os.remove('isErr.dat')
        if os.path.isfile('sErr.dat'):
            os.remove('sErr.dat')
        StatAnalysis(-1,0.005)                                                      #sErr.dat contains about 1000 iterations
        os.rename('MCMC','MCMC_'+str(iBurn*2)+'_500')

        if liMolgroups<>[]:
            ReflPar.fnPullMolgroup(liMolgroups,0)

        if bMCMCexists:
            shutil.rmtree('MCMC_'+str(iBurn)+'_500')

    return



#-------------------------------------------------------------------------------

def fMCPBS(iIterations, iConvergence=0.01, fConfidence=0.9546, sMode='is',iNodes=1,
    iIterationsPerCall='10', sQueue='default',iNumberOfPBSJobsInQueue=1):       #Monte Carlo for PBS batch system submission
                                                                                #fConfidence not yet used


    def fCleanUpDirectory(sJobID,iNotOutputFiles=0):
        subprocess.call(['rm','-f','run'+sJobID+".sh"])                         #delete run.sh files
        subprocess.call(['rm','-f',sJobID])                                     #delete file pointer for sucessfull PBS run

        if iNotOutputFiles==0:                                                  #and output files
            sOutName=sJobID+'.out'
            sErrName=sJobID+'.err'
            iSleepCounter=0
            while 1:
                if os.path.isfile(sOutName) and os.path.isfile(sErrName):       #wait for output files
                    subprocess.call(['rm','-f',sOutName])
                    subprocess.call(['rm','-f',sErrName])
                    break
                retcode=subprocess.call('ls '+os.path.expanduser('~/')+'*.OU',
                    shell=True, stdout=open(os.devnull,"w"))
                if retcode==0:                                                  #on some systems only
                    break                                                       #standardfiles are written
                time.sleep(1)
                iSleepCounter=iSleepCounter+1
                if iSleepCounter > 20:
                    print 'Waited 20s for output files to be written ... giving up.'
                    break

    def fCreateMultiCoreJobID(lSubProcesses):
        while 1:                                                                #create random name directory
            sRandomName=str(int(random.random()*1000000000000))                 #and check if the name already
            iNameExists=0                                                       #as a PBS Job ID
            for element in lSubProcesses:
                if sRandomName==element[1]:
                    iNameExists=1
            if not iNameExists==1:
                break
        return sRandomName                                                      #randomname will be job id


    def fKillAllProcesses():                                                    #kills all running processes
        for item in lSubProcesses:
            subprocess.call(['qdel',item[0]])                                   #delete submitted job
            time.sleep(2)                                                       #give system time
            print 'Delete directories ...'
            fCleanUpDirectory(item[1])


    lSubProcesses=[]
    iDoneIterations=0
    iChange=0
    fTimeAverage=0.
    iExitFlag=0


    try:

        while (iDoneIterations < iIterations) and (iExitFlag==0):
                                                                                #start new subprocesses
            while (iDoneIterations<(iIterations-len(lSubProcesses)*iIterationsPerCall)) and (len(lSubProcesses)<iNumberOfPBSJobsInQueue):
                sJobID = fCreateMultiCoreJobID(lSubProcesses)                   #sDirName is name of directory for Multicore architecture
                                                                                #and a Filepointer for the PBS architecture
                data=[]                                                         #create run batchfile
                #data.append('#PBS -l nodes=1:ppn=1\n')
                data.append('#PBS -l ncpus=1\n')
                data.append('#PBS -l cput=48:00:00\n')
                data.append('#PBS -l walltime=72:00:00\n')
                data.append('#PBS -e '+os.getcwd()+'/'+sJobID+'.err\n')
                data.append('#PBS -o '+os.getcwd()+'/'+sJobID+'.out\n')
                data.append('cd $PBS_O_WORKDIR\n')
                data.append('./rs.py -'+sMode+' '+str(iIterationsPerCall)+' '
                    +'-c'+' '+str(iConvergence)+' '+'-m'+' '+str(iNodes)+' '
                    +'-ipc'+' '+str(iIterationsPerCall/iNodes)+' '+'-id'+' '+sJobID+'\n')
                data.append('#end')
                runfile='run'+sDirName+".sh"
                file=open(runfile,"w")
                file.writelines(data)
                file.close()
                pid = subprocess.Popen(['qsub','-q',sQueue,runfile],
                    stdout==open(os.devnull,"w")).communicate()[0]              #ged pid from queue systm
                pid = pid.split()[0]                                            #remove newline at end of output

                lSubProcesses.append((pid,sJobID, time.time()))
                iChange=1
                time.sleep(2)                                                   #wait for ga_refl random number generator


            if iChange==1:                                                      #is there any change to report?
                iChange=0
                fActualTime=time.time()
                fProjectedTime=0
                fInProcessTime=0
                for element in lSubProcesses:
                    fInProcessTime=fInProcessTime+(fActualTime-element[2])      #add up time already spent in not finished proc.

                fProjectedTime=(iIterations-iDoneIterations)*fTimeAverage       #total cpu time for remaining iterations
                fProjectedTime=fProjectedTime-fInProcessTime                    #credit for total cpu time alrady done
                fProjectedTime=fProjectedTime/len(lSubProcesses)                #get real time

                if fProjectedTime<0:
                    fProjectedTime=0
                lTD=time.gmtime(fProjectedTime)
                lTA=time.gmtime(fTimeAverage)
                print ''
                print time.ctime(fActualTime)
                print '-%s Monte Carlo Error analysis using %i PBS jobs in queue.' % (sMode, iNumberOfPBSJobsInQueue)
                print 'Computing %i iterations per call' % iIterationsPerCall
                print '%i of %i iterations done.' % (iDoneIterations, iIterations)
                print ''
                print 'Process list:'
                for i in range(len(lSubProcesses)):
                    lTS=time.gmtime(fActualTime-lSubProcesses[i][2])
                    print '%i: PID %s rs.py ID %s running for %id %ih %imin %is' % (i+1, lSubProcesses[i][0],
                        lSubProcesses[i][1], lTS[2]-1, lTS[3], lTS[4], lTS[5])
                if fTimeAverage>0:
                    print ''
                    print 'Average time per iteration: %id %ih %imin %is' % (lTA[2]-1, lTA[3], lTA[4], lTA[5])
                    print 'Projected finish in %id %ih %imin %is' % (lTD[2]-1, lTD[3], lTD[4], lTD[5])
                    print 'on %s' % (time.ctime(fActualTime+fProjectedTime))
                print ''
                print ''

            time.sleep(30)                                                      #wait before checking for finished subprocesses

            while 1:

                iFinishedProcessFound=0

                for i in range(len(lSubProcesses)):                             #check for finished sub processes
                    sJobID=lSubProcesses[i][1]
                    pid=lSubProcesses[i][0]

                    if os.path.isfile(sJobID):                                  #look up if process is finished
                        fCleanUpDirectory(sJobID)

                        iDoneIterations=iDoneIterations+iIterationsPerCall
                        iChange=1
                        iFinishedProcessFound=1
                        fDeltaTime=time.time()-lSubProcesses[i][2]
                        fTimeAverage=fTimeAverage*(float(iDoneIterations-
                            iIterationsPerCall))/float(iDoneIterations)+fDeltaTime/float(iDoneIterations)
                        del lSubProcesses[i]                                    #remove entry from list
                        break                                                   #because we changed len(lSubProcesses)

                                                                                #check for incorrectly finished subprocesses
                    retcode=subprocess.call(['qstat',pid],
                        stdout=open(os.devnull,"w"))
                    if retcode!=0:                                              #check for orphaned processes
                        fCleanUpDirectory(sJobID,1)
                        if (time.time()-lSubProcesses[i][2])<180:
                            iExitFlag=1
                            print '============PBS Error============'
                            print 'Process termination within 3 min.'
                            print '================================='
                        del lSubProcesses[i]                                    #remove entry from list
                        iChange=1
                        iFinishedProcessFound=1
                        break                                                   #because we changed len(lSubProcesses)


                if (iFinishedProcessFound==0) or (len(lSubProcesses)==0) or (iExitFlag==1):
                    break

    finally:
        fKillAllProcesses()
        subprocess.call('rm -f '+os.path.expanduser('~/')+'*.OU',
            shell=True, stdout=open(os.devnull,"w"))                            #delete standard files
        subprocess.call('rm -f '+os.path.expanduser('~/')+'*.ER',
            shell=True, stdout=open(os.devnull,"w"))                            #(on some systems)
                                                                                #this has to go through shell



#-------------------------------------------------------------------------------
def DisplayMolgroups(plotname):                                                 #loads parameters and covariance matrix and prints all
    ReflPar.fnPlotMolgroups(plotname)

#-------------------------------------------------------------------------------
def DisplayFit(plotname):
    ReflPar.fnPlotFit(plotname)

#-------------------------------------------------------------------------------
def Monitor():                                                                  #monitors all profile.dat files
                                                                                #in a given directory and displays
                                                                                #the nSLD profiles
    import Tkinter as Tk
    import Image            #PIL
    import ImageTk          #PIL

    class MyApp:

        def __init__(self,root):
            self.root=root                                                      #stores the root pointer

            self.frame=Frame(root)                                              #create Frame in GUI
            self.frame.pack()
            self.modtime=os.path.getmtime('profile0.dat')                       #get initial modification time
                                                                                #for profile0.dat
            print 'Monitoring ... Stop with Ctrl-C'

            self.fnMainBody(1)                                                  #draw initial nSLD profile

        def fnMainBody(self,iInitial=0):

            acttime=os.path.getmtime('profile0.dat')
            if (acttime!=self.modtime) or (iInitial==1):                        #check if profile has been updated or
                                                                                #if this is initial run

                self.modtime=acttime                                            #update latest modification time

                pylab.title('nSLD profile')

                iCounter=0
                while 1:
                    sFileName='profile'+str(iCounter)+'.dat'
                    if os.path.isfile(sFileName):                               #load profile if existent
                        file=open(sFileName,'r')
                        data=file.readlines()
                        file.close()
                        data = data[1:]                                         #delete header

                        k=0; l=0
                        zlist=[]; rholist=[]
                        for line in data:                                       #extract nSLD profile data line by line
                            splitline=line.split()
                            zlist.append(float(splitline[0]))
                            rholist.append(float(splitline[1])*1e6)

                        pylab.plot(zlist,rholist,label='profile'+str(iCounter)) #plot

                        iCounter += 1

                    else:
                        break

                pylab.legend(loc=0)
                pylab.xlabel('z / Angstroem')
                pylab.ylabel('rho / 1e-6 Angstroem^-2')
                pylab.savefig('imgdata.png',format='png')                       #save plot as file
                pylab.close()


                im=Image.open('imgdata.png')                                    #display file with TkInter
                self.photo = ImageTk.PhotoImage(im)

                if iInitial==1:                                                 #check if previous image exists
                                                                                #no: create canvas
                    self.canvas = Canvas(self.frame, height=im.size[1]+20, width=im.size[0]+20)
                    self.canvas.pack(side=LEFT,fill=BOTH,expand=1)
                else:                                                           #yes: delete previous image
                    self.canvas.delete(self.item)

                                                                                #attach image
                self.item = self.canvas.create_image(10,10,anchor=NW, image=self.photo)
                iInitial=0

            root.update()                                                       #update plot
            root.after(500,self.fnMainBody)                                     #Call this function again after
                                                                                #500 ms







    root=Tk.Tk()
    myapp=MyApp(root)
    root.mainloop()
#-------------------------------------------------------------------------------


def Result():                                                                   # loads parameters and covariance matrix and prints all
    ReflPar.fnLoadAndPrintPar()


#-------------------------------------------------------------------------------
def SErr(iiterations, convergence=0.01, sStatMode='is', fConfidence=0.9546):
    try:
        ReflPar.fnLoadAndPrintPar()
        ReflPar.fnBackup()
        filelist=ReflPar.fnLoadFileListAndChangeToLocal()
        ReflPar.fnMake()
    except:
        file=open('StatAborted','w')                                            #indicate that all iterations are done
        file.write('statistical analysis aborted \n')
        file.close()
        return ()
    try:
        if not os.path.isfile(sStatMode+'Err.dat'):
            file=open(sStatMode+'Err.dat','w')
            file.write("Chisq "+" ".join(ReflPar.fnGetSortedParNames())+'\n')
            file.close()
        for i in range(iiterations):
            print 'Iteration #', i
            ReflPar.fnModifyAndCopyFiles(filelist)
            if sStatMode=="is":
                Auto(convergence)                                               #automatic fitting, independent start
            elif sStatMode=='s':
                AutoFinish2(convergence)                                        #automatic fitting dependent start
            if os.path.isfile('pop_bak.dat'):
                subprocess.call(['cp','pop_bak.dat','pop.dat'])                 #copy newest population into
                                                                                #pop.dat, this is done here because
                                                                                #sometimes AutoFinish does not find
                                                                                #a better fit with the Levenberg-M.
                                                                                #and then by default, the generated
                                                                                #population is not used, but we want
                                                                                #to have it here for the statistical
                                                                                #analysis
            ReflPar.fnLoadParameters()
            chisq=ReflPar.fnGetChiSq()
            file = open(sStatMode+'Err.dat','a')
            file.write(str(chisq)+" "+" ".join(ReflPar.fnGetSortedParValues())+'\n')
            file.close()

        file=open('StatFinished','w')                                           #indicate that all iterations are done
        file.write('statistical analysis finished \n')
        file.close()
    except:
        file=open('StatAborted','w')                                            #indicate that all iterations are done
        file.write('statistical analysis aborted \n')
        file.close()
    finally:
        ReflPar.fnRemoveBackup()                                                # always restore working dir,
        subprocess.call(['rm','-f','*.mce'])                                    # also in case of crash
        ReflPar.fnTruncateParDat()

#-------------------------------------------------------------------------------
def StatAnalysis(fConfidence=0.9546,sparse=0):                                  # Analyzes MC output
    ReflPar.fnAnalyzeStatFile(fConfidence,sparse)

#-------------------------------------------------------------------------------
def StatTable(sTableName,fConfidence=0.9546):                                   # Produces Latex formatted table
    ReflPar.fnStatTable(sTableName,fConfidence)                                 # and the actual statistical data

#-------------------------------------------------------------------------------






# main programm from command line
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv)==1:
        print ''
        print 'Reflscript usage:'
        print '-----------------------------------------------------------------'
        print '-a [c]               Auto fit'
        print '-aprof               Creates average profiles from ContArray'
        print '                     ContDimZ and ContDimRho files'
        print '-conf cf             Confidence for statistical calculatioins'
        print '                     0<=cf<=1, if cf < 0 than cf is interpreted'
        print '                     in units of (negative) sigma'
        print '-cont [z rho]        Create a contour/image plot from sErr.dat'
        print '                     or isErr.dat'
        print '-d [c]               Displacement error analysis'
        print '-env [zGrid sigma]   Calculates nSLD envelopes from stat file'
        print '-f [c]               Finish fit'
        print '-fit filename        Save fit, rho, res in a png file'
        print '                     called filename_fit.png'
        print '-ipc [ipc]           iterations per call/job for -pbs and -m'
        print '-is n [c]            Independent statistical error analysis'
        print '-l [c]               Auto fit finish with correct roughness'
        print '-m [m]               use m parallel processes on one node'
        print '-mol filename        Save molgroups plot in a png file'
        print '                     called filename_mol.png'
        print '-mon                 Monitors a running fit'
        print '-pbs [j workq]       Use pbs to submit j jobs to workq queue'
        print '-pull molgroups      Creates stat profile for a number of molgroups'
        print '-r                   Print results'
        print '-s n [c]             Statistical error analysis with n iterations'
        print '-stat                Analysis of previously calculated stat. data'
        print '                     For n<1 a maximum of 1000 iterations is'
        print '                     carried out until the MC converges and all'
        print '                     parameters do not show a relative change of'
        print '                     more than n compared to the fit interval'
        print '-sparse n            Uses only a subset of n iterations from stat'
        print '                     file. If n<1 than a random chance of f is '
        print '                     applied to each line that it is used.'
        print '-t                   Truncate par.dat'
        print ''
        print 'Displacement Error Analysis:'
        print 'For each parameter that should be analyzed, the script expects'
        print 'a tag with the syntax"!rstag min max s !" in the line where the'
        print 'parameter is initialized with the pars_add command. The tag'
        print 'should be embedded in a comment section starting with //. The '
        print 'tag parameters min, max, and s give the relative(!) range and'
        print 'stepsize of the displacement of the parameters. For each '
        print 'parameter, an output file is created that contains the fixed '
        print 'parameter value and the achieved chi squared'
        print ''
        print 'Statistical Error Analysis'
        print 'The script creates n data sets by varying the measured'
        print 'reflectivity using random-generated normal deviates'
        print 'applied to the uncertainty of the measured data points'
        print 'An automated fit determines the fit parameters and they'
        print 'are stored in a file called "sErr.dat". A histogram analysis'
        print 'should be carried out afterwards with a software like Igor.'
        print ''
        print 'Independent Statistical Error Analysis'
        print 'This option is equal to "-s" but does not load the stored'
        print 'population of a previous fit at start-up. The fit parameters'
        print 'are stored in isErr.dat.'
        print ''
        print 'Multiprocessor support'
        print 'Designed for workstations with multicore architecture. The'
        print 'attribute m defines the number of parallel processes used.'
        print ''
        print 'PBS batch support'
        print 'The j attribute defines the number of jobs submitted once'
        print 'at a time. The workq argument is the queue name. PBS and the'
        print 'Multi-process support -m can be combined. In this case the'
        print '-ipc option defines the number of MC iterations per process'
        print 'The number of iterations per PBS job is then ips * m with'
        print 'm being the number of parallel processes.'
        print ''
        print 'Contour/Image Plots'
        print 'Using SErr.dat or iSErr.dat an array with 3D-data usable'
        print 'for contour or image plots is created. The output file is'
        print 'ContArr#model.dat. 3D data is generated for all models'
        print 'specified in setup.c. For each array two files with axis'
        print 'label information are written intended to be used with Igor.'
        print 'The attributes z and rho determine the bin size of the'
        print 'grid used for mapping of the nSLD profiles.'
        print ''
        print 'All fit and error analysis calls may have an additional'
        print 'attribute c that sets the convergence condition of the fit'
        print 'The fit is viewed at as converged if the change of chi2 after'
        print 'a complete alternation over genetic algorigthm, amoeba and'
        print 'Levenberg Marquardt is less than c. A c of 0.001 is the'
        print 'default'
        print 'Example: ./rs.py -iS 1000 0.01 decreases the convergence'
        print 'condition by one order of magnitude'
        print 'Envelopes'
        print 'The option -short startes a less precise but significantly'
        print 'faster algorithm. The precision of the envelopes is'
        print '+/- 0.1 sigma. -env takes two arguments, the bin size'
        print 'for z and a sigma parameter. Presets are 0.5 in both cases'
        print 'The sigma parameter defines the spacing in units of '
        print 'sigma for which envelopes are calculated. A sigma parameter'
        print 'of 0 saves all calculates envelopes'
    else:

        sStatMode='none'
        iSummarizeStatistics=0
        sArchitecture=''
        sJobSubmission=''
        sWorkqueue='default'
        sJobID=''
        sTableTemplate='tabletemplate.tex'
        sMolgroup=''
        bShortFlag=False
        iCm=0
        iContour=0
        iContrast=-1
        iEnv=0
        fZGrid=0.5
        fSigma=0.5
        fRhoGrid=1e-8
        iIterationsPerCall=0
        iNumberOfMCMCIterations=64000
        iNumberOfParallelThreads=1
        iNumberOfPBSJobsInQueue=1
        iNumberOfPbsJobs=10
        iNumberOfMonteCarloIterations=1000
        iStatTable=0
        fConvergence=0.01
        fConfidence=0.9546
        fMCConvergence=0
        fSparse=0
        iSummarizeStatistic=0
        iPullMolgroups=0
        liMolgroups=[]

        ReflPar=CReflPar(CGaReflInteractor(),CRefl1DInteractor())

        for i in range(1,len(sys.argv)):
            arg=sys.argv[i]
            if sys.argv[i] == '-a' :
                if len(sys.argv)>i+1:
                    fConvergence=float(sys.argv[i+1])
                Auto(fConvergence)
            elif sys.argv[i] == '-aprof':
                AvgProfile()
            elif sys.argv[i] == '-cm':
                iCm=1
                sMolgroup=sys.argv[i+1]
                i += 1
            elif sys.argv[i] == '-conf':
                fConfidence=float(sys.argv[i+1])
                i += 1
            elif sys.argv[i] == '-contrast':
                    iContrast=int(sys.argv[i+1])
                    i+=1
            elif sys.argv[i] == '-f' :
                if len(sys.argv)>i+1:
                    fConvergence=float(sys.argv[i+1])
                AutoFinish(fConvergence)
            elif sys.argv[i] == '-fit' :
                import numpy as np
                import matplotlib.pyplot as plt
                if len(sys.argv)>i+1:
                    plotname=sys.argv[i+1]
                else:
                    plotname=''
                DisplayFit(plotname)
            elif sys.argv[i] == '-l' :
                if len(sys.argv)>i+1:
                    fConvergence=float(sys.argv[i+1])
                AutoFinish2(fConvergence)
            elif sys.argv[i] == '-mol' :
                import numpy as np
                import matplotlib.pyplot as plt
                if len(sys.argv)>i+1:
                    plotname=sys.argv[i+1]
                else:
                    plotname=''
                DisplayMolgroups(plotname)
            elif sys.argv[i] == '-mon' :
                from Tkinter import *
                import matplotlib
                import pylab
                Monitor()
            elif sys.argv[i] == '-t' :
                ReflPar.fnTruncateParDat()
            elif sys.argv[i] == '-cont' :
                if len(sys.argv)>i+2:
                    fZGrid=float(sys.argv[i+1])
                    fRhoGrid=float(sys.argv[i+2])
                    i += 2
                iContour=1
            elif sys.argv[i] == '-env' :
                iEnv=1
                try:
                    if len(sys.argv)>i+1:
                        fZGrid=float(sys.argv[i+1])
                        i += 1
                    if len(sys.argv)>i+1:
                        fSigma=float(sys.argv[i+1])
                        i += 1
                except:
                    pass
            elif sys.argv[i] == '-short':
                bShortFlag=True
            elif sys.argv[i] == '-sparse':
                fSparse=float(sys.argv[i+1])
                i+=1
            elif sys.argv[i] == '-d' :
                if len(sys.argv)>i+1:
                    fConvergence=float(sys.argv[i+1])
                DErr(fConvergence)
            elif sys.argv[i] == '-is':
                sStatMode='is'
                if len(sys.argv)>i+1:
                    temp=float(sys.argv[i+1])
                    if temp>=1:
                        iNumberOfMonteCarloIterations=int(sys.argv[i+1])
                        fMCConvergence=0
                    else:
                        iNumberOfMonteCarloIterations=1000
                        fMCConvergence=temp
                        iIterationsPerCall=1
                    i += 1
                if len(sys.argv)>i+1:
                    try:
                        fConvergence=float(sys.argv[i+1])
                        i += 1
                    except:
                        i=i
            elif sys.argv[i] == '-s':
                sStatMode='s'
                if len(sys.argv)>i+1:
                    temp=float(sys.argv[i+1])
                    if temp>=1:
                        iNumberOfMonteCarloIterations=int(sys.argv[i+1])
                        fMCConvergence=0
                    else:
                        iNumberOfMonteCarloIterations=1000
                        fMCConvergence=temp
                        iIterationsPerCall=1
                    i += 1
                if len(sys.argv)>i+1:
                    try:
                        fConvergence=float(sys.argv[i+1])
                        i += 1
                    except:
                        i=i
            elif sys.argv[i] == '-m':
                sArchitecture='Multicore'
                if len(sys.argv)>i+1:
                    iNumberOfParallelThreads=int(sys.argv[i+1])
                    i += 1
            elif sys.argv[i] == '-MCMC':
                if len(sys.argv)>i+1:
                    iNumberOfMCMCIterations=int(sys.argv[i+1])
                    i += 1
                fMCMC(iNumberOfMCMCIterations)
            elif sys.argv[i] == '-ipc':
                if len(sys.argv)>i+1:
                    iIterationsPerCall=int(sys.argv[i+1])
                    i += 1
            elif sys.argv[i] == '-id':
                if len(sys.argv)>i+1:
                    sJobID=sys.argv[i+1]
                    i += 1
            elif sys.argv[i] == '-c':
                if len(sys.argv)>i+1:
                    fConvergence=float(sys.argv[i+1])
                    i += 1
            elif sys.argv[i] == '-pbs':
                sJobSubmission='PBS'
                if len(sys.argv)>i+1:
                    iNumberOfPBSJobsInQueue=int(sys.argv[i+1])
                    i += 1
                if len(sys.argv)>i+1:
                    sWorkqueue=sys.argv[i+1]
                    i += 1
            elif sys.argv[i] == '-pull':
                iPullMolgroups=1
                liMolgroups=sys.argv[i+1:]
            elif sys.argv[i] == '-r' :
                Result()
            elif sys.argv[i] == '-stat':
                iSummarizeStatistic=1
            elif sys.argv[i] == '-stattable':
                iStatTable=1
                sTableTemplate = sys.argv[i+1]



        if iContour==1:
            if sStatMode=='none':                                               #Contour mode needs is based on a MC statistics
                sStatMode='is'                                                  #choose default
            fContour(fZGrid,fRhoGrid,sStatMode)
        elif iCm==1:
            fCalculateMolgroups(fConfidence)
        elif iEnv==1:
            if sStatMode=='none':                                               #Contour mode needs is based on a MC statistics
                sStatMode='is'                                                  #choose default
            fEnvelope(fZGrid,fSigma,sStatMode,bShortFlag,iContrast)
        elif sStatMode=='is' or sStatMode=='s':
            if sJobSubmission=='':
                if sArchitecture=='':
                    SErr(iNumberOfMonteCarloIterations,fConvergence,sStatMode,fConfidence)
                if sArchitecture=='Multicore':
                    if not iIterationsPerCall:
                        iIterationsPerCall=10
                    fMCMultiCore(iNumberOfMonteCarloIterations,fMCConvergence,fConvergence,fConfidence,sStatMode,
                    iNumberOfParallelThreads,iIterationsPerCall,sJobID)
            elif sJobSubmission=='PBS':
                if not iIterationsPerCall:
                    iIterationsPerCall=1
                if sArchitecture=='':
                    fMCPBS(iNumberOfMonteCarloIterations,fConvergence,fConfidence,sStatMode,
                        1,iIterationsPerCall,sWorkqueue,iNumberOfPBSJobsInQueue)
                if sArchitecture=='Multicore':
                    fMCPBS(iNumberOfMonteCarloIterations,fConvergence,fConfidence,sStatMode,
                        iNumberOfParallelThreads,iIterationsPerCall*iNumberOfParallelThreads,
                        sWorkqueue,iNumberOfPBSJobsInQueue)
        elif iSummarizeStatistic==1:
            StatAnalysis(fConfidence,fSparse)
        elif iStatTable==1:
            StatTable(sTableTemplate,fConfidence)
        elif iPullMolgroups:
            ReflPar.fnPullMolgroup(liMolgroups,fSparse)


