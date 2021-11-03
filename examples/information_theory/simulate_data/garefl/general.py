#!/usr/bin/env python
# -*- coding: utf-8 -*-
# general purpose for 16-Sep-09 F.H.

import math, os, string
import numpy
import scipy.special, scipy.optimize, scipy.stats
import matplotlib.pyplot as plt
import random
import copy


#-------------------------------------------------------------------------------
#General usefull routines

#appends column to a container and adds reference axis as well if not present
def fnAppend(data,container,columnname,referenceaxis):

    if referenceaxis not in container:
        container[referenceaxis]=data[referenceaxis][:]

    containercolumnname=columnname+str(len(container.keys()))
    container[containercolumnname]=data[columnname][:]


#Deletes leading and trailing zeros from entries in data while leaving the
#referenceaxis untouched. This will lead to non-uniform lengths of the lists
#contained in the dictionary
def fnCutZeros(data,referenceaxis):

    for element in data:
        if element != referenceaxis:
            while True:
                if not data[element][0]:
                    del data[element][0]
                else:
                    break
            while True:
                if not data[element][-1]:
                    del data[element][-1]
                else:
                    break



#Integrates columnname in data and returns the integral in data2
def fnIntegrate(data,columnname,referenceaxis):

    data2= {referenceaxis: data[referenceaxis][:], columnname: []}

    sum=0
    for i in range(len(data[columnname])):
        sum+=data[columnname][i]
        data2[columnname].append(sum)

    return data2



#Loads single column data into a dictionary {"columnname",[data]}
#it appends the data found in the file to the provided dictionary 'data'
#it will skip columns with names which are either in the exception list
#or appends data with name extension "_2nd" that are already present in the data list

def fnLoadSingleColumns(sFileName, data={}, exceptions=[], header=True, headerline=[], LoadList=[]):

    file  = open(sFileName,"r")
    content=file.readlines()
    file.close()

    if header==True:                                                            #if headerline in file, read it from there
        splitheaderline=content[0].split()
        content=content[1:]
    else:                                                                       #else use the one given as an attribute
        splitheaderline=headerline

    for i,columnname in enumerate(splitheaderline):
        if LoadList==[] or (columname in LoadList):
            data[columnname]=[]

    for line in content:
        splitline=line.split()
        if splitline!=[]:
            for i,column in enumerate(splitline):
                if LoadList==[] or (splitheaderline[i] in LoadList):
                    data[splitheaderline[i]].append(float(splitline[i]))

    return data

#Loads single row data into a dictionary ["rowname",[data]]
#it appends the data found in the file to the provided list 'data'
#it will skip rows with names which are either in the exception list
#or appends data with name extension "_2nd" that are already present in the data list

def fnLoadSingleRows(sFileName, data={}, exceptions=[]):

    file  = open(sFileName,"r")
    content=file.readlines()
    file.close()

    for line in content:
        splitline=line.split()
        if splitline[0] not in exceptions:
            data[splitline[0]]=[]
            for entry in range(1, len(splitline)):
                data[splitline[0]].append(splitline[entry])

    return data

#Normalizes integral of columnname in data to value, reference axis does not have to
#be equally spaced

def fnNormalizeIntegral(data,columnname,referenceaxis,value):

    sum=0
    for i in range(1,len(data[columnname])):
        sum+=(data[columnname][i]+data[columnname][i-1])/2*(data[referenceaxis][i]-data[referenceaxis][i-1])

    fCorr=value/sum

    for i in range(len(data[columnname])):
        data[columnname][i]*=fCorr

#normalizes columnname such that the last value becomes 'value'

def fnNormalizeLast(data, columnname, value):

    fCorr=value/data[columnname][-1]
    for i in range(len(data[columnname])):
        data[columnname][i]*=fCorr


#normalizes columnname such that the maximum value becomes 'value'

def fnNormalizeMax(data, columnname, value):

    maxvalue=data[columnname][0]
    for d in data[columnname]:
        if d>maxvalue:
            maxvalue=d

    fCorr=value/maxvalue
    for i in range(len(data[columnname])):
        data[columnname][i]*=fCorr


def fnSave2DArray(sFilename, array):

    file=open(sFilename,"w")

    for element1 in array:
        for element2 in element1:
            file.write(str(element2)+" ")
        file.write("\n")

    file.close()

#saves all data out to a file
def fnSaveSingleColumns(sFilename, data):

    file=open(sFilename,"w")

    for element in data:
        file.write(element+" ")
    file.write("\n")

    for i in range(len(data[data.keys()[0]])):
        for element in data:
            file.write(str(data[element][i])+" ")
        file.write("\n")

    file.close()

def fnStat(data,referenceaxis):

    statdata={referenceaxis : [], 'msigma' : [], 'median' : [], 'psigma' : []}

    for i in range(len(data[referenceaxis])):
        dl=[]
        for element in data.keys():
            if element != referenceaxis:
                dl.append(data[element][i])

        statdata[referenceaxis].append(data[referenceaxis][i])
        statdata['msigma'].append(scipy.stats.scoreatpercentile(dl,15.85))
        statdata['median'].append(scipy.stats.scoreatpercentile(dl,50.00))
        statdata['psigma'].append(scipy.stats.scoreatpercentile(dl,84.15))

    return statdata




def fnSkew(array):
    return scipy.stats.moment(array,moment=3)/scipy.stats.std(array)**3

def fnKurtosis(array):
    return scipy.stats.moment(array,moment=4)/scipy.stats.std(array)**4-3

def fnWeightedSdv(array,weightsarray):
    average=numpy.average(array,weights=weightsarray)
    sdv=0
    for i in range(len(array)):
        sdv=sdv+weightsarray[i]*(array[i]-average)*(array[i]-average)
        #print sdv, weightsarray[i], array[i], average
    sdv=sdv/numpy.sum(weightsarray)
    return numpy.sqrt(sdv)



#-------------------------------------------------------------------------------
#special routines
def fnASKS():                                                                   #computes the weighted average, weighted sdv,
                                                                                #kurtosis and skew
                                                                                #from a file containing columns of data for
                                                                                #each column the first column is treated as
                                                                                #the x-axis, the other columns are the weights
    data=[]
    data=fnLoadSingleColumns("Data.txt",data)
    result=[]
    for dataset in data:
        result.append([])
        result[-1].append(dataset[0])                                           #append file name first

        print dataset[0]

        result[-1].append(numpy.average(data[0][1],weights=dataset[1]))         #average column over z-axis data[0][1]
        result[-1].append(fnWeightedSdv(data[0][1],dataset[1]))
        result[-1].append(fnKurtosis(dataset[1]))
        result[-1].append(fnSkew(dataset[1]))

    file=open("result.txt",'w')
    for line in result:
        for element in line:
            file.write(str(element)+' ')
        file.write('\n')
    file.close()


#Loads all files in Directory and Adds a header
def fnAddHeader():

    data={}
    exceptions=[]
    headerline=["z","protnSL","deutnSL","area"]
    header=False

    for files in os.listdir("."):
        if files.endswith(".txt"):
            fnLoadSingleColumns(files,data,exceptions,header,headerline)
            fnSaveSingleColumns(files,data)

#Loads columns from file, rescales them by a factor fScaling, and saves them into a result file
#adds two data sets with exception of the parameters in exceptionadding and
#exceptionloading
def fnAdd(exceptionadding=["z","tot","tot_2nd"], exceptionloading=[]):

    data1={}; data2={}
    data=fnLoadSingleColumns("dmpc_liq_135ns_byatom_inner_wat.txt",data1, exceptionloading)
    data=fnLoadSingleColumns("dmpc_liq_135ns_byatom_outter_wat.txt",data2, exceptionloading)

    for element1 in data1:
        if element1 not in exceptionadding:
            for element2 in data2:
                if element2==element1:
                    for i in range(len(data[element1])):
                        data1[element1][i]=data1[element1][i]+data2[element2][i]
                    break

    fnSaveSingleColumns("dmpc_liq_135ns_byatom_added.txt", data1)

def fnBin2D(filename,list1,list2,start1,stop1,step1,start2,stop2,step2,Polar=False,Sigma=0,Deg=True, AvgVariable=[], Avg=False, fileext='LogL',NewTiltOrigin=0.0, NewOrientOrigin=0.0):

    print('Initialize Binning Array...')
    arr1=[]; arr2=[]                                          #initialize array with zeros
    avgarr1=[]; avgarr2=[]

    if Polar==False:
        for i in range(start1,stop1,step1):
            arr1.append([])
            for j in range(start2,stop2,step2):
                arr1[i].append(0.)
    else:                                           #polar coordinates come with two plots
        d=-90; i=0
        while d<=90:
            arr1.append([])
            arr2.append([])
            avgarr1.append([])
            avgarr2.append([])
            dd=-90; ii=0
            while dd<=90:
                arr1[i].append(0.)
                arr2[i].append(0.)
                avgarr1[i].append(0.)
                avgarr2[i].append(0.)
                dd+=step1
                ii+=1
            d+=step1
            i+=1
            #print(d,i)


    print('Binning ...')

    for i in range(len(list1)):
        f1=list1[i]                                                             #get floating point entries
        f2=list2[i]
                                                                               #apply correction
        f2=numpy.radians(f2)
        f20=numpy.radians(NewOrientOrigin)
        f1=numpy.radians(f1)
        f10=numpy.radians(NewTiltOrigin)

        RotZ=numpy.array([[numpy.cos(f2), -1*numpy.sin(f2),0],[numpy.sin(f2), numpy.cos(f2),0],[0,0,1]])
        RotX=numpy.array([[1,0,0],[0,numpy.cos(f1),-1*numpy.sin(f1)],[0,numpy.sin(f1),numpy.cos(f1)]])
        RotZ0=numpy.array([[numpy.cos(f20), -1*numpy.sin(f20),0],[numpy.sin(f20), numpy.cos(f20),0],[0,0,1]])
        RotX0=numpy.array([[1,0,0],[0,numpy.cos(f10),-1*numpy.sin(f10)],[0,numpy.sin(f10),numpy.cos(f10)]])

        RotProduct=numpy.dot(RotX,RotZ)
        RotProduct=numpy.dot(RotZ0,RotProduct)
        RotProduct=numpy.dot(RotX0,RotProduct)

        f2=numpy.arctan(RotProduct[2][0]/RotProduct[2][1])
        if f2<0:
            f2+=numpy.pi
        f2_2=f2+numpy.pi

        f1=numpy.arctan(RotProduct[2][1]/RotProduct[2][2]/numpy.cos(f2))
        if f1<0:
            f1+=numpy.pi
        f1_2=numpy.arctan(RotProduct[2][1]/RotProduct[2][2]/numpy.cos(f2_2))
        if f1_2<0:
            f1_2+=numpy.pi

        if round(numpy.cos(f1),10)==round(RotProduct[2][2],10):
            f2=numpy.degrees(f2)
            f1=numpy.degrees(f1)
        elif round(numpy.cos(f1_2),10)==round(RotProduct[2][2],10):
            f2=numpy.degrees(f2_2)
            f1=numpy.degrees(f1_2)
        else:
            print 'This is not good.'
            print('f1 %e f1_2 %e cos(f1) %e cos(f1_2) %e RotProduct[2][2]',f1,f1_2,numpy.cos(f1), numpy.cos(f1_2),RotProduct[2][2])
            f2=numpy.degrees(f2)
            f1=numpy.degrees(f1)



        if f2>=360:
            f2=f2-360
        if f2<0:
            f2=f2+360
        if f1<0:
            f1=f1*(-1)
            f2=f2+180
        if f1>=180:
            f1=180-(f1-180)
            f2=f2+180
        if f2>=360:
            f2=f2-360
        if f2<0:
            f2=f2+360

        bNegativeBeta=False
        if f1>90:
            f1=180-f1
            bNegativeBeta=True

        if Polar==False:
            n1=int((f1-start1)/step1+0.5)                                       #do the binning
            n2=int((f2-start2)/step2+0.5)
        else:
            fAngCorr=1
            if Deg==True:
                fAngCorr=2*numpy.pi/360
            n1=int((f1*numpy.cos(f2*fAngCorr))/step1+0.5)+int((90.)/step1)
            n2=int((f1*numpy.sin(f2*fAngCorr))/step2+0.5)+int((90.)/step1)

        #print(n1,n2)


        if Polar==True and bNegativeBeta==True:
            arr2[n1][n2]=arr2[n1][n2]+1.
        elif Polar==False or bNegativeBeta==False:
            arr1[n1][n2]=arr1[n1][n2]+1.

        if Avg:
            if Polar==True and bNegativeBeta==True:
                avgarr2[n1][n2]=((avgarr2[n1][n2]*(arr2[n1][n2]-1))+AvgVariable[i])/arr2[n1][n2]
            elif Polar==False or bNegativeBeta==False:
                avgarr1[n1][n2]=((avgarr1[n1][n2]*(arr1[n1][n2]-1))+AvgVariable[i])/arr1[n1][n2]





    fnSave2DArray(filename+"_betaneg_binned.dat",arr2)                                   #write out to file
    fnSave2DArray(filename+"_betapos_binned.dat",arr1)                                   #write out to file

    if Sigma!=0:
        print('Smoothing...')
        arr1=blur_image(arr1, Sigma)
        arr2=blur_image(arr2, Sigma)
        avgarr1=blur_image(avgarr1, Sigma)
        avgarr2=blur_image(avgarr2, Sigma)


    print('Writing axis files...')

    if NewTiltOrigin!=0 or NewOrientOrigin!=0:
        filename=filename+'_NewTiltOrigin'+str(NewTiltOrigin)
        filename=filename+'_NewOrientOrigin'+str(NewOrientOrigin)


    if Polar==False:                                                            #write axises
        file=open(filename+".ax1","w")
        for i in range(start1,stop1+step1,step1):
            file.write(str(float(i)-0.5*float(step1))+"\n")
        file.close()
        file=open(filename+".ax2","w")
        for i in range(start2,stop2+step2,step2):
            file.write(str(float(i)-0.5*float(step1))+"\n")
        file.close()
    else:
        file=open(filename+".axp","w")
        for i in range(-1*(stop1-step1),stop1+step1,step1):
            file.write(str(float(i)-0.5*float(step1))+"\n")                     #Igor puts labels at edges of bins
        file.close()
        file=open(filename+".axc","w")
        for i in range(-1*(stop1-step1),stop1,step1):
            file.write(str(float(i))+"\n")                                      #but not for contour plot
        file.close()


    if not Avg:
        print('Calculating Probability Plot...')

                                                                                    #create probability array
        total=0.                                                                    #statistics on array, total sum
        maxval=0.                                                                   #maximum value
        for i in range(len(arr1)):
            for j in range(len(arr1[i])):
                total=total+arr1[i][j]+arr2[i][j]
                if arr1[i][j]>maxval:
                    maxval=arr1[i][j]
                if arr2[i][j]>maxval:
                    maxval=arr2[i][j]
        probarr1=[]; probarr2=[]                                                    #initialize probability array with zeros
        for i in range(len(arr1)):
            probarr1.append([])
            probarr2.append([])
            for j in range(len(arr1[i])):
                probarr1[i].append(1.)
                probarr2[i].append(1.)

        increment=0                                                                 #create probability array
        iSpacing=(int(maxval)+1)*2
        fSpacing=1/float(iSpacing)
        for dec in range(0,iSpacing):
            for i in range(len(arr1)):                                               #how often does dec value appear?
                for j in range(len(arr1[i])):
                    if arr1[i][j]>float(dec)*fSpacing*maxval and arr1[i][j]<=float(dec+1)*fSpacing*maxval:
                        increment=increment+arr1[i][j]                               #increase sum increment
                    if arr2[i][j]>float(dec)*fSpacing*maxval and arr2[i][j]<=float(dec+1)*fSpacing*maxval:
                        increment=increment+arr2[i][j]                               #increase sum increment

            for i in range(len(arr1)):                                               #write out into probability array
                for j in range(len(arr1[i])):
                    if arr1[i][j]>float(dec)*fSpacing*maxval and arr1[i][j]<=float(dec+1)*fSpacing*maxval:
                        probarr1[i][j]=1-float(increment)/total
                    if arr2[i][j]>float(dec)*fSpacing*maxval and arr2[i][j]<=float(dec+1)*fSpacing*maxval:
                        probarr2[i][j]=1-float(increment)/total

        print(total,maxval,increment)



        print('Saving Probability Plot...')

        fnSave2DArray(filename+"_betapos_prob.dat",probarr1)
        fnSave2DArray(filename+"_betaneg_prob.dat",probarr2)

    else:
        fnSave2DArray(filename+"_betapos_"+fileext+".dat",avgarr1)
        fnSave2DArray(filename+"_betaneg_"+fileext+".dat",avgarr2)




#Loads columns from file, rescales them by a factor fScaling, and saves them into a result file
def fnConvolute(fSigma=2.5):

    import scipy.signal

    exceptionloading=[]
    exceptionprocessing=["z"]

    data={}
    for i in range(0,1855,5):
        for j in range(0,360,5):

            filename="protein_tilt"
            filename+=str(i)
            filename+="_orie"
            filename+=str(j)
            filename+=".txt"
            data=fnLoadSingleColumns(filename,data,exceptionloading)

            #print data.keys()
            print filename

            lGaussian=[]


            for k in range(len(data[data.keys()[0]])):
                d=float(k)*0.5
                lGaussian.append(0.5/(fSigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-1*(d-float(len(data[data.keys()[0]])*0.5*0.5))*(d-float(len(data[data.keys()[0]])*0.5*0.5))/2/fSigma/fSigma))
            d=0
            for dd in lGaussian:
                d+=dd


            for element in data:
                    if element not in exceptionprocessing:
                        data[element]=scipy.signal.convolve(data[element],lGaussian,"same")

            filename= 'protein_tilt'
            filename+=str(i)
            filename+="_orie"
            filename+=str(j)
            filename+=".txt"

            fnSaveSingleColumns(filename, data)

def fnConvoluteSingle(filename='',fSigma=2.5):

    import scipy.signal

    exceptionloading=[]
    exceptionprocessing=["z"]

    data={}
    data=fnLoadSingleColumns(filename,data,exceptionloading)

    #print data.keys()
    print filename

    lGaussian=[]


    for k in range(len(data[data.keys()[0]])):
        d=float(k)*0.5
        lGaussian.append(0.5/(fSigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-1*(d-float(len(data[data.keys()[0]])*0.5*0.5))*(d-float(len(data[data.keys()[0]])*0.5*0.5))/2/fSigma/fSigma))
    d=0
    for dd in lGaussian:
        d+=dd


    for element in data:
            if element not in exceptionprocessing:
                data[element]=scipy.signal.convolve(data[element],lGaussian,"same")

    filename+="._convoluted.txt"

    fnSaveSingleColumns(filename, data)

# create log likelihood plot
def fnCreatePenetrationPlot(NewTiltOrigin=0.0, NewOrientOrigin=0.0):
    fnCreateLogLPlot(NewTiltOrigin=NewTiltOrigin, NewOrientOrigin=NewOrientOrigin, filename='sErr.dat', parameter1='beta', parameter2='gamma', parameter3='penetration', start1=0, stop1=93, step1=3, start2=0, stop2=360, step2=5, Polar=True,Sigma=0,Deg=True,fileext='penetration')


def fnCreateLogLPlot(NewTiltOrigin=0.0, NewOrientOrigin=0.0, filename='sErr.dat', parameter1='beta', parameter2='gamma', parameter3='Chisq', start1=0, stop1=93, step1=3, start2=0, stop2=360, step2=5, Polar=True,Sigma=0,Deg=True, fileext='LogL'):

    data={}

    print('Load file...')
    fnLoadSingleColumns(filename, data)

    print('Analyze Data...')
    fnBin2D(filename,data[parameter1],data[parameter2],start1,stop1,step1,start2,stop2,step2,Polar,Sigma,Deg,data[parameter3],True, fileext, NewTiltOrigin=NewTiltOrigin, NewOrientOrigin=NewOrientOrigin)


def fnCreateProbabilityPlot(NewTiltOrigin=0.0, NewOrientOrigin=0.0, filename='sErr.dat', parameter1='beta', parameter2='gamma', start1=0, stop1=93, step1=3, start2=0, stop2=360, step2=5, Polar=True,Sigma=0,Deg=True):

    data={}

    print('Load file...')
    fnLoadSingleColumns(filename, data)

    print('Analyze Data...')
    fnBin2D(filename,data[parameter1],data[parameter2],start1,stop1,step1,start2,stop2,step2,Polar,Sigma,Deg, fileext='prob',NewTiltOrigin=NewTiltOrigin,NewOrientOrigin=NewOrientOrigin)

#uses the nSLD envelopes and the area envelopes to calculate the fraction of deuterated
#material and integrates it (running integral)
def fnnSLDFracIntEnvelope(nSLDmin=1.8e-6, nSLDmax=6.4e-6, normalization=1):

    random.seed()
    envelope_area={}
    envelope_nSLD={}
    envelope_frac1={}
    envelope_frac2={}
    exceptionloading=[]
    envelope_area=fnLoadSingleColumns("pulledmolgroups_area.dat",envelope_area, exceptionloading)
    envelope_frac1=copy.deepcopy(envelope_area)
    envelope_frac2=copy.deepcopy(envelope_area)
    envelope_nSLD=fnLoadSingleColumns("pulledmolgroups_nSLD.dat",envelope_nSLD, exceptionloading)
    intenvelope={}
    referenceaxis="zaxis"

    for key in envelope_area.keys():
        if key != referenceaxis:
            for i in range(len(envelope_area[key])):
                envelope_frac2[key][i]=envelope_area[key][i]*(envelope_nSLD[key][i]-nSLDmin)/(nSLDmax-nSLDmin)
                envelope_frac1[key][i]=envelope_area[key][i]-envelope_frac2[key][i]

    print 'hallo'

    for key in envelope_frac1.keys():
        if key != referenceaxis:
            data4=fnIntegrate(envelope_frac1,key,referenceaxis)
            fnNormalizeLast(data4,key,normalization)
            fnAppend(data4,intenvelope,key,referenceaxis)

    fnSaveSingleColumns('envelopeintfrac1.dat',intenvelope)
    fnSaveSingleColumns('envelopeintfracstat1.dat',fnStat(intenvelope,referenceaxis))

    intenvelope={}

    for key in envelope_frac2.keys():
        if key != referenceaxis:
            data4=fnIntegrate(envelope_frac2,key,referenceaxis)
            print key
            fnNormalizeLast(data4,key,normalization)
            fnAppend(data4,intenvelope,key,referenceaxis)

    fnSaveSingleColumns('envelopeintfrac2.dat',intenvelope)
    fnSaveSingleColumns('envelopeintfracstat2.dat',fnStat(intenvelope,referenceaxis))


#uses the nSLD envelopes and the area envelopes to calculate the fraction of deuterated
#material
def fnnSLDFracEnvelope(nSLDmin=1.8e-6, nSLDmax=6.4e-6):

    random.seed()
    envelope_area={}
    envelope_nSLD={}
    envelope_frac1={}
    envelope_frac2={}
    exceptionloading=[]
    envelope_area=fnLoadSingleColumns("pulledmolgroups_area.dat",envelope_area, exceptionloading)
    envelope_frac1=copy.deepcopy(envelope_area)
    envelope_frac2=copy.deepcopy(envelope_area)
    envelope_nSLD=fnLoadSingleColumns("pulledmolgroups_nSLD.dat",envelope_nSLD, exceptionloading)
    intenvelope={}
    referenceaxis="zaxis"

    for key in envelope_area.keys():
        if key != referenceaxis:
            for i in range(len(envelope_area[key])):
                envelope_frac2[key][i]=envelope_area[key][i]*(envelope_nSLD[key][i]-nSLDmin)/(nSLDmax-nSLDmin)
                envelope_frac1[key][i]=envelope_area[key][i]-envelope_frac2[key][i]


    fnSaveSingleColumns('envelopefrac1.dat',envelope_frac1)
    fnSaveSingleColumns('envelopefracstat1.dat',fnStat(envelope_frac1,referenceaxis))
    fnSaveSingleColumns('envelopefrac2.dat',envelope_frac2)
    fnSaveSingleColumns('envelopefracstat2.dat',fnStat(envelope_frac2,referenceaxis))



#calculates the continuous integral of an envelope
def fnIntEnvelope():

    random.seed()
    envelope={}
    exceptionloading=[]
    envelope=fnLoadSingleColumns("pulledmolgroups_area.dat",envelope, exceptionloading)
    intenvelope={}
    referenceaxis="zaxis"

    for key in envelope.keys():
        if key != referenceaxis:
            data4=fnIntegrate(envelope,key,referenceaxis)
            fnNormalizeLast(data4,key,1)
            fnAppend(data4,intenvelope,key,referenceaxis)

    fnSaveSingleColumns('envelopeint.dat',intenvelope)
    fnSaveSingleColumns('envelopeintstat.dat',fnStat(intenvelope,referenceaxis))

#calculates the difference and stats for a pair of sets of envelopes
def fnDiffEnvelope():

    random.seed()
    envelope1={}
    envelope2={}
    exceptionloading=[]
    envelope1=fnLoadSingleColumns("envelope1.dat",envelope1, exceptionloading)
    envelope2=fnLoadSingleColumns("envelope2.dat",envelope2, exceptionloading)
    diffenvelope={}
    intenvelope={}
    referenceaxis="zaxis"
    columnname='envelope'

    for _ in range(5000):

            randomcolumn1=referenceaxis
            randomcolumn2=referenceaxis

            while randomcolumn1==referenceaxis:
                j=int(random.random()*len(envelope1[referenceaxis]))
                randomcolumn1=envelope1.keys()[j]

            while randomcolumn2==referenceaxis:
                k=int(random.random()*len(envelope2[referenceaxis]))
                randomcolumn2=envelope2.keys()[k]

            data1={referenceaxis: envelope1[referenceaxis][:], columnname : envelope1[randomcolumn1][:] }
            data2={referenceaxis: envelope2[referenceaxis][:], columnname : envelope2[randomcolumn2][:] }

            fnCutZeros(data1,referenceaxis)
            fnCutZeros(data2,referenceaxis)

            fnUniformLength([data1,data2])

            fnNormalizeIntegral(data1,columnname,referenceaxis,85045.0)
            fnNormalizeIntegral(data2,columnname,referenceaxis,69636.0)

            data3=fnSubtract(data1,data2,columnname,referenceaxis)
            fnAppend(data3,diffenvelope,columnname,referenceaxis)

            data4=fnIntegrate(data3,columnname, referenceaxis)
            fnNormalizeLast(data4,columnname,140)
            fnAppend(data4,intenvelope,columnname,referenceaxis)

    fnSaveSingleColumns('envelopedif.dat',diffenvelope)
    fnSaveSingleColumns('envelopeint.dat',intenvelope)
    fnSaveSingleColumns('envelopedifstat.dat',fnStat(diffenvelope,referenceaxis))
    fnSaveSingleColumns('envelopeintstat.dat',fnStat(intenvelope,referenceaxis))


#reads a number of files depending on two variables xx and yy and creating
#a c-style output array containing the data of those files
#the array will have the format [xx][yy][columns] with the reference column
#being the first one
#xxrange, and yyrange have the format [start, end, stepsize]

def fnCreateCArray(filename1, filename2, filename3, xxrange, yyrange, datacolumns=[]):

    file=open("result.txt","w");

    file.write("double arr ["+str((xxrange[1]-xxrange[0])/xxrange[2])+"] ")
    file.write("["+str((yyrange[1]-yyrange[0])/yyrange[2])+"] ")
    file.write("["+str(len(datacolumns)+1)+"] {")
    for xx in range(xxrange[0],xxrange[1],xxrange[2]):
        file.write("{")
        for yy in range(yyrange[0], yyrange[1], yyrange[2]):
            file.write("{")

            if xx==0:                                                           #custom adaption to files
                yyp=0
            else:
                yyp=yy

            data={}
            data=fnLoadSingleColumns(filename1+str(xx)+filename2+str(yyp)+filename3, data, [], False, datacolumns)
            for columnname in datacolumns:
                column=data[columnname]
                file.write("{")

                i=0
                for element in column:
                    file.write(str(element))
                    i=i+1
                    if element!=column[len(column)-1]:
                        file.write(", ")
                        if i==10:
                            file.write("\n")
                            i=0
                file.write("}")
                if columnname!=datacolumns[len(datacolumns)-1]:
                    file.write(",\n")

            file.write("}")
            if yy!=(yyrange[1]-yyrange[2]):
                file.write(",\n")

        file.write("}")
        if xx!=(xxrange[1]-xxrange[2]):
            file.write(",\n")

    file.write("};\n")
    file.close()


#reads a referencefile of the format like N. Kucerka's gamma3b .cmp/.cmp file
#and adds the there grouped data columns and saves the groups under their group
#names in an output file
#originally designed for conveniently adding up sub-molecular groups from indi-
#vidual atoms
def fnGroupByFile(referenceaxis='z'):

    reference={}
    exceptionloading=["#", "##"]
    reference=fnLoadSingleRows("NandaDMPC.cmd",reference,exceptionloading)

    data={}
    data=fnLoadSingleColumns("dmpc_liq_135ns_byatom_wrapped.txt",data, exceptionloading)

    results={}

    results[referenceaxis]=data[referenceaxis][:]

    for category in reference:
        results[category]=[]
        for dataitem in reference[category]:
            try:
                if results[category]==[]:
                    results[category]=data[dataitem][:]
                else:
                    for i in range(len(data[dataitem])):
                        results[category][i]=results[category][i]+data[dataitem][i]
            except:
                print("Could not interprete "+category+" "+dataitem)

    fnSaveSingleColumns("dmpc_liq_135ns_byatom_grouped.txt", results)

def fnPeakAdjust(iNumberOfPeaks=1):

    def fnPA(value1, value2, iPeaked, iMaxPeak, iLastTrendUp):

        i2Return=0
        tolerance=0.98;

        if abs(value2)>=abs(value1)*tolerance and abs(value2)<=abs(value1)/tolerance:
            i2Return=iLastTrendUp
        elif abs(value1)<abs(value2):
            i2Return=1
        else:
            i2Return=0

        if iPeaked<iMaxPeak:
            if iLastTrendUp==1 and abs(value2)<abs(value1)*tolerance:
                iPeaked+=1
        else:
            if abs(value2)>abs(value1)/tolerance:
                value2=abs(value1)/tolerance

        return i2Return,iPeaked,abs(value2)


    exceptionloading=[]

    data={}
    data=fnLoadSingleColumns("sErr.dat",data,exceptionloading)

    keys=data.keys()
    print keys
    for i in range(len(data[keys[0]])):

        j=1
        while 1:
            if 'vf_on'+str(j) in keys:
                print '%(st)s ' % {'st':data['vf_on'+str(j)][i]}
                j+=1
            else:
                print('TO ')
                break

        j=1
        iPeaked=0; iMaxPeak=1; i2=0; vf_old=0;
        while 1:
            if 'vf_on'+str(j) in keys:
                i2,iPeaked,dCorrected=fnPA(vf_old,data['vf_on'+str(j)][i],iPeaked,iMaxPeak,i2);
                data['vf_on'+str(j)][i]=dCorrected
                vf_old=dCorrected
                j+=1
            else:
                break

        j=1
        while 1:
            if 'vf_on'+str(j) in keys:
                print '%(st)s ' % {'st':data['vf_on'+str(j)][i]}
                j+=1
            else:
                print('\n')
                break

    fnSaveSingleColumns('sErr.dat',data)


def fnRescale(fScaling=5):

    exceptionloading=[]
    exceptionscaling=["z","tot","tot_2nd"]

    data={}
    data=fnLoadSingleColumns("gammaoutput.vol",data, exceptionloading)

    for element in data:
        if element not in exceptionscaling:
            for i in range(len(data[element])):
                data[element][i]=data[element][i]*fScaling

    fnSaveSingleColumns("voloutput.txt", data)




#wraps data set around limit disregarding the fact if the data lies within limits or not
#this happens by simply adding leading to an increased amount of the wrapped quantity
#data must have sufficient zero in the center to avoid overlapping
#z-coordinate is expected at the zero column
#stepsize is determined from data
def fnWrap(limit=20.905,refaxis='z'):

    exceptionloading=[]
    wrapping=["H1","H2","OH2"]

    data={}
    data=fnLoadSingleColumns("dmpc_liq_135ns_byatom_scaled.txt",data, exceptionloading)

    datalength=len(data[refaxis])
    datastart=data[refaxis][0]
    datastep=data[refaxis][1]-data[refaxis][0]
    lowerlimit=int((limit*(-1)-datastart)/datastep+0.5)
    upperlimit=int((limit-datastart)/datastep+0.5)

    for element in wrapping:
        boxlength=int(2*limit/datastep+0.5)

        temp=[]
        for i in range(len(data[element])):
            temp.append(0.0)

        for i in range(len(data[element])):
            if data[refaxis][i]<0:
                pos=i-lowerlimit+upperlimit
            else:
                pos=i-upperlimit+lowerlimit
            if (pos>0) and (pos<datalength):
                temp[pos]=temp[pos]+data[element][i]
            temp[i]=temp[i]+data[element][i]

        for i in range(236):                                                    #custom cut-off
            temp[i]=0.
        for i in range(509,len(data[element])):
            temp[i]=0.

        data[element]=temp[:]


    fnSaveSingleColumns("dmpc_liq_135ns_byatom_wrapped.txt", data)




