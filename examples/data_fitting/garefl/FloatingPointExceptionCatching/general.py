#!/usr/bin/env python
# -*- coding: utf-8 -*-
# general purpose for 16-Sep-09 F.H.

import math, os, string
import numpy
import scipy.special, scipy.optimize, scipy.stats
import matplotlib.pyplot as plt


#-------------------------------------------------------------------------------
#General usefull routines

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


#Loads columns from file, rescales them by a factor fScaling, and saves them into a result file
def fnConvolute(fSigma=2.5):

    import scipy.signal

    exceptionloading=[]
    exceptionprocessing=["z"]

    data={}
    for i in range(-90,95,5):
        for j in range(0,360,5):

            filename="chain1_tilt"
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

            filename= 'chain1conv_tilt'
            filename+=str(i)
            filename+="_orie"
            filename+=str(j)
            filename+=".txt"

            fnSaveSingleColumns(filename, data)




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

        if value2>=value1*tolerance and value2<=value1/tolerance:
            i2Return=iLastTrendUp
        elif value1<value2:
            i2Return=1
        else:
            i2Return=0

        if iPeaked<iMaxPeak:
            if iLastTrendUp==1 and value2<value1*tolerance:
                iPeaked+=1
        else:
            if value2>value1/tolerance:
                value2=value1/tolerance

        return i2Return,iPeaked,value2


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




