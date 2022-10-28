from __future__ import print_function
from os import path
from random import seed, normalvariate, random
from sys import exit, stdout

import os
import pandas
import pathlib

from molgroups.support import general


class CBaseAPI:
    def __init__(self, spath='.', mcmcpath='.', runfile=''):
        # script path, MCMC path, runfile (script file)
        self.spath = spath
        self.mcmcpath = mcmcpath
        self.runfile = runfile
        self.diParameters = {}

    def fnBackup(self):
        raise NotImplementedError()

    def fnBackupSimdat(self):
        raise NotImplementedError()

    @staticmethod
    def fnExtendQRange(liData, qmin=None, qmax=None, conserve_dq_q=False):
        for i in range(len(liData)):
            if qmax is not None:
                # first cut data at qrange
                liData[i][1] = liData[i][1][(liData[i][1].Q <= qmax)]
                # now add data points in case the q-range is too short
                while liData[i][1]['Q'].iloc[-1] < qmax:
                    newframe = pandas.DataFrame(liData[i][1][-1:], columns=liData[i][1].columns, copy=True)
                    step = liData[i][1]['Q'].iloc[-1] - liData[i][1]['Q'].iloc[-2]
                    step *= 1.05
                    newframe['Q'].iloc[-1] = liData[i][1]['Q'].iloc[-1] + step
                    if conserve_dq_q:
                        newframe['dQ'].iloc[-1] = liData[i][1]['dQ'].iloc[-1] / liData[i][1]['Q'].iloc[-1] * \
                                                  newframe['Q'].iloc[-1]
                    if newframe['Q'].iloc[-1] < qmax:
                        liData[i][1] = liData[i][1].append(newframe, ignore_index=True)
                    else:
                        break
            if qmin is not None:
                liData[i][1] = liData[i][1][(liData[i][1].Q >= qmin)]
                while liData[i][1]['Q'].iloc[0] > qmin:
                    newframe = pandas.DataFrame(liData[i][1][:1], columns=liData[i][1].columns, copy=True)
                    newframe['Q'].iloc[0] = 2 * liData[i][1]['Q'].iloc[0] - liData[i][1]['Q'].iloc[1]
                    if conserve_dq_q:
                        newframe['dQ'].iloc[0] = liData[i][1]['dQ'].iloc[0] / liData[i][1]['Q'].iloc[0] * \
                                                 newframe['Q'].iloc[0]
                    if newframe['Q'].iloc[0] > qmin:
                        liData[i][1] = newframe.append(liData[i][1], ignore_index=True)
                    else:
                        break
        return liData

    def fnLoadMolgroups(self, problem=None):
        diMolgroups = {}
        moldict = problem.moldat

        for group in moldict:
            tdata = (moldict[group]['header']).split()  # read header that contains molgroup data
            diMolgroups[tdata[1]] = {}
            diMolgroups[tdata[1]].update({'headerdata': {}})
            diMolgroups[tdata[1]]['headerdata'].update({'Type': tdata[0]})
            diMolgroups[tdata[1]]['headerdata'].update({'ID': tdata[1]})
            for j in range(2, len(tdata), 2):
                diMolgroups[tdata[1]]['headerdata'].update({tdata[j]: tdata[j + 1]})

            zax = moldict[group]['zaxis']
            areaax = moldict[group]['area']
            nslax = moldict[group]['nsl']
            diMolgroups[tdata[1]].update({'zaxis': zax, 'areaaxis': areaax, 'nslaxis': nslax})

        return diMolgroups

    def fnLoadsErr(self):
        sFileName = self.mcmcpath + '/isErr.dat'
        try:
            if path.isfile(self.mcmcpath+'/isErr.dat') and path.isfile(self.mcmcpath+'/sErr.dat'):
                print('-------------------------------')
                print('Found isErr.dat and sErr.dat ?!')
                print('Load isErr.dat as default.')
                print('-------------------------------')
            elif path.isfile(self.mcmcpath+'/isErr.dat'):  # check which type of MC output present
                print('Found isErr.dat\n')
            elif path.isfile(self.mcmcpath+'/sErr.dat'):
                print('Found sErr.dat\n')
                sFileName = self.mcmcpath+'/sErr.dat'

            diStatRawData = {'Parameters': self.fnLoadSingleColumnsIntoStatDict(sFileName)}
            return diStatRawData

        except IOError:
            print('Could not load ' + sFileName + '. \n')
            exit(1)

    # The sparse parameter loads only a fraction of the data, if sparse is larger or equal than 1 this equates to the
    # number of lines loaded. If sparse is smaller than one this equates to a probability that any line will be stored
    @staticmethod
    def fnLoadSingleColumns(sFileName, data=None, exceptions=None, header=True, headerline=None, LoadList=None,
                            sparse=0):

        File = open(sFileName, "r")
        content = File.readlines()
        File.close()

        if data is None:
            data = {}

        if header:  # if headerline in file, read it from there
            splitheaderline = content[0].split()
            content = content[1:]
        else:  # else use the one given as an attribute
            splitheaderline = headerline

        for i, columnname in enumerate(splitheaderline):
            if LoadList is None or (columnname in LoadList):
                data[columnname] = []

        j = 0
        seed()
        for line in content:
            if sparse == 0 or (sparse >= 1 and j < sparse) or (1 > sparse > random()):
                splitline = line.split()
                if splitline:
                    for i, column in enumerate(splitline):
                        if LoadList is None or (splitheaderline[i] in LoadList):
                            try:
                                data[splitheaderline[i]].append(float(splitline[i]))
                            except:
                                data[splitheaderline[i]].append(splitline[i])
                j += 1

        return data

    # Loads single row data into a dictionary ["rowname",'Values',[data]]
    # it appends the data found in the file to the provided list 'data'
    # it will skip rows with names which are either in the exception list
    # or appends data with name extension "_2nd" that are already present in the data list
    # The sparse parameter loads only a fraction of the data, if sparse is larger or equal than 1 this equates to the
    # number of lines loaded. If sparse is smaller than one this equates to a probability that any line will be stored
    @staticmethod
    def fnLoadSingleColumnsIntoStatDict(sFileName, data=None, exceptions=None, header=True, headerline=None,
                                        LoadList=None, sparse=0):

        if not data: data = {}
        file = open(sFileName, "r")
        content = file.readlines()
        file.close()

        if data is None:
            data = {}

        if header:  # if headerline in file, read it from there
            splitheaderline = content[0].split()
            content = content[1:]
        else:  # else use the one given as an attribute
            splitheaderline = headerline

        for i, columnname in enumerate(splitheaderline):
            if LoadList is None or (columnname in LoadList):
                data[columnname] = {}
                data[columnname]['Values'] = []

        j = 0
        seed()
        for line in content:
            if sparse == 0 or (sparse >= 1 and j < sparse) or (sparse < 1 and random() < sparse):
                splitline = line.split()
                if splitline:
                    for i, column in enumerate(splitline):
                        if LoadList is None or (splitheaderline[i] in LoadList):
                            data[splitheaderline[i]]['Values'].append(float(splitline[i]))
                j += 1

        return data

    @staticmethod
    def fnLoadSingleRows(sFileName, data=None, exceptions=None):
        file = open(sFileName, "r")
        content = file.readlines()
        file.close()

        if data is None:
            data = {}

        for line in content:
            splitline = line.split()
            if splitline[0] not in exceptions:
                data[splitline[0]] = []
                for entry in range(1, len(splitline)):
                    data[splitline[0]].append(splitline[entry])
        return data

    @staticmethod
    def fnMCModifyFile(filelist):
        """
        reads original reflectivity files and modifies them with normal deviates
        works also with any other 3 or 4 column file types
        """

        seed()  # initialize random number generator

        for reflfile in filelist:  # iterate over all refl files
            file = open(reflfile)
            data = file.readlines()
            file.close()
            newdata = []
            iFileType = 0
            for line in data:
                if '#' in line:  # do not modify comment lines
                    newdata.append(line)
                else:
                    columns = line.split()  # access columns
                    if not iFileType:  # Stuart's Mod for 3 Column i.e. NIST
                        if len(columns) == 3:
                            iFileType = 1
                        elif len(columns) == 4:
                            iFileType = 2  # Stuart's Mod for 4 Column i.e. ISIS
                        else:
                            iFileType = 99  # unrecognized format
                    if not len(columns):
                        print('Found empty line in data file.')
                    if len(columns) == 3 and iFileType == 1:
                        try:
                            fvalue = float(columns[1])  # reflectivity
                            ferror = float(columns[2])  # error
                            columns[1] = str(fvalue + normalvariate(0, 1) * ferror)  # modify reflectivity
                            newline = ''
                            for column in columns:  # glue columns together to one line
                                newline = newline + column + ' '
                            newline = newline[:-1] + '\n'  # add newline to end of line
                            newdata.append(newline)  # add line to new data file
                        except:
                            print('-----------------------------------')
                            print('Data file %s corrupt.' % reflfile)
                            print('File was identified being NIST type')
                            print('-----------------------------------')
                            raise RuntimeError('Corrupt Data file.')
                    elif len(columns) == 4 and iFileType == 2:
                        try:
                            fvalue = float(columns[2])  # reflectivity
                            ferror = float(columns[3])  # error
                            columns[2] = str(fvalue + normalvariate(0, 1) * ferror)  # modify reflectivity
                            newline = ''
                            for column in columns:  # glue columns together to one line
                                newline = newline + column + ' '
                            newline = newline[:-1] + '\n'  # add newline to end of line
                            newdata.append(newline)  # add line to new data file
                        except:
                            print('-----------------------------------')
                            print('Data file %s corrupt.' % reflfile)
                            print('File was identified being ISIS type.')
                            print('-----------------------------------')
                            raise RuntimeError('Corrupt Data file.')

                    else:
                        print('-----------------------------------------------------')
                        print('Filetype not recognized or contains errors: %s' % reflfile)
                        print('-----------------------------------------------------')
                        raise RuntimeError('Corrupt Data file.')

            file = open(path.split(reflfile)[-1] + '.mce', "w")  # write modified file into .mce file in
            file.writelines(newdata)  # working directory
            file.close()

    def fnRestoreFit(self):
        raise NotImplementedError()

    def fnRemoveBackup(self):  # deletes the backup directory
        raise NotImplementedError()

    @staticmethod
    def fnSaveSingleColumns(sFilename, data):
        # saves all data out to a file

        file = open(sFilename, "w")

        for element in data:
            file.write(element + " ")
        file.write("\n")

        for i in range(len(data[list(data.keys())[0]])):
            for element in data:
                file.write(str(data[element][i]) + " ")
            file.write("\n")

        file.close()
        # saves all data out to a file

    @staticmethod
    def fnSaveSingleColumnsFromStatDict(sFilename, data, skipentries=None):

        if skipentries is None:
            skipentries = []
        File = open(sFilename, "w")

        for element in data:
            if element not in skipentries:
                File.write(element + " ")
        File.write("\n")

        for i in range(len(data[list(data.keys())[0]]['Values'])):
            for element in data:
                if element not in skipentries:
                    File.write(str(data[element]['Values'][i]) + " ")
            File.write("\n")

        File.close()

    def fnSimulateData(self, diExpression, liData):
        raise NotImplementedError()

    def fnWriteConstraint2Runfile(self, liExpression):
        raise NotImplementedError()
