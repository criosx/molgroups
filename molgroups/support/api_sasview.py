from __future__ import print_function
from math import fabs, pow, sqrt
from os import path
from random import seed, normalvariate, random
from re import VERBOSE, IGNORECASE, compile
from sys import exit, stdout
from subprocess import call, Popen
from time import sleep
import numpy
import pandas
import shutil
import glob
import os

from molgroups.support import api_bumps


class CSASViewAPI(api_bumps.CBumpsAPI):
    def __init__(self, spath='.', mcmcpath='.', runfile='', load_state=True):
        super().__init__(spath, mcmcpath, runfile, load_state=load_state)
