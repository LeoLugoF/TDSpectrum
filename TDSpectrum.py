import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from math import exp

#########################################################################
#                              TDSpectrum                               #
# Made by:      LeoLugoF                                                #                 
# Date:         30/May/2019                                             #
# Requirements: Mathplotlib, numpy.                                     #
#                                                                       #
# Description:                                                          #
# This program reads gaussian output files with .log or .out termation. #
# It reads the oscilator forces and wave lenght for each excitation,    #
# and shows a UV-Vis spectrum, but the data can also be saved.          #
# Data is stored with name of the old file + ".txt".                    #
# How the calculations are done can be consulted at the gaussian page:  #
# https://gaussian.com/uvvisplot/                                       #
# Note: Y values are in absorbance e(L mol-1 cm-1)                      #
#                                                                       #
# Command line:[python3] [*.py] [*.log/*.out] [Sigma] [MinWave] ...     #
#              [MaxWave] [NONE/-s]                                      #
# Arguments:                                                            #
#     [Sigma]   = Value of sigma; Gaussian recommends 0.4 (eV).         #
#     [MinWave] = Minimum Wavelenght (nm)                               #
#     [MaxWave] = Maximum Wavelenght (nm)                               #
#     [-s]      = Saves only the data; doesn't show the graph           #
#     [NONE]    = Shows the graph only without saving data              #
#                                                                       #
# Examples:                                                             #
# Example 1:    python TDSpectrum.py YourFile.log 0.4 300 500           #
# Example 2:    python TDSpectrum.py YourFile.out 0.2 200 600 -s        #
#                                                                       #
# The first example will show only the UV-Vis plot.                     #
# The second example will save only the data without showing the plot.  #
#########################################################################

class Global:
    """Global variables; Stores information."""
    WaveLenghts = np.array([])
    Forces = np.array([])
    XValues = np.array([])
    YValues = np.array([])
    ShowPlot = True

def ReadFile(FilePath):
    """Reads the output file and stores the information."""
    fstream = open(FilePath,"r")
    lines = fstream.readlines()
    fstream.close()
    for line in lines:
        if "Excited State" in line and "<S**2>=" in line:
            i = 0
            for subsentence in line.split(" "):
                if(len(subsentence) > 1):
                    if i == 6:
                        # This element always corresponds to the Wavelenght (nm)
                        Global.WaveLenghts = np.append(Global.WaveLenghts,float(subsentence))
                        i += 1
                    elif i == 8:
                        # This element always corresponds to the oscilator force (F)
                        Global.Forces = np.append(Global.Forces,float(subsentence.split("=")[1]))
                        break
                    else:
                        i += 1                  
    return

def DoCalcs(Sigma,MinWave,MaxWave):
    """Calculates the Y values from the MinWave and MaxWave giving with the sigma value."""
    CMm1 = Sigma*(10**7)*0.000806556
    NMm1 = 0.000806556*Sigma
    Global.XValues = np.arange(MinWave,MaxWave,1)
    Global.YValues = np.zeros(len(Global.XValues))
    Matrix = np.zeros((len(Global.XValues),len(Global.WaveLenghts)))
    #Row number
    i = 0
    for row in Matrix:
        #Column Number
        j = 0
        for cell in row:
            Constant = 130629740*(Global.Forces[j]/CMm1)
            Matrix[i,j] = Constant*exp(-((((1/Global.XValues[i])-(1/Global.WaveLenghts[j]))/NMm1)**2))      
            j += 1            
        i += 1
    #Sum columns
    i = 0
    for Row in Matrix:
        Summatory = 0
        for Cell in Row:
            Summatory += Cell
        Global.YValues[i] = Summatory
        i += 1
    return

def ShowGraph():
    """Shows the plot,"""
    fig, ax = plt.subplots()
    ax.plot(Global.XValues,Global.YValues)
    plt.xlabel("λ(nm)")
    plt.ylabel("e(L mol-1 cm-1)")
    plt.show()

def SaveFile():
    """Stores the x and y data into a text file."""
    SaveFilePath =  FilePath.split(".")[0] + ".txt"
    f = open(SaveFilePath,"a")
    i = 0
    for XValue in Global.XValues:
        f.write(str(XValue) + "\t" + str(Global.YValues[i]) + "\n")
        i += 1
    f.close()


FilePath = ""
i = 0
#Reads the extra comment arguments giving
for arg in sys.argv:
    if ".out" in arg or ".log" in arg or ".OUT" in arg or ".LOG" in arg:
        FilePath = os.getcwd() + "\\" + arg
    elif "-s" in arg:
        Global.ShowPlot = False
    else:
        try:
            Number = float(arg)
            if i == 0:
                Sigma = Number
            if i == 1:
                MinWave = Number
            if i == 2:
                MaxWave = Number
            i += 1
        except:
            pass
        
#If no comment arguments are giving it will ask for it.
if FilePath == "":
    FilePath = input("Please Insert the file path: ")
    ReadFile(FilePath)
    Sigma = input("Sigma Value: ")
    MinWave = input("Min WaveLenght (nm): ")
    MaxWave = input("Max WaveLenght (nm): ")

ReadFile(FilePath)
if(len(Global.WaveLenghts) == 0):
    print("No excited states found.")
else:
    DoCalcs(float(Sigma),float(MinWave),float(MaxWave))
    if Global.ShowPlot is True:
        ShowGraph()
    else:
        SaveFile()
