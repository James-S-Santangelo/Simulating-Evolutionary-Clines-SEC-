#!/usr/local/bin/python2.7

# Import necessary modules
import sys as sys
import os
import glob as glob
from datetime import datetime
import shutil


def merge_csv():
    '''Function merges csvs in directory.

    Merges .csv files with a user specified parameter value
    in it's name. Useful for merging only csv's with data
    corresponding to specific simulations.

    Returns:
    1. merged csv
    '''
    # Specify the parameter value for which datasets should be merged
    # param = sys.argv[1]
    # value = sys.argv[2]
    # parameter = param + value

    # Change directory to path containing .csv files to merge
    os.chdir(sys.argv[1])

    # Current date as local variable
    datestring = datetime.strftime(datetime.now(), '%Y%m%d')

    # Append all .csv files to python list
    filelist = []
    for file in glob.glob("*.csv"):
        # if parameter in file:
        filelist.append(file)

    # Create file that will be the merged .csv file.
    merged = open(datestring + "_" + sys.argv[2] + ".csv", "a")

    # For first .csv file in filelist, read all lines
    for line in open(filelist[0]):
        merged.write(line)

    # Read all other .csv files, skipping the header.
    for file in filelist[1:len(filelist)]:
        f = open(file)
        f.next()  # skip the header
        for line in f:
            merged.write(line)
        f.close()
    merged.close()

    for file in filelist:
        shutil.move(file, "/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/To_delete")


if __name__ == '__main__':
    merge_csv()  # Call function if script is run 'main'
