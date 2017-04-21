#!/usr/local/bin/python2.7

# Import necessary modules
import csv
import sys as sys
import os
import glob as glob
from datetime import datetime


def merge_csv():
	'''Function merges all csv''s in directory.

	Directory of csv's and name of exported csv
	taken from command line.
	'''
	# Change directory to path containing .csv files to merge
	os.chdir(sys.argv[1])

	# Current date as local variable
	datestring = datetime.strftime(datetime.now(), '%Y%m%d')

	# Append all .csv files to python list
	filelist = []
	for csv in glob.glob("*.csv"):
		filelist.append(csv)

	# Create file that will be the merged .csv file.
	merged = open(datestring + "_" + sys.argv[2],"a")

	# For first .csv file in filelist, read all lines
	for line in open(filelist[0]):
		merged.write(line)

	# Read all other .csv files, skipping the header.
	for csv in filelist[1:len(filelist)]:
		f = open(csv)
		f.next() # skip the header
		for line in f:
			merged.write(line)
		f.close()
	merged.close()

if __name__ == '__main__':
	merge_csv() # Call function if script is run 'main'