#!/usr/bin/env python

import argparse, os
import itertools as it
import numpy as np
import fastatools as ft	#Available at https://github.com/jtladner/Modules
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict

#Example command:

parser = argparse.ArgumentParser(description='''...''')

#parser.add_argument("-o", "--outputName", default="XXX", metavar='\b', help="... [XXX]")
#parser.add_argument("-v", "--visualization", default="XXX", metavar='\b', help="... [XXX]")
parser.add_argument("-w", "--windowSize", type=int, metavar='\b', help="Window size. If not provided, amplicon size is used as default.")
parser.add_argument("-s", "--stepSize", default=400, type=int, metavar='\b', help="Step size (i.e. number of nucleotides to move between each window). [XXX]")
parser.add_argument("-f", "--alignment", metavar='\b', help="Alignment file")
parser.add_argument("-n", "--referenceName", metavar='\b', help="Name of reference sequence, source of provided **unaligned** coords")
parser.add_argument("-a", "--start", type=int, metavar='\b', help="Start coord of amplicon of interest in **unaligned** reference seq")
parser.add_argument("-b", "--stop", type=int, metavar='\b', help="Stop coord of amplicon of interest in **unaligned** reference seq")
parser.add_argument("-m", "--batchMode", default=None, metavar='\b', help="You can use this flag to run the script in batch mode. If used, it should be followed by the path to a tsv file with four columns and one row per alignment. Columns, respectively, should correspond to --alignment, --referenceName, --start, --stop. In this mode, the output filenames will be generated based on the input file names. [default: None]")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-t", "--divThreshold", required=True, type=int, metavar='\b', help="Percent divergent threshold")

args = parser.parse_args()


def seqWindows(s1, s2, windowSize, stepSize):
	#Remove dashes if present in same position of both sequences
	ns1=""
	ns2=""
	for i in range(len(s1)):
		if not (s1[i]== "-" and s2[i]== "-"):
			ns1+= s1[i]
			ns2+= s2[i]
	#Create a list of windows for each sequence
	ns1windowsL=[]
	ns2windowsL=[]
	for i in range(0, len(ns1)-windowSize+1, stepSize):
		ns1window= ns1[i:i+windowSize]
		ns2window= ns2[i:i+windowSize]
		ns1windowsL.append(ns1window)
		ns2windowsL.append(ns2window)
	if i<len(ns1)-1:
		ns1windowsL.append(ns1[-windowSize:])
		ns2windowsL.append(ns2[-windowSize:])

	return ns1windowsL, ns2windowsL

def calcDivergence(s1, s2):
	diff=0
	tot=0
	for i, n1 in enumerate(s1):
		n2= s2[i]
		if n1 != "-" and n2 != "-":
			if n1 == n2:
				tot+=1
			else: #n1 != n2
				diff+=1
				tot+=1
	if tot!= 0:
		div=(diff/tot)*100
		return div
	else:
		return None


#Prep for batch mode
inputcolumnNames=["Alignment","referenceSeq","ampliconStart", "ampliconStop"]
if args.batchMode:
	inputDF = pd.read_csv(args.batchMode, sep='\t', names=inputcolumnNames)
else:
	inputDF = pd.DataFrame([args.alignment, args.referenceName, args.start, args.stop], columns=inputcolumnNames)

data = []
# for i, row in inputDF.iterrows():
for row in range(len(inputDF)):

	#Reading in alignment file, returns dictionary containing name:aligned sequence.
	alignmentD= ft.read_fasta_dict_upper(inputDF["Alignment"][row])

	#Determine amplicon coords within alignment, using provided reference seq and **unaligned** coords
	cnt=0
	for i, n in enumerate(alignmentD[inputDF["referenceSeq"][row]], 1):
		if n != "-":
			cnt+=1
			if cnt == inputDF["ampliconStart"][row]:
				alignedStart=i
			elif cnt == inputDF["ampliconStop"][row]:
				alignedStop=i
				break
	
	#Opening output file
	#fout= open(args.outputName, "w")
	hd=0
	for s1, s2 in it.combinations(alignmentD.values(),2):
		#Calculate divergence across amplicon
		s1amplicon=s1[alignedStart:alignedStop]
		s2amplicon=s2[alignedStart:alignedStop]
		ampDiv= calcDivergence(s1amplicon, s2amplicon)
		if not ampDiv:
			continue
	
		#Create a list of windows for each sequence
		if not args.windowSize:
			args.windowSize= inputDF["ampliconStop"][row] - inputDF["ampliconStart"][row]
		s1windowsL, s2windowsL= seqWindows(s1, s2, args.windowSize, args.stepSize)
	
		#Calculate divergence at each window, across whole genome
		windowsaboveThresh=0
		windowDivL=[] #List w/ divergences for each window
		for i in range(len(s1windowsL)):
			windowDiv= calcDivergence(s1windowsL[i], s2windowsL[i])
			if windowDiv:
				windowDivL.append(windowDiv)
				if windowDiv > args.divThreshold:
					windowsaboveThresh+=1
		wgSWDiv= windowsaboveThresh/len(s1windowsL)*100	# % of windows > args.divThreshold% divergent
		
		#Calculate divergence across whole genome (i.e. not using SW)
		wgoverallDiv= calcDivergence(s1, s2)
		if not wgoverallDiv:
			continue
		
		data.append([os.path.basename(inputDF["Alignment"][row]), ampDiv, wgSWDiv, wgoverallDiv])
		print(data)
	outputcolumnNames=["Alignment","ampDiv","wgDiv (SW)","wgDiv (overall)"]
	
outputDF = pd.DataFrame(data, columns=outputcolumnNames)

# 	#Descriptive statistics
# 	if hd == 0:
# 		header= "Sequence #1\tSequence #2\tMaximum*\tQ3*\tMedian*\tQ1*\tMinimum*\tIQR*\tMean*\t\t*%% of windows >%f%% divergent" % args.divThreshold
# 		fout.write(header)
# 	maximum= max(windowDivL)
# 	q3= np.percentile(windowDivL, 75, interpolation = 'midpoint')
# 	median= np.percentile(windowDivL, 50, interpolation = 'midpoint')
# 	q1= np.percentile(windowDivL, 25, interpolation = 'midpoint')
# 	minimum= min(windowDivL)
# 	IQR= q3-q1
# 	mean= np.mean(windowDivL)
# 
# 	line= "\n%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (s1NAME,s2NAME,maximum,q3,median,q1,minimum,IQR,mean)		
# 	fout.write(line)
# 
# 	hd+=1
# 
# fout.close()

# print("ampDivL:",ampDivL)
# print("wgSWDivL:",wgSWDivL)

#Generate scatterplot with pairwise divergence comparisons (% of windows > args.divThreshold% divergent vs. amplicon divergence)
sns.lmplot(x="ampDiv", y="wgDiv (SW)", data=outputDF, fit_reg=False, hue="Alignment", legend=False)
plt.legend(loc='lower right')
plt.show()

# fig, ax = plt.subplots(figsize=(8,6))
# ax.scatter(ampDivL,wgSWDivL, s=50, alpha=0.50, c="#c51b8a")
# ylabelName= "%% of windows >%d%% divergent" % args.divThreshold
# ax.set_ylabel(ylabelName, fontsize=20)
# ax.set_xlabel("% divergence at amplicon", fontsize=20)
# ax.set_ylim(0,100)
# ax.set_xlim(0,100)
# ax.set_title("Pairwise comparison of divergence across\nwhole genome vs amplicon", fontsize=18)
# ax.tick_params(axis="both", which="major", labelsize=15)
# fig.savefig("scatter.pdf")
