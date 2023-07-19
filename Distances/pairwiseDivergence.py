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

parser.add_argument("-o", "--outputName", default="divergenceSummary.tsv", metavar='\b', help="... [XXX]")
parser.add_argument("-v", "--visualization", default="scatterplot-pairwiseDiv.png", metavar='\b', help="... [XXX]")
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

#Opening output file
fout= open(args.outputName, "w")
hd=0
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

	for nseq1, nseq2 in it.combinations(alignmentD.items(),2):
		s1name,s1= nseq1[0],nseq1[1]
		s2name,s2= nseq2[0],nseq2[1]
		
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
		
		#Descriptive statistics for output TSV
		if hd == 0:
			header= "Alignment\tSequence #1\tSequence #2\tampliconDiv (%%)\twholegenomeDiv (%%)\t%% of windows >%.2f%% divergent\tMaximum*\tQ3*\tMedian*\tQ1*\tMinimum*\tIQR*\tMean*\t\t*Of all divergences, from each window" % args.divThreshold
			fout.write(header)
		maximum= max(windowDivL)
		q3= np.percentile(windowDivL, 75, interpolation = 'midpoint')
		median= np.percentile(windowDivL, 50, interpolation = 'midpoint')
		q1= np.percentile(windowDivL, 25, interpolation = 'midpoint')
		minimum= min(windowDivL)
		IQR= q3-q1
		mean= np.mean(windowDivL)
		
		line= "\n%s\t%s\t%s\t%f\t%f\t%f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (inputDF["Alignment"][row],s1name,s2name,ampDiv,wgoverallDiv,wgSWDiv,maximum,q3,median,q1,minimum,IQR,mean)
		fout.write(line)
		
		hd+=1
		
	outputcolumnNames=["Alignment","ampDiv","wgDiv (SW)","wgDiv (overall)"]
outputDF = pd.DataFrame(data, columns=outputcolumnNames)
fout.close()


####Visualizations
fig, ax = plt.subplots(2,1, figsize=(7,6))

#Generate scatterplot with pairwise divergence comparisons (% of windows > args.divThreshold% divergent vs. amplicon divergence)
#legend_map = 
v0 = sns.scatterplot(x="ampDiv", y="wgDiv (SW)", ax=ax[0], data=outputDF, hue="Alignment", s=50, alpha=0.75, legend=False)
v0.set_ylabel(("%% of windows\n>%d%% divergent" % args.divThreshold), fontsize=14)
v0.set_xlabel("Divergence across amplicon (%)", fontsize=14)
v0.tick_params(axis="both", which="major", labelsize=11)
v0.set_ylim(0,100)
v0.set_xlim(0,100)
v0.set_title("Pairwise comparison of whole genome divergence (via SW)\nvs. amplicon divergence", fontsize=14)

#Generate scatterplot with pairwise divergence comparisons (overall whole genome divergence vs. amplicon divergence)
v1 = sns.scatterplot(x="ampDiv", y="wgDiv (overall)", ax=ax[1], data=outputDF, hue="Alignment", s=50, alpha=0.75, legend=False)
v1.axline((0, 0), linestyle='dotted', color='black', slope=1)
v1.set_ylabel("Divergence across\nwhole genome (%)", fontsize=14)
v1.set_xlabel("Divergence across amplicon (%)", fontsize=14)
v1.tick_params(axis="both", which="major", labelsize=11)
v1.set_ylim(0,100)
v1.set_xlim(0,100)
v1.set_title("Pairwise comparison of whole genome divergence (overall)\nvs. amplicon divergence", fontsize=14)

fig.tight_layout(pad=2.5)
#fig.legend(bbox_to_anchor=(1.02, 1), loc='center right', title='Alignment')
fig.savefig(args.visualization, bbox_inches='tight', dpi=200)
