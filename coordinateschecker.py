#! /usr/bin/python3

# -*- coding: utf-8 -*-

## Importing python libraries
from __future__ import print_function
import argparse, os.path, re, sys, warnings
from BCBio import GFF
from Bio import SeqIO
from time import strftime

#import BCBio.GFF
#import csv
#import fractions
#import glob
#import multiprocessing
#import numpy
#import os
#import subprocess
#import time
#from Bio import SeqFeature
#try:
#	from Bio.Alphabet import IUPAC
#except ImportError:
#	IUPAC = None
#from Bio.Seq import Seq
#from Bio.SeqFeature import FeatureLocation
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqUtils.ProtParam import ProteinAnalysis
#from collections import OrderedDict, defaultdict
#from itertools import product
#from pathlib import Path
#from scipy import signal

warnings.filterwarnings("ignore")

## Processing the arguments
parser = argparse.ArgumentParser(description='coordinateschecker.py is a script that checks if the genes are found in the same coordinates than a reference and retrieves a table with the IDs and the numnber of TP, FP and FN')
basic_group = parser.add_argument_group('Basic options for coordinateschecker [REQUIRED]')
basic_group.add_argument("--reference", dest="referencefile", type=str, required=True, help='Reference file as a GenBank file', metavar="GENBANKFILE")
basic_group.add_argument("--input", dest="inputfile", type=str, required=True, help='Input file from the ORF caller program', metavar="INPUTFILE")

advanced_general_group = parser.add_argument_group('Advanced general options for genecounter [OPTIONAL]')
advanced_general_group.add_argument("--out", dest="rootoutput", type=str, help='Name of the outputs files (without extension)', metavar="OUTPUTNAME")
advanced_general_group.add_argument("--typefile", dest="typefile", type=str, default='gff_cds', help='Type of inputfile. At this moment, it can process GenBank, GFF3, MGA output files and outfiles from other ORF callers (Default: gff)', metavar="TYPEFILE")

args = parser.parse_args()

referencefile = args.referencefile

root_output = args.rootoutput
if not root_output:
	root_output = '{}_counted'.format(os.path.splitext(args.inputfile)[0])

typefile = args.typefile
if not typefile:
	typefile = "gff_cds"

## Printing the header of the program 
print("This is Coordinateschecker.py")
print("Written by Enrique Gonzalez Tortuero")
print("Local time: ", strftime("%a, %d %b %Y %H:%M:%S"))

## Processing the reference file
reference = {}
with open(referencefile, "r") as reffile:
	for gb_record in SeqIO.parse(reffile, "genbank"):
		reference[gb_record.name] = []
		for feat in gb_record.features:
			if feat.type == 'CDS':
				coords = '{}..{}'.format(feat.location.start, feat.location.end)
				reference[gb_record.name].append(coords)

## Processing the inputfile
observations = {}
if typefile == "genbank":
	with open(args.inputfile, "r") as infile:
		for gb_record in SeqIO.parse(infile, typefile):
			observations[gb_record.name] = []
			for feat in gb_record.features:
				if feat.type == 'CDS':
					coords = '{}..{}'.format(feat.location.start, feat.location.end)
					observations[gb_record.name].append(coords)
elif typefile == "gff_cds":
	# For PRODIGAL, METAPRODIGAL and FRAGGENESCAN
	with open(args.inputfile, "r") as infile:
		for rec in GFF.parse(infile):
			observations[rec.id] = []
			for feat in rec.features:
				if feat.type == 'CDS':
					coords = '{}..{}'.format(feat.location.start, feat.location.end)
					observations[rec.id].append(coords)					
elif typefile == "gff_gene":
	# For AUGUSTUS and GLIMMER
	with open(args.inputfile, "r") as infile:
		for rec in GFF.parse(infile):
			observations[rec.id] = []
			for feat in rec.features:
				if feat.type == 'gene':
					coords = '{}..{}'.format(feat.location.start, feat.location.end)
					observations[rec.id].append(coords)					
elif typefile == "mga":
	with open(args.inputfile, "r") as infile:
		for line in infile:
 			if line.startswith("#"):
 				if not line.startswith("# gc") and not line.startswith("# self"):
 					seqid = re.split("\s", line)[1]
 					observations[seqid] = []
 			else:
 				linesplit = re.sub("\s+", ",", line)
 				linesplit2 = re.split(",", linesplit)
 				coords = '{}..{}'.format(int(linesplit2[1])-1, int(linesplit2[2])) # Start coordinates had added 1 bp in the MGA output.
 				observations[seqid].append(coords)
elif typefile == "lst":
	with open(args.inputfile, "r") as infile:
		lines = infile.readlines()[7:]
		for line in lines:
			if line.startswith("FASTA definition line:"):
				seqid = re.match(re.compile('FASTA\sdefinition\sline:\s(\S+)\s'), line).groups()[0]
				observations[seqid] = []
			elif line.startswith("Predicted genes"):
				continue
			else:
				line = line.lstrip()
				linesplit = re.sub("\s+", ",", line)
				linesplit2 = re.split(",", linesplit)
				if not linesplit2[0] == "Gene" and not linesplit2[0] == "#" and not linesplit2[0] == "":
					startcoord = linesplit2[2]
					endcoord = linesplit2[3]
					if startcoord.startswith("<"):
						startcoord = startcoord.replace("<", "")
					if endcoord.startswith(">"):
						endcoord = endcoord.replace(">", "")
					coords = '{}..{}'.format(int(startcoord)-1, int(endcoord)) # Start coordinates had added 1 bp in the GeneMarkS output.
					observations[seqid].append(coords)
else:
	sys.exit('Inputfile is not interpretable')

## Calculating TP, FP, FN
dicttable = {}
for genomeid in reference:
	for genomeid2 in observations:
		if genomeid == genomeid2:
			dicttable[genomeid] = {}
			dicttable[genomeid]["TP"] = len(set(reference[genomeid]).intersection(observations[genomeid]))
			dicttable[genomeid]["FP"] = len(set(observations[genomeid])-set(observations[genomeid]).intersection(reference[genomeid]))
			dicttable[genomeid]["FN"] = len(set(reference[genomeid])-set(reference[genomeid]).intersection(observations[genomeid]))

# Printing the output
with open("%s.tsv" % root_output, "w") as tablefile:
	tablefile.write("Genome\tTP\tFP\tFN\n")
	for genome in dicttable:
		tablefile.write("%s\t%s\t%s\t%s\n" % (genome, dicttable[genome]["TP"], dicttable[genome]["FP"], dicttable[genome]["FN"]))
