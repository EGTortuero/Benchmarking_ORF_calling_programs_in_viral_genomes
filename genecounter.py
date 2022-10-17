#! /usr/bin/python3

# -*- coding: utf-8 -*-

## Importing python libraries
from __future__ import print_function
import argparse, os.path, re, sys
from BCBio import GFF
from Bio import SeqIO
from time import strftime

## Processing the arguments
parser = argparse.ArgumentParser(description='genecounter.py is a script that interprets output files from ORF callers and retrieve a table with the IDs and the numnber of genes')
basic_group = parser.add_argument_group('Basic options for genecounter [REQUIRED]')
basic_group.add_argument("--input", dest="inputfile", type=str, required=True, help='Input file from the ORF caller program', metavar="INPUTFILE")

advanced_general_group = parser.add_argument_group('Advanced general options for genecounter [OPTIONAL]')
advanced_general_group.add_argument("--out", dest="rootoutput", type=str, help='Name of the outputs files (without extension)', metavar="OUTPUTNAME")
advanced_general_group.add_argument("--typefile", dest="typefile", type=str, default='gff', help='Type of inputfile. At this moment, it can process GenBank, GFF3, MGA output files and outfiles from other ORF callers (Default: gff)', metavar="TYPEFILE")

args = parser.parse_args()

inputfile = args.inputfile

typefile = args.typefile
if not typefile:
	typefile = "gff"

root_output = args.rootoutput
if not root_output:
	root_output = '{}_counted'.format(os.path.splitext(args.inputfile)[0])

## Printing the header of the program 
print("This is Genecounter.py")
print("Written by Enrique Gonzalez Tortuero")
print("Local time: ", strftime("%a, %d %b %Y %H:%M:%S"))

## Retrieving the inputfile for the keyIDs
if typefile == "genbank":
	with open("%s.tsv" % root_output, "w") as tablefile, open(inputfile, "r") as infile:
		tablefile.write("ID\tCDSs\n")
		for gb_record in SeqIO.parse(infile, typefile):
			counter = 0
			for feat in gb_record.features:
				if feat.type == 'CDS':
					counter += 1
			tablefile.write("%s\t%i\n" % (gb_record.name, counter))
elif typefile == "gff":
	with open("%s.tsv" % root_output, "w") as tablefile, open(inputfile, "r") as infile:
		tablefile.write("ID\tCDSs\n")
		for rec in GFF.parse(infile):
			counter = 0
			for feat in rec.features:
				if feat.type == 'gene': # Change to "CDS" for PRODIGAL, METAPRODIGAL and FRAGGENESCAN outputs; "gene" for AUGUSTUS and GLIMMER
					counter += 1
			tablefile.write("%s\t%i\n" % (rec.id, counter))
elif typefile == "mga":
	infotable = {}
	with open("%s.tsv" % root_output, "w") as tablefile, open(inputfile, "r") as infile:
		tablefile.write("ID\tCDSs\n")
		for line in infile:
			if line.startswith("#"):
				if not line.startswith("# gc") and not line.startswith("# self"):
					seqid = re.split("\s", line)[1]
					counter = 0
			else:
				counter += 1
				infotable[seqid] = counter
		for seq, nocds in infotable.items():
			tablefile.write("%s\t%i\n" % (seq, nocds))
elif typefile == "lst":
	infotable = {}
	with open("%s.tsv" % root_output, "w") as tablefile, open(inputfile, "r") as infile:
		tablefile.write("ID\tCDSs\n")
		lines = infile.readlines()[7:]
		for line in lines:
			if line.startswith("FASTA definition line:"):
				seqid = re.match(re.compile('FASTA\sdefinition\sline:\s(\S+)\s'), line).groups()[0]
				counter = 0
			elif line.startswith("Predicted genes"):
				continue
			else:
				linesplit = re.sub("\s+", ",", line)
				linesplit2 = re.split(",", linesplit)
				firstelementtoignore = linesplit2.pop(0)
				if not linesplit2[0] == "Gene" and not linesplit2[0] == "#":
					counter += 1
					infotable[seqid] = counter
		for seq, nocds in infotable.items():
			tablefile.write("%s\t%i\n" % (seq, nocds))
else:
	sys.exit('Inputfile is not interpretable')
