#!/usr/bin/env python
# -*- coding: utf-8 -*-
"This script reads SNP density information from multiple samples and calculates and outputs average SNP density for each chromosome block"
import argparse
import sys
import os
from collections import defaultdict
import numpy as np

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str, nargs="+", help="A folder with individual SNP density files")
parser.add_argument("outfile", type=str, help="")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Function to list all files in the folder
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".txt")]

# Function to extract SNP density from the input file
def extract_snpdensity(source, snpdict):
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            line = line.split()
            # Adjust to accept lines with at least 9 columns (since your data has 9 columns)
            if len(line) >= 9:
                try:
                    # Extract the SNP density from the last column (index 8)
                    snpdensity = float(line[8])  # SNP density is in the 9th column
                    scaffold = f"{line[1]}_{line[2]}_{line[3]}"  # Use columns 2, 3, and 4 as key
                    snpdict[scaffold].append(snpdensity)
                except ValueError:
                    print(f"Error: Non-numeric value in line {line}")
            else:
                print(f"Warning: Skipping line with unexpected format: {line}")
    return snpdict

# Function to calculate average SNP density
def average_snpdensity(snpdensity):
    averagedict = defaultdict(list)
    for sequence in snpdensity:
        snplist = snpdensity[sequence]
        if len(snplist) == 0:
            continue
        sum_snplist = sum(snplist)
        average = (sum_snplist / len(snplist)) if len(snplist) > 0 else 0
        snp_values = '|'.join(map(str, snplist))
        logaverage = np.log2(average + 1)
        sumsnp = sum_snplist
        averagedict[sequence].append(snp_values)
        averagedict[sequence].append(sumsnp)
        averagedict[sequence].append(average)
        averagedict[sequence].append(logaverage)
    return averagedict

# Main function
def main():
    if len(args.infolder) == 1:
        if os.path.isdir(args.infolder[0]):
            infiles = list_folder(args.infolder[0])
        else:
            infiles = args.infolder
    else:
        infiles = args.infolder

    if not infiles:
        print("Error: No files found in the input directory.")
        sys.exit(1)

    print(f"Files to process: {infiles}")

    snpdensitydict = defaultdict(list)
    for infile in infiles:
        print(f"Processing file: {infile}")
        snpdensity = extract_snpdensity(infile, snpdensitydict)
        print(f"Number of chromosome blocks with SNP density = {len(snpdensity)}")

    average = average_snpdensity(snpdensity)
    print(f"Number of chromosome blocks with average SNP density = {len(average)}")

    with open(args.outfile, "w") as outfile:
        print(f"Writing output to {args.outfile}")
        header = "Chromosome,WindowStart,WindowEnd,SnpDensityList,Sum,Average,Logaverage"
        outfile.write(header)
        outfile.write("\n")
        for sequence in average:
            # Split sequence only by the last two underscores to avoid problems with extra underscores in chromosome names
            parts = sequence.rsplit("_", 2)
            if len(parts) == 3:
                chromosome, windowStart, windowEnd = parts
                outfile.write(f"{chromosome},{windowStart},{windowEnd},{average[sequence][0]},{average[sequence][1]},{average[sequence][2]},{average[sequence][3]}\n")
            else:
                print(f"Error: Unexpected sequence format: {sequence}")

if __name__ == '__main__':
    main()
