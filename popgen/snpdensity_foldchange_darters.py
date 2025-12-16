#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"This script reads in female and male SNP density data and calculates and outputs M:F SNP density fold change for each chromosome block"
import argparse
import sys
import os
from collections import defaultdict

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("females", type=str, help="A file with average SNP density data from females")
parser.add_argument("males", type=str, help="A file with average SNP density data from males")
parser.add_argument("outfile", type=str, help="Output file containing fold change values")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Function to extract SNP density from the input files
def extract_snps(source):
    snpsdict = defaultdict(list)
    with open(source, "r") as infile:
        next(infile)  # Skip header
        for line in infile:
            line = line.rstrip().split(",")
            # Chromosome names may have underscores, so use the first 3 fields carefully
            chromosome = line[0]
            windowStart = line[1]
            windowEnd = line[2]
            chromosome_block = chromosome + "_" + windowStart + "_" + windowEnd
            sum_snpdensity = line[4]
            average_snpdensity = line[5]
            logaverage_snpdensity = line[6]
            snpsdict[chromosome_block].append(sum_snpdensity)
            snpsdict[chromosome_block].append(average_snpdensity)
            snpsdict[chromosome_block].append(logaverage_snpdensity)
    return snpsdict

# Function to combine SNP density from females and males and compute fold change
def combine_density(females_snps, males_snps):
    combined_dict = defaultdict(list)
    for block in females_snps:
        fem_sum = float(females_snps[block][0])
        fem_average = float(females_snps[block][1])
        fem_logaverage = float(females_snps[block][2])
        mal_sum = float(males_snps[block][0])
        mal_average = float(males_snps[block][1])
        mal_logaverage = float(males_snps[block][2])
        combined_logaverage = mal_logaverage - fem_logaverage
        combined_dict[block].append(combined_logaverage)
        combined_dict[block].append(mal_average)
        combined_dict[block].append(mal_logaverage)
        combined_dict[block].append(fem_average)
        combined_dict[block].append(fem_logaverage)
    return combined_dict

# Main function
def main():
    # Extract SNP data from females and males files
    females_snps = extract_snps(args.females)
    print("Number of female chromosome blocks =", len(females_snps))
    males_snps = extract_snps(args.males)
    print("Number of male chromosome blocks =", len(males_snps))

    # Check if the chromosome blocks match between females and males
    for block in females_snps:
        if block not in males_snps:
            print("Error: Block missing in males data:", block)
    for block in males_snps:
        if block not in females_snps:
            print("Error: Block missing in females data:", block)

    # Combine SNP density for matched blocks
    combined_density = combine_density(females_snps, males_snps)
    print("Number of filtered chromosome blocks =", len(combined_density))

    # Write the combined data to the output file
    with open(args.outfile, "w") as outfile:
        header = "ChromosomeBlock,Chromosome,WindowStart,WindowEnd,MFLogaverage,Maverage,Mlogaverage,Faverage,Flogaverage"
        outfile.write(header)
        outfile.write("\n")
        for block in combined_density:
            # Since chromosome names can contain underscores, we use slicing for the first 3 fields
            fields = block.split("_")
            chromosome = "_".join(fields[:-2])  # Join all fields except the last two for the chromosome
            windowStart = fields[-2]
            windowEnd = fields[-1]
            outfile.write(block)
            outfile.write(",")
            outfile.write(chromosome)
            outfile.write(",")
            outfile.write(windowStart)
            outfile.write(",")
            outfile.write(windowEnd)
            outfile.write(",")
            outfile.write(str(combined_density[block][0]))
            outfile.write(",")
            outfile.write(str(combined_density[block][1]))
            outfile.write(",")
            outfile.write(str(combined_density[block][2]))
            outfile.write(",")
            outfile.write(str(combined_density[block][3]))
            outfile.write(",")
            outfile.write(str(combined_density[block][4]))
            outfile.write("\n")

if __name__ == '__main__':
    main()
