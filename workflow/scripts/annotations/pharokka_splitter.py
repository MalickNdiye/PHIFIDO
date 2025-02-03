# This script opens the output of Pharokka and splits it into separate files for each vOTU

import os
import sys
import pandas as pd
import argparse as arg

# set arguments
parser = arg.ArgumentParser(description='Split Pharokka output into separate files for each vOTU')
parser.add_argument('-i', '--input', help='Pharokka output directory', required=True)
parser.add_argument('-d', '--drep', help='dereplication data', required=True)
parser.add_argument('-v', '--vOTU', help='target vOTU', required=True)
parser.add_argument('-o', '--output', help='Output directory', required=True)

args = parser.parse_args()

print("Starting Pharokka splitter")
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# read dereplication data, header is first row
drep = pd.read_csv(args.drep, sep='\t', header=0)

# filter rows with column "vOTU" equal target vOTU
drep_filt=drep[drep["vOTU"]==args.vOTU]
if drep_filt.empty:
    print("No vOTU " + args.vOTU + " found in dereplication data")
    sys.exit()

# get list of genomes
genomes = drep_filt["genome"].tolist()

# create directories "gbk", "gff", "faa", "fna" in output directory
os.mkdir(args.output + "/gbk")
os.mkdir(args.output + "/gff")
os.mkdir(args.output + "/faa")
os.mkdir(args.output + "/fna")

# in the pharokka output directory there are the directories "single_gffs", "single_gbks", "single_fastas", "single_faas"
# copy the files corresponding to the genomes in the list to the corresponding directories in the output directory
print("Copying annotation files for vOTU " + args.vOTU)
for genome in genomes:
    os.system("cp " + args.input + "/single_gffs/" + genome + ".gff " + args.output + "/gff/")
    os.system("cp " + args.input + "/single_gbks/" + genome + ".gbk " + args.output + "/gbk/")
    os.system("cp " + args.input + "/single_fastas/" + genome + ".fasta " + args.output + "/fna/")
    os.system("cp " + args.input + "/single_faas/" + genome + ".faa " + args.output + "/faa/")

# concatenate all files in each directory in a single file namesd "all_genomes" + extension
os.system("cat " + args.output + "/gff/* > " + args.output + "/gff/all_genomes.gff")
os.system("cat " + args.output + "/gbk/* > " + args.output + "/gbk/all_genomes.gbk")
os.system("cat " + args.output + "/fna/* > " + args.output + "/fna/all_genomes.fna")
os.system("cat " + args.output + "/faa/* > " + args.output + "/faa/all_genomes.faa")


print("Done, all files copired in " + args.output)

