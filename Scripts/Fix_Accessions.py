#!/bin/python3
##Last Modifed Apr19 by ALR
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Input", help="Input BED file with unclustered PACs e.g. hg19_PACs.bed. Must be sorted by count e.g. '$ sort -k4 -rn In.bed > In.sorted.bed'")
parser.add_argument("Output", help="Unmasked Output BED file for clustered annotated PACs")
parser.add_argument("Names", help="AccessionNames fasta")
args = parser.parse_args()

InFile = str(args.Input)
NameFile = str(args.Names)
Output = str(args.Output)

Names = {}
with open(NameFile, 'r') as In:
    line = In.readline()
    while line:
        Data = line.split()
        Names[Data[0]] = Data[1]
        line = In.readline()
    
Out = open(Output, 'w')
with open(InFile, 'r') as In:
    line = In.readline()
    while line:
        Data = line.split()
        Exon = Data[3].split("_")
        Accession = Exon[0] + "_" + Exon[1]
        Name = Names[Accession]
        Data[3] = Name + "_" + "_".join(Exon[2:])
        Out.write('\t'.join(Data) + '\n') 
        line = In.readline()