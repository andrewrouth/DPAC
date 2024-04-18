#!/bin/python3
##Last Modifed Apr19 by ALR
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Input", help="Input BED file with unclustered PACs e.g. hg19_PACs.bed. Must be sorted by count e.g. '$ sort -k4 -rn In.bed > In.sorted.bed'")
parser.add_argument("Output", help="Unmasked Output BED file for clustered annotated PACs")
parser.add_argument("Names", help="AccessionNames fasta")
parser.add_argument("--Description", help="Provide Description/Comment for BED9+ annotation")
args = parser.parse_args()

InFile = str(args.Input)
NameFile = str(args.Names)
if args.Description:
    Description = str(args.Description)
else:
    Description = 'PAC_Clusters'
    
Output = str(args.Output)

AllPACs = {}
m=0
Header = ''
with open(InFile, 'r') as In:
    line = In.readline()
    while line:
        if 'track' not in line:
            Data = line.split()
            if Data[3] == '.':
                m+=1
                Exon = Data[0] + ':' + Data[1]
                Name = "Unannotated_intergenic" + "_" + Data[0] + ':' + Data[1] + "_" + "PAC-" + str(m)
            else:
                Exon = Data[3]
                n=1
                Name = Exon + ":" + Exon[6] + "_PAC-" + str(n)
                while Name in AllPACs:
                    n+=1
                    Name = Root + "_" + Exon[2] + "_" + Data[0]  + ":" + Exon[-2] + "_PAC-" + str(n)
            AllPACs[Name] = Data
        else:
            Header = line
        line = In.readline()

with open(Output, 'w') as Out:
    Out.write(Header)
    for i in AllPACs:
        Out.write('\t'.join([AllPACs[i][0], AllPACs[i][1], AllPACs[i][2], i,  AllPACs[i][4], AllPACs[i][5], '\n']))