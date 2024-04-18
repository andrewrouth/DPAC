#!/bin/python3
##Last Modifed Dec19 by ALR
from subprocess import check_output
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Input", help="Input BED file with unclustered PASs e.g. hg19_PASs.bed. Must be sorted by count e.g. '$ sort -k4 -rn In.bed > In.sorted.bed'")
parser.add_argument("Output", help="Unmasked Output BED file for clustered annotated PACs")
parser.add_argument("Strand", help="State Strand, F or R")
parser.add_argument("Genome", help="Genome_Path fasta")
parser.add_argument("--ACount", help="Allowed ACount. Default = 12")
parser.add_argument("--Description", help= "Description in string format")
args = parser.parse_args()

InFile = str(args.Input)
Genome = str(args.Genome)

if str(args.Strand) == 'F':
    Strand = '+'
elif str(args.Strand) == 'R':
    Strand = '-'
else:
    print("ERROR: incorrect strandedness entered, must be 'F' of 'R'.")
    Strand = 'STRAND'
if args.Description:
    Description = str(args.Description)
else:
    Description = 'DPAC_PAC'    
if args.ACount:
    ACount = int(args.ACount)
else:
    ACount = 12 

###################
def Rev_Comp(Seq):
        Seq = Seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(Seq)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)[::-1]

print("Masking internally primed poly(A) sites")
Output = open('UNNAMED_' + str(args.Output), 'w')
Output2 = open('MASKED_' + str(args.Output), 'w')
Output2.write('track name="Masked_PASs" description="Masked_PASs_' + Description + '" visibility=pack color=0,0,255\n')
with open(InFile, 'r') as In:
        line = In.readline().rstrip()
        while line:
                Data = line.split('\t')
                Loc = int(Data[1])
                if Strand == "+":
                        cmd = Data[0] + ":" + str(Loc + 2) + "-" + str(Loc + 21)
                elif Strand == "-":
                        cmd = Data[0] + ":" + str(Loc - 20) + "-" + str(Loc - 1)
                else:
                        print("Enter either F or R for strandedness")
                        break
                try: 
                    Seq = check_output(['samtools', 'faidx', Genome, cmd], universal_newlines=True).split()[1]
                except:
                        cmd = 'chr' + cmd
                        #line = 'chr' + line
                        try:
                            Seq = check_output(['samtools', 'faidx', Genome, cmd], universal_newlines=True)
                        except:
                                Seq = ''
                if Seq:
                        Seq = Seq.upper()
                        if Strand == "-":
                                Seq = Rev_Comp(Seq)
                        else:
                                pass
                        Count = [Seq.count("A"), Seq.count("T"), Seq.count("G"), Seq.count("C")]
                        if Count[0] > ACount:
                                Output2.write(line + '\n')
                        else:
                            Output.write(line + '\n')
                else:
                        print("PAS_Clustering.py:  Failed locus in index: ", cmd[3:])
                        Output2.write(line + '\n')
                line = In.readline().rstrip()
Output.close()
Output2.close()