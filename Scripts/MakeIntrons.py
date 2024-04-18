#!/bin/python3
##Last Modifed Dec19 by ALR
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Exons", help="Input BED file clustered sorted fixed exons")
parser.add_argument("Output", help="Output BED file for clustered annotated PACs")
args = parser.parse_args()

InFile = str(args.Exons)
Output = str(args.Output)

PrevFStop = ''
PrevFGene = ''
PrevFChr = ''
PrevRStop = ''
PrevRGene = ''
PrevRChr = ''

Out = open(Output, 'w')
with open(InFile, 'r') as In:
    line = In.readline().rstrip()
    while line:
        Data = line.split()
        Chr = Data[0]
        Start = int(Data[1])
        Stop = int(Data[2])
        Dir = Data[5]
        Gene = Data[3].split("_")[0]
        if Dir == '+':
            Name = Data[3].split("_")[0] + '_exon_0_0_' + Chr + '_' + str(Start) + '_f'
        else:
            Name = Data[3].split("_")[0] + '_exon_0_0_' + Chr + '_' + str(Stop) + '_r'
        Score = Data[4]
        if Dir == '+':
            if Gene == PrevFGene and Chr == PrevFChr:
                IntronStart = PrevFStop
                IntronStop = Start
                IntronName = Gene + '_intron_0_0_' + Chr + '_' + str(IntronStart) + '_f'
                Out.write(Chr + '\t' + str(IntronStart) + '\t' + str(IntronStop) + '\t' + IntronName + '\t1\t' + Dir + '\n')
                PrevFStop = Stop
            else:
                PrevFStop = Stop
                PrevFGene = Gene
                PrevFChr = Chr
        else:
            if Gene == PrevRGene and Chr == PrevRChr:
                IntronStart = PrevRStop
                IntronStop = Start
                IntronName = Gene + '_intron_0_0_' + Chr + '_' + str(IntronStop) + '_r'
                Out.write(Chr + '\t' + str(IntronStart) + '\t' + str(IntronStop) + '\t' + IntronName + '\t1\t' + Dir + '\n')
                PrevRStop = Stop
            else:
                PrevRStop = Stop
                PrevRGene = Gene
                PrevRChr = Chr
        Out.write(Chr + '\t' + str(Start) + '\t' + str(Stop) + '\t' + Name + '\t1\t' + Dir + '\n')
        line = In.readline().rstrip()
Out.close()