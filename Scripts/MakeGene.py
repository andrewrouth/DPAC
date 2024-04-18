#!/bin/python3
##Last Modifed Dec19 by ALR
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Exons", help="Input BED file clustered sorted fixed exons")
parser.add_argument("Output", help="Output BED file for clustered annotated PACs")
args = parser.parse_args()

InFile = str(args.Exons)
Output = str(args.Output)

PrevFStart = 0
PrevFStop = 0
PrevFGene = ''
PrevFChr = ''

PrevRStart = 0
PrevRStop = 0
PrevRGene = ''
PrevRChr = ''

F=False
R=False

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
            if Gene == PrevFGene and Chr == PrevFChr:
                PrevFStop = Stop
            else:
                if F:
                    Out.write(PrevFChr + '\t' + str(PrevFStart) + '\t' + str(PrevFStop) + '\t' + PrevFGene + '\t1\t' + Dir + '\n')
                else:
                    F=True
                PrevFStart = Start
                PrevFStop = Stop
                PrevFGene = Gene
                PrevFChr = Chr
        else:
            if Gene == PrevRGene and Chr == PrevRChr:
                PrevRStop = Stop
            else:
                if R:    
                    Out.write(PrevRChr + '\t' + str(PrevRStart) + '\t' + str(PrevRStop) + '\t' + PrevRGene + '\t1\t' + Dir + '\n')
                else:
                    R=True
                PrevRStart = Start
                PrevRStop = Stop
                PrevRGene = Gene
                PrevRChr = Chr
        line = In.readline().rstrip()
Out.close()