# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:22:59 2020

@author: AndrewRouth
"""

Counttable = 'HSPOutput_Gene-counts.txt'
Output = 'HSPOutput_Gene-counts_0-120_forclustering.txt'

InputFiles = ['S2_HS_0-vs-S2_HS_120_HSPOutput_Gene_DESeq2-results.csv']

Genes = set()

for i in InputFiles:
    with open(i, 'r') as In:
        Dheader = In.readline()
        Data = In.readline()
        while Data:
            Data = Data.split(',')
            Name = Data[1].replace('"','')
            try:
                padj = float(Data[7])
                fc = float(Data[3])
                if padj < 0.05:
                    if fc > 1 or fc < -1:
                        Genes.add(Name)
                else:
                    pass
            except:
                pass
            Data = In.readline()
            
n=0
Out = open(Output, 'w')
with open(Counttable, 'r') as In:
    header = In.readline()
    Out.write(header)
    Data = In.readline()
    while Data:
        temp = Data.split()
        if temp[0] in Genes:
            n+=1
            Out.write('\t'.join(temp[0:4] + temp[10:])+ '\n')
        else:
            pass
        Data = In.readline()
Out.close()