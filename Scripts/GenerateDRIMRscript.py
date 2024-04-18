#!/bin/python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("COMP", help="Meta Data file")
parser.add_argument("Output", help="R Script Name")
parser.add_argument("-CountAdj", action='store_true', help="Include Intronic PACs")
#parser.add_argument("MinCount", help="R Script Name")
args = parser.parse_args()

COMP = str(args.COMP)
Output = str(args.Output)
#MinCount = str(args.MinCount)
Conditions = COMP.split('-vs-')
Conditions[0] = '1_' + Conditions[0]
Conditions[1] = '2_' + Conditions[1]
if args.CountAdj:
    CountAdj = True
else:
    CountAdj = False

def MakeOpen(Group):
    with open(Output + '_' + Group + '-counts_drimseq.txt', 'r') as In:
        File = Output + '_' + Group + '-counts_drimseq_fix.txt'
        Out = open(File, 'w')
        header = In.readline()
        Out.write('gene_id\tfeature_id\t' + header)
        Samples = header.split()
        Replicates = int(len(Samples)/2)
        ##Make Sample list from header
        Samples = '", "'.join(Samples)
        ##Make Conditions list from header
        Conds = []
        for i in range(Replicates):
            Conds.append('"' + Conditions[0] + '"')        
        for i in range(Replicates):
            Conds.append('"' + Conditions[1] + '"')
        Conds = ', '.join(Conds)
        ##Write count data into DRIMseq compatible file
        data = In.readline()
        while data:
            data = data.split()
            if Group == 'Exon':
                Name = data[0].split('_')[0]
                if CountAdj == True:
                    Counts = data[1:]
                    if Counts[:Replicates] == ['0']*Replicates or Counts[Replicates:] == ['0']*Replicates:
                        Counts = [str(float(i) + 1) for i in data[1:]]
                        Out.write('\t'.join([Name] + [data[0]] + Counts))
                    else:
                        Out.write('\t'.join([Name] + data))
                else:
                    Out.write('\t'.join([Name] + data))
                Out.write('\n')
            elif Group == 'PAC':
                Name = data[0].split('_PAC')[0]
                if CountAdj == True:
                    Counts = [str(float(i) + 1) for i in data[1:]]
                    Out.write('\t'.join([Name] + [data[0]] + Counts))
                else:
                    Out.write('\t'.join([Name] + data))
                Out.write('\n')
            data = In.readline()
        Out.close()
        return ''.join(['counts <- read.delim("' + File + '", header=TRUE)\n'
                        'samples <- data.frame(sample_id = c("' + Samples + '"), group = c(' + Conds + '))\n'
                        'd <- dmDSdata(counts = counts, samples = samples)\n'
                        'd <- dmFilter(d, min_samps_gene_expr = ' + str(Replicates * 2) + ', min_samps_feature_expr = ' + str(Replicates) +
                        ', min_gene_expr = 10, min_feature_prop = 0.05)\n'
                        'design_full <- model.matrix(~ group, data = samples(d))\n'
                        'set.seed(343434)\n'
                        'd <- dmPrecision(d, design = design_full)\n'
                        'd <- dmFit(d, design = design_full, verbose = 1)\n'
                        'd <- dmTest(d, coef = "group' + Conditions[1] +'", verbose = 1)\n'
                        'Outfile = paste("' + Output + '_' + Group + '_DRIMseq-results.csv", sep="")\n'
                        'write.csv(results(d, level = "feature"), file=Outfile)\n'
                        'Outfile = paste("' + Output + '_' + Group + '_DRIMseq-proportions.csv", sep="")\n'
                        'write.csv(proportions(d), file=Outfile)\n\n'])

Out = open(Output + '_DRIMRscript.txt', 'w')
Out.write(''.join(['#!/usr/bin/env Rscript\n' +
                    'args = commandArgs(trailingOnly=TRUE)\n' +
                    'library("DRIMSeq")\n' +
                    'library("ggplot2")\n\n']))

Groups = ['Exon', 'PAC']
for j in Groups:
    Out.write(MakeOpen(j))
Out.close()


