#!/bin/python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("Metadata", help="Meta Data file")
parser.add_argument("COMP", help="Meta Data file")
parser.add_argument("Output", help="R Script Name")
args = parser.parse_args()

Metadata = str(args.Metadata)
COMP = str(args.COMP)
Output = str(args.Output)

Conditions = {}
CondOrder = []
CondOrder1 = []
CondOrder2 = []

class ClassCondition(object):
    def __init__(self, Name):
        self.Name = Name
        self.Samples = []
        self.Indices = []

with open(Metadata, 'r') as In:
    header = In.readline()
    data = In.readline()
    n=1
    while data:
        [Name, Condition] = data.split()[0:2]
        if Condition not in Conditions:
            Conditions[Condition] = ClassCondition(Condition)
            CondOrder.append(Condition)
        else:
            pass
        CondOrder1.append(Name)
        CondOrder2.append(Condition)
        Conditions[Condition].Samples.append(Name)
        Conditions[Condition].Indices.append(n)
        n+=1
        data = In.readline()

for i in Conditions:
   # Conditions[i].Indices = ",".join([str(i) for i in Conditions[i].Samples.values()])
    Conditions[i].Indices = ",".join([str(i) for i in Conditions[i].Indices])

def MakeOpen(Group):
    Conds1 = []
    for i in CondOrder1:
        Conds1.append('"' + i + '"')
    Conds1 = ', '.join(Conds1)
    Conds2 = []
    for i in CondOrder2:
        Conds2.append('"' + i + '"')
    Conds2 = ', '.join(Conds2)
    File = Output + '_' + Group + '-counts_drimseq_fix.txt'
    Out = open(File, 'w')
    with open(Output + '_' + Group + '-counts_drimseq.txt', 'r') as In:
        header = In.readline()
        Out.write('gene_id\tfeature_id\t' + header)
        data = In.readline()
        while data:
            data = data.split()
            if Group == 'Exon':
                Name = data[0].split('_')[0]
                Out.write('\t'.join([Name] + data))
                Out.write('\n')
            elif Group == 'PAC':
                Name = data[0].split('_PAC')[0]
                Out.write('\t'.join([Name] + data))
                Out.write('\n')
            data = In.readline()
    Out.close()
    return ''.join(['counts <- read.delim("' + File + '", header=TRUE)\n'
                    'samples <- data.frame(sample_id = c(' + Conds1 + '), group = c(' + Conds2 + '))\n'
                    'd <- dmDSdata(counts = counts, samples = samples)\n'
                    'd <- dmFilter(d, min_samps_gene_expr = ' + str(len(CondOrder1)) + ', min_samps_feature_expr = ' + str(len(CondOrder1)/2) +
                    ', min_gene_expr = 10, min_feature_expr = 10)\n'
                    'design_full <- model.matrix(~ group, data = samples(d))\n'
                    'set.seed(343434)\n'
                    'd <- dmPrecision(d, design = design_full)\n'
                    'd <- dmFit(d, design = design_full, verbose = 1)\n'
                    'd <- dmTest(d, coef = "group' + CondOrder2[0] +'", verbose = 1)\n'
                    'Outfile = paste("' + Output + '_' + Group + '_DRIMseq-results.csv", sep="")\n'
                    'write.csv(results(d, level = "feature"), file=Outfile)\n'
                    'Outfile = paste("' + Output + '_' + Group + '_DRIMseq-proportions.csv", sep="")\n'
                    'write.csv(proportions(d), file=Outfile)\n\n'])

NumConds = len(Conditions)
Comparisons = [(i,j) for i in range(NumConds) for j in range(NumConds) if i < j]
Groups = ['Exon', 'PAC']

Out = open(Output + '_DRIMRscript.txt', 'w')
Out.write(''.join(['#!/usr/bin/env Rscript\n' +
                    'args = commandArgs(trailingOnly=TRUE)\n' +
                    'library("DRIMSeq")\n' +
                    'library("ggplot2")\n\n']))
for j in Groups:
    Out.write(MakeOpen(j))
Out.close()


