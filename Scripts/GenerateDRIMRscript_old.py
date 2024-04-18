#!/bin/python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("Metadata", help="Meta Data file")
parser.add_argument("Output", help="R Script Name")
args = parser.parse_args()

Metadata = str(args.Metadata)
Output = str(args.Output)

Conditions = {}
CondOrder = []
CondOrders = []
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
        CondOrders.append(Name)
        Conditions[Condition].Samples.append(Name)
        Conditions[Condition].Indices.append(n)
        n+=1
        data = In.readline()

for i in Conditions:
    Conditions[i].Indices = ",".join([str(i) for i in Conditions[i].Indices])

def MakeOpen():
    Conds = []
    for i in CondOrder:
        Conds.append('"' + i + '"')
    Conds = ', '.join(Conds)
    Conds2 = []
    for i in CondOrders:
        Conds2.append('"' + i + '"')
    Conds2 = ', '.join(Conds2)
    return ''.join(['counts <- read.delim(' + Output + '_PAC-counts_drimseq.txt, header=TRUE)\n'
                    'samples <- data.frame(sample_id = c(' + Conds2 + '), group = c(' + Conds + '))\n'
                    'd <- dmDSdata(counts = counts, samples = samples)\n'
                    'd <- dmFilter(d, min_samps_gene_expr = ' + str(len(Conds2)) + ', min_samps_feature_expr = 3, min_gene_expr = 10, min_feature_expr = 10)\n'
                    'design_full <- model.matrix(~ group, data = samples(d))\n'
                    'd <- dmPrecision(d, design = design_full)\n'
                    'd <- dmFit(d, design = design_full, verbose = 1)\n'
                    'd <- dmTest(d, coef = "groupKD", verbose = 1)\n'
                    'Outfile = paste("' + Output + '_DRIMseq-results.csv", sep="")\n'
                    'write.csv(results, file=Outfile)'])

NumConds = len(Conditions)
Out = open(Output + '_DRIMRscript.txt', 'w')
Out.write(''.join(['#!/usr/bin/env Rscript\n' +
                    'args = commandArgs(trailingOnly=TRUE)\n' +
                    'library("DRIMSeq")\n' +
                    'library("ggplot2")\n']))
Out.write(MakeOpen())
Out.close()


