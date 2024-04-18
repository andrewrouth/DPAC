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
        Conditions[Condition].Samples.append(Name)
        Conditions[Condition].Indices.append(n)
        n+=1
        data = In.readline()

for i in Conditions:
   # Conditions[i].Indices = ",".join([str(i) for i in Conditions[i].Samples.values()])
    Conditions[i].Indices = ",".join([str(i) for i in Conditions[i].Indices])

def MakeSingle(Cond1, Cond2, Group):
    return ''.join(['res <- results(dds, contrast = c("condition","', Conditions[Cond2].Name, '","', Conditions[Cond1].Name,'"), filterFun=ihw, parallel=TRUE)\n',
                   'resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE))[c(', Conditions[Cond1].Indices, ',', Conditions[Cond2].Indices,')], by="row.names", sort=FALSE)\n',
                   'Outfile = paste("', Conditions[Cond1].Name,'-vs-', Conditions[Cond2].Name, '_', Output, '_', Group, '-counts_drimseq.txt")\n'
                   'write.table(data.frame(counts(dds, normalized=TRUE))[c(', Conditions[Cond1].Indices, ',', Conditions[Cond2].Indices,')], file=Outfile, sep="\t", quote = FALSE)\n'
                   'names(resdata)[1] <- "Gene"\n'
                   'Outfile = paste("', Conditions[Cond1].Name,'-vs-', Conditions[Cond2].Name, '_', Output, '_', Group, '_DESeq2-results.csv", sep="")\n'
                   'write.csv(resdata, file=Outfile)\n'])

def MakeOpen(Group):
    reps = []
    for i in CondOrder:
        reps.append('rep("' + i + '", ' + str(len(Conditions[i].Samples)) + ')')
    reps = ', '.join(reps)
    return ''.join(['Infile = paste(args[1], "_' + Group + '-counts.txt", sep="")\n' +
                    'counts <- read.delim(Infile, row.names = 1, header=TRUE)\n' +
                    'counts <- as.matrix(counts)\n' +
                    'condition <- factor(c(' + reps + '))\n' +
                    'coldata <- data.frame(row.names=colnames(counts), condition)\n' +
                    'dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)\n' +
                    'dds <- DESeq(dds, fitType="local", parallel=TRUE)\n'
                    'vsd <- vst(dds, blind=FALSE)\n'
                    'outcount = paste(args[1], "_' + Group + '_normcounts_Genes.txt", sep="")\n'
                    'normcounts <- counts(dds, normalized = TRUE)\n'
                    'write.table(normcounts, file=outcount, sep="\t", quote = FALSE)\n'
                    'outpdf = paste(args[1], "_' + Group + '_PCA_plot_Genes.pdf", sep="")\n'
                    'pdf(outpdf)\n'
                    'plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name),vjust=2)\n'
                    'dev.off()\n'])

NumConds = len(Conditions)
Comparisons = [(i,j) for i in range(NumConds) for j in range(NumConds) if i < j]
Out = open(Output + '_Comparisons.txt', 'w')
for i in Comparisons:
    Out.write('-'.join([CondOrder[i[0]], 'vs', CondOrder[i[1]]]))
    Out.write('\n')
Out.close()

Groups = ['Gene', 'Exon', 'PAC']

Out = open(Output + '_Rscript.txt', 'w')
Out.write(''.join(['#!/usr/bin/env Rscript\n' +
                    'args = commandArgs(trailingOnly=TRUE)\n' +
                    'library("DESeq2")\n' +
                    'library("BiocParallel")\n' +
                    'register(MulticoreParam(4))\n' +
                    'library("ggplot2")\n' +
                    'library("IHW")\n']))
for j in Groups:
    Out.write(MakeOpen(j))
    for i in Comparisons:
        print(j, CondOrder[i[0]], 'vs', CondOrder[i[1]])
        Out.write(MakeSingle(CondOrder[i[0]], CondOrder[i[1]], j))
Out.close()


