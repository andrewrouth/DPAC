import pickle as pickle
import argparse
#import numpy as np
#from APA_Terminal_Exon import Gene

parser = argparse.ArgumentParser()
parser.add_argument("Root", help="Root of Datafiles")
parser.add_argument("Output", help="Report Name")
args = parser.parse_args()
Root = str(args.Root)

GeneDict = pickle.load(open(Root + '.GeneDict.py.pi','rb'))


##Class for unannotated PACs
class Unannotated(object):
    def __init__(self, Name, baseMean, log2FoldChange, pvalue, padj, weight, counts):
        self.Name = Name
        self.Counts = counts
        self.baseMean = round(float(baseMean), 4)
        self.Locus = PAC_library[Name]
        try:
            self.FoldChange = round(float(log2FoldChange), 4)
            self.pvalue = round(float(pvalue), 4)
            self.padj = round(float(padj), 4)
        except:
            ##For cases were e.g. padj == 'NA'
            self.log2FoldChange = log2FoldChange
            self.pvalue = pvalue
            self.padj = padj

    def __str__(self):
        return  str([self.Name,
                self.baseMean,
                self.FoldChange,
                self.pvalue,
                self.padj])

##Class for Gene DE info
class Gene(object):
    def __init__(self, Name, baseMean, log2FoldChange, pvalue, padj, weight, counts):
        self.Exons = {}
        self.Introns = {}
        self.Name = Name
        self.Counts = counts
        self.MeanCountsCond1 = str(round((sum(np.array(self.Counts[:Replicates], dtype=float)) / Replicates), 4))
        self.MeanCountsCond2 = str(round((sum(np.array(self.Counts[Replicates:], dtype=float)) / Replicates), 4))

        #GeneResult = [DE, TE, APA]
        self.GeneResult = ['','No','No']
        
        #ExonResult = [#Exons wPAC, ExonFracs Cond1, ExonFracs Cond2, Fraction Change, Changers, Exon Names, httpslink]  
        self.ExonResult = [0, [], [], [], [], [], []]
        
        ##PACResult = [#PACs per Exon, PD Cond1, PD Cond2, PDUI, Names]
        self.PACResult = [[], [], [], [], [], []]
        
        self.baseMean = round(float(baseMean), 4)
        try:
            self.FoldChange = round(float(log2FoldChange), 4)
            self.pvalue = round(float(pvalue), 4)
            self.padj = round(float(padj), 4)
        except:
            self.FoldChange = log2FoldChange
            self.pvalue = pvalue
            self.padj = padj

    def FinishExons(self):
        #Conditions = 2
        TotalCounts = [0] * Conditions# * Replicates
        for i in self.Exons:
            TotalCounts[0] += sum(np.array(self.Exons[i].Counts[:Replicates], dtype=float))
            TotalCounts[1] += sum(np.array(self.Exons[i].Counts[Replicates:], dtype=float))
        for i in self.Exons:
            self.Cond1Count = sum(np.array(self.Exons[i].Counts[:Replicates], dtype=float))
            self.Cond2Count = sum(np.array(self.Exons[i].Counts[Replicates:], dtype=float))
            self.Exons[i].NormCond1Count = self.Cond1Count/TotalCounts[0]
            #self.ExonResult[1].append(NormCond1Count)
            self.Exons[i].NormCond2Count = self.Cond2Count/TotalCounts[1]
            self.Exons[i].FractionChange = round((self.Exons[i].NormCond2Count - self.Exons[i].NormCond1Count), 4)
            #self.ExonResult[2].append(NormCond2Count)
            #self.ExonResult[3].append(self.Exons[i].FractionChange)
            #self.ExonResult[0] = len(self.Exons) 

    def __str__(self):
        return  str([self.Name,
                self.baseMean,
                self.FoldChange,
                self.pvalue,
                self.padj,
                [i for i in self.Exons]])

##Class for Intron DE info, invoked unless IncIntron == True, then treated as exons
class Intron(object):
    def __init__(self, Name, baseMean, log2FoldChange, pvalue, padj, weight, counts):
        self.Name = Name
        self.Counts = counts
        self.PACs = {}
        self.GeneName = Name.split("_")[0]
        self.baseMean = round(float(baseMean), 4)
        self.Fraction = self.baseMean / GeneDict[self.GeneName].baseMean
        self.FractionChange = 0.0
        try:
            self.FoldChange = round(float(log2FoldChange), 4)
            self.pvalue = round(float(pvalue), 4)
            self.padj = round(float(padj), 4)
        except:
            self.FoldChange = log2FoldChange
            self.pvalue = pvalue
            self.padj = padj

    def FinishPACs(self):
        #Conditions = 2 ##Change this later or learn from metadata
        TotalCounts = [0] * Conditions# * Replicates
        for i in self.PACs:
            TotalCounts[0] += sum(np.array(self.PACs[i].Counts[:Replicates], dtype=float))
            TotalCounts[1] += sum(np.array(self.PACs[i].Counts[Replicates:], dtype=float))
        for i in self.PACs:
            Cond1Count = sum(np.array(self.PACs[i].Counts[:Replicates], dtype=float))
            Cond2Count = sum(np.array(self.PACs[i].Counts[Replicates:], dtype=float))
            NormCond1Count = Cond1Count/TotalCounts[0]
            NormCond2Count = Cond2Count/TotalCounts[1]
            self.PACs[i].FractionChange = round((NormCond2Count - NormCond1Count), 4)

    def __str__(self):
        return  str([self.Name,
                self.baseMean,
                self.FoldChange,
                self.pvalue,
                self.padj,
                self.Fraction,
                [i for i in self.PACs]])

##Class for Exon DE info
class Exon(object):
    def __init__(self, Name, baseMean, log2FoldChange, pvalue, padj, weight, counts):
        self.Name = Name
        self.Counts = counts
        self.Cond1Count = sum(np.array(self.Counts[:Replicates], dtype=float))
        self.Cond2Count = sum(np.array(self.Counts[Replicates:], dtype=float))
        self.PACs = {}
        self.GeneName = Name.split("_")[0]
        self.baseMean = round(float(baseMean), 4)
        self.Fraction = self.baseMean / GeneDict[self.GeneName].baseMean
        self.FractionChange = 0.0
        self.Cond1BED = []
        self.Cond2BED = []
        #self.PD1 = 0
        #self.PD2 = 0
        self.PDs = np.array([0.0] * Replicates * Conditions)
        self.PDUI = 'NA'
        try:
            self.FoldChange = round(float(log2FoldChange), 4)
            self.pvalue = round(float(pvalue), 4)
            self.padj = round(float(padj), 4)
        except:
            self.FoldChange = log2FoldChange
            self.pvalue = pvalue
            self.padj = padj
        if Links:
            try:
                Link = Exon_library[Name][0] + ":" + Exon_library[Name][1] + "-" + Exon_library[Name][2]
                self.Link = ''.join(['=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?org=',
                        org,
                        '&db=',
                        db,
                        '&position=',
                        Link,
                        '", "',
                        Name,
                        '")'])
            except:
                self.Link = Name
        else:
            self.Link = Name

    def FinishPACs(self):
        #Conditions = 2 ##Change this later or learn from metadata
        TotalCounts = [0] * Conditions# * Replicates
        TotalCountsArray = np.array([0.0] * Replicates * Conditions)
        RunningCounts = np.array([0.0] * Replicates * Conditions)
        for i in self.PACs:
            if self.PACs[i].Fraction > 0.05:
                TotalCounts[0] += sum(np.array(self.PACs[i].Counts[:Replicates], dtype=float))
                TotalCounts[1] += sum(np.array(self.PACs[i].Counts[Replicates:], dtype=float))
                TotalCountsArray += np.array(self.PACs[i].Counts, dtype=float)
        TEMP = []
        for i in self.PACs:
            if self.PACs[i].Fraction > 0.05:
                TEMP.append(self.PACs[i].Name)
        TEMP.sort(key=lambda x: int(x.split('-')[-1]))
        #Start = i.split(':')[-1].split("_")[0]
        try:
            if self.PACs[i].Locus[3] == '+':
                Start = self.PACs[TEMP[0]].Locus[1]
                #Start = i.split(':')[-1].split("_")[0]
                Stop = self.PACs[TEMP[-1]].Locus[2]
                Length = int(Stop) - int(Start)
                RunningCount1 = TotalCounts[0]
                RunningCount2 = TotalCounts[1]
                RunningCounts += TotalCountsArray
                LastCoord = Start
                for i in TEMP:
                    Cond1Count = sum(np.array(self.PACs[i].Counts[:Replicates], dtype=float))
                    Cond2Count = sum(np.array(self.PACs[i].Counts[Replicates:], dtype=float))
                    CurrentCounts = np.array(self.PACs[i].Counts, dtype=float)
                    ##Fraction Calc
                    NormCond1Count = Cond1Count/TotalCounts[0]
                    NormCond2Count = Cond2Count/TotalCounts[1]
                    self.PACs[i].FractionChange = round((NormCond2Count - NormCond1Count), 4)
                    ##Make annotation for bedgraph
                    BEDEntry = self.PACs[i].Locus[0] + '\t' + LastCoord + '\t' + self.PACs[i].Locus[2] + '\t'
                    self.Cond1BED.append(BEDEntry + str(round(RunningCount1,2)) + '\n')
                    self.Cond2BED.append(BEDEntry + str(round(RunningCount2,2)) + '\n')
                    ##PDUI Calc
                    SegmentFraction = (int(self.PACs[i].Locus[1]) - int(LastCoord))/ float(Length)
                    #self.PD1 += (SegmentFraction * (RunningCount1/float(TotalCounts[0])))
                    #self.PD2 += (SegmentFraction * (RunningCount2/float(TotalCounts[1])))
                    self.PDs += (SegmentFraction * (RunningCounts/TotalCountsArray))
#                    if 'TSHZ2' in self.Name:
#                        print(self.Name, self.PD1, self.PD2, self.PDs)
#                        print(RunningCount1, RunningCount2, RunningCounts)
#                        print(TotalCountsArray, TotalCounts)
#                        print(SegmentFraction, Length, Start, Stop, LastCoord)
                    RunningCount1 -= Cond1Count
                    RunningCount2 -= Cond2Count
                    RunningCounts -= CurrentCounts
                    LastCoord = self.PACs[i].Locus[2]
            else:
                TEMP = TEMP[::-1]
                Start = self.PACs[TEMP[0]].Locus[2]
                Stop = self.PACs[TEMP[-1]].Locus[1]     ##LH of PAC
                Length = int(Start) - int(Stop)         ##Max Distance
                RunningCount1 = TotalCounts[0]
                RunningCount2 = TotalCounts[1]
                RunningCounts += TotalCountsArray
                LastCoord = Start
                for i in TEMP:
                    Cond1Count = sum(np.array(self.PACs[i].Counts[:Replicates], dtype=float))
                    Cond2Count = sum(np.array(self.PACs[i].Counts[Replicates:], dtype=float))
                    CurrentCounts = np.array(self.PACs[i].Counts, dtype=float)
                    ##Fraction Calc
                    NormCond1Count = Cond1Count/TotalCounts[0]
                    NormCond2Count = Cond2Count/TotalCounts[1]
                    self.PACs[i].FractionChange = round((NormCond2Count - NormCond1Count), 4)
                    ##Make annotation for bedgraph
                    BEDEntry = self.PACs[i].Locus[0] + '\t' + self.PACs[i].Locus[1] + '\t' + LastCoord + '\t'
                    self.Cond1BED.append(BEDEntry + str(round(RunningCount1,2)) + '\n') 
                    self.Cond2BED.append(BEDEntry + str(round(RunningCount2,2)) + '\n')
                    ##PDUI Calc
                    SegmentFraction = (int(LastCoord) - int(self.PACs[i].Locus[1]))/ float(Length)
                    #self.PD1 += (SegmentFraction * (RunningCount1/float(TotalCounts[0])))
                    #self.PD2 += (SegmentFraction * (RunningCount2/float(TotalCounts[1])))
                    self.PDs += (SegmentFraction * (RunningCounts/TotalCountsArray))
#                    if 'DPYD' in self.Name:
#                        print(self.Name, self.PD1, self.PD2, self.PDs)
#                        print(RunningCount1, RunningCount2, RunningCounts)
#                        print(TotalCountsArray, TotalCounts)
#                        print(SegmentFraction, Length, Start, Stop, LastCoord)
                    RunningCount1 -= Cond1Count
                    RunningCount2 -= Cond2Count
                    RunningCounts -= CurrentCounts
                    LastCoord = self.PACs[i].Locus[1]
            #self.PDUI = round(np.nanmean(self.PDs[Replicates:]) - np.nanmean(self.PDs[:Replicates]), 4)
            #self.PDUI = self.PDs[Replicates:] - self.PDs[:Replicates]
            if any([float(k)<MinCount for k in TotalCountsArray]):
                self.PDUI = 'NA'
            else:
                self.PDUI = round((np.nansum(self.PDs[Replicates:])/float(Replicates)) - (np.nansum(self.PDs[:Replicates])/float(Replicates)), 4)
            #self.PDUI = round((self.PD1 - self.PD2), 4)
        except:
            print('Error: ', self.Name, " has error in Coord")

    def __str__(self):
        return  str([self.Name,
                self.baseMean,
                self.FoldChange,
                self.pvalue,
                self.padj,
                self.Fraction,
                [i for i in self.PACs]])

##Class for annotated PACs DE info
class PAC(object):
    def __init__(self, Name, baseMean, log2FoldChange, pvalue, padj, weight, counts):
        self.Name = Name
        self.Counts = counts
        self.baseMean = round(float(baseMean), 4)
        self.GeneName = Name.split("_")[0]
        self.ExonName = Name.split("_PAC")[0]
        self.Locus = PAC_library[Name]
        self.Direction = Gene_library[self.GeneName]
        self.Start = self.Locus[1]
        if 'exon' in Name or 'end' in Name or IncIntrons == True:
            self.Fraction = self.baseMean / GeneDict[self.GeneName].Exons[self.ExonName].baseMean
        else:
            self.Fraction = self.baseMean / GeneDict[self.GeneName].Introns[self.ExonName].baseMean
        self.FractionChange = 0.0
        try:
            self.FoldChange = round(float(log2FoldChange), 4)
            self.pvalue = round(float(pvalue), 4)
            self.padj = round(float(padj), 4)
        except:
            self.log2FoldChange = log2FoldChange
            self.pvalue = pvalue
            self.padj = padj

    def __str__(self):
        return  str([self.Name,
                self.baseMean,
                self.FoldChange,
                self.padj])


xs, ys, cs = [],[],[]

m=0
for i in GeneDict:
        if (int(GeneDict[i].ExonResult[0]) > 0):
            j = 0
            if GeneDict[i].PACResult[1][j] != 'NA' and GeneDict[i].PACResult[2][j] != 'NA':
                m+=1
                xs.append(GeneDict[i].PACResult[1][j])
                ys.append(GeneDict[i].PACResult[2][j])
                if GeneDict[i].PACResult[4][j] == 'Shortening':
                    cs.append('c')
                elif GeneDict[i].PACResult[4][j] == 'Lengthening':
                    cs.append('y')
                elif GeneDict[i].PACResult[4][j] == 'Both':
                    cs.append('magenta')
                else:
                    cs.append('k')
            j+=1
            while (j < int(GeneDict[i].ExonResult[0])):
                if GeneDict[i].PACResult[1][j] != 'NA' and GeneDict[i].PACResult[2][j] != 'NA':
                    m+=1
                    xs.append(GeneDict[i].PACResult[1][j])
                    ys.append(GeneDict[i].PACResult[2][j])
                    if GeneDict[i].PACResult[4][j] == 'Shortening':
                        cs.append('c')
                    elif GeneDict[i].PACResult[4][j] == 'Lengthening':
                        cs.append('y')
                    elif GeneDict[i].PACResult[4][j] == 'Both':
                        cs.append('magenta')
                    else:
                        cs.append('k')
                j+=1
