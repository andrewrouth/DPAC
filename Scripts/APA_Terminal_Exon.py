#!/bin/python3
import argparse
import numpy as np
#import pickle as pickle
np.seterr(divide='ignore', invalid='ignore')

parser = argparse.ArgumentParser()
parser.add_argument("Root", help="Root of Datafiles")
parser.add_argument("Output", help="Report Name")
parser.add_argument("PACs", help="Report Name")
parser.add_argument("--Replicates", help="Report Name")
parser.add_argument("--Conditions", help="Report Name")
parser.add_argument("--Fraction", help="Report Name")
parser.add_argument("--MinCount", help="Report Name")
parser.add_argument("-Introns", action='store_true', help="Include Intronic PACs")
parser.add_argument("-Links", action='store_true', help="Make Exon Links to USCS browsers")
parser.add_argument("--Annos", help="If Links: enter Annotation file name; e.g.: mm10_Annos_23Dec19.bed")
parser.add_argument("--db", help="If Links: enter database; e.g.: mouse,mm10")
parser.add_argument("-NEFormat", action='store_true', help="Instead of carriage return use in excel output, write a numbered line per exon")
parser.add_argument("-DRIMSeq", action='store_true', help="Toggle whether to include DRIMSeq")
args = parser.parse_args()

Root = str(args.Root)
##DESeq2 Results
Input_Gene = Root + '_Gene_DESeq2-results.csv'
Input_Exon = Root + '_Exon_DESeq2-results.csv'
Input_PAC = Root + '_PAC_DESeq2-results.csv'

if args.Introns:
    IncIntrons = True
else:
    IncIntrons = False
if args.Replicates:
    Replicates = int(args.Replicates)
else:
    Replicates = 3
if args.Conditions:
    Conditions = int(args.Conditions)
else:
    Conditions = 2
if args.Fraction:
    ReqFractionChange = float(args.Fraction)
else:
    ReqFractionChange = 0.1
if args.MinCount:
    MinCount = int(args.MinCount)
else:
    MinCount = 10
if args.Links:
    org, db = str(args.db).split(',')
    Links = True
else:
    Links = False
if args.NEFormat:
    NEFormat = True
else:
    NEFormat = False
if args.DRIMSeq:
    DRIMSeq = True
    ##DRIMseq Results
    DRIMproportions_Exon = Root + '_Exon_DRIMseq-proportions.csv'
    DRIMproportions_PAC = Root + '_PAC_DRIMseq-proportions.csv'
    DRIMresults_Exon = Root + '_Exon_DRIMseq-results.csv'
    DRIMresults_PAC = Root + '_PAC_DRIMseq-results.csv'
else:
    DRIMSeq = False

##Read in PAC BED file into dict
Gene_library = {}
PAC_library = {}
with open(str(args.PACs),'r') as In:
    line = In.readline()
    while line:
        data = line.split()
        PAC_library[data[3]] = (data[0], data[1], data[2], data[5])
        Gene = data[3].split("_")[0]
        Gene_library[Gene] = data[5]
        line = In.readline()
    
if Links:
    Exon_library = {}
    with open(str(args.Annos),'r') as In:
        line = In.readline()
        while line:
            data = line.split()
            if data[5] == '+':
                Name = "_".join(data[3].split("_")[:2]) + "_" + data[0] + ':' + data[1]
                Exon_library[Name] = [data[0], data[1], data[2]]
            else:
                Name = "_".join(data[3].split("_")[:2]) + "_" + data[0] + ':' + data[2]
                Exon_library[Name] = [data[0], data[2], data[1]]
            line = In.readline()    
    
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
        
        if args.DRIMSeq:
            self.DRIMGeneResult = ['','No','No']
            self.DRIMExonResult = [0, [], [], [], [], [], []]
            self.DRIMPACResult = [[], [], [], [], [], []]
            
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
            self.Exons[i].NormCond2Count = self.Cond2Count/TotalCounts[1]
            self.Exons[i].FractionChange = round((self.Exons[i].NormCond2Count - self.Exons[i].NormCond1Count), 4)

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
        if args.DRIMSeq:
            self.DRIMpadj = 'NA'
            self.DRIMproportionCond1 = ''
            self.DRIMproportionCond2 = ''
            self.DRIMFractionChange = 0.0
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
                    self.PDs += (SegmentFraction * (RunningCounts/TotalCountsArray))
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
                    self.PDs += (SegmentFraction * (RunningCounts/TotalCountsArray))
                    RunningCount1 -= Cond1Count
                    RunningCount2 -= Cond2Count
                    RunningCounts -= CurrentCounts
                    LastCoord = self.PACs[i].Locus[1]
            if any([float(k)<MinCount for k in TotalCountsArray]):
                self.PDUI = 'NA'
            else:
                self.PDUI = round((np.nansum(self.PDs[Replicates:])/float(Replicates)) - (np.nansum(self.PDs[:Replicates])/float(Replicates)), 4)
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
        if args.DRIMSeq:
            self.DRIMpadj = 'NA'
            self.DRIMproportionCond1 = ''
            self.DRIMproportionCond2 = ''
            self.DRIMFractionChange = 0.0

    def __str__(self):
        return  str([self.Name,
                self.baseMean,
                self.FoldChange,
                self.padj])

GeneDict = {}
##Reads in all Gene DE info
with open(Input_Gene, 'r') as In:
        Header = In.readline()
        Headings = Header.rstrip()[1:-1].split('","')[1:]
        Data = In.readline().rstrip()
        while Data:
            Data = Data.split(',')
            Name = Data[1].replace('"', '')
            GeneDict[Name] = Gene(Name, Data[2], Data[3], Data[6], Data[7], Data[8], Data[9:])
            Data = In.readline().rstrip()

##Reads in all Exon DE info
with open(Input_Exon, 'r') as In:
        Header = In.readline()
        Headings = Header.rstrip()[1:-1].split('","')[1:]
        Data = In.readline().rstrip()
        while Data:
            Data = Data.split(',')
            Name = Data[1].replace('"', '')
            if 'Unannotated' in Name:
                pass
            else:
                GeneName = Name.split('_')[0]
                if 'exon' in Name or 'end' in Name or IncIntrons == True:
                    GeneDict[GeneName].Exons[Name] = Exon(Name, Data[2], Data[3], Data[6], Data[7], Data[8], Data[9:])
                else:
                    GeneDict[GeneName].Introns[Name] = Intron(Name, Data[2], Data[3], Data[6], Data[7], Data[8], Data[9:])
            Data = In.readline().rstrip()

##Reads in all PAC DE info
UnannotatedDict = {}
with open(Input_PAC, 'r') as In:
        Header = In.readline()
        Headings = Header.rstrip()[1:-1].split('","')[1:]
        Data = In.readline().rstrip()
        while Data:
            Data = Data.split(',')
            Name = Data[1].replace('"', '')
            if 'Unannotated' in Name:
                UnannotatedDict[Name] = Unannotated(Name, Data[2], Data[3], Data[6], Data[7], Data[8], Data[9:])
            else:
                GeneName = Name.split('_')[0]
                if 'exon' in Name or 'end' in Name or IncIntrons == True:
                    ExonName = Data[1].replace('"', '').split('_PAC')[0]
                    GeneDict[GeneName].Exons[ExonName].PACs[Name] = PAC(Name, Data[2], Data[3], Data[6], Data[7], Data[8], Data[9:])
                else:
                    IntronName = Data[1].replace('"', '').split('_PAC')[0]
                    GeneDict[GeneName].Introns[IntronName].PACs[Name] = PAC(Name, Data[2], Data[3], Data[6], Data[7], Data[8], Data[9:])
            Data = In.readline().rstrip()


##DRIMSeq Data
if args.DRIMSeq:
    with open(DRIMproportions_Exon, 'r') as In:
            Header = In.readline()
            Data = In.readline().rstrip()
            while Data:
                Data = Data.split(',')
                Gene = Data[1].replace('"', '')
                Exon = Data[2].replace('"', '')
                try:
                    ProportionCond1 = float(Data[3])
                    ProportionCond2 = float(Data[-1])
                    Change = ProportionCond1 - ProportionCond2
                except:
                    ProportionCond1 = 'NA'
                    ProportionCond2 = 'NA'
                    Change = 0.0
                GeneDict[Gene].Exons[Exon].DRIMproportionCond1 = ProportionCond1
                GeneDict[Gene].Exons[Exon].DRIMproportionCond2 = ProportionCond2
                GeneDict[Gene].Exons[Exon].DRIMFractionChange = Change
                Data = In.readline().rstrip()
    ##DRIMSeq Data        
    with open(DRIMresults_Exon, 'r') as In:
            Header = In.readline()
            Data = In.readline().rstrip()
            while Data:
                Data = Data.split(',')
                Gene = Data[1].replace('"', '')
                Exon = Data[2].replace('"', '')
                try:
                    DRIMpadj = float(Data[-1])
                except:
                    DRIMpadj = 'NA'
                GeneDict[Gene].Exons[Exon].DRIMpadj = DRIMpadj
                Data = In.readline().rstrip()
    ##DRIMSeq Data
    with open(DRIMproportions_PAC, 'r') as In:
            Header = In.readline()
            Data = In.readline().rstrip()
            while Data:
                Data = Data.split(',')
                Gene = Data[1].split("_")[0].replace('"', '')
                Exon = Data[1].split("_PAC")[0].replace('"', '')
                PAC = Data[2].replace('"', '')
                try:
                    ProportionCond1 = float(Data[3])
                    ProportionCond2 = float(Data[-1])
                    Change = ProportionCond1 - ProportionCond2
                except:
                    ProportionCond1 = 'NA'
                    ProportionCond2 = 'NA'
                    Change = 0.0
                GeneDict[Gene].Exons[Exon].PACs[PAC].DRIMproportionCond1 = ProportionCond1
                GeneDict[Gene].Exons[Exon].PACs[PAC].DRIMproportionCond2 = ProportionCond2
                GeneDict[Gene].Exons[Exon].PACs[PAC].DRIMFractionChange = Change
                Data = In.readline().rstrip()
    ##DRIMSeq Data
    with open(DRIMresults_PAC, 'r') as In:
            Header = In.readline()
            Data = In.readline().rstrip()
            while Data:
                Data = Data.split(',')
                Gene = Data[1].split("_")[0].replace('"', '')
                Exon = Data[1].split("_PAC")[0].replace('"', '')
                PAC = Data[2].replace('"', '')
                try:
                    DRIMpadj = float(Data[-1])
                except:
                    DRIMpadj = 'NA'            
                GeneDict[Gene].Exons[Exon].PACs[PAC].DRIMpadj = DRIMpadj
                Data = In.readline().rstrip()

#After reading in, calculate adbundance of PACs within genes and within exons/introns
for i in GeneDict:
    GeneDict[i].FinishExons()
    for j in GeneDict[i].Exons:
        GeneDict[i].Exons[j].FinishPACs()
    for j in GeneDict[i].Introns:
        GeneDict[i].Introns[j].FinishPACs()
#pickle.dump(GeneDict, open(Root + '.GeneDict.py.pi','wb'))

def FindSizeChange(Results):
    #Start with Doubles.
    global APAC, APAB, APAL
    if len(Results) == 2:  
        if Results == ['UP', 'nc'] or Results == ['UP', 'DOWN'] or Results == ['nc', 'DOWN']:
            Result = '_Shortening'
            APAC += 1
        elif Results == ['nc', 'UP'] or Results == ['DOWN', 'UP'] or Results == ['DOWN', 'nc']:
            Result = '_Lengthening'
            APAL += 1   
        return Result
    elif len(Results) > 2:
        ##Single changer:
        if Results.count('DOWN') + Results.count('UP') == 1:
            if Results[0] == 'UP' or Results[-1] == 'DOWN':
                Result = '_Shortening'
                APAC += 1
            elif Results[0] == 'DOWN' or Results[-1] == 'UP':
                Result = '_Lengthening'
                APAL += 1
            else:
                Result = '_Both'
                APAB += 1
        else:
            ##MultiChange
            ##Changes must be in order, else are 'both'
            ##Allowed combinations are:
                #Cond1 [UP [ncx] DOWN]
                #Cond2 [DOWN [ncx] UP]
            ##Collapse adjacent results
            x = Results
            y = ['']
            for i in x:
                if i != y[-1]:
                    y.append(i)
            y = y[1:]
            if y == ['UP', 'nc', 'DOWN']:
                Result = '_Shortening'
                APAC += 1
            elif y == ['DOWN', 'nc', 'UP']:
                Result = '_Lengthening'
                APAL += 1
            else:
                Result = '_Both'
                APAB += 1
    return Result

##MAKE REPORT
DEn, TEn, APAn, APATEn, APAL, APAC, APAB = 0,0,0,0,0,0,0

with open(str(args.Output),'w') as Output:
    Heading = ['Gene', 'baseMean', 'FoldChange', 'padj', 'GeneResult',
               'Exons', 'baseMean', 'FractionChange', 'padj', 'ExonResult', 'TE_Result',
               'PACs', 'baseMean', 'FractionChange', 'padj', 'PACResult', 'APA_Result']
    Output.write('\t'.join(Heading))
    Output.write ('\n')
    for i in GeneDict:
        DE = False
        TE = False
        APA = False
        ##GENE RESULTS
        if GeneDict[i].padj != 'NA' and GeneDict[i].padj < 0.1:
            if GeneDict[i].FoldChange > 0.585:
                Result = 'UP'
                DE = True
                DEn += 1
            elif GeneDict[i].FoldChange < -0.585:
                Result = 'DOWN'
                DE = True
                DEn += 1
            else:
                Result = 'NC'
        else:
            Result = 'NC'
        GeneEntry = [i, str(GeneDict[i].baseMean), str(GeneDict[i].FoldChange), str(GeneDict[i].padj), Result, '']
        Output.write('\t'.join(GeneEntry))
        GeneDict[i].GeneResult[0] = Result
        Direction = Gene_library[i]
        ##EXONS RESULTS
        TEMP = []
        for j in GeneDict[i].Exons:
            if GeneDict[i].Exons[j].Fraction > 0.05:
                TEMP.append(GeneDict[i].Exons[j])
        ##Make sure Exons are sequential
        TEMP.sort(key=lambda x: x.Name)
        NumSites = len(TEMP)
        if NumSites == 1:
            GeneDict[i].ExonResult[4].append('NA')
            GeneDict[i].ExonResult[0] += 1
            GeneDict[i].ExonResult[1].append(1.0)
            GeneDict[i].ExonResult[2].append(1.0)
            GeneDict[i].ExonResult[3].append(0.0)
            GeneDict[i].ExonResult[5].append(TEMP[0].Name)
            GeneDict[i].ExonResult[6].append(TEMP[0].Link)     
        else:
            Results = []
            Hits = 0
            ExonEntry = [ [], [],[],[], 'NC', 'NA', '']
            n=0
            for j in TEMP:
                n+=1
                ExonEntry[0].append(j.Name)
                ExonEntry[1].append(j.baseMean)
                ExonEntry[2].append(j.FractionChange)
                ExonEntry[3].append(j.padj)
                GeneDict[i].ExonResult[0] += 1
                GeneDict[i].ExonResult[1].append(round(j.NormCond1Count, 4))
                GeneDict[i].ExonResult[2].append(round(j.NormCond2Count, 4))
                GeneDict[i].ExonResult[3].append(round(j.FractionChange, 4))
                GeneDict[i].ExonResult[5].append(j.Name)
                GeneDict[i].ExonResult[6].append(j.Link)
                if j.padj != 'NA' and j.padj < 0.1:
                    if j.FractionChange > ReqFractionChange:
                        Results.append('UP')
                        GeneDict[i].ExonResult[4].append('Up')
                        Hits += 1
                    elif j.FractionChange < -ReqFractionChange:
                        Results.append('DOWN')
                        GeneDict[i].ExonResult[4].append('Down')
                        Hits += 1
                    else:
                        Results.append('nc')
                        GeneDict[i].ExonResult[4].append('NA')
                else:
                    Results.append('nc')
                    GeneDict[i].ExonResult[4].append('NA')
            if Hits == 0:
                ExonEntry[-2] = 'NotSig'           
            elif Hits == 1:
                ##Single Changer
                if GeneDict[i].Cond1Count >= MinCount and GeneDict[i].Cond2Count >= MinCount:
                    ExonEntry[-2] = 'TE'
                    n=0
                    for j in Results:
                        n+=1
                        if 'UP' in j or 'DOWN' in j:
                            ExonEntry[-3] = 'Exon#' + str(n) + j
                else:
                    ExonEntry[-2] = 'NotSig'
            else:
                if GeneDict[i].Cond1Count >= MinCount and GeneDict[i].Cond2Count >= MinCount:
                    ExonEntry[-2] = 'TE'
                    ExonEntry[-3] = ''
                    n=0
                    for j in Results:
                        n+=1
                        if 'UP' in j or 'DOWN' in j:
                            String = 'Exon#' + str(n) + j
                            ExonEntry[-3] += String
                else:
                    ExonEntry[-2] = 'NotSig'
            if ExonEntry[-2] == 'TE':
                TE = True
                GeneDict[i].GeneResult[1] = 'Yes'
                TEn += 1
            Output.write('\t'.join([str(k) for k in ExonEntry]))
        ##PAC RESULTS
        for j in GeneDict[i].Exons:
            if GeneDict[i].Exons[j].Fraction > 0.05:
                TEMP = []
                for k in GeneDict[i].Exons[j].PACs:
                    if GeneDict[i].Exons[j].PACs[k].Fraction > 0.05:
                        TEMP.append(GeneDict[i].Exons[j].PACs[k])
                ##Make sure PACs are sequential to determine lengthening vs shortening
                TEMP.sort(key=lambda x: int(x.Name.split('-')[-1]))
                NumSites = len(TEMP)
                GeneDict[i].PACResult[5].append(GeneDict[i].Exons[j].Name)
                if NumSites == 1:
                    GeneDict[i].PACResult[0].append(1)      #Number of PACs
                    GeneDict[i].PACResult[1].append('NA')    #PD1
                    GeneDict[i].PACResult[2].append('NA')    #PD2
                    GeneDict[i].PACResult[3].append('NA')    #PDUI
                    GeneDict[i].PACResult[4].append('NA')   #DPAC Results
                    #GeneDict[i].PACResult[5].append(GeneDict[i].Exons[j].Name)   #Name
                else:
                    GeneDict[i].PACResult[0].append(len(GeneDict[i].Exons[j].PACs))
                    if any([float(k)<MinCount for k in GeneDict[i].Exons[j].Counts]):
                        GeneDict[i].PACResult[1].append('NA')
                        GeneDict[i].PACResult[2].append('NA')
                    else:
                        GeneDict[i].PACResult[1].append(round(np.nansum(GeneDict[i].Exons[j].PDs[:Replicates])/float(Replicates), 4))
                        GeneDict[i].PACResult[2].append(round(np.nansum(GeneDict[i].Exons[j].PDs[Replicates:])/float(Replicates), 4))
                    GeneDict[i].PACResult[3].append(GeneDict[i].Exons[j].PDUI)
                    Results = []
                    Hits = 0
                    PACEntry = [ [], [],[],[], 'NC', 'NA', '']
                    n=0
                    for k in TEMP:
                        n+=1
                        PACEntry[0].append(k.Name)
                        PACEntry[1].append(k.baseMean)
                        PACEntry[2].append(k.FractionChange)
                        PACEntry[3].append(k.padj)
                        #GeneDict[i].PACResult[5].append(k.Name)
                        if k.padj != 'NA' and k.padj < 0.1:
                            if k.FractionChange > ReqFractionChange:
                                Results.append('UP')
                                Hits += 1
                            elif k.FractionChange < -ReqFractionChange:
                                Results.append('DOWN')
                                Hits += 1
                            else:
                                Results.append('nc')
                        else:
                            Results.append('nc')
                    if Hits == 0:
                        PACEntry[-2] = 'NotSig'
                    elif Hits == 1:
                        ##Single Changer
                        if GeneDict[i].Exons[j].Cond1Count >= MinCount and GeneDict[i].Exons[j].Cond2Count >= MinCount:
                            PACEntry[-2] = 'APA'
                            n=0
                            for k in Results:
                                n+=1
                                if 'UP' in k or 'DOWN' in k:
                                    PACEntry[-3] = 'PAC#' + str(n) + '-' + k
                        else:
                            PACEntry[-2] = 'NotSig'
                    else:
                        if GeneDict[i].Exons[j].Cond1Count >= MinCount and GeneDict[i].Exons[j].Cond2Count >= MinCount:
                            PACEntry[-2] = 'APA'
                            PACEntry[-3] = ''
                            n=0
                            for k in Results:
                                n+=1
                                if 'UP' in k or 'DOWN' in k:
                                    String = 'PAC#' + str(n) + '-' + k
                                    PACEntry[-3] += String
                        else:
                            PACEntry[-2] = 'NotSig'
                    if PACEntry[-2][:3] == 'APA':
                        APAn += 1
                        APA = True
                        GeneDict[i].GeneResult[2] = 'Yes'
                        ##Find Shortening vs Lengthening vs both
                        if Direction == '-':
                            Results = Results[::-1]
                        else:
                            pass
                        PACEntry[-2] += FindSizeChange(Results)
                        GeneDict[i].PACResult[4].append(PACEntry[-2][4:])
                        padj = []
                        for l in PACEntry[-4]:
                            if l != 'NA':
                                padj.append(max(l, 0.00001))
                            else:
                                padj.append(1.0)
                    else:
                        GeneDict[i].PACResult[4].append('NA')
                    Output.write('\t'.join([str(k) for k in PACEntry]))
        Output.write('\n')
        if TE and APA:
            APATEn += 1                    

##Process using DRIM padj
if args.DRIMSeq:
    for i in GeneDict:
        TE = False
        APA = False
        Direction = Gene_library[i]
        ##EXONS RESULTS
        TEMP = []
        for j in GeneDict[i].Exons:
            if GeneDict[i].Exons[j].Fraction > 0.05:
                TEMP.append(GeneDict[i].Exons[j])
        ##Make sure Exons are sequential
        TEMP.sort(key=lambda x: x.Name)
        NumSites = len(TEMP)
        if NumSites == 1:
            GeneDict[i].DRIMExonResult[4].append('NA')
            GeneDict[i].DRIMExonResult[0] += 1
            GeneDict[i].DRIMExonResult[1].append(1.0)
            GeneDict[i].DRIMExonResult[2].append(1.0)
            GeneDict[i].DRIMExonResult[3].append(0.0)
            GeneDict[i].DRIMExonResult[5].append(TEMP[0].Name)
            GeneDict[i].DRIMExonResult[6].append(TEMP[0].Link)     
        else:
            Results = []
            Hits = 0
            ExonEntry = [ [], [],[],[], 'NC', 'NA', '']
            n=0
            for j in TEMP:
                n+=1
                ExonEntry[0].append(j.Name)
                ExonEntry[1].append(j.baseMean)
                ExonEntry[2].append(j.DRIMFractionChange)
                ExonEntry[3].append(j.DRIMpadj)
                GeneDict[i].DRIMExonResult[0] += 1
                GeneDict[i].DRIMExonResult[1].append(j.NormCond1Count)
                GeneDict[i].DRIMExonResult[2].append(j.NormCond2Count)
                GeneDict[i].DRIMExonResult[3].append(j.FractionChange)
                GeneDict[i].DRIMExonResult[5].append(j.Name)
                GeneDict[i].DRIMExonResult[6].append(j.Link)
                if j.DRIMpadj != 'NA' and j.DRIMpadj < 0.1:
                    if j.FractionChange > ReqFractionChange:
                        Results.append('UP')
                        GeneDict[i].DRIMExonResult[4].append('Up')
                        Hits += 1
                    elif j.FractionChange < -ReqFractionChange:
                        Results.append('DOWN')
                        GeneDict[i].DRIMExonResult[4].append('Down')
                        Hits += 1
                    else:
                        Results.append('nc')
                        GeneDict[i].DRIMExonResult[4].append('NA')
                else:
                    Results.append('nc')
                    GeneDict[i].DRIMExonResult[4].append('NA')
            if Hits == 0:
                ExonEntry[-2] = 'NotSig'           
            elif Hits == 1:
                ##Single Changer
                if GeneDict[i].Cond1Count >= MinCount and GeneDict[i].Cond2Count >= MinCount:
                    ExonEntry[-2] = 'TE'
                    n=0
                    for j in Results:
                        n+=1
                        if 'UP' in j or 'DOWN' in j:
                            ExonEntry[-3] = 'Exon#' + str(n) + j
                else:
                    ExonEntry[-2] = 'NotSig'
            else:
                if GeneDict[i].Cond1Count >= MinCount and GeneDict[i].Cond2Count >= MinCount:
                    ExonEntry[-2] = 'TE'
                    ExonEntry[-3] = ''
                    n=0
                    for j in Results:
                        n+=1
                        if 'UP' in j or 'DOWN' in j:
                            String = 'Exon#' + str(n) + j
                            ExonEntry[-3] += String
                else:
                    ExonEntry[-2] = 'NotSig'
            if ExonEntry[-2] == 'TE':
                TE = True
                GeneDict[i].DRIMGeneResult[1] = 'Yes'
        ##PAC RESULTS
        for j in GeneDict[i].Exons:
            if GeneDict[i].Exons[j].Fraction > 0.05:
                TEMP = []
                for k in GeneDict[i].Exons[j].PACs:
                    if GeneDict[i].Exons[j].PACs[k].Fraction > 0.05:
                        TEMP.append(GeneDict[i].Exons[j].PACs[k])
                ##Make sure PACs are sequential to determine lengthening vs shortening
                TEMP.sort(key=lambda x: int(x.Name.split('-')[-1]))
                NumSites = len(TEMP)
                GeneDict[i].DRIMPACResult[5].append(GeneDict[i].Exons[j].Name)
                if NumSites == 1:
                    GeneDict[i].DRIMPACResult[0].append(1)      #Number of PACs
                    GeneDict[i].DRIMPACResult[1].append('NA')    #PD1
                    GeneDict[i].DRIMPACResult[2].append('NA')    #PD2
                    GeneDict[i].DRIMPACResult[3].append('NA')    #PDUI
                    GeneDict[i].DRIMPACResult[4].append('NA')   #DPAC Results
                    #GeneDict[i].DRIMPACResult[5].append(TEMP[0].Name)   #Name
                else:
                    GeneDict[i].DRIMPACResult[0].append(len(GeneDict[i].Exons[j].PACs))
                    if any([float(k)<MinCount for k in GeneDict[i].Exons[j].Counts]):
                        GeneDict[i].DRIMPACResult[1].append('NA')
                        GeneDict[i].DRIMPACResult[2].append('NA')
                    else:
                        GeneDict[i].DRIMPACResult[1].append(round(np.nansum(GeneDict[i].Exons[j].PDs[:Replicates])/float(Replicates), 4))
                        GeneDict[i].DRIMPACResult[2].append(round(np.nansum(GeneDict[i].Exons[j].PDs[Replicates:])/float(Replicates), 4))
                    GeneDict[i].DRIMPACResult[3].append(GeneDict[i].Exons[j].PDUI)
                    Results = []
                    Hits = 0
                    PACEntry = [ [], [],[],[], 'NC', 'NA', '']
                    n=0
                    for k in TEMP:
                        n+=1
                        PACEntry[0].append(k.Name)
                        PACEntry[1].append(k.baseMean)
                        PACEntry[2].append(k.FractionChange)
                        PACEntry[3].append(k.DRIMpadj)
                        #GeneDict[i].DRIMPACResult[5].append(k.Name)
                        if k.DRIMpadj != 'NA' and k.DRIMpadj < 0.1:
                            if k.FractionChange > ReqFractionChange:
                                Results.append('UP')
                                Hits += 1
                            elif k.FractionChange < -ReqFractionChange:
                                Results.append('DOWN')
                                Hits += 1
                            else:
                                Results.append('nc')
                        else:
                            Results.append('nc')
                    if Hits == 0:
                        PACEntry[-2] = 'NotSig'
                    elif Hits == 1:
                        ##Single Changer
                        if GeneDict[i].Exons[j].Cond1Count >= MinCount and GeneDict[i].Exons[j].Cond2Count >= MinCount:
                            PACEntry[-2] = 'APA'
                            n=0
                            for k in Results:
                                n+=1
                                if 'UP' in k or 'DOWN' in k:
                                    PACEntry[-3] = 'PAC#' + str(n) + '-' + k
                        else:
                            PACEntry[-2] = 'NotSig'
                    else:
                        if GeneDict[i].Exons[j].Cond1Count >= MinCount and GeneDict[i].Exons[j].Cond2Count >= MinCount:
                            PACEntry[-2] = 'APA'
                            PACEntry[-3] = ''
                            n=0
                            for k in Results:
                                n+=1
                                if 'UP' in k or 'DOWN' in k:
                                    String = 'PAC#' + str(n) + '-' + k
                                    PACEntry[-3] += String
                        else:
                            PACEntry[-2] = 'NotSig'
                    if PACEntry[-2][:3] == 'APA':
                        APA = True
                        GeneDict[i].DRIMGeneResult[2] = 'Yes'
                        ##Find Shortening vs Lengthening vs both
                        if Direction == '-':
                            Results = Results[::-1]
                        else:
                            pass
                        PACEntry[-2] += FindSizeChange(Results)
                        GeneDict[i].DRIMPACResult[4].append(PACEntry[-2][4:])
                    else:
                        GeneDict[i].DRIMPACResult[4].append('NA')

##Make BEDGraphs
Out1 = open(Root + 'Cond1.bedgraph.txt', 'w')
Out2 = open(Root + 'Cond2.bedgraph.txt', 'w')
            
Out1.write('track type=bedGraph\tname="Cond1"\tdescription="S2"\talwaysZero=on visibility=full\twindowingFunction=maximum\tcolor=255,0,0\taltColor=0,100,200\n')
Out2.write('track type=bedGraph\tname="Cond2"\tdescription="S2"\talwaysZero=on visibility=full\twindowingFunction=maximum\tcolor=255,155,0\taltColor=0,100,200\n')

for Gene in GeneDict:
    for Exon in GeneDict[Gene].Exons:
        Out1.write(''.join(GeneDict[Gene].Exons[Exon].Cond1BED))
        Out2.write(''.join(GeneDict[Gene].Exons[Exon].Cond2BED))
Out1.close()
Out2.close()

print("Total genes detected: %s" % len(GeneDict))
n = 0
for i in GeneDict:
    if len(GeneDict[i].Exons) > 1:
        n+=1
print("Number of genes with multiple terminal exons:%s" % n)
n=0
for i in GeneDict:
    for j in GeneDict[i].Exons:
        if len(GeneDict[i].Exons[j].PACs) > 1:
            n+=1
print("Number of exons with multiple poly(A) sites: %s" % n)
print("Number of differentially expressed genes (padj<0.1): %s" % DEn)
print("Number of APA events: %s" % APAn)
print("\t of which %s are Lengthening, %s are Shortening: and %s are both." % (APAL, APAC, APAB))
print("Number of alternative terminal exons (TE): %s" % TEn)
print("Number genes with both APA and TE: %s" % APATEn)

##Make New Report
Out = open(Root + '_DPAC2_Report.tsv', 'w')
Header = ['Gene Name',
          'Mean Counts Cond1',
          'Mean Counts Cond2',
          'Fold Change', 
          'pAdj',
          'DE Result',
          'Any Splicing APA? (DPAC)',
          'Any Tandem APA? (DPAC)',          
          'Any Splicing APA? (DRIM)',
          'Any Tandem APA? (DRIM)',
          '#Exons with a PAC',
          'Exon Names',
          'Exon fraction  Cond1',
          'Exon fraction  Cond2',
          'Exon Fraction Change',
          '"Exon Results\r\n(NA = Not significant)"',
          '#PACs per Exon',
          'PDU per Exon Cond1',
          'PDU per Exon Cond2',
          'delta-PDUI per Exon',
          'DPAC Result',
          'DRIM Result',
          'Links to USCS Gene Browser (one link per cell):']
if args.DRIMSeq:
    pass
else:
    Header  = Header[:8] + Header[10:21] + Header[22:]

Header2 = ['gene_id', 'feature_id', 'Cond1_1', 'Cond1_2', 'Cond1_3', 'Cond2_1', 'Cond2_2', 'Cond2_3']

if NEFormat == True:
    Out.write('\t' + '\t'.join(Header) + '\n')
else:
    Out.write('\t'.join(Header) + '\n')

m=0
for i in GeneDict:
    if NEFormat == True:
        Refresh = 1
        Out.write('\t'.join([str(Refresh),
                            GeneDict[i].Name, 
                            str(GeneDict[i].MeanCountsCond1), 
                            str(GeneDict[i].MeanCountsCond2),
                            str(GeneDict[i].FoldChange),
                            str(GeneDict[i].padj),
                            GeneDict[i].GeneResult[0], #DE Result
                            GeneDict[i].GeneResult[1], #Splicing
                            GeneDict[i].GeneResult[2], #Tandem                            
                            GeneDict[i].DRIMGeneResult[1], #Splicing
                            GeneDict[i].DRIMGeneResult[2], #Tandem
                            str(GeneDict[i].ExonResult[0]),
                            '']))
        if (int(GeneDict[i].ExonResult[0]) > 0):
            j = 0
            Out.write('\t'.join([str(GeneDict[i].ExonResult[5][j]),
                                str(GeneDict[i].ExonResult[1][j]),
                                str(GeneDict[i].ExonResult[2][j]),
                                str(GeneDict[i].ExonResult[3][j]),
                                (GeneDict[i].ExonResult[4][j]),
                                str(GeneDict[i].PACResult[0][j]), #Number of PACs
                                str(GeneDict[i].PACResult[1][j]), #'PDU per Exon Cond1',
                                str(GeneDict[i].PACResult[2][j]), #'PDU per Exon Cond2',
                                str(GeneDict[i].PACResult[3][j]), #'PDUI,
                                str(GeneDict[i].PACResult[4][j]), #DPAC Results
                                str(GeneDict[i].DRIMPACResult[4][j]), #DRIM Results
                                str(GeneDict[i].ExonResult[6][j]),
                                '\n']))
            j+=1
            Refresh+=1
            while (j < int(GeneDict[i].ExonResult[0])):
                spacer=[str(Refresh),GeneDict[i].Name,'','','','','','','','','','','']
                Out.write('\t'.join(spacer))
                Out.write('\t'.join([str(GeneDict[i].ExonResult[5][j]),
                                    str(GeneDict[i].ExonResult[1][j]),
                                    str(GeneDict[i].ExonResult[2][j]),
                                    str(GeneDict[i].ExonResult[3][j]),
                                    (GeneDict[i].ExonResult[4][j]),
                                    str(GeneDict[i].PACResult[0][j]),
                                    str(GeneDict[i].PACResult[1][j]),
                                    str(GeneDict[i].PACResult[2][j]),
                                    str(GeneDict[i].PACResult[3][j]),
                                    str(GeneDict[i].PACResult[4][j]),
                                    str(GeneDict[i].DRIMPACResult[4][j]),
                                    str(GeneDict[i].ExonResult[6][j]),
                                    '\n']))
                j+=1
                Refresh+=1
        else:
            Out.write('\n')
    else:
        if Links:
            AnyLinks = '\t'.join([str(n) for n in GeneDict[i].ExonResult[6]])
        else:
            AnyLinks = ''
        if args.DRIMSeq:
            Out.write('\t'.join([GeneDict[i].Name, 
                        str(GeneDict[i].MeanCountsCond1), 
                        str(GeneDict[i].MeanCountsCond2),
                        str(GeneDict[i].FoldChange),
                        str(GeneDict[i].padj),
                        GeneDict[i].GeneResult[0],
                        GeneDict[i].GeneResult[1],
                        GeneDict[i].GeneResult[2],                        
                        GeneDict[i].DRIMGeneResult[1],
                        GeneDict[i].DRIMGeneResult[2],
                        str(GeneDict[i].ExonResult[0]),
                        '"' + '\r\n'.join(GeneDict[i].ExonResult[5]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].ExonResult[1]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].ExonResult[2]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].ExonResult[3]]) + '"',
                        '"' + '\r\n'.join(GeneDict[i].ExonResult[4]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[0]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[1]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[2]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[3]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[4]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].DRIMPACResult[4]]) + '"',
                        AnyLinks,
                        '\n']))
        else:
            Out.write('\t'.join([GeneDict[i].Name, 
                        str(GeneDict[i].MeanCountsCond1), 
                        str(GeneDict[i].MeanCountsCond2),
                        str(GeneDict[i].FoldChange),
                        str(GeneDict[i].padj),
                        GeneDict[i].GeneResult[0],
                        GeneDict[i].GeneResult[1],
                        GeneDict[i].GeneResult[2],                        
                        str(GeneDict[i].ExonResult[0]),
                        '"' + '\r\n'.join(GeneDict[i].ExonResult[5]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].ExonResult[1]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].ExonResult[2]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].ExonResult[3]]) + '"',
                        '"' + '\r\n'.join(GeneDict[i].ExonResult[4]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[0]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[1]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[2]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[3]]) + '"',
                        '"' + '\r\n'.join([str(n) for n in GeneDict[i].PACResult[4]]) + '"',
                        AnyLinks,
                        '\n']))
Out.close()


xs, ys, cs = [],[],[]
z = 0
if NEFormat == True:
    import matplotlib.pyplot as plt
    Out2 = open(Root + '_DPAC2_PDs.txt', 'w')
    Out2.write('#1.2\n')
    Out2.write('rows\tcolumns\n')
    Out2.write('\t'.join(Header2) + '\n')
    for i in GeneDict:
        for j in range(len(GeneDict[i].PACResult[0])):
            ExonName = GeneDict[i].PACResult[5][j]
            if GeneDict[i].PACResult[1][j] != 'NA' and GeneDict[i].PACResult[2][j] != 'NA':
              #if GeneDict[i].PACResult[4][j] != 'NA':
                z+=1
                Out2.write('\t'.join([GeneDict[i].Name] + ['na'] + [str(x) for x in GeneDict[i].Exons[ExonName].PDs]) + '\n')
                xs.append(GeneDict[i].PACResult[1][j])
                ys.append(GeneDict[i].PACResult[2][j])
                if GeneDict[i].PACResult[4][j] == 'Shortening':
                    cs.append('c')
                elif GeneDict[i].PACResult[4][j] == 'Lengthening':
                    cs.append('y')
                elif GeneDict[i].PACResult[4][j] == 'Both':
                    cs.append('magenta')
                else:
                    cs.append('grey')
    print(z)
    Out2.close()
    #PAPER FIGURES
    plt.plot(np.arange(0,11)/10.0,np.arange(0,11)/10.0, color='r')
    plt.plot(np.arange(0,10)/10.0,np.arange(1,11)/10.0, color='r', ls='--')
    plt.plot(np.arange(1,11)/10.0,np.arange(0,10)/10.0, color='r', ls='--')
    plt.scatter(xs,ys,color=cs, s=10, alpha=0.75)
    
    
xs, ys, cs = [],[],[]
z = 0
if NEFormat == True:
    for i in GeneDict:
        for j in GeneDict[i].Exons:
            for k in GeneDict[i].Exons[j].PACs:
                if type(GeneDict[i].Exons[j].PACs[k].padj) == float and GeneDict[i].Exons[j].PACs[k].padj < 0.1 and GeneDict[i].Exons[j].PACs[k].DRIMpadj != 'NA':
                    z+=1
                    xs.append(GeneDict[i].Exons[j].PACs[k].padj)
                    ys.append(float(GeneDict[i].Exons[j].PACs[k].DRIMpadj))
#                if GeneDict[i].PACResult[4][j] == 'Shortening':
#                    cs.append('c')
#                elif GeneDict[i].PACResult[4][j] == 'Lengthening':
#                    cs.append('y')
#                elif GeneDict[i].PACResult[4][j] == 'Both':
#                    cs.append('magenta')
#                else:
#                    cs.append('grey')
    print(z)
    #PAPER FIGURES
    plt.scatter(xs,ys, s=10, alpha=0.75)
    
