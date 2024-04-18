#!/usr/bin/python3
import argparse
import gzip
from re import findall
try:
    import pysam
    PYSAM = True
except:
    print("Pysam not installed; no support for .BAM files.")
    PYSAM = False
parser = argparse.ArgumentParser()
parser.add_argument("MetaFile", help="Meta Data File")
parser.add_argument("--Description", help= "Description in string format")
parser.add_argument("--MinScore", help= "Description in string format")
parser.add_argument("--MinCount", help= "Description in string format")
parser.add_argument("--ReqAs", help= "Description in string format, Default = 10")
args = parser.parse_args()

if args.Description:
    Description = str(args.Description)
else:
    Description = 'PACSeq'
    
if args.MinScore:
    MinScore = int(args.MinScore)
else:
    MinScore = 1
    
if args.ReqAs:
    ReqAs = int(args.ReqAs)
else:
    ReqAs = 25
    
if args.MinCount:
    MinCount = int(args.MinCount)
else:
    MinCount = 1

Dict = {}   
##OLD SCORING SCHEME FROM FIRST SUBMISSION
#Scoring_Scheme = {20:0,
#                  21:1,
#                  22:2,
#                  23:4,
#                  24:9,
#                  25:16,
#                  26:25,
#                  27:36,
#                  28:49,
#                  29:64,
#                  30:81,
#                  31:100}

##NEW SCORING SCHEME
##Includes all As above 10. 0-9 = 0; 10=10, 11=11 ... 30=30, 31=31, >31=31.

def Score_pA(pAs):
        Score = 0
        for i in range(ReqAs,31):
            Score += pAs[i] * (i)
        Score += sum(pAs[31:]) * 31
        return Score
    
###MAIN SCRIPT        
Files = open(str(args.MetaFile),'r')
File = Files.readline()
File = Files.readline()
while File:
    File = File.split()
    if PYSAM:
        FileName = File[0] + '_Hisat2-mapping.bam'
        In = pysam.AlignmentFile(FileName,'rb')
    else:
        FileName = File[0] + '_Hisat2-mapping.sam.gz'
        In = gzip.open(FileName,'rt')
    print(FileName)
    for data in In:
        Data = str(data).split()
        if Data[0][0] != '@':
            #Not a Header
            if Data[1] == '0' or Data[1] == '16':
                #Only primary alignments (SE reads only)
                if Data[1] == '0':
                    Dir = '+'
                elif Data[1] == '16':
                    Dir = '-'
                else:
                    pass
                Cigar = findall(r"[^\W\d_]+|\d+", Data[5])
                ## Breaks Cigar string in to numbers and characters
                Seq = Data[9]
                Name = Data[0]
                #pALength = int(Name.split('_pA=')[1])
                pACoords = Name.split('_')
                n=0
                for i in pACoords:
                    if 'pA=' in i:
                        pALength = int(i.split('=')[1])
                Ref = Data[2]
                if PYSAM:
                    Ref = data.reference_name
                if Dir == '+':
                    StartPos = int(Data[3]) - 1
                    #SAM format in 1=leftmost nt
                    n=0
                    for i in Cigar:
                        if i == 'D' or i == 'N' or i == 'M':
                            ## Deletion, Recombination and Mapping
                            StartPos += int(Cigar[n-1])
                        elif i == 'I':
                            ##Insertion in Read relative to reference
                            StartPos -= int(Cigar[n-1])
                        else:
                            ##Other mappings types?
                            ##Pads here 'S'
                            pass
                        n+=1
                    ##Correct for BEDFILE to represent poly(A) Cleavage site rather than first nuc of poly(A)-tail
                    if not PYSAM:
                        StartPos -= 1
                    EndPos = StartPos + 1
                    ##Giveswidth to BEDGraph file
                else:
                    if PYSAM:
                        EndPos = int(Data[3]) + 1
                        StartPos = int(Data[3])
                    else:
                        EndPos = int(Data[3])
                        StartPos = int(Data[3]) - 1
                DataPoint = (Ref + '\t' + str(StartPos) + '\t' + str(EndPos) + '\t' + Dir)
                if DataPoint not in Dict:
                    Dict[DataPoint] = [1,[0]*300]
                    Dict[DataPoint][1][pALength] += 1
                else:
                    Dict[DataPoint][0] += 1
                    Dict[DataPoint][1][pALength] += 1 
    File = Files.readline()
In.close()

OutF = open(Description + '_combined-PASs.f.bedgraph','w')
OutR = open(Description + '_combined-PASs.r.bedgraph','w')
for i in Dict:
    Score = Score_pA(Dict[i][1])
    Count = sum(Dict[i][1])
    if Score >= MinScore and Count >= MinCount:
        if i[-1] == '+':
            OutF.write(i[:-1] + str(Count) + '\t' + str(Score) + '\t+\t' + '\n')
        else:
            OutR.write(i[:-1] + str(Count) + '\t' + str(Score) + '\t-\t' + '\n')
OutF.close()
OutR.close()
    
