import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Meta_Data", help="String of metadata files")
parser.add_argument("Output", help="Merged Output Root name")
parser.add_argument("-Introns", action='store_true', help="Include Intronic PACs")
args = parser.parse_args()

Inputs = str(args.Meta_Data)
Output = str(args.Output)
if args.Introns:
    IncIntrons = True
else:
    IncIntrons = False

PACs = {}
Exons = {}
Genes = {}
Roots = []

##Compile metadata
with open(Inputs,'r') as In:
        header = In.readline()
        sample = In.readline()
        while sample:
                Root = sample.split()[0]
                Roots.append(Root)
                sample = In.readline()

##Populate Dictionary
for Root in Roots:
    Names = set()
    print(Root)
    with open(Root + '_PAC_overlaps.txt','r') as In:
        data = In.readline().split()
        while data:
            if 'exon' in data[4] or 'end' in data[4] or IncIntrons == True:
                PAC = data[4]
                Exon = data[4].split('_PAC')[0]
                Gene = data[4].split("_")[0]
                Count = int(float(data[3]))
                ##PACs
                if PAC not in PACs:
                    PACs[PAC] = {}
                else:
                    pass
                if Root not in PACs[PAC]:
                    PACs[PAC][Root] = Count
                else:
                    PACs[PAC][Root] += Count
                ##Exons
                if Exon not in Exons:
                    Exons[Exon] = {}
                else:
                    pass
                if Root not in Exons[Exon]:
                    Exons[Exon][Root] = Count
                else:
                    Exons[Exon][Root] += Count
                ##Gene
                if Gene not in Genes:
                    Genes[Gene] = {}
                else:
                    pass
                if Root not in Genes[Gene]:
                    Genes[Gene][Root] = Count
                else:
                    Genes[Gene][Root] += Count
            else:
                pass
            data = In.readline().split()

PAC_Events = {}
for i in PACs:
        event = i + '\t'
        temp = []
        for j in Roots:
                if j in PACs[i]:
                        temp.append(str(PACs[i][j]))
                else:
                        temp.append('0')
        counts = '\t'.join(temp)
        counts += '\n'
        PAC_Events[event] = counts
        
with open(Output + '_PAC-counts.txt','w') as Out:
        Out.write('\t')
        Out.write('\t'.join(Roots))
        Out.write('\n')
        for i in PAC_Events:
            Out.write(i)  
            Out.write(PAC_Events[i])
            
#with open(Output + '_PAC-counts_drimseq.txt','w') as Out:
#        Out.write('gene_id\tfeature_id\t')
#        Out.write('\t'.join(Roots))
#        Out.write('\n')
#        for i in PAC_Events:
#            Out.write(i.split('_PAC')[0])
#            Out.write('\t' + i)
#            Out.write(PAC_Events[i])  

Exon_Events = {}
for i in Exons:
        event = i + '\t'
        temp = []
        for j in Roots:
                if j in Exons[i]:
                        temp.append(str(Exons[i][j]))
                else:
                        temp.append('0')
        counts = '\t'.join(temp)
        counts += '\n'
        Exon_Events[event] = counts
        
with open(Output + '_Exon-counts.txt','w') as Out:
        Out.write('\t')
        Out.write('\t'.join(Roots))
        Out.write('\n')
        for i in Exon_Events:
            Out.write(i)  
            Out.write(Exon_Events[i])

#with open(Output + '_Exon-counts_drimseq.txt','w') as Out:
#        Out.write('gene_id\tfeature_id\t')
#        Out.write('\t'.join(Roots))
#        Out.write('\n')
#        for i in Exon_Events:
#            Out.write(i.split("_")[0])
#            Out.write('\t' + i)
#            Out.write(Exon_Events[i])  

Gene_Events = {}
for i in Genes:
        event = i + '\t'
        temp = []
        for j in Roots:
                if j in Genes[i]:
                        temp.append(str(Genes[i][j]))
                else:
                        temp.append('0')
        counts = '\t'.join(temp)
        counts += '\n'
        Gene_Events[event] = counts

with open(Output + '_Gene-counts.txt','w') as Out:
        Out.write('\t')
        Out.write('\t'.join(Roots))
        Out.write('\n')
        for i in Gene_Events:
            Out.write(i)  
            Out.write(Gene_Events[i])
                
