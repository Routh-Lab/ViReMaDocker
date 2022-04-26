##      Copyright (c) 2013-2021 Andrew Laurence Routh
##      
##      Permission is hereby granted, free of charge, to any person obtaining a copy
##      of this software and associated documentation files (the "Software"), to deal
##      in the Software without restriction, including without limitation the rights
##      to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##      copies of the Software, and to permit persons to whom the Software is
##      furnished to do so, subject to the following conditions:
##      
##      The above copyright notice and this permission notice shall be included in
##      all copies or substantial portions of the Software.
##      
##      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##      IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##      FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##      AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##      LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##      OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
##      THE SOFTWARE.
##
##      -------------------------------------------------------------------------------------------
from os.path import exists
import gzip
from os import makedirs
import argparse
import ConfigViReMa as cfg
from math import fabs
from re import finditer, findall
import time
import numpy as np
from subprocess import check_output
SamHeaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']
cfg.ScrutSAM = False
start = time.time()

##      ------------------------------------------------------------------------------------------------------------
##      Compound_Handling_Script will determine whether trimmed nucleotides found between recombination
##      sites are present between the donor and accpetor sites for these two recombination events.
##      If so, they are assumed to have arisen due to multiple recombination events occuring within close proximity.
##      This is only a good assumption in viral genomes (where the number if possible matches is low) and
##      when there are sufficient nucleotides in the trimmed sequence. This number is set at the command line
##      wih the option: --Compound_Handling.  A default value of 10 is recommended.  This value must be
##      larger than the MicroInDel number.
##      ------------------------------------------------------------------------------------------------------------

def Compound_Handling_Script(Donor, DonorSite, Insertion, AcceptorSite, uDelDicts, RecDicts, ReadName):
    if "_RevStrand" in Donor:
            DonorA = ">" + Donor[:-10]
            Insertion = Rev_Comp(Insertion)
            DonorSite, AcceptorSite = AcceptorSite, DonorSite
    else:
            DonorA = ">" + Donor
    Frag = cfg.Genes[DonorA][int(DonorSite):int(AcceptorSite)]
    if Insertion in Frag:
        Hits = [m.start() for m in finditer(Insertion, Frag)]
        if len(Hits) == 1:
            if "_RevStrand" in Donor:
                #Unique Compound Recombination Site Found
                NewAcceptorSite = str(int(Hits[0]) + int(DonorSite) + 1)
                if int(NewAcceptorSite) - int(DonorSite) - 1 <= cfg.MicroInDel_Length:
                        AddToDict(Donor, Donor, NewAcceptorSite, DonorSite, uDelDicts, ReadName)
                else:
                        AddToDict(Donor, Donor, NewAcceptorSite, DonorSite, RecDicts, ReadName)
                NewDonorSite = str(int(Hits[0]) + int(DonorSite) + len(Insertion) + 1) 
                if int(AcceptorSite) - int(NewDonorSite) - 1 <= cfg.MicroInDel_Length:              
                        AddToDict(Donor, Donor, AcceptorSite, NewDonorSite, uDelDicts, ReadName)
                else:
                        AddToDict(Donor, Donor, AcceptorSite, NewDonorSite, RecDicts, ReadName)
            else:
                #Unique Compound Recombination Site Found
                NewAcceptorSite = str(int(Hits[0]) + int(DonorSite) + 1)
                if int(NewAcceptorSite) - int(DonorSite) - 1 <= cfg.MicroInDel_Length:
                        AddToDict(Donor, Donor, DonorSite, NewAcceptorSite, uDelDicts, ReadName)
                else:
                        AddToDict(Donor, Donor, DonorSite, NewAcceptorSite, RecDicts, ReadName)
                NewDonorSite = str(int(Hits[0]) + int(DonorSite) + len(Insertion) + 1) 
                if int(AcceptorSite) - int(NewDonorSite) - 1 <= cfg.MicroInDel_Length:              
                        AddToDict(Donor, Donor, NewDonorSite, AcceptorSite, uDelDicts, ReadName)
                else:
                        AddToDict(Donor, Donor, NewDonorSite, AcceptorSite, RecDicts, ReadName)
            return "HIT"
                
##      ----------------------------------------------------------------------------------------------------------
##      UniquifyReport() removes identical results.  Reads giving identical results may be PCR duplicates.  
##      Similarly, finding unique multiple unique reads over single recombination junctions validates recombinant
##      ----------------------------------------------------------------------------------------------------------

def UniquifyReport(FileIn, FileOut):
    print("Removing potential PCR duplicates...")
    TempSet = set()
    n = 0
    m = 0
    if cfg.FileOut[-3:] == '.gz':
        DeDupedData = gzip.open(cfg.FileOut,'wb')
    else:
        DeDupedData = open(cfg.FileOut,'w')
    if cfg.FileOut[-3:] == '.gz':
        OPEN = gzip.open(cfg.FileIn,'rt')
    else:
        OPEN = open(cfg.FileIn,'r')
    with OPEN as InputData:  
        line = InputData.readline()
        while line:
            if line[:3] in SamHeaders:
                DeDupedData.write(line)
            else:
                n+=1
                UniqueData = line.split()[1:9]
                if 'TC:i:' in line:
                    NSegs = int(line[line.index('TC:i:') + 5])
                    CurrentSeg = int(line[line.index('FI:i:') + 5])
                    while CurrentSeg < NSegs:
                        newline = InputData.readline()
                        line += newline
                        UniqueData += newline[1:9]
                        CurrentSeg = int(newline[newline.index('FI:i:') + 5])
                UniqueData = '\t'.join(UniqueData)
                if UniqueData not in TempSet:
                    DeDupedData.write(line)
                    m+=1
                    TempSet.add(UniqueData) 
            line = InputData.readline()
    DeDupedData.close()
    print("Total of %s reads/segments in original dataset" % n)
    print("%s reads remaining after removing potential PCR duplicates" % m)

##      ----------------------------------------------------------------------------------------
##      Function BedGraph_Plot() will find regions deleted or duplicated due to recombination and return
##      a string string of frequencies and nucleotide positions. 
##      ----------------------------------------------------------------------------------------

def BEDGraph_Plot():
        DelCover = {}
        InsCover = {}
        with open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_Recombination_Results.bed", "r") as Input_Data:
            TrackName = Input_Data.readline()
            for i in cfg.RefsLib1:
                if "_RevStrand" not in i:
                    DelCover[i] = [0] * len(cfg.Genes[i])
                    InsCover[i] = [0] * len(cfg.Genes[i])
            for line in Input_Data:
                line = line.split()
                Gene_Name = line[0]
                Donorsite = int(line[1])
                Acceptorsite = int(line[2])
                Count = int(line[4])
                Strand = line[5]
                if Strand == '+':
                    if Donorsite < Acceptorsite:
                        while Donorsite < Acceptorsite:
                            DelCover[Gene_Name][Donorsite] += Count
                            Donorsite +=1
                    elif Acceptorsite < Donorsite:
                        while Acceptorsite < Donorsite:
                            InsCover[Gene_Name][Acceptorsite] += Count
                            Acceptorsite +=1
        OutputFile = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_Conservation.bedgraph","w")
        OutputFile.write('track type=bedGraph name="Virus_Conservation" description="Virus_Conservation"\n')
        for i in DelCover:
                n = 0
                lastline = [i,0,0,0]
                for j in DelCover[i]:
                    CurrLine = [i, n, n+1, j]
                    if CurrLine[-1] != lastline[-1]:
                        OutputFile.write('\t'.join([str(k) for k in lastline]) + '\n')
                        lastline = CurrLine
                    else:
                        lastline[2] += 1
                    n+=1
                OutputFile.write('\t'.join([str(k) for k in lastline]) + '\n')
        OutputFile.close()
        OutputFile = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_Duplications.bedgraph","w")
        OutputFile.write('track type=bedGraph name="Virus_Duplications" description="Virus_Duplications"\n')
        for i in InsCover:
                n = 0
                lastline = [i,0,0,0]
                for j in InsCover[i]:
                    CurrLine = [i, n, n+1, j]
                    if CurrLine[-1] != lastline[-1]:
                        OutputFile.write('\t'.join([str(k) for k in lastline]) + '\n')
                        lastline = CurrLine
                    else:
                        lastline[2] += 1
                    n+=1
                OutputFile.write('\t'.join([str(k) for k in lastline]) + '\n')
        OutputFile.close()
        ##CuttingSite BEDGRAPH
#        OutF = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_cuttingsites.f.bedgraph","w")
#        OutF.write('track type=bedGraph name=Virus_CuttingSiteCoverage description="Virus_CuttingSiteCoverage"\n')
#        OutF.close()
#        OutR = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_cuttingsites.r.bedgraph","w")
#        OutR.write('track type=bedGraph name=Virus_RevStrand_CuttingSiteCoverage description="Virus_RevStrand_CuttingSiteCoverage"\n')
#        OutR.close()
#        for i in cfg.RefsLib1_CuttingSites:
#            if 'RevStrand' in i:
#                with open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_cuttingsites.r.bedgraph","a") as OutR:
#                    Start = 1
#                    Stop = 2
#                    Last = cfg.RefsLib1_CuttingSites[i][1]
#                    for j in range(1, len(cfg.RefsLib1_CuttingSites[i])):
#                        if cfg.RefsLib1_CuttingSites[i][j] != Last:
#                            OutR.write(i[:-10] + '\t' + str(Start) + '\t' + str(Stop) + '\t' + str(Last) + '\n')
#                            Start = j
#                            Stop = j+1
#                            Last = cfg.RefsLib1_CuttingSites[i][j]
#                        else:
#                            Stop += 1
#            else:
#                with open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_cuttingsites.f.bedgraph","a") as OutF:
#                    Start = 1
#                    Stop = 2
#                    Last = cfg.RefsLib1_CuttingSites[i][1]
#                    for j in range(1, len(cfg.RefsLib1_CuttingSites[i])):
#                        if cfg.RefsLib1_CuttingSites[i][j] != Last:
#                            OutF.write(i + '\t' + str(Start) + '\t' + str(Stop) + '\t' + str(Last) + '\n')
#                            Start = j
#                            Stop = j+1
#                            Last = cfg.RefsLib1_CuttingSites[i][j]
#                        else:
#                            Stop += 1
                            
##      ----------------------------------------------------------------------------------------
##      Function Rev_Comp() will return the Reverse Complement of a given DNA string
##      ----------------------------------------------------------------------------------------

def Rev_Comp(Seq):
        Seq = Seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(Seq)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)[::-1]

##      ----------------------------------------------------------------------------------------------------
##      Indices will find the locations of pertinent information in the results file as directed by the Code
##      ----------------------------------------------------------------------------------------------------
    
def Indices(List):
    n = 1
    Ms = []
    Xs = []
    for i in List:
        if i == "M":
            Ms.append(n)
            n+=2
        else:
            Xs.append(n)
            n+=1
    return [Ms, Xs]

##      -------------------------------------------------------------------------------------------------------
##      ExtractRefData() will find the names of the genes used in the virus or host genome references.
##      Bowtie-inspect must be in $PATH.
##      -------------------------------------------------------------------------------------------------------
            
def ExtractRefData():
        cfg.RefsLib1 = set()
        cfg.RefsLib2 = set()
        print("Extracting Virus Genes...")
        z = check_output([cfg.Aligner_Directory + 'bowtie-inspect', '-a', '1000000', cfg.Lib1], universal_newlines=True).splitlines()
        cfg.Genes = {}
        #Temp = []
        for i in z:
            if i[0] == '>':
                Name = i.split()[0][1:]
                Name = Name.rstrip()
                cfg.RefsLib1.add(Name)
                cfg.RefsLib1.add(Name + "_RevStrand")
                print(Name)
                Temp = []
            else:
                Temp.append(i.rstrip())
                Gene = "".join(Temp)
                cfg.Genes[Name] = Gene
        if cfg.Lib2:
            print("Extracting Host Genes...")
            z = check_output([cfg.Aligner_Directory + 'bowtie-inspect', '-a', '1000000', cfg.Lib2], universal_newlines=True).splitlines()
            for i in z:
                if i[0] == '>':
                    Name = i.split()[0][1:]
                    Name = Name.rstrip()
                    cfg.RefsLib2.add(Name)
                    cfg.RefsLib2.add(Name + "_RevStrand")
                    #print(Name)
                    Temp = []
                else:
                    if cfg.Compound_Handling or cfg.Defuzz or cfg.MicroInDel_Length or cfg.BED:
                        Temp.append(i.rstrip())
                        Gene = "".join(Temp)
                        cfg.Genes[Name] = Gene
                    else:
                        pass
        else:
            pass
        for Name in list(cfg.Genes.keys()):
            cfg.Genes[Name + "_RevStrand"] = Rev_Comp(cfg.Genes[Name])
        print("Finished extracting gene data")
        return cfg.RefsLib1, cfg.RefsLib2, cfg.Genes

##      -------------------------------------------------------------------------------------------------------
##      ExtractRefDataBWA() will find the names of the genes used in the virus or host genome references.
##      Input is typical FASTA file, as per the BWA command line.
##      -------------------------------------------------------------------------------------------------------

def ExtractRefDataBWA():
        cfg.RefsLib1 = set()
        cfg.RefsLib2 = set()
        cfg.Genes = {}
        print("Extracting Virus Gene Data...")
        with open(cfg.Lib1,'r') as FASTAIN:
            for line in FASTAIN:
                if line[0] == '>':
                    Name = line.split()[0][1:]
                    Name = Name.rstrip()
                    cfg.RefsLib1.add(Name)
                    cfg.RefsLib1.add(Name + "_RevStrand")
                    print(Name)
                    Temp = []
                else:
                    Temp.append(line.rstrip())
                    Gene = "".join(Temp)
                    cfg.Genes[Name] = Gene
        if cfg.Lib2:
            print("Extracting Host Gene Data...")
            with open(cfg.Lib2,'r') as FASTAIN:
                for line in FASTAIN:
                    if line[0] == '>':
                        Name = line.split()[0][1:]
                        Name = Name.rstrip()                        
                        cfg.RefsLib2.add(Name)
                        cfg.RefsLib2.add(Name + "_RevStrand")
                        print(Name)
                        Temp = []
                    else:
                        if cfg.Compound_Handling or cfg.Defuzz or cfg.MicroInDel_Length or cfg.BED:
                            Temp.append(line.rstrip())
                            Gene = "".join(Temp)
                            cfg.Genes[Name] = Gene
                        else:
                            pass
        for Name in list(cfg.Genes.keys()):
            cfg.Genes[Name + "_RevStrand"] = Rev_Comp(cfg.Genes[Name])
        print("Finished extracting gene data")
        return cfg.RefsLib1, cfg.RefsLib2, cfg.Genes
    
##      -------------------------------------------------------------------------------------------
##      AddToDict() takes the Donor and Acceptor sites and references for a given recombination event,
##      and collates them into a Dictionary which will later be written to a results file.
##      -------------------------------------------------------------------------------------------

def AddToDict(Donor, Acceptor, DonorSite, AcceptorSite, Dict, ReadName):
        LeftFuzz, RightFuzz, uHomology = 0,0,''
        if cfg.Defuzz:
                if "_RevStrand" in Donor:
                        if int(AcceptorSite) != int(DonorSite) - 1:
                                LeftFuzz, RightFuzz, uHomology = FindFuzz(Donor, DonorSite, Acceptor, AcceptorSite, (cfg.MaxFuzz+1))
                                if LeftFuzz != 0 or RightFuzz != 0:
                                        if cfg.Defuzz == 'Centre':
                                                Fuzz = RightFuzz - LeftFuzz
                                                DonorSite = str(int(DonorSite) + (int(Fuzz/2)))
                                                AcceptorSite = str(int(AcceptorSite) + (int(Fuzz/2)))
                                        elif cfg.Defuzz == 'Left':
                                                DonorSite = str(int(DonorSite) - LeftFuzz)
                                                AcceptorSite = str(int(AcceptorSite) - LeftFuzz)
                                        else: ##Right
                                                DonorSite = str(int(DonorSite) + RightFuzz)
                                                AcceptorSite = str(int(AcceptorSite) + RightFuzz)
                        else:
                                pass
                else:
                        if int(AcceptorSite) != int(DonorSite) + 1:
                                LeftFuzz, RightFuzz, uHomology = FindFuzz(Donor, DonorSite, Acceptor, AcceptorSite, (cfg.MaxFuzz+1))
                                if LeftFuzz != 0 or RightFuzz != 0:
                                        if cfg.Defuzz == 'Centre':
                                                Fuzz = RightFuzz - LeftFuzz
                                                DonorSite = str(int(DonorSite) + (int(Fuzz/2)))
                                                AcceptorSite = str(int(AcceptorSite) + (int(Fuzz/2)))
                                        elif cfg.Defuzz == 'Left':
                                                DonorSite = str(int(DonorSite) - LeftFuzz)
                                                AcceptorSite = str(int(AcceptorSite) - LeftFuzz)
                                        else: ##Right
                                                DonorSite = str(int(DonorSite) + RightFuzz)
                                                AcceptorSite = str(int(AcceptorSite) + RightFuzz)
                        else:
                                pass
        else:
                pass
        if RightFuzz + LeftFuzz <= cfg.MaxFuzz:
                if (Donor + "_to_" + Acceptor) not in Dict:
                        Dict[Donor + "_to_" + Acceptor] = {}
                else:
                        pass
                if cfg.ReadNamesEntry:
                        if (DonorSite + "_to_" + AcceptorSite) not in Dict[Donor + "_to_" + Acceptor]:
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_to_" + AcceptorSite] = [1, [ReadName + '_Fuzz=' + str(uHomology)]]
                        else:
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_to_" + AcceptorSite][0] += 1
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_to_" + AcceptorSite][1].append(ReadName + '_Fuzz=' + str(uHomology))
                elif cfg.FuzzEntry:
                        if (DonorSite +  "_" + str(uHomology)+ "_" + AcceptorSite) not in Dict[Donor + "_to_" + Acceptor]:
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_" + str(uHomology)+ "_" + AcceptorSite] = 1
                        else:
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_" + str(uHomology)+ "_" + AcceptorSite] += 1
                else:
                        if (DonorSite + "_to_" + AcceptorSite) not in Dict[Donor + "_to_" + Acceptor]:
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_to_" + AcceptorSite] = 1
                        else:
                                Dict[Donor + "_to_" + Acceptor][DonorSite + "_to_" + AcceptorSite] += 1
        else:
                pass
        return Donor + "_to_" + Acceptor, DonorSite + "_to_" + AcceptorSite

##      ----------------------------------------------------------------------------------------------------------------
##      FindFuzz() will determine the number of overlapping nucleotides in the donor region and the nucleotides
##      preceeding the acceptor site.
##      -----------------------------------------------------------------------------------------------------------------

def FindFuzz(Donor, DonorSite, Acceptor, AcceptorSite, MaxFuzz):
        LeftFuzz = 0
        RightFuzz = 0
        uHomology = ''
        LeftHomology = ''
        RightHomology = ''
        if "_RevStrand" in Donor:
                Gene = Donor[:-10]
                UpstreamDonorSeq = Rev_Comp(cfg.Genes[Gene][int(DonorSite) - 1:int(DonorSite) - 1 + MaxFuzz])
                DownstreamDonorSeq = Rev_Comp(cfg.Genes[Gene][int(DonorSite) - 1 + MaxFuzz:int(DonorSite)])
        else:
                Gene = Donor
                UpstreamDonorSeq = cfg.Genes[Gene][int(DonorSite) - MaxFuzz:int(DonorSite)]
                DownstreamDonorSeq = cfg.Genes[Gene][int(DonorSite):int(DonorSite) + MaxFuzz]
        if "_RevStrand" in Acceptor:
                Gene = Acceptor[:-10]
                UpstreamAcceptorSeq = Rev_Comp(cfg.Genes[Gene][int(AcceptorSite): int(AcceptorSite) + MaxFuzz])
                DownstreamAcceptorSeq = Rev_Comp(cfg.Genes[Gene][int(AcceptorSite) + MaxFuzz:int(AcceptorSite)])
        else:
                Gene = Acceptor
                UpstreamAcceptorSeq = cfg.Genes[Gene][int(AcceptorSite) - 1 - MaxFuzz: int(AcceptorSite) - 1]
                DownstreamAcceptorSeq = cfg.Genes[Gene][int(AcceptorSite) - 1:int(AcceptorSite) - 1 + MaxFuzz]
        for i in range(len(UpstreamAcceptorSeq)):
                try:
                    if UpstreamAcceptorSeq[-i-1] == UpstreamDonorSeq[-i-1]:
                            LeftFuzz += 1
                            LeftHomology += UpstreamAcceptorSeq[-i-1]
                    else:
                            break
                except:
                        break
        LeftHomology = LeftHomology[::-1]
        for i in range(len(DownstreamDonorSeq)):
                try:
                    if DownstreamAcceptorSeq[i] == DownstreamDonorSeq[i]:
                            RightFuzz += 1
                            RightHomology += DownstreamAcceptorSeq[i]
                    else:
                            break
                except:
                        break
        uHomology = LeftHomology + RightHomology
        return LeftFuzz, RightFuzz, uHomology

##      -----------------------------------------------------------------------------------------------------------------
##      AddInsToDict() takes the Donor and Acceptor sites, references and trimmed nucleotides for a given Insertion event
##      and collate them into a Dictionary which will later be written to a results file.
##      -----------------------------------------------------------------------------------------------------------------

def AddInsToDict(Donor, DonorSite, AcceptorSite, Insertion, Dict, ReadName):
    if Donor not in Dict:
        Dict[Donor] = {}
    else:
        pass
    if cfg.ReadNamesEntry:
        if (str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)) not in Dict[Donor]:
            Dict[Donor][str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)] = [1, [ReadName]]
        else:
            Dict[Donor][str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)][0] += 1
            Dict[Donor][str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)][1].append(ReadName)
    else:
        if (str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)) not in Dict[Donor]:
            Dict[Donor][str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)] = 1
        else:
            Dict[Donor][str(DonorSite) + "_" + Insertion + "_" + str(AcceptorSite)] += 1

def ContractX(x):
    while 'Mismatch' in x:
        y = x.index('Mismatch')
        newx= x[:y-2]
        newx[-1] += x[y-1].split("_")[1]
        newx[-1]+= x[y+3]
        newx += x[y+4:]
        x = newx
    while 'Sub' in x:
        y = x.index('Sub')
        newx= x[:y-2]
        newx[-1] += x[y-1].split("_")[1]
        newx[-1]+= x[y+3]
        newx += x[y+4:]
        x = newx
    return x

##      -----------------------------------------------------------------------
##      RecreateOldFormatfromSAM() Turns information from a SAM file into form
##      originally used in ViReMa to compile recombination results. 
##      -----------------------------------------------------------------------

def RecreateOldFormatfromSAM(lines):
    Totalline = []
    TotalCode = []
    Name = lines[0][0]
    Totalline.append(Name)
    for Seg in lines:
        line = []
        Code = []
        FLAG = Seg[1]
        if FLAG == '4':
            ##UnMapped is format: Name, Seq
            line.append(Seg[9])
        else:
            CIGAR = findall(r"[^\W\d_]+|\d+", Seg[5])
            Read = Seg[9]
            Start = int(Seg[3])
            Ref = Seg[2]
            PrevRec = False
            if FLAG == '0' or FLAG == '2048' or FLAG == '256':
                if Ref in cfg.RefsLib1:
                    ##Only make cutting site arrays for virus
                    FindCuttingSitesfromCIGAR(Seg[5], Start, cfg.Seed, Ref)
                else:
                    pass
                LenMapped = int(CIGAR[0])
                if CIGAR[1] == 'H':
                    CIGAR = CIGAR[2:]
                else:
                    pass
                if CIGAR[1] == 'S':
                    line.append(Read[:int(CIGAR[0])])
                    Read = Read[int(CIGAR[0]):]
                    if Code and Code[-1] == 'X':
                        Code[-2] = str(int(Code[-2]) + int(CIGAR[0]))
                    else:
                        Code += (CIGAR[:1]) + ['X']
                    CIGAR = CIGAR[2:]
                elif CIGAR[1] == 'M' and cfg.ScrutSAM and LenMapped <= cfg.Seed and len(CIGAR) > 2: 
                    if CIGAR[3] == 'N':
                        ## First Mapped Segment before Rec and was shorter than seed allowed in ViReMa 
                        ## (probably due to non-virema mapper)
                        line.append(Read[:int(CIGAR[0])])
                        Read = Read[int(CIGAR[0]):]
                        if Code and Code[-1] == 'X':
                            Code[-2] = str(int(Code[-2]) + int(CIGAR[0]))
                        else:
                            Code += (CIGAR[:1]) + ['X']
                        CIGAR = CIGAR[2:]
                    else:
                        pass
                else:
                    pass
                CurrentNt = int(Seg[3])
                while CIGAR:
                    if CIGAR[1] == 'M':
                        LenMapped = int(CIGAR[0])
                        if cfg.ScrutSAM and LenMapped <= cfg.Seed and PrevRec == True and len(CIGAR) < 3: 
                            ## Last Mapped Segment was shorter than seed allowed in ViReMa 
                            ## (probably due to non-virema mapper)
                            Xs = Read[:int(CIGAR[0])]
                            Read = Read[int(CIGAR[0]):]
                            Code += (CIGAR[:1]) + ['X']
                            line.append(Xs)
                            CurrentNt += int(CIGAR[0])
                            CIGAR = CIGAR[2:]
                            PrevRec = False
                        else:
                            line.append(Seg[2].split()[0])
                            line.append(str(CurrentNt) + '_' + str(CurrentNt + int(CIGAR[0]) - 1))
                            Read = Read[int(CIGAR[0]):]
                            Code += CIGAR[:2]
                            CurrentNt += int(CIGAR[0])
                            CIGAR = CIGAR[2:]
                            PrevRec = False                            
                    elif CIGAR[1] == 'X' or CIGAR[1] == 'S':
                        Xs = Read[:int(CIGAR[0])]
                        Read = Read[int(CIGAR[0]):]
                        Code += (CIGAR[:1]) + ['X']
                        line.append(Xs)
                        CurrentNt += int(CIGAR[0])
                        CIGAR = CIGAR[2:]
                        PrevRec = False
                    elif CIGAR[1] == 'I':
                        Is = Read[:int(CIGAR[0])]
                        Read = Read[int(CIGAR[0]):]
                        Code += (CIGAR[:1]) + ['U']
                        CIGAR = CIGAR[2:]
                        line.append(Is)
                        PrevRec = False
                    elif CIGAR[1] == 'N':
                        CurrentNt += int(CIGAR[0])
                        CIGAR = CIGAR[2:]
                        PrevRec = True
                    elif CIGAR[1] == 'D':
                        CurrentNt += int(CIGAR[0])
                        CIGAR = CIGAR[2:]
                        PrevRec = False
                    else: #H
                        CIGAR = CIGAR[2:]
                        PrevRec = False
            elif FLAG == '16' or FLAG == '2064' or FLAG == '272':
                if Ref in cfg.RefsLib1:
                    ##Only make cutting site arrays for virus
                    FindCuttingSitesfromCIGAR(Seg[5], Start, cfg.Seed, Ref + '_RevStrand')
                else:
                    pass
                if CIGAR[1] == 'H':
                    CIGAR = CIGAR[2:]
                else:
                    pass
                if CIGAR[1] == 'S':
                    line.append(Read[:int(CIGAR[0])])
                    Read = Read[int(CIGAR[0]):]
                    if Code and Code[-1] == 'X':
                        Code[-2] = str(int(Code[-2]) + int(CIGAR[0]))
                    else:
                        Code += (CIGAR[:1]) + ['X']
                    CIGAR = CIGAR[2:]
                else:
                    pass
                CurrentNt = int(Seg[3])
                while CIGAR:
                    if CIGAR[1] == 'M':
                        line.insert(0, str(CurrentNt + int(CIGAR[0]) - 1) + '_RevStrand_' + str(CurrentNt))
                        line.insert(0, Seg[2].split()[0] + '_RevStrand')
                        Read = Read[int(CIGAR[0]):]
                        Code = CIGAR[:2] + Code
                        CurrentNt += int(CIGAR[0])   
                        CIGAR = CIGAR[2:]
                    elif CIGAR[1] == 'X' or CIGAR[1] == 'S':
                        Xs = Read[:int(CIGAR[0])]
                        Read = Read[int(CIGAR[0]):]
                        Code = (CIGAR[:1]) + ['X'] + Code
                        line.insert(0, Xs)
                        CurrentNt += int(CIGAR[0])
                        CIGAR = CIGAR[2:]
                    elif CIGAR[1] == 'I':
                        Is = Read[:int(CIGAR[0])]
                        Read = Read[int(CIGAR[0]):]
                        Code = CIGAR[:1] + ['U'] + Code
                        CIGAR = CIGAR[2:]
                        line.insert(0, Is)
                    elif CIGAR[1] == 'N' or CIGAR[1] == 'D':
                        CurrentNt += int(CIGAR[0])
                        CIGAR = CIGAR[2:]
                    else: #H
                        CIGAR = CIGAR[2:]
        Totalline += line
        TotalCode += Code
    TotalCode = ''.join(TotalCode)
    Totalline.append(TotalCode)
    return Totalline

def ReverseEvents(line, Entries):
    newline = []
    newline.append(line[0])
    x = line[1:]
    for i in range(int(Entries/3))[::-1]:
        name = x[i*3]
        name = ''.join(name.split('_RevStrand'))
        newline.append(name)
        y = x[i*3+1]
        if "_" in y:
            if '_to_' in y:
                y = y.split("_")
                y = y[::-1]
                y = "_".join(y)
            else:
                y = y.split("_")
                y = y[::-1]
                y[1] = Rev_Comp(y[1])
                y = "_".join(y)
            newline.append(y)
            newline.append(x[i*3+2])
        else:
            z = x[i*3+2]
            z = Rev_Comp(z)
            y = str(int(y) - len(z) + 1)
            newline.append(y)
            newline.append(z)
    return newline
        
##      -------------------------------------------------------------------------------------------
##      FindCuttingSitesfromCIGAR() will read all the information collated in each Dictionary and write out the
##      -------------------------------------------------------------------------------------------

def FindCuttingSitesfromCIGAR(Cigar, Start, MinSegmentLength, Ref):
    ##Examples of CIGAR:
    ##120M all mapped:  Count all cutting sites in mapped regions
    ##60M1X59M all mapped:  Count all cutting sites in mapped regions
    ##60M60N60M deletion:  Count all cutting sites in mapped and deleted regions
    ##60M4I60M insertion:  Count all cutting sites in mapped regions, ignore insertion (cannot count twice)
    ##3S57M80H soft pad and hard pad: Count all cutting sites in mapped regions, ignore pads    n=0
    CIGAR = findall(r"[^\W\d_]+|\d+", Cigar)
    ##FIND length of read    
    #Adds = ['M', 'X', 'S', 'H', 'I']
    m=0        
    while CIGAR:
        #if CIGAR[1] in Adds:
        m+= int(CIGAR[0])
        CIGAR = CIGAR[2:]
    ##Find CSs
    Mask = np.array([1]*m)
    #Mask[:MinSegmentLength] = 0
    #Mask[-MinSegmentLength:] = 0
    Sites = np.array([0]*m)
    n=0
    CIGAR = findall(r"[^\W\d_]+|\d+", Cigar) 
    while CIGAR:
        if CIGAR[1] == 'M' or CIGAR[1] == 'X' or CIGAR[1] == 'D': #or CIGAR[1] == 'N' 
            ##Add ns
            ##Add cutting sites
            From = Start
            To = Start + int(CIGAR[0])
            for i in range(From, To):
                Sites[n] = i
                n+=1
                Start += 1
        elif CIGAR[1] == 'N':# or CIGAR[1] == 'D':
            ##Add cutting sites
            From = Start
            To = Start + int(CIGAR[0])
            for i in range(From, To):
                #Sites[n] = i
                Start += 1
        elif CIGAR[1] == 'S' or CIGAR[1] == 'H' or CIGAR[1] == 'I':# or CIGAR[1] == 'N':
            ##Add ns
            From = Start
            To = Start + int(CIGAR[0])
            for i in range(From, To):
                n+=1
                #Start += 1
        CIGAR = CIGAR[2:]
    Sites = Sites * Mask
    for i in Sites:
        cfg.RefsLib1_CuttingSites[Ref][i] += 1


##      -------------------------------------------------------------------------------------------
##      ResultsSort() is the backbone of the results compilation script.  It will take each line from the results
##      output from ViReMa and determine the nature of every discovered event.  This includes recombinations, insertions, subsitutions, etc.
##      The Code appended to the end of each line in the ViReMa output allows ResultsSort() to group particular events.
##      The finer details are then established based upon parameters set in the command-line (e.g. MicroInDel length), and upon the
##      information provided in the ViReMa output.
##      -------------------------------------------------------------------------------------------

def ResultsSort(File1):
        SamHeaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']

        ##      -------------------------------------------------------------------------------------------
        ##      WriteFinalDict() will read all the information collated in each Dictionary and write out the
        ##      results to the results files.
        ##      -------------------------------------------------------------------------------------------

        def WriteFinalDict(DictName, Mod):

            ##      ---------------------------------------------------------------
            ##      WritetoBedFile() will append the entry to the optional BED File
            ##      Sloppy handling of item names is why function is defined within
            ##      another function. Needs Re-write, but works for now.
            ##      ---------------------------------------------------------------
        
            def WritetoBEDFile(Genes, Entry, TargetFile):
                    # Genes is in format: 'Donor_to_Acceptor'
                    # Entry is in format: [DonorSite, 'to', AcceptorSite, '#', Count]
                    Genes = Genes.split("_to_")
                    ##Set allowed targetfiles
                    if cfg.MicroInDel_Length > 0:
                        x = [VirusuIns, VirusInsertions, VirusuDels, VirusRecs]
                        if cfg.Lib2:
                            y = [HostuIns, HostuDels, HostRecs]
                        else:
                            y = []
                    else:
                        x = [VirusRecs, VirusInsertions]
                        if cfg.Lib2:
                            y = [HostRecs]
                        else:
                            y = []
                    if len(Genes) == 2 and Genes[0] != Genes[1]:
                        ##Copy-Back or gene fusion
                        ##Use BEDPE format(SMC-RNA)
                        if Genes[0][-10:] == "_RevStrand":
                            Genes[0] = Genes[0][:-10]
                            Dir1 = '-'
                            DonorRightSeq = Rev_Comp(cfg.Genes[Genes[0]][int(Entry[0]) - cfg.Seed:int(Entry[0])])
                            DonorLeftSeq = Rev_Comp(cfg.Genes[Genes[0]][int(Entry[0]):int(Entry[0]) + cfg.Seed])
                        else:
                            Dir1 = '+'
                            DonorLeftSeq = cfg.Genes[Genes[0]][int(Entry[0]) - cfg.Seed:int(Entry[0])]
                            DonorRightSeq = cfg.Genes[Genes[0]][int(Entry[0]):int(Entry[0]) + cfg.Seed]
                        if Genes[1][-10:] == "_RevStrand":
                                Genes[1] = Genes[1][:-10]
                                Dir2 = '-'                            
                                AcceptorRightSeq = Rev_Comp(cfg.Genes[Genes[1]][int(Entry[2]) - cfg.Seed -1:int(Entry[2]) - 1])
                                AcceptorLeftSeq = Rev_Comp(cfg.Genes[Genes[1]][int(Entry[2]) - 1:int(Entry[2]) + cfg.Seed - 1])
                        else:
                                Dir2 = '+'
                                AcceptorLeftSeq = cfg.Genes[Genes[1]][int(Entry[2]) - cfg.Seed - 1:int(Entry[2]) - 1]
                                AcceptorRightSeq = cfg.Genes[Genes[1]][int(Entry[2]) - 1:int(Entry[2]) + cfg.Seed - 1]
                        if Genes[0] != Genes[1]:
                            if TargetFile in x:
                                NAME = 'Intergenic-Fusion'
                                BEDFILE = VirusFusions_BED
                            elif TargetFile in y:
                                NAME = 'Gene-Fusion'
                                BEDFILE = HostFusions_BED
                            elif TargetFile == VirustoHostRecs:
                                NAME = 'Virus-Host-Fusion'
                                BEDFILE = VirustoHostRecs_BED
                            else:
                                pass
                        else:
                            if TargetFile in x:
                                NAME = 'Copy/Snap-Back'
                                BEDFILE = VirusFusions_BED
                            elif TargetFile in y:
                                NAME = 'Gene-Fusion' 
                                BEDFILE = HostFusions_BED
                            else:
                                pass
                        BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], str(int(Entry[0]) - 1), Entry[0],
                                                                                   Genes[1], Entry[2], str(int(Entry[2]) + 1),
                                                                                   NAME, Entry[4], Dir1, Dir2)
                        BEDFILE.write(BED_OUTPUT)
                    else:
                        Start = int(Entry[0])
                        Stop = int(Entry[2])
                        try:
                            LeftCount = str(cfg.RefsLib1_CuttingSites[Genes[0]][Start])
                            RightCount = str(cfg.RefsLib1_CuttingSites[Genes[0]][Stop])
                        except:
                            LeftCount, RightCount = '-','-'
                        if Genes[0][-10:] == "_RevStrand":
                            Genes[0] = Genes[0][:-10]
                            #Genes[1] = Genes[1][:-10]
                            Dir = '-'
                            DonorRightSeq = Rev_Comp(cfg.Genes[Genes[0]][int(Entry[0]) - cfg.Seed - 1:int(Entry[0]) - 1])
                            DonorLeftSeq = Rev_Comp(cfg.Genes[Genes[0]][int(Entry[0]) - 1:int(Entry[0]) + cfg.Seed - 1])
                            AcceptorRightSeq = Rev_Comp(cfg.Genes[Genes[0]][int(Entry[2]) - cfg.Seed:int(Entry[2])])
                            AcceptorLeftSeq = Rev_Comp(cfg.Genes[Genes[0]][int(Entry[2]):int(Entry[2]) + cfg.Seed])
                        else:
                            Dir = '+'
                            DonorLeftSeq = cfg.Genes[Genes[0]][int(Entry[0]) - cfg.Seed:int(Entry[0])]
                            DonorRightSeq = cfg.Genes[Genes[0]][int(Entry[0]):int(Entry[0]) + cfg.Seed]
                            AcceptorLeftSeq = cfg.Genes[Genes[0]][int(Entry[2]) - cfg.Seed - 1:int(Entry[2]) - 1]
                            AcceptorRightSeq = cfg.Genes[Genes[0]][int(Entry[2]) - 1:int(Entry[2]) + cfg.Seed - 1]
                        if TargetFile in x:
                            if Dir == '+':
                                if Stop > Start + 1:
                                    NAME = 'Deletion'
                                elif Stop == Start + 1:
                                    NAME = 'Ins:' + Entry[1]
                                else:
                                    NAME = 'Duplication'
                            elif Dir == '-':
                                if Stop < Start - 1:
                                    NAME = 'Deletion'
                                elif Stop == Start - 1:
                                    NAME = 'Ins:' + Entry[1]
                                else:
                                    NAME = 'Duplication'
                        elif TargetFile in y:
                            if Dir == '+':
                                if Stop > Start + 1:
                                    NAME = 'Splice'
                                elif Stop == Start + 1:
                                    NAME = 'Ins:' + Entry[1]
                                else:
                                    NAME = 'Duplication'
                            elif Dir == '-':
                                if Stop < Start - 1:
                                    NAME = 'Splice'
                                elif Stop == Start - 1:
                                    NAME = 'Ins:' + Entry[1]
                                else:
                                    NAME = 'Back-Splice'                            
                        else:
                            pass
                        #BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], Entry[0], Entry[2], NAME, Entry[4], Dir, LeftCount, RightCount, DonorLeftSeq + "|" + DonorRightSeq, AcceptorLeftSeq + "|" + AcceptorRightSeq)
                        if TargetFile == VirusRecs:# and Genes[0] == Genes[1]:
                            BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], Entry[0], Entry[2], NAME, Entry[4], Dir, LeftCount, RightCount, DonorLeftSeq + "|" + DonorRightSeq, AcceptorLeftSeq + "|" + AcceptorRightSeq)
                            VirusRecs_BED.write(BED_OUTPUT)
                        elif TargetFile == VirusInsertions:
                            BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], Entry[0], Entry[2], NAME, Entry[4], Dir, LeftCount, RightCount, DonorLeftSeq + "|" + DonorRightSeq, AcceptorLeftSeq + "|" + AcceptorRightSeq)
                            VirusRecs_BED.write(BED_OUTPUT)
                        elif TargetFile == VirusuIns or TargetFile == VirusuDels:
                            BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], Entry[0], Entry[2], NAME, Entry[4], Dir, LeftCount, RightCount, DonorLeftSeq + "|" + DonorRightSeq, AcceptorLeftSeq + "|" + AcceptorRightSeq)
                            VirusuRecs_BED.write(BED_OUTPUT)
                        elif TargetFile == HostRecs:
                            BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], Entry[0], Entry[2], NAME, Entry[4], Dir)
                            HostRecs_BED.write(BED_OUTPUT)
                        elif TargetFile == HostuIns or TargetFile == HostuDels:
                            BED_OUTPUT = "%s\t%s\t%s\t%s\t%s\t%s\n" % (Genes[0], Entry[0], Entry[2], NAME, Entry[4], Dir)
                            HostuRecs_BED.write(BED_OUTPUT)
                        else:
                            pass
                ##      ---------------------------------------------------------------

            for k in DictName:
                    Libs = k.split("_to_")
                    n = 0
                    for i in Libs:
                            if i[-10:] == "_RevStrand":
                                    Libs[n] = Libs[n][:-10]
                            n+=1
                    if Mod == 'Recs':
                            if Libs[0] in cfg.RefsLib1 and Libs[1] in cfg.RefsLib1:
                                    TargetFile = VirusRecs
                            elif Libs[0] in cfg.RefsLib2 and Libs[1] in cfg.RefsLib2:
                                    TargetFile = HostRecs
                            else:
                                    TargetFile = VirustoHostRecs
                    elif Mod == 'uDel':
                            if Libs[0] in cfg.RefsLib1 and Libs[1] in cfg.RefsLib1:
                                    TargetFile = VirusuDels
                            elif Libs[0] in cfg.RefsLib2 and Libs[1] in cfg.RefsLib2:
                                    TargetFile = HostuDels
                            else:
                                    pass                        
                    elif Mod == 'uIns':
                            if Libs[0] in cfg.RefsLib1:
                                    TargetFile = VirusuIns
                            else:
                                    TargetFile = HostuIns
                    elif Mod == 'Ins':
                            if Libs[0] in cfg.RefsLib1:
                                    TargetFile = VirusInsertions
                            else:
                                    TargetFile = HostInsertions
                    elif Mod == 'Sub':
                            if Libs[0] in cfg.RefsLib1:
                                    TargetFile = VirusSubstitutions
                            else:
                                    TargetFile = HostSubstitutions
                    Temp = []
                    if cfg.ReadNamesEntry:
                            for i in DictName[k]:
                                    x = [(str(i) + "_#_" + str(DictName[k][i][0])).split("_"), DictName[k][i][1]]
                                    Temp.append(x)
                            Temp.sort(key=lambda a:int(a[0][4]), reverse=True)
                            TargetFile.write("@NewLibrary: " + str(k) + "\n")
                            for i in Temp:
                                    j = '_'.join(i[0])
                                    TargetFile.write(str(j) + "\n")
                                    for Names in i[1]:
                                            TargetFile.write(str(Names) + '\t')
                                    TargetFile.write('\n')
                            TargetFile.write("\n@EndofLibrary\n")
                    else:
                            for i in DictName[k]:
                                    x = (str(i) + "_#_" + str(DictName[k][i])).split("_")
                                    Temp.append(x)
                            Temp.sort(key=lambda a:int(a[4]), reverse=True)
                            TargetFile.write("@NewLibrary: " + str(k) + "\n")
                            for i in Temp:
                                    if cfg.BED:
                                            if cfg.Lib2:
                                                    if cfg.MicroInDel_Length > 0:
                                                            BEDableTargetFiles = [VirusRecs, VirusInsertions, VirusuDels, HostRecs, HostuDels, VirusuIns, HostuIns, VirustoHostRecs]
                                                    else:
                                                            BEDableTargetFiles = [VirusRecs, VirusInsertions, HostRecs, VirustoHostRecs]
                                            else:
                                                    if cfg.MicroInDel_Length > 0:
                                                            BEDableTargetFiles = [VirusRecs, VirusInsertions, VirusuDels, VirusuIns]
                                                    else:
                                                            BEDableTargetFiles = [VirusRecs, VirusInsertions]
                                            if TargetFile in BEDableTargetFiles:
                                                    WritetoBEDFile(k, i, TargetFile)
                                            else:
                                                    pass
                                    else:
                                            pass
                                    j = '_'.join(i)
                                    TargetFile.write(str(j) + "\t")
                            TargetFile.write("\n@EndofLibrary\n")

        #Dictionaries 
        InsDicts = {}
        SubDicts = {}
        uDelDicts = {}
        uInsDicts = {}
        RecDicts = {}

        #Counts for events used for printed summary
        Padcount = 0
        LongPadcount = 0
        Totalcount = 0
        uCount = 0
        InsCount = 0
        SubCount = 0
        CompoundCount = 0
        RecombCount = 0
        ViralRecombinationCount = 0
        HostRecombinationCount = 0
        ViraltoHostRecombinationCount = 0
        UnknownRecombinationCount = 0
        UnmappedReadsCount = 0
       # UnknownRecombinations = open(cfg.Output_Dir + cfg.FileTag + "Unknown_Recombinations.txt", "w")

        ##Allow Gzipped samfile readin
        print(File1)
        if File1[-3:] == '.gz':
            OPEN = gzip.open(File1,'rt')
        else:
            OPEN = open(File1,'r')
        with OPEN as InRecombs:
            wholeline = InRecombs.readline()
            while wholeline[:3] in SamHeaders:
                wholeline = InRecombs.readline()
            #Meat of script, reads each line from ViReMa.py output and stores details into dictionaries
            while wholeline:
                lines = [wholeline.split("\t")]
                Totalcount += 1
                if 'TC:i:' in wholeline:
                    NSegs = int(wholeline[wholeline.index('TC:i:') + 5])
                    CurrentSeg = int(wholeline[wholeline.index('FI:i:') + 5])
                    while CurrentSeg < NSegs:
                        wholeline = InRecombs.readline()
                        CurrentSeg = int(wholeline[wholeline.index('FI:i:') + 5])
                        lines.append(wholeline.split("\t"))
                else:
                    pass
                line = RecreateOldFormatfromSAM(lines)
                ReadName = line[0]
                if cfg.UMI:
                    UMI = ReadName.split(cfg.UMI)[-1]
                    if UMI in cfg.UMIs:
                        UMIHalt = True
                    else:
                        UMIHalt = False
                        cfg.UMIs.add(UMI)
                else:
                    UMIHalt = False
                if UMIHalt:
                    ##Skipping read as is a PCR duplicate. 
                    wholeline = InRecombs.readline()
                else:
                    Code = ''.join(findall(r"\D", line[-1]))
                    Index = Indices(Code)
                    MCount = Code.count("M")
                    if "M" not in Code:
                        #UnMapped Read
                        UnmappedReadsCount += 1
                    elif MCount == 1:
                        #Singly mapped Read
                        PadLongerThanSeed = False
                        for i in Index[1]:
                            if len(line[i]) > int(cfg.Seed):
                                PadLongerThanSeed = True
                            else:
                                pass
                        if PadLongerThanSeed:
                            #Mapable region in Seed
                            LongPadcount += 1
                            Padcount += 1
                        else:
                            #Single non-padded Alignment
                            Padcount += 1
                    else:
                        #Multiple mappings, means either recombination, insertion, or substitution.
                        n=0
                        UnRec = ''
                        for i in Index[0][:-1]:
                                Donor = line[i]
                                Ref = Donor
                                if "RevStrand" in line[i+1]:
                                    DonorSite = line[i+1].split("_")[2]
                                else:
                                    DonorSite = line[i+1].split("_")[1]
                                MappingStartPos = line[i+1].split("_")[0]
                                if "RevStrand" in line[i+1]:
                                    MappedReadData = cfg.Genes[str(Ref)][int(DonorSite)-1:int(MappingStartPos)]
                                    MappedReadData = Rev_Comp(MappedReadData)
                                else:   
                                    MappedReadData = cfg.Genes[str(Ref)][int(MappingStartPos)-1:int(DonorSite)]
                                if Index[0][n+1] == i + 2:
                                        #Recombination Event
                                        Acceptor = line[i+2]
                                        if "RevStrand" in line[i+3]:
                                                AcceptorSite = line[i+3].split("_")[0]
                                        else:
                                                AcceptorSite = line[i+3].split("_")[0]
                                        if Donor == Acceptor and "_RevStrand" in Donor and fabs(int(DonorSite) - int(AcceptorSite) - 1) <= cfg.MicroInDel_Length:
                                            if int(DonorSite) - int(AcceptorSite) - 1 < 0:
                                                    #MicroInsertion on negative strand
                                                    uCount += 1
                                                    DonorA = Donor[:-10]
                                                    Insertion = cfg.Genes[DonorA][int(DonorSite) - 1:int(AcceptorSite)]
                                                    Insertion = Rev_Comp(Insertion)
                                                    NewAcceptorSite = str(int(DonorSite) - 1)
                                                    AddInsToDict(Donor, DonorSite, NewAcceptorSite, Insertion, uInsDicts, ReadName)
                                            elif int(DonorSite) - int(AcceptorSite) - 1 > 0:
                                                    #MicroDeletion on negative strand
                                                    uCount += 1
                                                    x = AddToDict(Donor, Acceptor, DonorSite, AcceptorSite, uDelDicts, ReadName)
                                        elif Donor == Acceptor and "_RevStrand" not in Donor and fabs(int(DonorSite) - int(AcceptorSite) + 1) <= cfg.MicroInDel_Length:
                                            if int(DonorSite) - int(AcceptorSite) + 1 > 0:
                                                    #MicroInsertion
                                                    uCount += 1
                                                    DonorA = Donor
                                                    Insertion = cfg.Genes[DonorA][int(AcceptorSite) - 1:int(DonorSite)]
                                                    NewAcceptorSite = str(int(DonorSite) + 1)
                                                    AddInsToDict(Donor, DonorSite, NewAcceptorSite, Insertion, uInsDicts, ReadName)
                                            elif int(DonorSite) - int(AcceptorSite) + 1 < 0:
                                                    #MicroDeletion
                                                    uCount += 1
                                                    x = AddToDict(Donor, Acceptor, DonorSite, AcceptorSite, uDelDicts, ReadName)
                                        else:
                                                #Direct Recombination Event
                                                RecombCount += 1
                                                x = AddToDict(Donor, Acceptor, DonorSite, AcceptorSite, RecDicts, ReadName)
                                                if Donor in cfg.RefsLib1 and Acceptor in cfg.RefsLib1:
                                                        ViralRecombinationCount +=1
                                                elif Donor in cfg.RefsLib2 and Acceptor in cfg.RefsLib2:
                                                        HostRecombinationCount += 1
                                                else:
                                                        ViraltoHostRecombinationCount += 1
                                else:
                                        #Insertion between mapped Segments
                                        Acceptor = line[i+3]
                                        Insertion = line[i+2]
                                        if "RevStrand" in line[i+4]:
                                                AcceptorSite = line[i+4].split("_")[0]
                                                #Acceptor += "_RevStrand"
                                        else:
                                                AcceptorSite = line[i+4].split("_")[0]
                                        if Acceptor == Donor and "_RevStrand" in Donor and int(DonorSite) == (int(AcceptorSite) + 1):
                                            Insertion = Rev_Comp(Insertion)
                                            #Simple Insertion Event in negative strand
                                            if len(Insertion) >= cfg.MicroInDel_Length:
                                                InsCount += 1
                                                AddInsToDict(Donor, DonorSite, AcceptorSite, Insertion, InsDicts, ReadName)
                                            else:
                                                #MicroInDel on negative strand 
                                                uCount += 1
                                                AddInsToDict(Donor, DonorSite, AcceptorSite, Insertion, uInsDicts, ReadName)
                                        elif Acceptor == Donor and "_RevStrand" not in Donor and int(AcceptorSite) == (int(DonorSite) + 1):
                                            #Simple Insertion Event
                                            if len(Insertion) >= cfg.MicroInDel_Length:
                                                InsCount += 1
                                                AddInsToDict(Donor, DonorSite, AcceptorSite, Insertion, InsDicts, ReadName)
                                            else:
                                                uCount += 1
                                                AddInsToDict(Donor, DonorSite, AcceptorSite, Insertion, uInsDicts, ReadName)
                                        elif int(AcceptorSite) == (int(DonorSite) + len(Insertion) + 1) and Acceptor == Donor and "_RevStrand" not in Donor:
                                            #Direct Substitution
                                            if len(Insertion) <= cfg.Mismatches:
                                                #Mismatch, not Substitution
                                                MCount -= 1
                                                if MCount == 1:
                                                    #Singly mapped Read
                                                    PadLongerThanSeed = ''
                                                    for i in Index[1]:
                                                        if len(line[i]) >= int(cfg.Seed) or len(line[i]) >= int(cfg.Host_Seed):
                                                            PadLongerThanSeed += 'X'
                                                        else:
                                                            pass
                                                    if PadLongerThanSeed:
                                                        #Mapable region in pad
                                                        LongPadcount += 1
                                                        Padcount += 1
                                                    else:
                                                        #Single non-padded Alignment
                                                        Padcount += 1
                                                else:
                                                        pass
                                                    
                                            else:
                                                SubCount += 1
                                                AddInsToDict(Donor, DonorSite, AcceptorSite, Insertion, SubDicts, ReadName)
                                        elif int(DonorSite) == (int(AcceptorSite) + len(Insertion) + 1) and Acceptor == Donor and "_RevStrand" in Donor:
                                            #Direct Substitution on negative strand
                                            #Insertion = Rev_Comp(Insertion)
                                            if len(Insertion) <= cfg.Mismatches:
                                                #Mismatch, not Substitution
                                                MCount -= 1
                                                if MCount == 1:
                                                    #Singly mapped Read
                                                    PadLongerThanSeed = ''
                                                    for i in Index[1]:
                                                        if len(line[i]) >= int(cfg.Seed) or len(line[i]) >= int(cfg.Host_Seed):
                                                            PadLongerThanSeed += 'X'
                                                        else:
                                                            pass
                                                    if PadLongerThanSeed:
                                                        #Mapable region in pad
                                                        LongPadcount += 1
                                                        Padcount += 1
                                                    else:
                                                        #Single non-padded Alignment
                                                        Padcount += 1
                                                else:
                                                        pass
                                            else:
                                                    SubCount += 1
                                                    AddInsToDict(Donor, AcceptorSite, DonorSite, Insertion, SubDicts, ReadName)
                                        else:
                                                if len(Insertion) >= int(cfg.Seed) or len(Insertion) >= int(cfg.Host_Seed):
                                                    #Mapable Insertion/Recombination
                                                    UnknownRecombinationCount += 1
                                                    UnRec = 'Y'
                                                else:
                                                    if cfg.Compound_Handling and len(Insertion) > int(cfg.Compound_Handling) and Donor == Acceptor and Donor in cfg.RefsLib1:
                                                        #Compound Recombination
                                                        CompoundCount += 1
                                                        CompTest = Compound_Handling_Script(Donor, DonorSite, Insertion, AcceptorSite, uDelDicts, RecDicts, ReadName)
                                                        if CompTest == "HIT":
                                                            RecombCount += 2
                                                            ViralRecombinationCount += 2
                                                            CompoundCount += 1
                                                        else:
                                                            #Unknown Compound
                                                            UnRec = 'Y'
                                                            UnknownRecombinationCount += 1
                                                    else:
                                                        #Unknown Insertion.
                                                        UnknownRecombinationCount += 1
                                                        UnRec = 'Y'
                                n+=1 
                        if UnRec:
                                pass#UnknownRecombinations.write(wholeline)
                        else:
                                pass
                    wholeline = InRecombs.readline()        
        #Output Files for each type of event
        VirusSubstitutions = open(cfg.Output_Dir + cfg.FileTag + "Virus_Substitutions.txt","w")
        VirusInsertions = open(cfg.Output_Dir + cfg.FileTag + "Virus_Insertions.txt","w")
        VirusRecs = open(cfg.Output_Dir + cfg.FileTag + "Virus_Recombination_Results.txt","w")
        if cfg.MicroInDel_Length > 0:
                VirusuDels = open(cfg.Output_Dir + cfg.FileTag + "Virus_MicroDeletions.txt","w")
                VirusuIns = open(cfg.Output_Dir + cfg.FileTag + "Virus_MicroInsertions.txt","w")
        else:
                VirusuDels = False
                VirusuIns = False
        if cfg.Lib2:
                HostSubstitutions = open(cfg.Output_Dir + cfg.FileTag + "Host_Substitutions.txt","w")
                HostInsertions = open(cfg.Output_Dir + cfg.FileTag + "Host_Insertions.txt","w")
                HostRecs = open(cfg.Output_Dir + cfg.FileTag + "Host_Recombination_Results.txt","w")
                if cfg.MicroInDel_Length > 0:
                        HostuDels = open(cfg.Output_Dir + cfg.FileTag + "Host_MicroDeletions.txt","w")
                        HostuIns = open(cfg.Output_Dir + cfg.FileTag + "Host_MicroInsertions.txt","w")
                else:
                        HostuDels = False
                        HostuIns = False
                VirustoHostRecs = open(cfg.Output_Dir + cfg.FileTag + "Virus-to-Host_Recombination_Results.txt","w")
        else:
                 HostSubstitutions = False
                 HostInsertions = False
                 HostRecs = False
                 HostuDels = False
                 HostuIns = False
                 VirustoHostRecs = False
        if cfg.BED:
                #Create optional BED files.
                VirusRecs_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_Recombination_Results.bed","w")
                VirusRecs_BED.write('track name=Virus_Recombinations description="Virus_Recombinations" graphType=junctions\n')                
                VirusFusions_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_Fusions.BEDPE","w")
                VirusFusions_BED.write('track name=Virus_Fusions description="Virus_Fusions" graphType=BEDPE\n')
                if cfg.MicroInDel_Length > 0:
                        VirusuRecs_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus_MicroRecombinations.bed","w")
                        VirusuRecs_BED.write('track name=Virus_MicroInDels description="Virus_MicroInDels" graphType=junctions\n')
                else:
                        pass
                if cfg.Lib2:
                        HostRecs_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Host_Recombination_Results.bed","w")
                        HostRecs_BED.write('track name=Host_Recombinations description="Host_Recombinations" graphType=junctions\n')
                        HostFusions_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Host_Fusions.BEDPE","w")
                        HostFusions_BED.write('track name=Host_Fusions description="Host_Fusions" graphType=BEDPE\n')
                        VirustoHostRecs_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Virus-to-Host_Recombinations.BEDPE","w")
                        VirustoHostRecs_BED.write('track name=Virus-to-Host_Fusions description="Virus-to-Host_Fusions" graphType=BEDPE\n')
                        if cfg.MicroInDel_Length > 0:
                                HostuRecs_BED = open(cfg.Output_Dir + 'BED_Files/' + cfg.FileTag + "Host_MicroRecombinations.bed","w")
                                HostuRecs_BED.write('track name=Host_MicroInDels description="Host_MicroInDels" graphType=junctions\n')
                        else:
                                pass
                else:
                        pass
        ##      Take final Dictionaries of recombination events are write out to files
        print("Writing sorted results to individual output files...")
        WriteFinalDict(RecDicts, 'Recs')
        if cfg.MicroInDel_Length > 0:
                WriteFinalDict(uDelDicts, 'uDel')
                WriteFinalDict(uInsDicts, 'uIns')
        else:
                pass
        WriteFinalDict(InsDicts, 'Ins')
        WriteFinalDict(SubDicts, 'Sub')
        
        ##      Print summary
        print("---------------------------------------------------------------------------------------------------------------------")
        print("Total of %s reads have been analysed:" % Totalcount)
        print("%s were single mapping reads with pads, %s of which were longer than the chosen seed (%s nts)." % (Padcount, LongPadcount, cfg.Seed))
        print("%s Straight-forward Recombination Events detected"% RecombCount)
        print("of which %s were Viral Recombinations, %s were Host Recombinations and %s were Virus-to-Host Recombinations" % (ViralRecombinationCount, HostRecombinationCount, ViraltoHostRecombinationCount))
        if cfg.MicroInDel_Length > 0:
                print("%s were MicroIndels below a threshold of less than or equal to %s nucleotides." % (uCount, cfg.MicroInDel_Length))
        else:
                pass
        print("%s UnIdentified Insertion Events." % InsCount)
        print("%s Nucleotide Subsitution events, including mismatches that preserve the gene length." % SubCount)
        if cfg.Compound_Handling:
                print("%s Compound Recombination Events detected." % CompoundCount)
        else:
                pass
        print("%s events were Unknown or Ambiguous Recombination Events." % UnknownRecombinationCount)
        print("%s reads were completely unmapped." % UnmappedReadsCount)

        #Close all output files and finish
        InRecombs.close()
        #UnknownRecombinations.close()
        VirusRecs.close()
        if cfg.MicroInDel_Length > 0:
                VirusuIns.close()
                VirusuDels.close()
        else:
                pass
        VirusSubstitutions.close()
        VirusInsertions.close()
        if cfg.Lib2:
                HostRecs.close()
                VirustoHostRecs.close()
                HostInsertions.close()
                if cfg.MicroInDel_Length > 0:
                        HostuDels.close()
                        HostuIns.close()
                else:
                        pass
                HostSubstitutions.close()
        else:
                pass
        if cfg.BED:
                VirusRecs_BED.close()
                if cfg.MicroInDel_Length > 0:
                        VirusuRecs_BED.close()
                        VirusFusions_BED.close()
                else:
                        pass
                if cfg.Lib2:
                        HostRecs_BED.close()
                        if cfg.MicroInDel_Length > 0:
                                HostuRecs_BED.close()
                                HostFusions_BED.close()
                                VirustoHostRecs_BED.close()
                        else:
                                pass
                else:
                        pass
                BEDGraph_Plot()
        else:
                pass

##      -------------------------------------------------------------------------------------------
##      This module can be run seperately from the main ViReMa script. This may be useful when tweeking variables such
##      as the MicroInDel length or Compound_Handling, or when combining the results from multiple instances of ViReMa.
##      Consequently, the following code takes arguments from command line, and sends them to the config file for cross-module access.
##      Results Compilation is then initiated as normal.
##      -------------------------------------------------------------------------------------------

if __name__ == '__main__':
        print('\n-------------------------------------------------------------------------------------------')
        print('ViReMa_0.25 - Viral Recombination Mapper - Compilation Module')
        print('Last modified Jun 2021')
        print('-------------------------------------------------------------------------------------------')
        parser = argparse.ArgumentParser()
        parser.add_argument("Input_Data", help= "UnCompiled Results file from ViReMa run or other SAM file. Autodetectcion of .gz")
        parser.add_argument("--Header", help= "SAM header file if absent from input SAM file")
        parser.add_argument("-Overwrite", action='store_true', help= "Allow overwrite of previous ViReMa output. Will crash if you don't have permissions. Default = Off.")
        parser.add_argument("--Output_Tag", help= "Enter a tag name that will be appended to end of each output file.")
        parser.add_argument("-DeDup", action='store_true', help="Remove potential PCR duplicates. Default is off.")
        parser.add_argument("--UMI", help="Enter string/delimiter used in read name to define UMI barcode. Default is off.")
        parser.add_argument("-ReadNamesEntry", action='store_true', help="Append Read Names contributing to each compiled result. Default is off.")
        parser.add_argument("-FuzzEntry", action='store_true', help="Append Fuzz present in each recombination result. Default is off.")        
        parser.add_argument("--Defuzz", help="Choose how to defuzz data:  '5' to report at 5' end of fuzzy region, '3' to report at 3' end, or '0' to report in centre of fuzzy region. Default is no fuzz handling (similar to choosing Right - see Routh et al).")
        parser.add_argument("--MaxFuzz", help="Select maximum allowed length of fuzzy region. Recombination events with longer fuzzy regions will not be reported. Default is Seed Length.")
        parser.add_argument("--MicroInDel_Length", help= "Size of MicroInDels - these are common artifacts of cDNA preparation.  See Routh et al JMB 2012. Default size is 0)")
        parser.add_argument("--BackSplice_limit", help= "Size of Back-Splice or Duplication that is reported. Default size is 0)")
        parser.add_argument("--Compound_Handling", help= "Select this option for compound recombination event mapping (see manual for details). Enter number of nucleotides to map (must be less than Seed, and greater than number of nts in MicroInDel). Default is off.")
        parser.add_argument("--Aligner_Directory", help= "Enter a directory with aligner software.")
        parser.add_argument("--Output_Dir", help= "Enter a directory name that all compiled output files will be saved in.")
        parser.add_argument("-BED", action='store_true', help= "Output recombination data into BED files.")
        #parser.add_argument("-CoVaMa", action='store_true', help= "Make CoVaMa output data.")
        parser.add_argument("-NoViReMa", action='store_true', help= "Select if used non-ViReMa package to generate SAM file.")
        parser.add_argument("-ScrutSAM", action='store_true', help= "Select if used non-ViReMa package to generate SAM file and want to validate mapped segments are longer than requested seed length.")
        parser.add_argument("Virus_Index", help="Virus genome reference index key. e.g. FHV_Genome.txt")
        parser.add_argument("--Host_Index", help="Host genome reference index key, e.g. d_melanogaster_fb5_22")
        parser.add_argument("--Seed", help="Number of nucleotides in the Seed region. Default is 25.")
        parser.add_argument("--N", help= "Number of mismatches tolerated in mapped seed and in mapped segments. Default is 1.")
        parser.add_argument("--Host_Seed", help="Number of nucleotides in the Seed region when mapping to the Host Genome. Default is same as Seed value.")
        args = parser.parse_args()
        ## Handle options
        File1 = str(args.Input_Data)
        if args.Overwrite:
                cfg.Overwrite = True
        else:
                cfg.Overwrite = False
        if args.NoViReMa:
            cfg.NoViReMa = True
        else:
            cfg.NoViReMa = False
        if args.ScrutSAM:
            cfg.ScrutSAM = True
        else:
            cfg.ScrutSAM = False
                    
        if not cfg.NoViReMa:
            if args.Header:
                cfg.HeaderFile = str(args.Header)
                cfg.HeadRead = cfg.HeaderFile
            else:
                cfg.HeaderFile = ''
                cfg.HeadRead = File1
            if cfg.HeadRead[-3:] == '.gz':
                OPEN = gzip.open(cfg.HeadRead,'rt')
            else:
                OPEN = open(cfg.HeadRead,'r')
            with OPEN as InRecombs:
                    #Find arguments used in Mapping Phase from ViReMa.py
                    line = InRecombs.readline().split()
                    while line[0] in SamHeaders: 
                        if line[0] == '@PG':    
                            line = line[5:]
                            print(line)
                            #Find arguments used in Mapping Phase from ViReMa.py
                            Raw_Data = line[1]
                            if "--Host_Index" in line:
                                    cfg.Lib2 = line[line.index("--Host_Index")+1]
                                    #cfg.Genome2 = cfg.Lib2 + '.fa'
                            else:
                                    cfg.Lib2 = None
                            if "--Seed" in line:
                                    cfg.Seed = int(line[line.index("--Seed")+1])
                            else:
                                    cfg.Seed = 25
                            if "--Host_Seed" in line:
                                    cfg.Host_Seed = int(line[line.index("--Host_Seed")+1])
                            else:
                                    cfg.Host_Seed = cfg.Seed
                            if "--N" in line:
                                    cfg.Mismatches = int(line[line.index("--N")+1])
                            else:
                                    cfg.Mismatches = 1
                            break
                        line = InRecombs.readline().split()
        else:
            ##Parameters a re-specified here from a non-ViReMa aligner
            if args.Virus_Index:
                cfg.Lib1 = str(args.Virus_Index) 
                #cfg.Genome1 = cfg.Lib1 + '.fa'
            else:
                print('Error! Virus Index must be specific when using non-ViReMa aligner')
                print('Error! (There is no default setting here)')
            #cfg.Genome1 = cfg.Lib1 + '.fa'
            if args.Host_Index:
                    cfg.Lib2 = str(args.Host_Index)
            else:
                    cfg.Lib2 = None            
            if args.Seed:
                cfg.Seed = int(args.Seed)
            else:
                cfg.Seed = 25
            if args.Host_Seed:
                if int(args.Host_Seed) < cfg.Seed:
                        cfg.Host_Seed = cfg.Seed
                else:
                        cfg.Host_Seed = int(args.Host_Seed)
            else:
                    cfg.Host_Seed = cfg.Seed
            if args.N:
                cfg.Mismatches = int(args.N)
            else:
                cfg.Mismatches = 1
        cfg.Lib1 = str(args.Virus_Index) 
        if args.Output_Tag:
                cfg.FileTag = str(args.Output_Tag)
        else:
                cfg.FileTag = ''
        if args.Defuzz == '3':
                cfg.Defuzz = 'Right'
        elif args.Defuzz == '5':
                cfg.Defuzz = 'Left'
        elif args.Defuzz == '0':
                cfg.Defuzz = 'Centre'
        else:
                cfg.Defuzz = False
        if args.DeDup:
                cfg.DeDup = True
        else:
                cfg.DeDup = False        
        if args.UMI:
                cfg.UMI = str(args.UMI)
        else:
                cfg.UMI = False
        if args.ReadNamesEntry:
                cfg.ReadNamesEntry = True
        else:
                cfg.ReadNamesEntry = False
        if args.FuzzEntry:
                cfg.FuzzEntry = True
        else:
                cfg.FuzzEntry = False
        if args.MaxFuzz:
                cfg.MaxFuzz = int(args.MaxFuzz)
        else:
                cfg.MaxFuzz = int(cfg.Seed)
        if args.Compound_Handling:
                cfg.Compound_Handling = str(args.Compound_Handling)
        else:
                cfg.Compound_Handling = ''
        if args.MicroInDel_Length:
                cfg.MicroInDel_Length = int(args.MicroInDel_Length)
        else:
                cfg.MicroInDel_Length = 0
        if args.BackSplice_limit:
                cfg.BackSplice_limit = int(args.BackSplice_limit)
        else:
                cfg.BackSplice_limit = 0
        if args.BED:
                cfg.BED = True
        else:
                cfg.BED = False
        if args.Output_Dir:
                if not exists(str(args.Output_Dir)) or cfg.Overwrite == True:
                        cfg.Output_Dir = str(args.Output_Dir) + '/'
                        makedirs(cfg.Output_Dir)
                else:
                        print("Output Directory already exists!  Appending time to directory name to prevent overwrite.")
                        cfg.Output_Dir = str(args.Output_Dir) + str(int(time.time())) + '/'
                        makedirs(cfg.Output_Dir)
        else:
                pass
        if args.BED:
                if not exists(cfg.Output_Dir + 'BED_Files/'):
                        makedirs(cfg.Output_Dir + 'BED_Files/')
                else:
                        print("WARNING: BED Folder already present in output directory!")
        else:
                pass
        if args.Aligner_Directory:
            if cfg.Windows:
                cfg.Aligner_Directory = str(args.Aligner_Directory) + '\\'
            else:                
                cfg.Aligner_Directory = str(args.Aligner_Directory) + '/'
        else:
                cfg.Aligner_Directory = ''
        print("Finding reference gene data using Bowtie-Inspect")
        if cfg.Aligner =='bwa':
            cfg.RefsLib1, cfg.RefsLib2, cfg.Genes = ExtractRefDataBWA()
        else:
            cfg.RefsLib1, cfg.RefsLib2, cfg.Genes = ExtractRefData()
        cfg.RefsLib1_CuttingSites = {}
        for Name in cfg.RefsLib1:
            cfg.RefsLib1_CuttingSites[Name] = np.array([0]*(len(cfg.Genes[Name])+1))
        #for Name in cfg.RefsLib2:
        #    cfg.RefsLib2_CuttingSites[Name] = np.array([0]*(len(cfg.Genes[Name])+1))
        if cfg.DeDup:
                UniquifyReport(File1, 'DeDuped_' + File1)
                File1 = 'DeDuped_' + File1
        else:
                pass
        cfg.UMIs = set()
        print("Sorting Results and saving into individual outputs")
        ResultsSort(File1)

##      -------------------------------------------------------------------------------------------
##      End
##      -------------------------------------------------------------------------------------------
