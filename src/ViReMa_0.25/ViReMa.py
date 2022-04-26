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
##      ----------------------------------------------------------------------------------------
print('\n-------------------------------------------------------------------------------------------')
print('ViReMa Version 0.25 - Command-Line version - written by Andrew Routh')
print('Last modified Jun 2021')
print('-------------------------------------------------------------------------------------------')
##      ----------------------------------------------------------------------------------------

import time
start = time.time()
from subprocess import call
import re
from Compiler_Module import *
import ConfigViReMa as cfg
from os import makedirs, remove
from os.path import exists
import sys
import argparse
import gzip

##      -------------------------------------------------------------------------------------------
##      Take arguments from command line, and send them to the config file for cross-module access
##      -------------------------------------------------------------------------------------------

def MainArgs():
        parser = argparse.ArgumentParser()
        Output_group = parser.add_argument_group("Output Handling", "Specify where output files are put and how they are named")
        Error_group = parser.add_argument_group("Error Handling", "Customize the search parameters and error handling - described in ViReMa manuscript")
        Fuzz_group = parser.add_argument_group("Fuzz Handling", "Customize the how 'fuzzy' recombination events are handled and reported.")
        Host_group = parser.add_argument_group("Host Options", "Include host parameters")
        Runtime_group = parser.add_argument_group("Runtime Options", "Customize the runtime options")
        Alignment_group = parser.add_argument_group("Alignment Options", "Set alignment conditions")
        Compilation_group = parser.add_argument_group("Compilation Options", "Set some compilation options and extra output files")

        ##Required Options        
        parser.add_argument("-Windows", action='store_true', help= "Select this option if running ViReMa from a Windows/Cygwin shell.")
        parser.add_argument("Virus_Index", help="Virus reference genome in FASTA format. ViReMa will call bowtie to generate index. e.g. FHV_Genome.fasta")
        parser.add_argument("Input_Data", help= "File containing single reads in FASTQ format")
        parser.add_argument("-F", action='store_true', help="Select if raw data is in FASTA format.")
        parser.add_argument("Output_SAM", help= "Name of Output SAM file")
        
        ##Output Handling
        Output_group.add_argument("--Output_Tag", help= "Enter a tag name that will be appended to end of each output file.")
        Output_group.add_argument("--Output_Dir", help= "Enter a directory name that all compiled output files will be saved in.")
        Output_group.add_argument("-Overwrite", action='store_true', help= "Allow overwrite of previous ViReMa output. Will crash if you don't have permissions. Default = Off.")
        
        ##Error Handling Options        
        Error_group.add_argument("--N", help= "Number of mismatches tolerated in mapped seed and in mapped segments. Default is 1.")
        Error_group.add_argument("--Seed", help="Number of nucleotides in the Seed region. Default is 25.")
        Error_group.add_argument("--X", help="Number of nucleotides not allowed to have mismatches at 3' end and 5' of segment. \nOverrides seperate ThreePad and FivePad settings. Default is 5.")
        Error_group.add_argument("--ThreePad", help="Number of nucleotides not allowed to have mismatches at 3' end of segment. Default is 5.")
        Error_group.add_argument("--FivePad", help="Number of nucleotides not allowed to have mismatches at 5' end of segment. Default is 5.")
        Error_group.add_argument("--ErrorDensity", help="Number of mismatches allowed per N nts. \ne.g. 2,20 = 2 error within any 20nt window.")
        Error_group.add_argument("--MicroInDel_Length", help= "Size of MicroInDels - these are common artifacts of cDNA preparation. \nDefault size is 0)")
        Error_group.add_argument("--BackSplice_limit", help= "Size of Back-Splice or Duplication that is reported. \nDefault size is 0)")
        Error_group.add_argument("--Compound_Handling",  help= "Select this option for compound recombination event mapping (see manual for details). \nEnter number of nucleotides to map (must be less than Seed, and greater than number of nts in MicroInDel). Default is off.")
        Error_group.add_argument("--Internal_Pad",  help= "Enter number of nucleotides to allow within recombinatin events without creating a second segment map (must be less than Seed). Default is MicroInDel Length.")
        
        ##Fuzz Handling
        Fuzz_group.add_argument("--Defuzz", help="Choose how to defuzz data: '5' to report at 5' end of fuzzy region, \n3' to report at 3' end, or '0' to report in centre of fuzzy region. \nDefault is no fuzz handling (similar to choosing Right - see Routh et al).")
        Fuzz_group.add_argument("--MaxFuzz", help="Select maximum allowed length of fuzzy region. Recombination events with longer fuzzy regions will not be reported. Default is Seed Length.")
        Fuzz_group.add_argument("-FuzzEntry", action='store_true', help="Append Fuzz present in each recombination result. Default is off.")        
        
        ##Host Options
        Host_group.add_argument("--Host_Index", help="Host genome reference index key, e.g. d_melanogaster_fb5_22. Do not add trailing .n.ebwt")
        Host_group.add_argument("--Host_Seed", help="Number of nucleotides in the Seed region when mapping to the Host Genome. Default is same as Seed value.")
        
        ##Runtime Options        
        Runtime_group.add_argument("--MaxIters", help= "Set maximum number of iterations. Default = Length of longest read.")
        Runtime_group.add_argument("--p", help= "Enter number of available processors. Default is 1.")
            
        ##Alignment Options
        Alignment_group.add_argument("--Chunk", help= "Enter number of reads to process together. Default is 1M reads")
        Alignment_group.add_argument("--Pad",  help="Enter number of A's to add to 3' end of viral genome. 'Pads' are required if recombination occurs as end of genome. Default is off")
        Alignment_group.add_argument("--Aligner", help="Enter Alignment Software: 'bwa', 'bowtie'. Default is bowtie.")
        Alignment_group.add_argument("--Aligner_Directory", help="Specify directory containing aligner if this error: 'The system cannot find the file specified' appears early in ViReMa run")
        
        #Compilation Options
        Compilation_group.add_argument("-DeDup", action='store_true', help="Remove potential PCR duplicates. Default is 'off'.")
        Compilation_group.add_argument("--UMI", help="Enter string/delimiter used in read name to define UMI barcode. Default is off.")
        Compilation_group.add_argument("-ReadNamesEntry", action='store_true', help="Append Read Names contributing to each compiled result. Default is off.")
        Compilation_group.add_argument("-No_Compile", action='store_true', help= "Select this option if you do not wish to compile the results file into.  Maybe useful when combining results from different datasets.")
        Compilation_group.add_argument("-Only_Compile", action='store_true', help= "Select this option if you have already perform ViReMa once and wish to re-compile with different parameters.")
        Compilation_group.add_argument("-BED", action='store_true', help= "Output recombination data into BED files.")
        #Compilation_group.add_argument("-CoVaMa", action='store_true', help= "Make CoVaMa output data.")
                
        args = parser.parse_args()
        
        if args.Windows:
                cfg.Windows = True
        else:
                cfg.Windows = False
        if args.F:
            cfg.ReadType = '-f'
        else:
            cfg.ReadType = '-q'
            
        cfg.Lib1 = str(args.Virus_Index)
        bwts = ['.1.ebwt','.2.ebwt','.3.ebwt','.4.ebwt','.rev.1.ebwt','.rev.1.ebwt']
        for i in bwts:
            if i in cfg.Lib1:
                cfg.Lib1 = cfg.Lib1.split(i)[0]
            else:
                pass

        if args.Host_Index:
                cfg.Lib2 = str(args.Host_Index)
        else:
                cfg.Lib2 = None
        cfg.File1 = str(args.Input_Data)
        if cfg.Lib2: 
            for i in bwts:
                if i in cfg.Lib1:
                    cfg.Lib1 = cfg.Lib1.split(i)[0]
                else:
                    pass
        else:
            pass
        cfg.Working_Directory = ''
        if args.Overwrite:
                cfg.Overwrite = True
        else:
                cfg.Overwrite = False
        if args.Output_Dir:
            if not exists(cfg.Working_Directory + str(args.Output_Dir)):
                cfg.Output_Dir = cfg.Working_Directory + str(args.Output_Dir) + '/'
            else:
                if cfg.Overwrite == True:
                    cfg.Output_Dir = cfg.Working_Directory + str(args.Output_Dir) + '/'
                else:
                    print("Output Directory already exists!  Appending time to directory name to prevent overwrite.")
                    cfg.Output_Dir = cfg.Working_Directory + str(args.Output_Dir) + str(int(time.time())) + '/'
            try:
                makedirs(cfg.Output_Dir)
            except:
                print("WARNING!", cfg.Output_Dir, "is already present and -Overwrite is selected. Files therein will be overwritten.")
        else:
                cfg.Output_Dir = cfg.Working_Directory  
        if args.BED:
                if not exists(cfg.Output_Dir + 'BED_Files/'):
                        makedirs(cfg.Output_Dir + 'BED_Files/')
                else:
                        print("WARNING: BED Folder already present in output directory!")
        else:
                pass
        if not exists(str(args.Output_SAM)) or cfg.Overwrite == True:
                cfg.File3 = str(args.Output_SAM)
        else:
                print("SAM File already exists!  Appending time to directory name to prevent overwrite.")
                cfg.File3 = str(args.Output_SAM) + str(int(time.time()))
        if args.Seed:
                cfg.Seed = int(args.Seed)
        else:
                cfg.Seed = 25
        if args.N:
                cfg.Mismatches = int(args.N)
        else:
                cfg.Mismatches = 1
        if args.ErrorDensity:
                cfg.ErrorDensity = eval(str(args.ErrorDensity))
                cfg.EDMode = True
        else:
                cfg.EDMode = False
        if args.Compound_Handling:
                cfg.Compound_Handling = str(args.Compound_Handling)
        else:
                cfg.Compound_Handling = ''
        if args.MicroInDel_Length:
                cfg.MicroInDel_Length = int(args.MicroInDel_Length)
        else:
                cfg.MicroInDel_Length = 0
        if args.Internal_Pad:
                cfg.Internal_Pad = str(args.Internal_Pad)
        else:
                cfg.Internal_Pad = 0#cfg.MicroInDel_Length
        if args.BackSplice_limit:
                cfg.BackSplice_limit = int(args.BackSplice_limit)
        else:
                cfg.BackSplice_limit = cfg.MicroInDel_Length
        if args.ThreePad:
                cfg.ThreePad = int(args.ThreePad)
        else:
                cfg.ThreePad = 5
        if args.FivePad:
                cfg.FivePad = int(args.FivePad)
        else:
                cfg.FivePad = 5
        if args.X:
                cfg.FivePad = int(args.X)
                cfg.ThreePad = int(args.X)
        else:
                pass
        if args.p:
                cfg.Threads = str(args.p)
        else:
                cfg.Threads = '1'
        if args.Output_Tag:
                cfg.FileTag = str(args.Output_Tag) + "_"
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
        if args.MaxFuzz:
                cfg.MaxFuzz = int(args.MaxFuzz)
        else:
                cfg.MaxFuzz = cfg.Seed
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
        if args.Aligner_Directory:
#            if cfg.Windows:
#                cfg.Aligner_Directory = str(args.Aligner_Directory) + '\\'
#            else:                
            cfg.Aligner_Directory = str(args.Aligner_Directory) + '/'
        else:
                cfg.Aligner_Directory = ''         
        if args.Aligner == 'bwa':
                cfg.Aligner = 'bwa'
        else:
                cfg.Aligner = 'bowtie'
        if args.Chunk:
                cfg.Chunk = str(args.Chunk)
        else:
                cfg.Chunk = '1000000' #False
        if args.MaxIters:
                cfg.MaxIters = int(args.MaxIters)
        else:
                cfg.MaxIters = 10000
        if args.Host_Seed:
                if int(args.Host_Seed) < cfg.Seed:
                        cfg.Host_Seed = cfg.Seed
                else:
                        cfg.Host_Seed = int(args.Host_Seed)
        else:
                cfg.Host_Seed = cfg.Seed
        if args.Only_Compile:
                cfg.Map = False
        else:
                cfg.Map = True
        if args.No_Compile:
                cfg.Compile = False
        else:
                cfg.Compile = True
        if args.BED:
                cfg.BED = True
        else:
                cfg.BED = False        
#        if args.CoVaMa:
#                cfg.CoVaMa = True
#        else:
#                cfg.CoVaMa = False
                
##      ----------------------------------------------------------------------------------------
##      Function Countreads will determine the number of complete reads in the given input file.
##      ----------------------------------------------------------------------------------------

def Countreads(File, ReadType):
    if File[-3:] == '.gz':
        OPEN = gzip.open(File, 'rt')
    else:
        OPEN = open(File, 'r')
    with OPEN as CountReadsIn:
        NumberofReads = 0
        for line in CountReadsIn:
            if ReadType == "Q":
                NumberofReads += 0.25
            elif ReadType == "F":
                NumberofReads += 0.5
    return int(NumberofReads)

def MakeReadDict(FILE):
    global ReadDict
    if FILE[-3:] == '.gz':
        OPEN = gzip.open(cfg.Output_Dir + FILE, 'rt')
    else:
        OPEN = open(cfg.Output_Dir + FILE, 'r')
    if cfg.ReadType == '-f':
        with OPEN as IN:
            Name = IN.readline().rstrip()
            while Name:
                Seq = IN.readline().rstrip()
                Quals = 'D' * (len(Seq))
                ReadDict[Name.split()[0][1:]] = [Seq, Quals]
                Name = IN.readline().rstrip()        
    else:
        with OPEN as IN:
            Name = IN.readline().rstrip()
            while Name:
                Seq = IN.readline().rstrip()
                IN.readline()
                Quals = IN.readline().rstrip()
                Name = Name.split()[0][1:]
                if Name in ReadDict:
                    print('!!ERROR!! Duplicate read name found. Will cause erroneous mapping/fusion events.\n', Name)
                else:
                    pass#print("HERE")
                ReadDict[Name] = [Seq, Quals]
                Name = IN.readline().rstrip()

##      ----------------------------------------------------------------------------------------
##      Function Rev_Comp() will return the Reverse Complement of a given DNA string
##      ----------------------------------------------------------------------------------------

def Rev_Comp(Seq):
        Seq = Seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(Seq)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)[::-1]
            
##      ----------------------------------------------------------------------------------------
##      Find error density determines if read mapping is legitimate.
##      ----------------------------------------------------------------------------------------

Bases = ['A','G','T','C']
def FindErrorDensity(Align):
    ##Align in format e.g.: ['0', 'G', '72', 'G', '0', 'G', '1', 'C', '3', 'C', '0', 'C']
    NErrors = [Align.count(i) for i in Bases]
    N, Allowed = cfg.ErrorDensity
    ##default is (2,25)
    if sum(NErrors) <= N:
        return Align
    else:
        n=0
        while n <= len(Align) - (N*2+1):
            Mapped = 0
            for j in Align[n:N*2+n+1]:
                try:
                    Mapped += int(j)
                except:
                    Mapped += 1
            if Mapped >= Allowed:
                n+=2
            else:
                break
        #Ends on Disqualifying Nuc
        NewAlign = Align[:n+1]
    return NewAlign
       
def MakeCode(Align):
    Code = ''
    MappedLength = 0
    for i in Align:
        try: 
            MappedLength += int(i)
            Code += i + 'M'
        except:
            Code += '1X'
            MappedLength += 1
    for i in range(2, cfg.Mismatches+1)[::-1]:
        Code = Code.replace('1X' + '0M1X'*(i-1), (str(i) +'X'))
    return Code, MappedLength
         
##      ------------------------------------------------------------------------------------------------------------
##      For each read aligned and output to the temporary SAM file, the function FindReadMapping() will extract
##      data for any succesfully aligning portions of the read, and then write to a new TEMPREADS file any remaining
##      nucleotides that will be mapped in a subsequent alignment iteration. The mapping data is then summarised and
##      returned to functions Alignment() for compiling.
##      ------------------------------------------------------------------------------------------------------------

def FindReadMapping(output, CurrentSeed, Seed):
    if output[2] != '*':
        if len(output[9]) < Seed:
            #output[9] is the query sequence
            return "TOOSMALL", "U", "*", "*", output[9], "N"
        else:
            MismatchTag = [i[5:] for i in output if 'MD:Z:' in i][0]
            Align = cigar_regex.findall(MismatchTag)
            #output[12] is the default mismatches field in bowtie
            flag = bin(int(output[1]))
            #This is the bitwise FLAG from the standard .SAM format.  A flag of '4' means the read is unmapped.
            if len(flag) > 6 and flag[-5] == '1':
            #A flag of '16' means the read mapped to the reference in the reverse direction and so needs to be reverse complented to regain to the original read
                Align = Align[::-1]
                output[9] = Rev_Comp(output[9])
                Direction = '_RevStrand_'
            else:
                Direction = '_'
            if int(Align[0]) <= cfg.FivePad:
                #Here, a mismatched nucleotide has occcurred too near the 5' end of the mapped read.
                #Consequently, this read will be trimmed as if it has not aligned.
                #If there is a good mapping, it will be found in a subsequent iteration.
                Code = '%sX' % (output[9][0])
                if len(output[9][1:]) >= Seed:
                    return "NONE", Code, "*", "*", str(len(output[9][1:])), "Y", output[9], output[10],
                else:
                    return "NONE", Code, "*", "*", output[9][1:], "N", output[9], output[10]
            else:
                #Here, we find the number of mapped nucleotides including allowed mismatches (note, mismatches are not allowed at the ends of a segment as determined by the 'ThreePad' and 'FivePad' variables
                ##Decide here whether to enter Density Mode or Old mode
                if cfg.EDMode:
                    ##Learn how many mismatches to ignore, before entering Normal mode
                    Align = FindErrorDensity(Align)
                    ##mismatches allready accounted for
                    ##only need to check for 3' ends
                    while int(Align[-1]) < cfg.ThreePad:
                        Align = Align[:-2]
                        if not Align:
                            print("Error, ThreePad too long")
                            break
                    Code, MappedLength = MakeCode(Align)
                    if int(MappedLength) < CurrentSeed:
                        #After accounting for disallowed mismatches/density, the remaining mapped nucleotides are now shorter than the required Seed Length
                        #Therefore, there is no confident mapping.
                        Code = '%sX' % (output[9][0])
                        if len(output[9][1:]) >= Seed:
                                return "NONE", Code, "*", "*", str(len(output[9][1:])), "Y", output[9], output[10],
                        else:
                                return "NONE", Code, "*", "*", output[9][1:], "N", output[9], output[10]
                    else:
                        if Direction == '_RevStrand_':
                            #output[3] is the 1-based leftmost position of the clipped alignment from the .SAM format
                            Alignment = str(int(output[3]) + len(output[9]) - 1) + "_RevStrand_" + str(int(output[3]) + len(output[9]) - MappedLength)
                            #Means leftmost position of -ve mapping, plus readlength to get rightmost, to rightmost minus mapped
                        else:
                            Alignment = output[3] + "_" + str(int(output[3]) + MappedLength)
                        if len(output[9][(MappedLength):]) >= Seed:
                            #Means there are still enough unmapped nucleotide remaining after mapped section to form a new seed
                            ##Mapping in format ('SOME', '60M', 'NC_004146.1_FHV_RNA1.seq', '896_955', '', 'N')
                            return "SOME", Code, output[2], Alignment, str(len(output[9][(MappedLength):])), "Y", output[9][(MappedLength):], output[10][(MappedLength):]
                        else:
                            return "SOME", Code, output[2], Alignment, output[9][(MappedLength):], "N"
                else:
                    ##Run normal mode
                    if cfg.Mismatches >= 2 and len(Align) > 3 and int(Align[4]) >= cfg.ThreePad:
                        #Means if two mismatches are allowed, and if two mismatches are found, and if none of these mismatches are disqualifying
                        if Align[2] == '0':
                            #Means there are two adjacent but allowed mismatches
                            Code = '%sM2X%sM' % (Align[0], Align[4])
                        else:
                            #Means there are two non-adjacent and allowed mismatches
                            Code = '%sM1X%sM1X%sM' % (Align[0], Align[2], Align[4])
                        MappedLength = int(Align[0]) + int(Align[2]) + int(Align[4]) + 2
                        #Length of the three mapped sections plus the mismatches
                    elif cfg.Mismatches >= 1 and len(Align) > 1 and int(Align[2]) >= cfg.ThreePad:
                        #Means if one mismatch is allowed, and if one mismatch is found and it is not disqualifying
                        MappedLength = int(Align[0]) + int(Align[2]) + 1
                        Code = '%sM1X%sM' % (Align[0], Align[2])
                    else:
                        #Means no mismatches were found
                        MappedLength = int(Align[0])
                        Code = '%sM' % (str(MappedLength))
                    if int(MappedLength) < CurrentSeed:
                        #After accounting for disallowed mismatches, the remaining mapped nucleotides are now shorter than the required Seed Length
                        #Therefore, there is no confident mapping.
                        Code = '%sX' % (output[9][0])
                        if len(output[9][1:]) >= Seed:
                                return "NONE", Code, "*", "*", str(len(output[9][1:])), "Y", output[9], output[10],
                        else:
                                return "NONE", Code, "*", "*", output[9][1:], "N", output[9], output[10]
                    else:
                        if Direction == '_RevStrand_':
                            #output[3] is the 1-based leftmost position of the clipped alignment from the .SAM format
                            output[3] = str(int(output[3]) + len(output[9]) - 1)
                            if cfg.Mismatches >= 2 and len(Align) > 3 and int(Align[4]) >= cfg.ThreePad:
                            #Means if two mismatches are allowed, and if two mismatches are found, and if none of these mismatches are disqualifying
                                if Align[2] == '0':
                                    #Means there are two adjacent but allowed mismatches
                                    Alignment = output[3] + Direction + str(int(output[3]) - int(Align[0]) + 1) + '\t' + output[9][int(Align[0]):int(Align[0]) + 2] + '\t' + output[2] + '\t' + str(int(output[3]) - int(Align[0]) - 2) + Direction + str(int(output[3]) - MappedLength + 1)
                                    #This is for Reverse Strand. So, this means: mapping of first section + \t + identity of two allowed and adjacent mismatching nucleotides + \t + mapping of last section.
                                else:
                                    #Means there are two non-adjacent and allowed mismatches
                                    Alignment = output[3] + Direction + str(int(output[3]) - int(Align[0]) + 1) + '\t' + output[9][int(Align[0])] + '\t' + output[2] + '\t' + str(int(output[3]) - int(Align[0]) - 1) + Direction + str(int(output[3]) - int(Align[0]) - int(Align[2])) + '\t' + output[9][(int(Align[0]) + int(Align[2]) + 1)] + "\t" + output[2] + '\t' + str(int(output[3]) - int(Align[0]) - int(Align[2]) - 2) + Direction + str(int(output[3]) - MappedLength + 1)
                                    #This is for Reverse Strand. This means: mapping of first section + \t + identity of allowed mismatching nucleotide + \t + mapping of second section + \t + identity of second allowed mismatching nucleotide + \t + mapping of third section
                            elif cfg.Mismatches >= 1 and len(Align) > 1 and int(Align[2]) >= cfg.ThreePad:
                                #Means if one mismatch is allowed, and if one mismatch is found and it is not disqualifying
                                Alignment = output[3] + Direction + str(int(output[3]) - int(Align[0]) + 1) + '\t' + output[9][int(Align[0])] + '\t' + output[2] + '\t' + str(int(output[3]) - int(Align[0]) - 1) + Direction + str(int(output[3]) - MappedLength + 1)
                                #This is for Reverse Strand. So, this means: mapping of first section + \t + identity of allowed mismatching nucleotides + \t + mapping of last section.
                            else:
                                #Means no mismatches were found
                                Alignment = output[3] + Direction + str(int(output[3]) - MappedLength + 1)
                                #This is Reverse Strand. So, this means: mapping of whole read.
                        else:
                            if cfg.Mismatches >= 2 and len(Align) > 3 and int(Align[4]) >= cfg.ThreePad:
                            #Means if two mismatches are allowed, and if two mismatches are found, and if none of these mismatches are disqualifying
                                if Align[2] == '0':
                                    #Means there are two adjacent but allowed mismatches
                                    Alignment = output[3] + Direction + str(int(output[3]) + int(Align[0]) - 1) + '\t' + output[9][int(Align[0]):int(Align[0]) + 2] + '\t' + output[2] + '\t' + str(int(output[3]) + int(Align[0]) + 2) + Direction + str(MappedLength + int(output[3]) - 1)
                                    #This means: mapping of first section + \t + identity of two allowed and adjacent mismatching nucleotides + \t + mapping of last section.
                                else:
                                    #Means there are two non-adjacent and allowed mismatches
                                    Alignment = output[3] + Direction + str(int(output[3]) + int(Align[0]) - 1) + '\t' + output[9][int(Align[0])] + '\t' + output[2] + '\t' + str(int(output[3]) + int(Align[0]) + 1) + Direction + str(int(Align[0]) + int(Align[2]) + int(output[3])) + '\t' + output[9][(int(Align[0]) + int(Align[2]) + 1)] + '\t' + output[2] + '\t' + str(int(output[3]) + int(Align[0]) + int(Align[2]) + 2) + Direction + str(MappedLength + int(output[3]) - 1)
                                    #This means: mapping of first section + \t + identity of allowed mismatching nucleotide + \t + mapping of second section + \t + identity of second allowed mismatching nucleotide + \t + mapping of third section
                            elif cfg.Mismatches >= 1 and len(Align) > 1 and int(Align[2]) >= cfg.ThreePad:
                                #Means if one mismatch is allowed, and if one mismatch is found and it is not disqualifying
                                Alignment = output[3] + Direction + str(int(output[3]) + int(Align[0]) - 1) + '\t' + output[9][int(Align[0])] + '\t' + output[2] + '\t' + str(int(output[3]) + int(Align[0]) + 1) + Direction + str(MappedLength + int(output[3]) - 1)
                                #This means: mapping of first section + \t + identity of allowed mismatching nucleotides + \t + mapping of last section.
                            else:
                                #Means no mismatches were found
                                Alignment = output[3] + Direction + str(MappedLength + int(output[3]) - 1)
                                #This means: mapping of whole read.
                        if len(output[9][(MappedLength):]) >= Seed:
                            #Means there are still enough unmapped nucleotide remaining after mapped section to form a new seed
                            return "SOME", Code, output[2], Alignment, str(len(output[9][(MappedLength):])), "Y", output[9][(MappedLength):], output[10][(MappedLength):]
                        else:
                            return "SOME", Code, output[2], Alignment, output[9][(MappedLength):], "N"
    else:
        #No mapping was found for this read during this alignment
        if len(output[9]) < Seed:
            #Read was too short.
            return "TOOSMALL", "U", "*", "*", output[9], "N",
        else:
            #Code = '%sX' % (output[9][0])
            if len(output[9][1:]) >= Seed:
                #Read is still long enough for next iteration
                return "NONE", '%sX' % (output[9][0]), "*", "*", str(len(output[9][1:])), "Y", output[9], output[10]
            else:
                #Read is now too short for subsequent iterations.
                return "NONE", '%sX' % (output[9][0]), "*", "*", output[9][1:], "N", output[9], output[10]

##      ----------------------------------------------------------------------------------------
##      SAM File Handling
##      ----------------------------------------------------------------------------------------

def FindLengthMapped(CIGAR, Allowed):
    x = cigar_regex.findall(CIGAR)
    Mapped = sum([int(x[i-1]) for i, y in enumerate(x) if y in Allowed])
    return Mapped

def FindStartNuc(Coords, Ref):
    Coords = Coords.split("_")
    if Coords[1] == 'RevStrand':
        Ref += '_RevStrand'
        StartNuc = int(Coords[-1])
    else:
        StartNuc = int(Coords[0])
    return StartNuc, Ref

##      ----------------------------------------------------------------------------------------
##      SAM File Handling
##      ----------------------------------------------------------------------------------------

class ReadReport(object):
    __slots__ = ['Name', 'Segments', 'SEQ', 'QUAL', 'TAGS']    
    
    def __init__(self, Name):
        self.Name = Name
        self.Segments = []
        self.SEQ, self.QUAL = ReadDict[Name][0], ReadDict[Name][1]
        self.TAGS = []
    
    def AddSegment(self, Mapping):
            ## Example Mapping = ('SOME', '60M', 'NC_004146.1_FHV_RNA1.seq', '896_955', '', 'N')
            #Xsegment = [Nuc, Count, Code]
            #Msegment = [Mapping, StartNuc, Code] e.g. ['NC_004146.1_FHV_RNA1.seq', 298, '48M1X5M2X']
            #if self.Segments:
            try:
                ##Identify last segment:
                OldCode = cigar_regex.findall(self.Segments[-1][-1])
                if Mapping[0] == 'NONE':
                    if self.Segments[-1][-1][-1] == 'S':
                        ##Add Mismatch/trim to Softpad
                        ##Softpad Segment (e.g.) = ['NNN', 3, 'S']
                        self.Segments[-1][0] += Mapping[1][0]
                        self.Segments[-1][1] += 1
                    elif self.Segments[-1][-1][-1] == 'X':
                        ##Add Mismatch to Mismatch
                        OldCode[-2] = str(int(OldCode[-2]) + 1)
                        NewCode = ''.join(OldCode)
                        self.Segments[-1][-1] = NewCode
                    elif self.Segments[-1][-1][-1] == 'M':
                        ##Add Mismatch to Segment
                        OldCode = self.Segments[-1][-1]
                        NewCode = OldCode + '1X'
                        self.Segments[-1][-1] = NewCode
                else:
                    ##SOME mapping
                    if self.Segments[-1][-1][-1] == 'S':
                        ##Add Segment to Softpad
                        Ref = Mapping[2]
                        StartNuc, Ref = FindStartNuc(Mapping[3], Ref)
                        OldCode = str(self.Segments[-1][1]) + self.Segments[-1][-1]
                        if len(self.Segments) > 1:
                            if self.Segments[-2][-1][-1] == 'M':
                                print("Does this ever happen?", self.Name)
                                if self.Segments[-2][0] == Ref:
                                    PrevMapped = FindLengthMapped(self.Segments[-2][-1], ['M','X','N','D','P'])
                                    PreviousEndNuc = self.Segments[-2][1] + PrevMapped - 1
                                    if PreviousEndNuc == StartNuc - 1:
                                        ##THIS IS WHERE TO ADD CANONICAL INSERTIONS
                                        OldMappedCode = self.Segments[-2][-1]
                                        PadCode = self.Segments[-1][-1]
                                        PadCode = PadCode[:-1] + 'I'
                                        NewCode = OldMappedCode + PadCode + Mapping[1]
                                        self.Segments[-2][-1] = NewCode
                                        ##NOW ADD TAG
                                        Insertion = self.Segments[-1][0]
                                        self.Segments = self.Segments[:-1]
                                        if self.TAGS[len(self.Segments)]:
                                            self.TAGS[len(self.Segments)].append('XI:Z:' + Insertion)
                                        else:
                                            self.TAGS[len(self.Segments)] = ['XI:Z:' + Insertion]
                                    else:
                                        ##THIS is  complex insertion, consider entering TAG
                                        self.Segments[-1] = ([Ref, StartNuc, OldCode + Mapping[1]])
                                else:
                                    self.Segments[-1] = ([Ref, StartNuc, OldCode + Mapping[1]])
                            else:
                                self.Segments[-1] = ([Ref, StartNuc, OldCode + Mapping[1]])
                        else:
                            self.Segments[-1] = ([Ref, StartNuc, OldCode + Mapping[1]])
                    elif self.Segments[-1][-1][-1] == 'X':
                        ##Add Segment to Mismatch
                        ##Could be insertion here.
                        ##New Segment
                        Ref = Mapping[2]
                        StartNuc, Ref = FindStartNuc(Mapping[3], Ref)
                        if self.Segments[-1][0] == Ref:
                            ##SameRef
                            OldCode = cigar_regex.findall(self.Segments[-1][-1])
                            if int(OldCode[-2]) <= cfg.MicroInDel_Length:
                                #allowed insertion
                                if 'RevStrand' not in Ref:
                                    PrevMapped = FindLengthMapped(''.join(OldCode[:-2]), ['M','X','N','D','P'])
                                    LastNuc = self.Segments[-1][1] + PrevMapped - 1
                                    if LastNuc == StartNuc - 1:
                                        ##straight-forward insertion
                                        NewCode = ''.join(OldCode[:-1]) + 'I' + Mapping[1]
                                        self.Segments[-1][-1] = NewCode
                                        Insertion = ReadDict[self.Name][0][-int(OldCode[-2]):]
                                        self.TAGS.append('XI:Z:Ins:' + Insertion)
                                    else:
                                        #complex insertion/pad
                                        self.Segments.append([Ref, StartNuc, Mapping[1]])
                                else:
                                    LastNuc = self.Segments[-1][1]
                                    PrevMapped = FindLengthMapped(Mapping[1], ['M','X','N','D','P'])
                                    if LastNuc == StartNuc + PrevMapped:
                                        ##straight-forward insertion
                                        NewCode = ''.join(OldCode[:-1]) + 'I' + Mapping[1]
                                        self.Segments[-1][-1] = NewCode
                                        self.Segments[-1][1] -= PrevMapped
                                        Insertion = ReadDict[self.Name][0][-int(OldCode[-2]):]
                                        self.TAGS.append('XI:Z:Ins:' + Insertion)
                                    else:
                                        #complex insertion/pad
                                        self.Segments.append([Ref, StartNuc, Mapping[1]])                            
                            else:
                                ##Large insertion disallowed in SAMfile (will be detected by compiler)
                                self.Segments.append([Ref, StartNuc, Mapping[1]])
                        else:
                            #complex insertion/pad
                            self.Segments.append([Ref, StartNuc, Mapping[1]])
                    elif self.Segments[-1][-1][-1] == 'M':
                        ##Add Segment to Segment
                        Ref = Mapping[2]
                        StartNuc, Ref = FindStartNuc(Mapping[3], Ref)
                        OldCode = self.Segments[-1][-1]
                        if self.Segments[-1][0] == Ref:
                            PrevMapped = FindLengthMapped(OldCode, ['M','X','N','D','P'])
                            if "RevStrand" not in Ref:
                                PreviousEndNuc = self.Segments[-1][1] + PrevMapped - 1
                                if PreviousEndNuc == StartNuc - 1:
                                    #No Gap
                                    NewCode = OldCode + Mapping[1]
                                    self.Segments[-1][-1] = NewCode
                                elif PreviousEndNuc < StartNuc - 1:
                                    ##DELETION or RECOMBINATION
                                    Gap = StartNuc - PreviousEndNuc - 1
                                    if Gap > cfg.MicroInDel_Length:
                                        #Recombination/Splice
                                        NewCode = OldCode + str(Gap) + 'N' + Mapping[1]
                                    else:
                                        #MicroDeletion
                                        NewCode = OldCode + str(Gap) + 'D' + Mapping[1]
                                    self.Segments[-1][-1] = NewCode
                                else:
                                    ##SEQ Duplications/Insertions Here.
                                    Overlap = PreviousEndNuc - StartNuc + 1
                                    if Overlap <= cfg.BackSplice_limit:
                                        CurrMapping = cigar_regex.findall(Mapping[1])
                                        if Overlap >= FindLengthMapped(Mapping[1], ['M','X']):
                                            ##Overlap in mapping, but not enough mapped to confirm insertion
                                            ##Also means BackSplice_limit is likely set too high
                                            ##Treat as new segment
                                            self.Segments.append([Ref, StartNuc, Mapping[1]])
                                        else:
                                            ##Mismatch in next mapping, means insertion, not duplication
                                            RemovedCIGAR = CurrMapping[0:2]
                                            CurrMapping = CurrMapping[2:]
                                            RemovedLength = sum([int(RemovedCIGAR[i]) for i in range(0,len(RemovedCIGAR),2)])
                                            while RemovedLength <= Overlap:
                                                if CurrMapping and CurrMapping[1] == 'X':
                                                    self.TAGS.append('XI:Z:ImperfectDuplication:' + str(Overlap))
                                                    RemovedCIGAR += CurrMapping[:2]
                                                    CurrMapping = CurrMapping[2:]
                                                    RemovedLength = sum([int(RemovedCIGAR[i]) for i in range(0,len(RemovedCIGAR),2)])
                                                else:
                                                    RemovedCIGAR += CurrMapping[:2]
                                                    CurrMapping = CurrMapping[2:]
                                                    RemovedLength = sum([int(RemovedCIGAR[i]) for i in range(0,len(RemovedCIGAR),2)])
                                            NewMap = RemovedLength - Overlap
                                            if CurrMapping:
                                                NewMapping = [str(NewMap)] + ['M'] + CurrMapping
                                            else:
                                                NewMapping = [str(NewMap)] + ['M']
                                            NewMapping = ''.join(NewMapping)
                                            NewCode = OldCode + str(Overlap) + 'I' + NewMapping
                                            self.TAGS.append('XI:Z:Duplication:' + str(Overlap))
                                            self.Segments[-1][-1] = NewCode
                                    else:
                                        self.Segments.append([Ref, StartNuc, Mapping[1]])
                            elif "RevStrand" in Ref:
                                ##PrevNt should be Startnt +1
                                PreviousEndNuc = self.Segments[-1][1]
                                StartNuc = int(Mapping[3].split("_")[0])
                                if PreviousEndNuc > StartNuc + 1:
                                    ##DELETION or RECOMBINATION
                                    Gap = PreviousEndNuc - StartNuc - 1
                                    if Gap > cfg.MicroInDel_Length:
                                        #Recombination/Splice
                                        NewCode = OldCode + str(Gap) + 'N' + Mapping[1]
                                    else:
                                        #MicroDeletion
                                        NewCode = OldCode + str(Gap) + 'D' + Mapping[1]
                                    self.Segments[-1][-1] = NewCode
                                    self.Segments[-1][1] -= (FindLengthMapped(Mapping[1], ['M','X']) + Gap)
                                elif PreviousEndNuc == StartNuc + 1:
                                    #No Gap
                                    NewCode = OldCode + Mapping[1]
                                    self.Segments[-1][-1] = NewCode
                                else:
                                    ##SEQ Duplications Here.
                                    Overlap = StartNuc - PreviousEndNuc + 1
                                    if Overlap <= cfg.BackSplice_limit:
                                        CurrMapping = cigar_regex.findall(Mapping[1])
                                        if Overlap >= FindLengthMapped(Mapping[1], ['M','X']):
                                            ##Overlap in mapping, but not enough mapped to confirm insertion
                                            ##Also means BackSplice_limit is likely set too high
                                            ##Treat as new segment
                                            StartNuc = int(Mapping[3].split("_")[-1])
                                            self.Segments.append([Ref, StartNuc, Mapping[1]])
                                        else:
                                            ##Mismatch in next mapping, means insertion, not duplication
                                            RemovedCIGAR = CurrMapping[0:2]
                                            CurrMapping = CurrMapping[2:]
                                            RemovedLength = sum([int(RemovedCIGAR[i]) for i in range(0,len(RemovedCIGAR),2)])
                                            while RemovedLength <= Overlap:
                                                if CurrMapping and CurrMapping[1] == 'X':
                                                    self.TAGS.append('XI:Z:ImperfectDuplication:' + str(Overlap))
                                                    RemovedCIGAR += CurrMapping[:2]
                                                    CurrMapping = CurrMapping[2:]
                                                    RemovedLength = sum([int(RemovedCIGAR[i]) for i in range(0,len(RemovedCIGAR),2)])
                                                else:
                                                    RemovedCIGAR += CurrMapping[:2]
                                                    CurrMapping = CurrMapping[2:]
                                                    RemovedLength = sum([int(RemovedCIGAR[i]) for i in range(0,len(RemovedCIGAR),2)])
                                            NewMap = RemovedLength - Overlap
                                            if CurrMapping:
                                                NewMapping = [str(NewMap)] + ['M'] + CurrMapping
                                                NewLength = NewMap + sum([int(CurrMapping[i]) for i in range(0,len(CurrMapping),2)])
                                            else:
                                                NewMapping = [str(NewMap)] + ['M']
                                                NewLength = NewMap
                                            NewMapping = ''.join(NewMapping)
                                            NewCode = OldCode + str(Overlap) + 'I' + NewMapping
                                            self.TAGS.append('XI:Z:Duplication:' + str(Overlap))
                                            self.Segments[-1][-1] = NewCode
                                            self.Segments[-1][1] -= NewLength
                                    else:
                                        StartNuc = int(Mapping[3].split("_")[-1])
                                        self.Segments.append([Ref, StartNuc, Mapping[1]])
                        else:
                            ##Different Reference, New Segment
                            self.Segments.append([Ref, StartNuc, Mapping[1]])
            except:
                ##New Segment
                if Mapping[0] == 'NONE':
                    self.Segments.append([Mapping[1][0], 1, 'S'])
                elif Mapping[0] == 'SOME':
                    Ref = Mapping[2]
                    StartNuc, Ref = FindStartNuc(Mapping[3], Ref)
                    self.Segments.append([Ref, StartNuc, Mapping[1]])
                else:
                    pass
            try:
                ##Learn if finished
                RemainingNucs = int(Mapping[4])
            except:
                try:
                    ##READ IS FINISHED
                    RemainingNucs = len(Mapping[4])
                    if self.Segments[-1][-1][-1] == 'S':
                        ##Add Final Softpad to unmapped read
                        self.Segments[-1][0] += Mapping[1][0]
                        self.Segments[-1][1] += RemainingNucs
                    elif self.Segments[-1][-1][-1] == 'X':
                        ##Add Final Softpad to mismatch
                        OldCode = cigar_regex.findall(self.Segments[-1][-1])
                        CurrentX = int(OldCode[-2])
                        OldCode[-2] = str(CurrentX + RemainingNucs)
                        OldCode[-1] = 'S'
                        NewCode = ''.join(OldCode)
                        self.Segments[-1][-1] = NewCode
                    elif self.Segments[-1][-1][-1] == 'M':
                        ##Add Final Softpad to Segment
                        if RemainingNucs:     
                            NewCode = self.Segments[-1][-1] + str(RemainingNucs) + 'S'
                            self.Segments[-1][-1] = NewCode
                        else:
                            pass
                except:
                    pass   
                
    def __str__(self):
        return str(self.Segments) + str(self.TAGS)

##      ----------------------------------------------------------------------------------------
##      SAM File Handling
##      ----------------------------------------------------------------------------------------

class SAM_Alignment(object):
    def __init__(self, Name):
        self.Segments = []
        self.QNAME = Name
        self.SEQ, self.QUAL = ReadDict[Name][0], ReadDict[Name][1]
        self.FLAG = '4'
        self.RNAME = '*'
        self.POS = '0'
        self.MAPQ = '0'
        self.CIGAR = '*'
        self.RNEXT = '*'
        self.PNEXT = '0'
        self.TLEN = '0'
        self.TAGS = []

    def RevCigar(self, x):
        x = cigar_regex.findall(x)[::-1]
        NewCigar = ''
        for i in range(0, len(x), 2):
            NewCigar += x[1+i]
            NewCigar += x[0+i]
        self.CIGAR = NewCigar

    def AddTag(self, Tag):
        self.TAGS.append(Tag)
    
    def Output(self):
        ##Make NM tag
        x = cigar_regex.findall(self.CIGAR)
        Count = sum([int(x[i-1]) for i, y in enumerate(x) if y == 'X'])
        self.TAGS.append("NM:i:" + str(Count))
        self.TAGS = '\t'.join(sorted(self.TAGS))
        return '\t'.join([self.QNAME, self.FLAG, self.RNAME, self.POS,
                          self.MAPQ, self.CIGAR, self.RNEXT, self.PNEXT,
                          self.TLEN, self.SEQ, self.QUAL, self.TAGS]) + '\n'

    def __str__(self):
        return ' '.join([self.QNAME, self.FLAG, self.RNAME, self.POS,
                        self.MAPQ, self.CIGAR, self.RNEXT, self.PNEXT,
                        self.TLEN, self.SEQ, self.QUAL, ' '.join(sorted(self.TAGS))])

def CompleteSAMRead(Name):
    global Report
    Segments = SAMDict[Name].Segments
    Tags = SAMDict[Name].TAGS
    NSegs = len(Segments)
    n=1
    LengthMapped = ['', 0]
    for Seg in Segments:
        Result = SAM_Alignment(Name)
        if n in Tags:
            [Result.AddTag(i) for i in Tags[n]]
        if n == 1 and 'M' not in Seg[-1]:
            #unmapped:
            Report.write(Result.Output())
        else:
            #Mapping:
            Ref = Seg[0]
            Pos = Seg[1]
            Result.MAPQ = '255'
            if NSegs > 1:
                Result.AddTag('TC:i:' + str(NSegs))
                Result.AddTag('FI:i:' + str(n))
            else:
                pass
            if NSegs == n:
                ##LastSegment
                CurrentMapping = FindLengthMapped(Seg[2], ['M','X','S','I','P'])
                Result.SEQ = Result.SEQ[LengthMapped[1]:LengthMapped[1] + CurrentMapping]
                Result.QUAL = Result.QUAL[LengthMapped[1]:LengthMapped[1] + CurrentMapping]
            else:
                if Seg[-1][-1] == 'X':
                    ##Here there are unmapped nucs between two mapped segments
                    ## This warrants closer attention, e.g. Compound handling                    
                    OldCode = cigar_regex.findall(Seg[-1])
                    LastPad = int(OldCode[-2])
                    if LastPad > cfg.Internal_Pad:
                        Seg[-1] = Seg[-1][:-1] + 'S'
#                        if 'RevStrand' in Seg[0]:
#                            ## Adjustment of startsite in negative sense if using pad.
#                            Seg[1] += LastPad
#                        else:
#                            pass
                    else:
                        pass
                else:
                    pass
                CurrentMapping = FindLengthMapped(Seg[2], ['M','X','S','I','P'])
                HardPadLength = len(Result.SEQ) - LengthMapped[1] - CurrentMapping
                Seg[-1] += str(HardPadLength) + 'H'
                Result.SEQ = Result.SEQ[LengthMapped[1]:LengthMapped[1] + CurrentMapping]
                Result.QUAL = Result.QUAL[LengthMapped[1]:LengthMapped[1] + CurrentMapping]
                if '_RevStrand' in Segments[n][0]:
                    Result.RNEXT = Segments[n][0][:-10]
                else:
                    Result.RNEXT = Segments[n][0]
                Result.PNEXT = str(Segments[n][1])
            CIGAR = LengthMapped[0] + Seg[2]
            if '_RevStrand' in Ref:
                Result.RNAME = Ref[:-10]
                if n == 1:
                    Result.FLAG = '16'
                else:
                    Result.FLAG = '2064'
                Result.SEQ = Rev_Comp(Result.SEQ)
                Result.QUAL = Result.QUAL[::-1]
                Result.RevCigar(CIGAR)
                Result.POS = str(Pos)
            else:
                Result.RNAME = Ref
                if n == 1:
                    Result.FLAG = '0'
                else:
                    Result.FLAG = '2048'
                Result.CIGAR = CIGAR
                Result.POS = str(Pos)
            LengthMapped[1] += CurrentMapping
            LengthMapped[0] = str(LengthMapped[1]) + 'H'
            Report.write(Result.Output())
            n+=1

##      ----------------------------------------------------------------------------------------
##      Function AddToReportDict() adds details of any alignments, mismatches or trimmed nucleotide to temporary dictionary.
##      ----------------------------------------------------------------------------------------

def AddToReportDict(Name, Mapping, Iter):
    ## Example Mapping = ('SOME', '60M', 'NC_004146.1_FHV_RNA1.seq', '896_955', '', 'N')
    global SAMDict
    if Name not in SAMDict:
            SAMDict[Name] = ReadReport(Name)
            SAMDict[Name].AddSegment(Mapping)
    else:
            SAMDict[Name].AddSegment(Mapping)
    if Mapping[5] == 'N' or Iter > cfg.MaxIters:
        CompleteSAMRead(Name)
        del SAMDict[Name]
        del ReadDict[Name]
    else:
        pass

##      ----------------------------------------------------------------------------------------
##      Function Alignment() will take read data and attempt to align it to the reference genomes (Virus first, Host second).
##      Bowtie must be in your $PATH.
##      If the Seed of the read successfully aligns to a reference genome, bowtie will continue to align the remaining nucleotides after the Seed.
##      Alignment() will extract all the successfully aligned nucleotides and the remaining unaligned nucleotides will be written to a new temporary read file.
##      If there is no succesful alignment, Alignment() will trim one nucleotide from the beginning of the read and report.
##      Again, the remaining nucleotides will be written to a new temporary file which will be used for subsequent alignmen.
##      ----------------------------------------------------------------------------------------

def Alignment(ReadsIn, ReadType, Seed, Iter):
        global SamHeaderSet
        ####  Entries are annotated accordingly to the sequence read name, therefore read names MUST be unique.  Take care when using paired-end reads.
        SamHeaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']
        #Run Bowtie/BWA using Virus Genome.   Bowtie/BWA must be in your PATH.  Remove the --mm options if operating in Windows or cygwin.
        if cfg.Aligner == 'bwa':
            with open(cfg.Output_Dir + 'TEMPSAI1', 'w') as outfilesai:
                call([cfg.Aligner_Directory + 'bwa', 'aln', '-k', str(cfg.Mismatches), '-l', str(Seed), '-n', '10000', '-o', '0', '-t', cfg.Threads, cfg.Lib1, ReadsIn], stdout = outfilesai)
            with open(cfg.Output_Dir + 'TEMPSAM1', 'w') as outfilesam:
                call([cfg.Aligner_Directory + 'bwa', 'samse', cfg.Lib1, cfg.Output_Dir + 'TEMPSAI1', ReadsIn], stdout = outfilesam)
        else:
            if cfg.Windows:
                call([cfg.Aligner_Directory + 'bowtie', ReadType, '-n', str(cfg.Mismatches), '-l', str(Seed), '-e', '100000', '--quiet', '-p', cfg.Threads, '--best', '-S', cfg.Lib1, ReadsIn, cfg.Output_Dir + 'TEMPSAM1'])
            else:
                call([cfg.Aligner_Directory + 'bowtie', ReadType, '-n', str(cfg.Mismatches), '-l', str(Seed), '-e', '100000', '--quiet', '--mm', '-p', cfg.Threads, '--best', '-S', cfg.Lib1, ReadsIn, cfg.Output_Dir + 'TEMPSAM1'])
        NumHostReads = 0
        with open(cfg.Output_Dir + 'TEMPSAM1','r') as SAMIN1:
            TempReads = open(cfg.Output_Dir + 'TEMPREADS', 'w')
            if cfg.Lib2:
                #If using two genomes, open a new file to write any reads that did not map to virus genome
                HostAttemptReads = open(cfg.Output_Dir + "TEMPREADS2", "w")
            else:
                pass
            for line in SAMIN1:
                line = line.split('\t')
                Name = line[0]
                if Name[:3] in SamHeaders:
                    if Name[:3] == '@SQ':
                        SamHeaderSet.add('\t'.join(line))
                    else:
                        pass
                else:
                    Mapping = FindReadMapping(line, Seed, Seed)
                    if Mapping[0] == "NONE":
                        #No mapping to virus genone was found.  If a host genome is provided, write read out to new tempread file for an extra alignment
                        if cfg.Lib2:
                            if Mapping[5] == "Y" and Iter < cfg.MaxIters:
                            #'Y' is just a tag to say that the read has enough nucleotides remaining to be used in subsequent iterations.
                                if int(Mapping[4]) >= cfg.Host_Seed:
                                    #Only write read to tempfile used for host alignment is that read is longer than the chosen Host_Seed length, otherwise, skip host alignment
                                    #HostAttemptReads.write("@" + str(Name) + "\n" + str(Mapping[6]) + "\n+\n" +str(Mapping[7])+ "\n" )
                                    HostAttemptReads.write("@%s\n%s\n+\n%s\n" % (str(Name), str(Mapping[6]), str(Mapping[7])) )
                                    NumHostReads += 1
                                else:
                                #Proceed to next iteration without host mapping and trim first nucleotide
                                    #TempReads.write("@" + str(Name) + "\n" + str(Mapping[6][1:]) + "\n+\n" + str(Mapping[7][1:])+ "\n" )
                                    TempReads.write("@%s\n%s\n+\n%s\n" % (str(Name), str(Mapping[6][1:]), str(Mapping[7][1:])) )
                                    AddToReportDict(Name, Mapping, Iter)
                            else:
                            #'N' is just a tag to sat that the read is too short for subsequent iterations and so will not be written to the temp read file.
                                AddToReportDict(Name, Mapping, Iter)
                        else:
                            if Mapping[5] == "Y" and Iter < cfg.MaxIters:
                                #Proceed to next iteration without host mapping and trim first nucleotide
                                #TempReads.write("@" + str(Name) + "\n" + str(Mapping[6][1:]) + "\n+\n" + str(Mapping[7][1:])+ "\n" )
                                TempReads.write("@%s\n%s\n+\n%s\n" % (str(Name), str(Mapping[6][1:]), str(Mapping[7][1:])) )
                            else:
                                pass
                            AddToReportDict(Name, Mapping, Iter)
                    else:
                        if Mapping[5] == "Y" and Iter < cfg.MaxIters:
                            #Proceed to next iteration without host mapping.
                            #TempReads.write("@" + str(Name) + "\n" + str(Mapping[6]) + "\n+\n" + str(Mapping[7])+ "\n" )
                            TempReads.write("@%s\n%s\n+\n%s\n" % (str(Name), str(Mapping[6]), str(Mapping[7])) )
                        else:
                            pass
                            #'N' is just a tag to sat that the read is too short for subsequent iterations and so will not be written to the temp read file.
                        AddToReportDict(Name, Mapping, Iter)
        if cfg.Lib2:
            HostAttemptReads.close()
            if NumHostReads > 0:
                #Run Bowtie/BWA using Host Genome.   Bowtie/BWA must be in your PATH.  Remove the --mm and -p options if operating in Windows or cygwin.  With a large host genome, this will dramatically increase runtime as the reference genome will have to be loaded into temporary memory with each iteration.
                if cfg.Aligner == 'bwa':
                    with open(cfg.Output_Dir + 'TEMPSAI2', 'w') as outfilesai:
                        call([cfg.Aligner_Directory + 'bwa', 'aln', '-k', str(cfg.Mismatches), '-l', str(cfg.Host_Seed), '-n', '10000', '-o', '0', '-t', cfg.Threads, cfg.Lib2, ReadsIn], stdout = outfilesai)
                    with open(cfg.Output_Dir + 'TEMPSAM2', 'w') as outfilesam:
                        call([cfg.Aligner_Directory + 'bwa', 'samse', cfg.Lib2, cfg.Output_Dir + 'TEMPSAI2', ReadsIn], stdout = outfilesam)
                else:
                    if cfg.Windows:
                        call([cfg.Aligner_Directory + 'bowtie', '-q', '-n', str(cfg.Mismatches), '-l', str(cfg.Host_Seed), '-e', '100000', '--quiet', '-p', cfg.Threads, '--best', '-S', cfg.Lib2, cfg.Output_Dir + 'TEMPREADS2', cfg.Output_Dir + 'TEMPSAM2'])
                    else:
                        call([cfg.Aligner_Directory + 'bowtie', '-q', '-n', str(cfg.Mismatches), '-l', str(cfg.Host_Seed), '-e', '100000', '--quiet', '--mm', '-p', cfg.Threads, '--best', '-S', cfg.Lib2, cfg.Output_Dir + 'TEMPREADS2', cfg.Output_Dir + 'TEMPSAM2'])
                with open(cfg.Output_Dir + 'TEMPSAM2', 'r') as SAMIN2:
                    for line in SAMIN2:
                        line = line.split('\t')
                        Name = line[0]
                        if Name[:3] in SamHeaders:
                            if Name[:3] == '@SQ':
                                SamHeaderSet.add('\t'.join(line))
                            else:
                                pass
                        else:
                            Mapping = FindReadMapping(line, cfg.Host_Seed, Seed)
                            if Mapping[0] != "NONE":
                                if Mapping[5] == "Y" and Iter < cfg.MaxIters:
                                    TempReads.write("@%s\n%s\n+\n%s\n" % (str(Name), str(Mapping[6]), str(Mapping[7])) )
                                else:
                                    pass
                                AddToReportDict(Name, Mapping, Iter)
                            else:
                                if Mapping[5] == "Y" and Iter < cfg.MaxIters:
                                    TempReads.write("@%s\n%s\n+\n%s\n" % (str(Name), str(Mapping[6][1:]), str(Mapping[7][1:])) )
                                else:
                                    pass
                                AddToReportDict(Name, Mapping, Iter)
            else:
                    pass
        TempReads.close()

##      ----------------------------------------------------------------------------------------
##      Function IterateAlignments() will call Bowtie and beging alignments starting with supplied input file
##      and continuing with the generated TEMPREAD files until there no reads left.
##      If providing a FASTA file, this get converted to a FASTQ file with uniform quality
##      scores of PHRED == 30 after the first Bowtie iteration
##      ----------------------------------------------------------------------------------------

def IterateAlignments(File):
        #MakeReadDict(File)
        if cfg.ReadType == '-f':
                print("%s reads in input file." % Countreads(cfg.Output_Dir + File, 'F'))
                Alignment(cfg.Output_Dir + File, '-f', cfg.Seed, 1)
        else:
                print("%s reads in input file." % Countreads(cfg.Output_Dir + File, 'Q'))
                Alignment(cfg.Output_Dir + File, '-q', cfg.Seed, 1)
        ReadsRemaining = Countreads(cfg.Output_Dir + 'TEMPREADS', 'Q')
        print("%s reads remaining to be aligned after first iteration." % (ReadsRemaining))
        Iter = 1
        while ReadsRemaining > 0 and Iter < cfg.MaxIters:
            Iter += 1
            if Iter == cfg.MaxIters:
                    print("ViReMa will quit early as max Iterations was set at", cfg.MaxIters)
            else:
                    pass
            Alignment(cfg.Output_Dir + 'TEMPREADS', '-q', cfg.Seed, Iter)
            ReadsRemaining = Countreads(cfg.Output_Dir + 'TEMPREADS', 'Q')
            print("%s reads remaining to be aligned after %s iterations." % (ReadsRemaining, Iter))

##      ----------------------------------------------------------------------------------------
##      RUN MAIN SCRIPT
##      ----------------------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    MainArgs()
    if cfg.Aligner_Directory:
        print("Aligner Directory =", cfg.Aligner_Directory)
    else:
        print("Aligner Directory must be in PATH/Registry or be assigned using --Aligner_Directory.")
    print("Output Directory =", cfg.Output_Dir)
    if not exists(cfg.Lib1 + '.1.ebwt') and cfg.Aligner == 'bowtie':
        call([cfg.Aligner_Directory + 'bowtie-build', cfg.Lib1, cfg.Lib1])
    elif not exists(cfg.Lib1 + '.amb') and cfg.Aligner == 'bwa':
        call([cfg.Aligner_Directory + 'bwa', 'index', cfg.Lib1])
    else:
        print("Virus Index Found")
    cigar_regex = re.compile(r"[^\W\d_]+|\d+")
    SamHeaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']
    SamHeaderSet = set()
    cfg.UMIs = set()
    if cfg.Map:
        if cfg.Chunk:
            Report = gzip.open(cfg.Output_Dir + cfg.File3 + '_temp', "wt")
            if cfg.File1[-3:] == '.gz':
                OPEN = gzip.open(cfg.File1, 'rt')
            else:
                OPEN = open(cfg.File1, 'r')
            with OPEN as MAINFILE:
                ChunkNum = 0
                Read = MAINFILE.readline()
                while Read:
                    ChunkedReads = open(cfg.Output_Dir + 'ChunkedReads', 'w')
                    ReadNum = 0
                    while ReadNum < int(cfg.Chunk):
                        if Read:
                            if cfg.ReadType == '-f':
                                ChunkedReads.write(Read)
                                ChunkedReads.write(MAINFILE.readline())
                            else:
                                ChunkedReads.write(Read)
                                ChunkedReads.write(MAINFILE.readline())
                                ChunkedReads.write(MAINFILE.readline())
                                ChunkedReads.write(MAINFILE.readline())
                            Read = MAINFILE.readline()
                            ReadNum += 1
                        else:
                            ReadNum += 1
                    ChunkedReads.close()
                    ChunkNum += 1
                    SAMDict = {}
                    ReadDict = {}
                    MakeReadDict('ChunkedReads')
                    print("Beginning alignments on Chunk Number %s" % ChunkNum)
                    IterateAlignments('ChunkedReads')
                    print("Appending Results from Chunk Number %s to: " % ChunkNum, str(cfg.File3))
                    if ChunkNum == 1:
                        if cfg.File3[-3:] == '.gz':
                            Header = gzip.open(cfg.Output_Dir + cfg.File3, "wt")
                        else:
                            Header = open(cfg.Output_Dir + cfg.File3, "w")
                       # Header = gzip.open(cfg.Output_Dir + cfg.File3, "w")
                        [Header.write(i) for i in SamHeaderSet]
                        Header.write('@PG\tID:ViReMa\tPN:ViReMa\tVN:0.14\tCL:' + str(' '.join(sys.argv)) + '\n')
            Report.close()
        else:
                SAMDict = {}
                ReadDict = {}
                MakeReadDict(cfg.File1)
                print("Beginning alignments")
                Report = gzip.open(cfg.Output_Dir + cfg.File3 + '_temp', "wt")
                IterateAlignments(cfg.File1)
                Report.close()
                print("Reporting Results to: ", str(cfg.Output_Dir + cfg.File3))
                if cfg.File3[-3:] == '.gz':
                    Header = gzip.open(cfg.Output_Dir + cfg.File3, "wt")
                else:
                    Header = open(cfg.Output_Dir + cfg.File3, "w")
                #Header = gzip.open(cfg.Output_Dir + cfg.File3, "w")
                [Header.write(i) for i in SamHeaderSet]
                Header.write('@PG\tID:ViReMa\tPN:ViReMa\tVN:0.14\tCL:' + str(' '.join(sys.argv)) + '\n')
        In = gzip.open(cfg.Output_Dir + cfg.File3 + '_temp', "rt")
        for line in In:
            Header.write(str(line))
        In.close()
        Header.close()
        remove(cfg.Output_Dir + cfg.File3 + '_temp')
        remove(cfg.Output_Dir + 'ChunkedReads')
        remove(cfg.Output_Dir + 'TEMPREADS')
        remove(cfg.Output_Dir + 'TEMPSAM1')
        if cfg.Aligner == 'bwa':
            remove(cfg.Output_Dir + 'TEMPSAI1')
        else:
            pass
        if cfg.Lib2:
            remove(cfg.Output_Dir + 'TEMPREADS2')
            remove(cfg.Output_Dir + 'TEMPSAM2')
            if cfg.Aligner == 'bwa':
                remove(cfg.Output_Dir + 'TEMPSAI2')
            else:
                pass
        else:
            pass
    else:
        print("No mapping - only compiling")
    if cfg.Compile:
            if cfg.Aligner =='bwa':
                cfg.RefsLib1, cfg.RefsLib2, cfg.Genes = ExtractRefDataBWA()
            else:
                cfg.RefsLib1, cfg.RefsLib2, cfg.Genes = ExtractRefData()
            cfg.RefsLib1_CuttingSites = {}
            #cfg.RefsLib2_CuttingSites = {}
            for Name in cfg.RefsLib1:
                cfg.RefsLib1_CuttingSites[Name] = np.array([0]*(len(cfg.Genes[Name])+1))
            #for Name in cfg.RefsLib2:
             #   cfg.RefsLib2_CuttingSites[Name] = np.array([0]*(len(cfg.Genes[Name])+1))
            if cfg.DeDup:
                UniquifyReport(cfg.Output_Dir + cfg.File3, 'DeDuped_' + cfg.File3)
                cfg.File3 = 'DeDuped_' + cfg.File3
            else:
                pass
            print("Compiling Results and saving into individual outputs")
#            if cfg.CoVaMa:
#                cfg.Named_Output = open(cfg.Output_Dir + cfg.FileTag + "CoVaMa_Output.txt", "w")
#            else:
#                pass
            ResultsSort(cfg.Output_Dir + cfg.File3)
            
finish = time.time()
print("Time to complete in seconds: ", int(finish - start))
##      ----------------------------------------------------------------------------------------
##      End
##      ----------------------------------------------------------------------------------------
