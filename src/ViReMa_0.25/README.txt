ViReMa Version 0.24
Last Modified: May-21

Test Data - FHV (8 files)
- FHV_10k.txt		Contains ten thousand reads from the Flock House Virus Dataset: SRP013296
- FHV_Genome_padded.txt	Contains reference genes for Flock House Virus with long 3' terminal A residues.
- FHV_Genome_padded.*.ebwt 	These are the built index sequences for the FHV padded genome using Bowtie-Build v0.12.9
- FHV_P7R2_rep2_100k.txt	Contains raw data from Jaworski et al PLoS Path paper
- FHV_P7R2_rep2_100k_virema.bam Contains example mapping of data using ViReMa2

Compiler_Module.py
Module or Stand-alone script used to compile output results from ViReMa.py.  Runs from Command-line.

ConfigViReMa.py
This script carries the global variables used by both ViReMa.py and Compiler_Module.py

README.txt
Includes instructions to run ViReMa.

ViReMa.py
Runs ViReMa (Viral-Recombination-Mapper) from commmand line.

ViReMa_GUI.py
Runs ViReMa (Viral-Recombination-Mapper) from GUI (requires GOOEY).





Before you Start:

ViReMa is a simple python script and so should not require any special installation.  

ViReMa requires python version 3.7 and Bowtie version 0.12.9. ViReMa only uses modules packaged as a standard with Python version 2.7. 

Bowtie and Bowtie-Inspect must be in your $PATH.

Indexes for reference genomes must be built with Bowtie-Build. For maximum sensitivity, please add a terminal pad using 'A' nucleotides to the end of your genome sequence before creating virus reference indexes using Bowtie-Build.  This pad must be longer than the length of the reads being aligned.  Without these pads, ViReMa will fail to detect recombination events occuring at the edges of the viral genome. 



ViReMa is run from the command line:

>python ViReMa.py Virus_Index Input_Data Output_Data [args]


Example using test data:

>home/ViReMa0.1/python ViReMa.py Test_Data/FHV_Genome_padded Test_Data/FHV_10k.txt FHV_recombinations.txt --Seed 20 --MicroInDel_Length 5 


ViReMa will take read data and attempt to align it to the reference genomes (Virus first, Host second). If the Seed of the read successfully aligns to a reference genome, bowtie will continue to align the remaining nucleotides after the Seed.Alignment() will extract all the successfully aligned nucleotides and the remaining unaligned nucleotides will be written to a new temporary read file. If there is no succesful alignment, Alignment() will trim one nucleotide from the beginning of the read and report. Again, the remaining nucleotides will be written to a new temporary file which will be used for subsequent alignment.



Required arguments:

  Virus_Index		

Virus genome reference. e.g. FHV_Genome.txt
Enter full path if the index is not in the current working directory, even when that index is stored in your Bowtie-0.12.9/indexes folder.  E.g.:  ../../Desktop/Bowtie-0.12.9/indexes/FHV_Genome


  Input_Data            

File containing single reads in FASTQ or FASTA format.


  Output_Data           

Destination file for results.  This is be saved in the current working directory.  
