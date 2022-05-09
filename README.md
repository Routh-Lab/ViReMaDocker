# **ViReMa Docker**
The above repository contains a Dockerfile for ViReMa version 0.25 with associated example data. The original ViReMa files can be downloaded [here](https://sourceforge.net/projects/virema/).

## **What is ViReMa?**
**ViReMa** is an algorithm developed by the [Routh and Johnson Labs](https://www.utmb.edu/routhlab/home) providing **"a versatile platform for rapid, sensitive and nucleotide-resolution detection of recombination junctions in viral genomes using next-generation sequencing data".**

## **Why Docker?**
**Docker** makes setting up and detecting viral recombination events with ViReMa easy. Having a Docker container handles installation of program dependencies and versioning. Additionally, containers can be run on Windows, Mac, and Linux operating systems.

## **Setup (~10 minutes)**
1. Download [Docker Desktop](https://www.docker.com/products/docker-desktop/).
2. Download this GitHub repository's files as a .zip file and unzip it as a new folder.
3. Open the command line on your operating system. The next three steps involve entering lines into the command line.
4. Change the current directory using the command line to the folder you just unzipped with the file "virema" and the folder "src". This command looks something like ```cd C:\Users\user\Documents\GitHub\ViReMaDocker``` replacing "C:\Users\user\Documents\GitHub\ViReMaDocker" with your folder path.
5. Now, build the image with ```docker build -t virema . -f virema```. If you want to add your own data first, add appropriate files to the "src/TestData" folder before building the image. Messages should indicate if everything builds smoothly. Here, the "virema" after "-t" can be replaced with whatever text you would like to name/tag the container with.
6. You can now make a container using ```docker run -d -t virema``` or by going to the "Images" tab in Docker Desktop and hitting "Run" when hovering over the "virema" image.
7. With Docker Desktop open, navigate to the "Containers/Apps" tab.
8. If the container is running, hit "CLI". If the container is stopped, hit the "Start" button first.
9. Type ```python ./ViReMa.py ./FHV_Genome_padded ./FHV_10k.txt FHV_recombinations.SAM --Seed 20 --MicroInDel_Length 5 -BED --Output_Dir FHV_Test``` into the "CLI". This command runs the ViReMa algorithm on flock house virus example data. If everything worked, it should produce a folder called "FHV_Test" with a file called "FHV_recombinations.SAM" inside. You can visualize this with the commands ```ls``` followed by ```cd FHV_Test``` and one more ```ls```.
10. Go back to the command line on your operating system. The next two steps involve entering lines into the command line.
11. With the container still running, use ```docker container ls``` and copy the container id for the "virema" container.
12. Export your files from the container to your computer. Use the command ```cd C:\Users\user\Documents\Results``` replacing "C:\Users\user\Documents\Results" with where you want the folder to show up on your local machine.  Run ```docker cp <container-id>:/FHV_Test .```, replacing "```<container-id>```" with the container id from Step 11. The '.' puts the files into the local directory we just used "cd" to get into.
13. When finished, stop the container in Docker Desktop or using the command line.

## **Analysis**
You can then create figures and visualizations with the newly created files by following tutorials [here](https://jayeung12.github.io/).

## **Details**
<details>
  <summary>Click to expand</summary>
  
## ViReMa Version 0.25
### Last Modified: Jun-21

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
</details>

<br>

**Read the updated ViReMa paper with vignettes:**<br>
Sotcheff S, Zhou Y, Sun Y, Johnson JE, Torbett B, Routh AL. ViReMa: A Virus Recombination Mapper of Next-Generation Sequencing data characterizes diverse recombinant viral nucleic acids. bioRxiv. 2022:2022.03.12.484090. doi:10.1101/2022.03.12.484090

**Read the original ViReMa paper:**<br>
Andrew Routh, John E. Johnson, Discovery of functional genomic motifs in viruses with ViReMa-a Virus Recombination Mapper-for analysis of next-generation sequencing data, Nucleic Acids Research, Volume 42, Issue 2, 1 January 2014, Page e11, https://doi.org/10.1093/nar/gkt916

<br>
