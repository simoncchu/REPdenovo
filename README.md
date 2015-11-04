# REPdenovo
@2015 by Chong Chu and Yufeng Wu. This software is provided ``as is” without warranty of any
kind. In no event shall the author be held responsible for any damage resulting from the
use of this software. The program package, including source codes, executables, and this
documentation, is distributed free of charge.
A manuscript for this software is submitted to PLoS One.
If you use this program in a publication, please cite the following reference:  
Chong Chu, Rasmus Nielsen and Yufeng Wu,
REPdenovo: Inferring de novo repeat motifs from short sequence reads, PLoS One, 2015, submitted.
Please check back this site for more up-to-date information on citing this software.

## **Functionalities and Usage of REPdenovo**
REPdenovo is designed for constructing repeats directly from sequence reads. It based on the idea of frequent k-mer assembly. REPdenovo provides many
functionalities, and can generate much longer repeats than existing tools. The overall pipeline is shown in the mannual file.
REPdenovo supports the following main functionalities.

1. Assembly. This step performs k-mer counting. Then we find frequent k-mers whose frequencies are over certain threshold. We then assemble these frequent k-mers into consensus repeats (in the form of contigs). Then we merge the constructed contigs to more completeness ones.

2. Scaffolding. We use paired-end reads to connect repeat contigs into scaffolds, also provide the average coverage (indicates the copy number) for each constructed repeats.

## **Dependencies**
The current released version of REPdenovo runs on Linux OS. And REPdenovo needs the following tools to be installed in the machine you are working on.

1. Python 2.7 or higher version is required to run REPdenovo. 

2. A k-mer counting tool. REPdenovo uses Jellyfish program for performing k-mer counting. Jellyfish can be downloaded from https://github.com/gmarcais/Jellyfish.

3. A reads assembler. REPdenovo uses Velvet at this point. In the future, we may support different assembler. Velvet can be downloaded from: https://www.ebi.ac.uk/~zerbino/velvet/. Caution: if you want to assemble k-mers that are longer than 30 bp, you need to recompile Velvet to let it work with longer sequence length. For example: make ’MAXKMERLENGTH=60’. This makes Velvet work for k-mer length up to 60. 

4. Reads mapping. REPdenovo uses bwa mem. BWA (version 0.7 or later) can be downloaded from https://github.com/lh3/bwa.

5. Sequence processing utilities. These include the commonly used samtools. Our code also uses bamtools (https://github.com/pezmaster31/bamtools), but bamtools is not required to be installed.

## **Download and Install**
First, download the whole folder from https://github.com/Reedwarbler/REPdenovo, including the subfolder TERefiner and ContigsMerger-v0.1.9.

By default, users can directly run the tool and there is no need to install if you have all the dependencies installed. However, on some machines users may fail to run the pre-compiled tools TERefiner_1 and ContigsMerger, then users need to compile by themselves and run the follow commands:

cd TERefiner  &&  make  &&  cd .. 

cd ContigsMerger-v0.1.9  &&  make  &&  cd .. 

cp ./TERefiner/TERefiner_1 ./  &&  cp ./ContigsMerger-v0.1.9/ContigsMerger ./ 

chmod +x ./TERefiner_1  &&  chmod +x ./ContigsMerger 

## **Preparing inputs**
REPdenovo takes sequence reads in the FASTQ format (uncompressed or compressed in .fastq.gz format). A raw reads file which list the path, mean and standard derivation of the insert size should be provided in the format:

**Read-file-path group mean-insert-size insert-size-standard-derivation**

For single end reads, group, mean-insert-size and insert-size-standard-derivation should be set to -1. 

For paired-end reads, left raw reads and right raw reads should be in separate files, and in two lines (one line for left raw reads, and the other for right raw reads). 
The "group" should be same for these two lines. 
**Users can find sample files (for both paired-end reads and single-end reads) from the same folder in this github cite.**

REPdenovo needs a configuration file, which tells REPdenovo the basic settings. Users can find one sample from the same folder in this github cite.

Here, we give an explanation on the parameters.
**In general, you can find all the entries in the sample configuration file.** For some parameters, 
the values shown in the example are perhaps those that you should use (especially those are said to not change below).

MIN_REPEAT_FREQ.  This is the cutoff of k-mers that are considered to be frequent for assembly. Note that this is the relative to the average coverage of the sequence reads. The average coverage of the sequence reads is calculated by the number of reads, reads length and the genome size. 

RANGE_ASM_FREQ_DEC and RANGE_ASM_FREQ_GAP: these are used for assembly. Usually you don't need to change these.

K_MIN, K_MAX and K_INC: the smallest value, maximum value and increment of K. REPdenovo can use different K. In the example shown in the figure, three K values will be used: 30, 40 and 50.

K_DFT: default value of K value. This is equivalent of setting K_MIN = K_MAX = K_DFT.

READ_LENGTH: length of reads.

READ_DEPTH: reads depth.

THREADS: how many threads to use to run Jellyfish, BWA and ContigsMerger.
GENOME_LENGTH: the length of the genome. Can only provide an approximate one. 

ASM_NODE_LENGTH_OFFSET: if set to -1, then require each k-mer in
the repeat be frequent. That is, all k-mers in a repeat is considered to be frequent.

MIN_CONTIG_LENGTH: the minimum length of the contigs output.
IS_DUPLICATE_REPEATS: ratio used to check whether two repeats are duplicate. If the similarity is over this threshold, then see the two repeats are duplicate.

COV_DIFF_CUTOFF, MIN_SUPPORT_PAIRS, MIN_FULLY_MAP_RATIO, : used by REPdenovo in improving quality of assembled repeats. You don't usually need to change these.

TR_SIMILARITY: REPdenovo merges two assembled repeats if their similarity is over this threshold.

JELLYFISH_PATH: set to the path of Jellyfish executable. If already added the PATH, then use GLOBAL in this field.

VELVET_PATH: set to the path of Velvet executable. If already added the PATH, then use GLOBAL in this field.

BWA_PATH: set to the path of BWA executable. If already added the PATH, then use GLOBAL in this field.

SAMTOOLS_PATH: set to the path of samtools executable. If already added the PATH, then use GLOBAL in this field.

REFINER_PATH: set to the path of TERefiner_1 executable. If put at the same folder as main.py, then use GLOBAL in this field.

CONTIGS_MERGER_PATH: set to the path of ContigsMerger executable. If put at the same folder as main.py, then use GLOBAL in this field.

OUTPUT_FOLDER : where to output the results. It is relative to the installation folder of REPdenovo.

VERBOSE. If set to be 1, output more information about the current running states of REPdenovo.


## **Basic usage**
Run the assembly part:

python ./main.py Assembly configuration-file-name raw-reads-file-name 
 
Then run the scaffolding and analysis part (not work for single-end reads):

python ./main.py Scaffolding configuration-file-name raw-reads-file-name

## **Ouptput**
Two main output:

contigs.fa which contains the constructed repeats.

X_contig_pairs_info.txt_cov_info_with_cutoff.txt contains the repeat coverage information.
