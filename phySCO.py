#!/bin/env python3

"""
#### **phySCO** is a python script to infer maximum likelihood _phy_logenomic tree using BUSCO _s_ingle-_c_opy _o_rthologous genes. ####

General description: given a directory containing busco results for a group of species, this script computes the maximum likelihood phylogenetic tree of all the species.
This script is currently able to process just AMINO ACID SEQUENCES of complete single-copy BUSCO genes.
In fact, amino acids sequences are always returned by any BUSCO analysis run mode (either genome, transcriptome or protein).

REQUIRED SOFTWARES: mafft, trimal, IQTREE2

To run this script, just copy the BUSCO output folders for each species with their default structure into a dedicated directory.
Then, make sure that each species folder specifies a unique identifier, which should be placed at the beginning of the name followed by an underscore.
Here you can find an example of input directory: 
input_dir/
├── uniqueID1_whateveryouWant/
|   └── {busco default outputs for species uniqueID1}
├── uniqueID2_whateveryouWant/
|   └── {busco default outputs for species uniqueID2}
├── uniqueID3_whateveryouWant/
|   └── {busco default outputs for species uniqueID3}
...
└── uniqueIDN_whateveryouWant/
    └── {busco default outputs for species uniqueIDN}

Written by: Filippo Nicolini
Last update: 26/11/2024

Aknowledgements. Many thanks to Niccolò Righetti (https://github.com/NiccoloRighetti), whose work has pushed me to write this code.

"""

#-------------------------------------------------------------------------------

import subprocess, argparse, sys, os, glob, math
from argparse import RawTextHelpFormatter


############################
#     Define functions     #
############################

# Type function for argparse: a float within predefined bounds
def range_limited_float_type(arg):
    try:
        f = float(arg)
    except ValueError:    
        raise argparse.ArgumentTypeError("Must be a floating point number.")
    if f < 0 or f > 1:
        raise argparse.ArgumentTypeError("Argument must be within " + str(0) + " and " + str(1) + "!")
    return f

# Function to create a new directory, after checking for it's existence
def create_new_directory(dir_name):
    if os.path.isdir(dir_name):
        subprocess.run(f"rm -rf {dir_name}",
                    shell = True)
        subprocess.run(f"mkdir -p {dir_name}",
                        shell = True)
    else:
        subprocess.run(f"mkdir -p {dir_name}",
                        shell = True)

# Function to align fasta files with mafft in auto mode
def align_fasta(input_fasta, output_fasta):
    try:
        mafft_process = subprocess.run(f"mafft --auto --thread 15 {input_fasta} > {output_fasta}",
                                       shell = True,
                                       capture_output = True,
                                       text = True)
        
        mafft_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err)

# Function to trim alignments with trimAl in auto mode
def trim_alignments(input_fasta, output_fasta):
    try:
        trimal_process = subprocess.run(f"trimal -in {input_fasta} -out {output_fasta} -automated1",
                                       shell = True,
                                       capture_output = True,
                                       text = True)
        
        trimal_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err)


##########################################
#     Define arguments of the script     #
##########################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "### WELCOME TO phySCO - *phy*logenetic trees from busco *s*ingle-*c*opy *o*rthologous genes ###\n"
                                 "This script computes the Maximum Likelihood phylogenetic tree of a group of species, using BUSCO complete single-copy genes.\n"
                                 "BUSCO is supposed to have been already run.",
                                 formatter_class = RawTextHelpFormatter)

# Define options/arguments/parameters
parser.add_argument("-i", "--input_dir",
                    help = "Directory containing the default BUSCO output for each species.",
                    required = True)

parser.add_argument("-o", "--output_dir", 
                    help = "Name of the output directory. (Default = ML_phylogey)",
                    default = "ML_phylogeny")

parser.add_argument("-t", "--occupancy_threshold",
                    type = range_limited_float_type,
                    help = "Range: 0-1. The minimum percentage of species required to keep a\n"
                           "gene. E.g., if this argument is set to 0.8, all the genes that are\n"
                           "present in less than the 80%% of species will be discarded from\n"
                           "the phylogenetic analysis. (Default = 0.85)",
                    default = 0.85)

parser.add_argument("-g", "--genes_to_keep",
                    type = int,
                    help = "The number of genes you want to use to infer the phylogenetic tree.\n"
                           "E.g., if you want to speed up the phylogenetic analysis, you can\n"
                           "use just 10 random genes. Do not include this flag if you want to use all\n"
                           "the available genes.")

parser.add_argument("-m", "--merge_partitions",
                    action = "store_false",
                    help = "If used, this flag enables partition merging of IQTREE\n"
                           "(--merge --rcluster-max 25). (Default = False)")

parser.add_argument("-s", "--sequence_type",
                    help = "Choose the type of fasta file to process, i.e., amino acids (AA)\n"
                           "or nucleotides (NT). (Default = AA)",
                    choices = ["AA","NT"],
                    default = "AA")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


#---------------------------------------------------------------------------------

########################################
#     Create the output directory      #
########################################

print()

print(f"Creating output directory in {args.output_dir}/")
create_new_directory(args.output_dir)

print()


##################################################################################
#     Aggregate BUSCO complete single-copy genes into separate fasta files       #
##################################################################################

if args.input_dir.endswith("/"):
    INPUT_DIR = args.input_dir[:-1]
else:
    INPUT_DIR = args.input_dir

# Get the list of unique identifiers per each species
speciesID_list = [os.path.basename(i).split("_")[0] for i in os.listdir(INPUT_DIR)]
print(f"{len(speciesID_list)} species identifier found in {INPUT_DIR}/:")
print("    " + ", ".join(speciesID_list) + ".")
print(f"    {int(round(len(speciesID_list)*args.occupancy_threshold))} is the minimum number of species to keep a gene for downstream analysis\n"
      f"    (occupancy threshold: {int(args.occupancy_threshold*100)}%).")
if args.genes_to_keep is not None:
    print(f"    {args.genes_to_keep} genes will be used for downstream analyses.")
else:
    print(f"    All the available genes will be used for downstream analyses.")
print()

# Get the list of all the complete BUSCO genes found in the input species
gene_file = args.output_dir + "/complete_busco_genes.ls"

print("Retrieving the list of all the complete BUSCO genes found in the input species...")
try:
    awk_process = subprocess.run(f"cat {INPUT_DIR}/*/run*/full_table.tsv | "
                                 "awk -F \"\t\" '{if ($2 == \"Complete\") print $1}' | "
                                 "sort -u | "
                                 f"shuf > {gene_file}",
                                 shell = True)
    
    awk_process.check_returncode()

except subprocess.CalledProcessError as err:
    print("An error occured:", err)

# Store the obtained gene IDs in a list
gene_list = []
with open(gene_file) as input_list:
    for line in input_list.readlines():
        gene_list.append(line.strip())

print(f"    {len(gene_list)} unique BUSCO gene identifiers found.")

# Create a directory to store fasta files
fasta_dir = f"{args.output_dir}/fasta_files"
create_new_directory(fasta_dir)

# Counters for discarded and kept genes
discarded_genes = 0
kept_genes = 0

# Check if fasta files exist
if args.sequence_type == "NT":
    if not glob.glob(f"{INPUT_DIR}/*/run*/busco_sequences/single_copy_busco_sequences/*.fna"):
        print()
        print(f"### Mmh, you selected the {args.sequence_type} fasta type (*.fna), but no such file has been found. ###\n"
              f"### You may want to check your input directory ({INPUT_DIR}/). ###\n"
              "Quitting the analysis...")
        print()
    
        quit()
    else:
        print(f"    The analyses will be performed on nucleotide sequences!")
        print()

elif args.sequence_type == "AA":
    if not f"{INPUT_DIR}/*/run*/busco_sequences/single_copy_busco_sequences/*.faa":
        print()
        print(f"### Mmh, you selected the {args.sequence_type} fasta type (*.faa), but no such file has been found. ###\n"
              f"### You may want to check your input directory ({INPUT_DIR}/). ###\n"
              "Quitting the analysis...")
        print()

        quit()
        
    else:
        print(f"    The analyses will be performed on amino acid sequences!")
        print()


# Generate separate fasta files with genes from species meeting inclusion threshold
print(f"Generating fasta files for each gene with all the annotated species...")
for gene in gene_list:

    if args.sequence_type == "NT":
        # Get the list of each genes in each species (nucleotides)
        gene_file_list = glob.glob(f"{INPUT_DIR}/*/run*/busco_sequences/single_copy_busco_sequences/{gene}.fna")

    elif args.sequence_type == "AA":
        # Get the list of each genes in each species (amino acids)
        gene_file_list = glob.glob(f"{INPUT_DIR}/*/run*/busco_sequences/single_copy_busco_sequences/{gene}.faa")

    # Get the occupancy thershold for each gene
    occupancy = len(gene_file_list)/len(os.listdir(INPUT_DIR))

    # print(f"{gene=} {occupancy=} {args.occupancy_threshold=}")

    # Check if occupancy requirements are met, in which case do not continue with the analysis
    if float(occupancy) < float(args.occupancy_threshold):
        discarded_genes += 1
        continue
    else:
        kept_genes += 1

    # For each gene in each species, change the header to include just the unique species identifier
    for id in speciesID_list:
        if args.sequence_type == "NT":
            input_fasta = glob.glob(f"{INPUT_DIR}/{id}*/run*/busco_sequences/single_copy_busco_sequences/{gene}.fna")
            output_fasta = f"{args.output_dir}/fasta_files/{gene}.fna"
        elif args.sequence_type == "AA":
            input_fasta = glob.glob(f"{INPUT_DIR}/{id}*/run*/busco_sequences/single_copy_busco_sequences/{gene}.faa")
            output_fasta = f"{args.output_dir}/fasta_files/{gene}.faa"

        if input_fasta:
            # Progressively add sequences to each fasta file
            subprocess.run(f"touch {output_fasta}",
                           shell = True,
                           text = True)

            subprocess.run(f"cat {input_fasta[0]} | sed -E 's/^>.+$/>{id}/' >> {output_fasta}",
                           shell = True,
                           text = True)
        else:
            continue
        
    # if a maximum number of genes to use has been set, check whether to break the cycle or not
    if args.genes_to_keep is not None and kept_genes == args.genes_to_keep:
        break

if discarded_genes == 1:
    print(f"    {discarded_genes} gene was discarded because didn't meet the occupancy threshold ({args.occupancy_threshold}).")
else:
    print(f"    {discarded_genes} genes were discarded because didn't meet the occupancy threshold ({args.occupancy_threshold}).")

print(f"    {kept_genes} genes were kept for downstream phylogenomic analysis.")
print("Done")


##################################################
#     Align BUSCO complete single-copy genes     #
##################################################

print()
print("Aligning generated files with mafft...")

# Create a directory to store alignments
alignments_dir = f"{args.output_dir}/fasta_files_aligned"
create_new_directory(alignments_dir)

fasta_file_list = os.listdir(fasta_dir)

print(f"    {len(fasta_file_list)} fasta files found in {fasta_dir}/...")

# Align fasta files with mafft (auto mode)
for file in fasta_file_list:
    INPUT_FASTA = fasta_dir + "/" + file
    if args.sequence_type == "NT":
        OUTPUT_FASTA = alignments_dir + "/" + file.replace(".fna","_aligned.fna")
    elif args.sequence_type == "AA":
        OUTPUT_FASTA = alignments_dir + "/" + file.replace(".faa","_aligned.faa")

    # print(f"{INPUT_FASTA=} {OUTPUT_FASTA=}")

    align_fasta(INPUT_FASTA, OUTPUT_FASTA)

print("Done")


###############################################################
#     Trim alignments of BUSCO complete single-copy genes     #
###############################################################

print()
print("Trimming generated alignments with trimAl...")

# Create a directory to store trimmed alignments
trim_dir = f"{args.output_dir}/fasta_files_aligned_trimmed"
create_new_directory(trim_dir)

alignments_file_list = os.listdir(alignments_dir)

print(f"    {len(alignments_file_list)} alignment files found in {alignments_dir}/...")

# Trim alignments with trimAl (automated1 mode)
for file in alignments_file_list:
    INPUT_FASTA = alignments_dir + "/" + file
    if args.sequence_type == "NT":
        OUTPUT_FASTA = trim_dir + "/" + file.replace(".fna","_trim.fna")
    elif args.sequence_type == "AA":
        OUTPUT_FASTA = trim_dir + "/" + file.replace(".faa","_trim.faa")
    # print(f"{INPUT_FASTA=} {OUTPUT_FASTA=}")

    trim_alignments(INPUT_FASTA, OUTPUT_FASTA)

print("Done")


################################################################
#     Infer ML phylogenetic tree out of trimmed alignments     #
################################################################

if args.merge_partitions == False:
    print()
    print("Inferring ML phylogenetic tree with IQ-TREE...\n"
          "    Partition merging is enabled.\n"
          "    This may take a while...")

    try:
        iqtree_process = subprocess.run(f"iqtree2 -p {trim_dir} -m MFP --merge --rcluster-max 25 -nstop 500 -T AUTO -bb 1000 --runs 3 --prefix {args.output_dir}/MLtree",
                                        shell = True,
                                        capture_output = True,
                                        text = True)
        
        iqtree_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err)
else:
    print()
    print("Inferring ML phylogenetic tree with IQ-TREE...\n"
          "    Partition merging is disabled.\n"
          "    This may take a while...")

    try:
        iqtree_process = subprocess.run(f"iqtree2 -p {trim_dir} -m MFP -nstop 500 -T AUTO -bb 1000 --runs 3 --prefix {args.output_dir}/MLtree",
                                        shell = True,
                                        capture_output = True,
                                        text = True)
        
        iqtree_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err)


print("Done")
