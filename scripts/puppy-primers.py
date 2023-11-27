#!/usr/bin/env python3

# Script: puppy-primers
# Part of PUPpy distribution.
# Designs taxon-specific primers for user-defined organisms.

# Author = Hans Ghezzi
# Credits = Hans Ghezzi, Juan Burckhardt
# Licence = GPL3
# Version = 1.0.0
# Maintainer = Hans Ghezzi
# Email = hans.ghezzi@gmail.com
# Status = Development

# Load libraries
import primer3
import pandas as pd
import os
from glob import glob
import argparse
from Bio import SeqIO
from colorama import Fore, Style
import subprocess
import time
import dask.distributed as Client
import dask.dataframe as dd
import re
import shutil

######################################################## flags ###########################################################################

parser = argparse.ArgumentParser(
    prog="PUPpy",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
    PUPpy: A Phylogenetically Unique Primer pipeline in PYthon for high resolution classification of defined microbial communities. Version 1.0
    """,
)

# Setting up flags in python
parser.add_argument(
    "-p",
    "--primers_type",
    help="Design unique or shared primers among the target bacterial group",
    default="unique",
    choices=["unique", "group"],
)
parser.add_argument(
    "-in",
    "--intended",
    help="Directory containing the CDS files for the species to design taxon-specific primers",
    default="",
    required=True,
)
parser.add_argument(
    "-i",
    "--input",
    help="Input file to generate primers. Either 'ResultDB.tsv' OR 'UniqueGenesList.tsv' file from `puppy-align` must be provided",
    default="",
    required=True,
)
parser.add_argument(
    "-o",
    "--outdir",
    help="Relative path to the output folder",
    default="Primer3_output",
)
parser.add_argument(
    "-ng",
    "--genes_number",
    help="Number of genes per species to design primers",
    default=5,
    type=int,
)
parser.add_argument(
    "-np",
    "--primers_number",
    help="Number of primer pairs to design for each gene",
    default=4,
    type=int,
)
parser.add_argument(
    "-ops", "--optimal_primer_size", help="Primer optimal size", default=20, type=int
)
parser.add_argument(
    "-mips", "--min_primer_size", help="Primer minimum size", default=18, type=int
)
parser.add_argument(
    "-maps", "--max_primer_size", help="Primer maximum size", default=22, type=int
)
parser.add_argument(
    "-optm",
    "--optimal_primer_Tm",
    help="Primer optimal melting temperature",
    default=60.0,
    type=float,
)
parser.add_argument(
    "-mitm",
    "--min_primer_Tm",
    help="Primer minimum melting temperature",
    default=58.0,
    type=float,
)
parser.add_argument(
    "-matm",
    "--max_primer_Tm",
    help="Primer maximum melting temperature",
    default=63.0,
    type=float,
)
parser.add_argument(
    "-tmd",
    "--max_Tm_diff",
    help="Maximum Tm difference between the primer pair",
    default=2.0,
    type=float,
)
parser.add_argument(
    "-migc",
    "--min_primer_gc",
    help="Primer minimum GC content",
    default=40.0,
    type=float,
)
parser.add_argument(
    "-magc",
    "--max_primer_gc",
    help="Primer maximum GC content",
    default=60.0,
    type=float,
)
parser.add_argument(
    "-s",
    "--product_size_range",
    help="Product size range",
    default=[75, 150],
    nargs="+",
    type=int,
)
parser.add_argument(
    "-mpolyx", "--max_poly_x", help="Maximum poly X allowed", default=3, type=int
)
parser.add_argument("-GCc", "--GC_clamp", help="Primer GC clamp", default=1, type=int)

args = parser.parse_args()

######################################################## variables ###########################################################################

# Path to folder with CDS files. Takes user-defined directory with target species (-in or --intended)
CDS = os.path.abspath(args.intended)
# Path to CDS filenames with .fna extension
cds_filenames = CDS + "/*.fna"

# Path to output directory. Takes default or user-defined outdir (-o or --outdir)
OUT = os.path.abspath(args.outdir)

# Check if output folder exists. If it does, exit the program, otherwise create it.
while True:
    if os.path.exists(OUT):
    # Get user input (yes or no)
        user_input = input("Do you want to OVERWRITE the data in the existing output directory? (yes/no): ").strip().lower()

        # Check the user's response
        if user_input == "yes":
            # Overwrite output directory
            print("Overwriting output directory...")
            directory_to_overwrite = OUT
            try:
                # Remove the existing directory and its contents
                shutil.rmtree(directory_to_overwrite)
            except FileNotFoundError:
                # Handle the case where the directory doesn't exist
                pass

            # Recreate the directory (an empty directory)
            try:
                os.makedirs(directory_to_overwrite)
            except OSError as e:
                # Handle any errors that occur while creating the directory
                print("Error creating directory: " + str(e))
            break #Exit loop if user provides valid input

        elif user_input == "no":
            # User chose not to perform the task
            print("Please choose a different output directory name. Exiting...")
            break #Exit loop if user provides valid input
        else:
            # Handle invalid input
            print("Invalid input. Please enter 'yes' or 'no'.")

    else:
        os.mkdir(OUT)
        break #Exit loop if output dir doesn't exist initially


# Path to the individual primer3 output files, with .tsv extension
primer3_files = os.path.join(OUT, r"primer3_files")
primer_filenames = primer3_files + "/*.tsv"

# Extract basename of alignment file from user-provided path. Takes user-defined input (-i or --input)
INPUT = os.path.basename(args.input)
# Path to the file with gene numbers.
geneNums = os.path.join(OUT, "geneNumbers.tsv")

######################################################## Functions ###########################################################################

def extract_seq_dict(names, bacteria_name):
    """
    Extract filenames based on the user-defined species (-t flag).
    Store FASTA headers and sequences of target species in a dictionary
    """
    # Extract filename to use in SeqIO parsing
    for f in glob(names):
        if bacteria_name in f:
            record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
            return record_dict

def extract_seq(gene_name, dictionary):
    """
    Store FASTA sequences as a variable to then design primers with primer3
    """
    for key in dictionary:
        if gene_name in key:
            sequence = str(
                dictionary[key].seq
            )  # Save the fasta sequence ONLY, as a string
            return sequence

def mainGroup(f, df, name, CDSes, found, missing):
    """
    From a pile output by primer3, parse the file and put it into a dictionary.
    Then, change the dictionary into a pandas dataframe.
    """
    i = fBase.index("-")
    gene = fBase[Index + 1 : i]
    Dictionary = extract_seq_dict(CDSes, name)
    seq = extract_seq(gene, Dictionary)
    Found = found
    Missing = missing

    # parse contents into dictionary
    with open(f, "r") as f1:
        data = f1.readlines()
        data_dic = lines_to_dicGroup(data, name, gene, seq, Found, Missing)
    # load dictionary into pandas dataframe
    df = pd.DataFrame.from_dict(data_dic)

    return df

def mainUnique(f, df, name, CDSes):
    """
    From a pile output by primer3, parse the file and put it into a dictionary.
    Then, change the dictionary into a pandas dataframe.
    """
    i = f.index(".tsv")
    gene = f[Index + 1 : i]
    Dictionary = extract_seq_dict(CDSes, name)
    seq = extract_seq(gene, Dictionary)

    # parse contents into dictionary
    with open(f, "r") as f1:
        data = f1.readlines()
        data_dic = lines_to_dicUnique(data, name, gene, seq)
    # load dictionary into pandas dataframe
    df = pd.DataFrame.from_dict(data_dic)

    return df

def lines_to_dicGroup(lines, bacteria_name, gene_name, gene_seq, targetsFound, targetsMissing):
    """
    Parse through the lines extracted from the primer3 output file to obtain
    different parameters. Store them into a dictionary
    """
    # initiate dictionary
    keys = [
        "species",
        "gene",
        "Full_sequence",
        "primer_option",
        "pair_penalty_score",
        "amplicon_size",
        "F_primer",
        "R_primer",
        "F_length",
        "R_length",
        "F_tm",
        "R_tm",
        "F_GC",
        "R_GC",
        "num_targets_found",
        "num_targets_missing"
    ]
    primer_dic = {key: [] for key in keys}
    # initialize counter
    count = 0
    for line in lines:
        # split line to obtain the name of the statistic and the value
        name, value = line.split("\t")

        # Dont store more than what was inputted by the user (-1 added since the notation starts at 0)
        if count > args.primers_number - 1:
            break
        # continue if the line does not have count in it
        if str(count) not in name:
            continue

        # check that it is a primer statistic
        if int(name.split("_")[2]) == count:
            # pair penalty score
            if "PRIMER_PAIR_" + str(count) + "_PENALTY" == name:
                primer_dic["pair_penalty_score"].append(float(value))
            # left primer seq and length
            elif "PRIMER_LEFT_" + str(count) + "_SEQUENCE" == name:
                primer_dic["F_primer"].append(value.strip())
                primer_dic["F_length"].append(len(value))
            # right primer seq and length
            elif "PRIMER_RIGHT_" + str(count) + "_SEQUENCE" == name:
                primer_dic["R_primer"].append(value.strip())
                primer_dic["R_length"].append(len(value))
            # left primer TM
            elif "PRIMER_LEFT_" + str(count) + "_TM" == name:
                primer_dic["F_tm"].append(float(value))
            # right primer TM
            elif "PRIMER_RIGHT_" + str(count) + "_TM" == name:
                primer_dic["R_tm"].append(float(value))
            # left primer GC%
            elif "PRIMER_LEFT_" + str(count) + "_GC_PERCENT" == name:
                primer_dic["F_GC"].append(float(value))
            # right primer GC%
            elif "PRIMER_RIGHT_" + str(count) + "_GC_PERCENT" == name:
                primer_dic["R_GC"].append(float(value))
            # amplicon size
            elif "PRIMER_PAIR_" + str(count) + "_PRODUCT_SIZE" == name:
                primer_dic["amplicon_size"].append(int(value))
                # add all other information since this is the last value for a primer pair
                primer_dic["species"].append(bacteria_name)
                primer_dic["gene"].append(gene_name)
                primer_dic["Full_sequence"].append(gene_seq)
                primer_dic["num_targets_found"].append(targetsFound)
                primer_dic["num_targets_missing"].append(targetsMissing)
                primer_dic["primer_option"].append(count)
                # add +1 to counter for next set of primers
                count += 1
    return primer_dic

def lines_to_dicUnique(lines, bacteria_name, gene_name, gene_seq):
    """
    Parse through the lines extracted from the primer3 output file to obtain
    different parameters. Store them into a dictionary
    """
    # initiate dictionary
    keys = [
        "species",
        "gene",
        "Full_sequence",
        "primer_option",
        "pair_penalty_score",
        "amplicon_size",
        "F_primer",
        "R_primer",
        "F_length",
        "R_length",
        "F_tm",
        "R_tm",
        "F_GC",
        "R_GC",
    ]
    primer_dic = {key: [] for key in keys}
    # initialize counter
    count = 0
    for line in lines:
        # split line to obtain the name of the statistic and the value
        name, value = line.split("\t")

        # Dont store more than what was inputted by the user (-1 added since the notation starts at 0)
        if count > args.primers_number - 1:
            break
        # continue if the line does not have count in it
        if str(count) not in name:
            continue

        # check that it is a primer statistic
        if int(name.split("_")[2]) == count:
            # pair penalty score
            if "PRIMER_PAIR_" + str(count) + "_PENALTY" == name:
                primer_dic["pair_penalty_score"].append(float(value))
            # left primer seq and length
            elif "PRIMER_LEFT_" + str(count) + "_SEQUENCE" == name:
                primer_dic["F_primer"].append(value.strip())
                primer_dic["F_length"].append(len(value))

            # right primer seq and length
            elif "PRIMER_RIGHT_" + str(count) + "_SEQUENCE" == name:
                primer_dic["R_primer"].append(value.strip())
                primer_dic["R_length"].append(len(value))
            # left primer TM
            elif "PRIMER_LEFT_" + str(count) + "_TM" == name:
                primer_dic["F_tm"].append(float(value))

            # right primer TM
            elif "PRIMER_RIGHT_" + str(count) + "_TM" == name:
                primer_dic["R_tm"].append(float(value))
            # left primer GC%
            elif "PRIMER_LEFT_" + str(count) + "_GC_PERCENT" == name:
                primer_dic["F_GC"].append(float(value))

            # right primer GC%
            elif "PRIMER_RIGHT_" + str(count) + "_GC_PERCENT" == name:
                primer_dic["R_GC"].append(float(value))
            # amplicon size
            elif "PRIMER_PAIR_" + str(count) + "_PRODUCT_SIZE" == name:
                primer_dic["amplicon_size"].append(int(value))
                # add all other information since this is the last value for a primer pair
                primer_dic["species"].append(bacteria_name)
                primer_dic["gene"].append(gene_name)
                primer_dic["Full_sequence"].append(gene_seq)
                primer_dic["primer_option"].append(count)
                # add +1 to counter for next set of primers
                count += 1
    return primer_dic

def parse_alignments(dataframe):
    """
    # Split dataframe query and target columns to isolate species and gene name in each row
    """
    queryNames = dataframe["query"].str.split(r"_cds_|.peg.", n=1, expand=True)
    dataframe.insert(loc=0, column="query_species", value=queryNames[0])
    dataframe.insert(loc=1, column="query_gene", value="cds_" + queryNames[1])
    dataframe.loc[dataframe["query"].str.contains("peg"), "query_gene"] = dataframe[
        "query_gene"
    ].str.replace("cds_", "peg.")

    # Change target column
    targetNames = dataframe["target"].str.split(r"_cds_|.peg.", n=1, expand=True)
    dataframe.insert(loc=2, column="intended_targets", value=targetNames[0])
    dataframe.insert(loc=3, column="target_gene", value="cds_" + targetNames[1])
    dataframe.loc[dataframe["target"].str.contains("peg"), "target_gene"] = dataframe[
        "target_gene"
    ].str.replace("cds_", "peg.")

def partialParsing(dataframe):
    """
    # Split dataframe query and target columns to isolate species and gene name in each row
    """
    queryNames = dataframe["query"].str.split(r"-|_cds_|.peg.", n=0, expand=True)
    dataframe.insert(loc=0, column="species_name", value=queryNames[0])
    dataframe.insert(loc=1, column="species_gene", value="cds_" + queryNames[2])
    dataframe.drop(columns=["query"], inplace=True)

def rename_fasta(f, species_name):
    """
    Function to modify all the FASTA headers in each input CDS file.
    """
    with open(f, "r") as f1:
        data = f1.readlines()
        new_lines = rename_fasta_headers(data, species_name)

    with open(f, "w") as f2:
        f2.writelines(new_lines)

def rename_fasta_headers(list_of_lines, species_name):
    """
    Function to select FASTA header lines to which apply name changes
    """
    # run for loop to parse every line and change the headers
    new_lines = []
    for line in list_of_lines:
        if line.startswith(">"):
            new_lines.append(change_header(line, species_name))
        else:
            new_lines.append(line)
    return new_lines

def change_header(line, name):
    """
    Function to rename FASTA headers to include species name.
    """
    # find the | in the string
    if "|" in line:
        pipe_index = line.index("|") + 1
    else:
        pipe_index = 1

    # add the species name to the string
    if name in line:
        return line
    else:
        header_name = line[pipe_index:]
        if ".peg." in line:
            separator_index = header_name.index(".peg.")
            New_line = (
            ">lcl|"
            + name
            + "-"
            + header_name[:separator_index]
            + "_cds_"
            + header_name[separator_index+5:]
            )
        else:
            separator_index = header_name.index("_")
            New_line = (
                ">lcl|"
                + name
                + "-"
                + header_name[:separator_index]
                + "_cds"
                + header_name[separator_index:]
            )
        return New_line

######################################################## run program ###########################################################################

# Execute function 'rename_fasta' to rename FASTA headers in input CDS files.
print("Renaming FASTA headers (if needed)...")
for f in glob(cds_filenames):
    name = os.path.basename(f)
    i = name.index("_cds")
    rename_fasta(f, name[:i])

# Load dataframe with all the alignments
if INPUT == "ResultDB.tsv":
    ## Load alignment file (ResultDB.tsv) as a dataframe
    # Define header line for all the dataframes
    Line = [
        "query",
        "target",
        "qlen",
        "tlen",
        "alnlen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "pident",
        "qcov",
        "tcov",
        "evalue",
    ]
    # Determine size of alignments file (ResultDB.tsv):
    ALN_Size = os.path.getsize(args.input) / 1000000000
    # Load ResultDB.tsv
    print(
        f"\nLoading the alignments file ({args.input}, size = {ALN_Size} GB) as a Dask Dataframe..."
    )

    # Load input file 'ResultDB.tsv' as dask dataframe
    daskDF = dd.read_csv(
        args.input, sep="\t", names=Line, lineterminator="\n", low_memory=False, blocksize="64MB"
    )
    print("Converting Dask Dataframe to Pandas dataframe...")
    ResultDB = daskDF.compute() # convert dask to pandas dataframe

    # Check if df ResultDB has duplicate index values. If yes, remove:
    if ResultDB.index.is_unique == False:   # Duplicated index values present
        ResultDB=ResultDB.loc[~ResultDB.index.duplicated(), :]
    elif ResultDB.index.is_unique == True:   # No duplicated index values present
        pass

    # Parse alignments
    parse_alignments(ResultDB)

# Workflow to design taxon-specific primers UNIQUE to each target member provided
if args.primers_type == "unique":
    # Calculate the number of genes in each target species
    print(
        "Calculating the number of genes in each target species. Output will be stored in file: {}.".format(
            geneNums
        )
    )
    # Open each CDS.fna file of the target organisms
    with open(geneNums, "w") as file:
        for f in glob(cds_filenames):
            name = os.path.basename(f)
            i = name.index("_cds")
            # Calculate the number of genes for each species in the defined community
            fileOpen = open(f)
            n = 0
            for line in fileOpen:
                if line.startswith(">"): # Count gene numbers by the number of FASTA headers found
                    n += 1
            fileOpen.close()
            file.write(name[:i] + "\t" + str(n) + "\n")
    # Sort the genes count file in alphabetical order
    lines = open(geneNums, "r").readlines()
    output = open(geneNums, "w")
    for LINE in sorted(lines, key=lambda LINE: LINE.split()[0]):
        output.write(LINE)
    output.close()
    
    # Worflow if users provide the file 'ResultDB.tsv' as --input
    if INPUT == "ResultDB.tsv":
        print(f"\nParsing alignments from MMseqs2 output file")
        # Parse Parsed_ResultDB.tsv file to find unique genes in queryDB
        df_uniqueQuery = ResultDB.drop_duplicates(
            subset=["query"], keep=False
        ).sort_values("query")
        # Error message if no unique genes in the query database are found
        if df_uniqueQuery.empty:
            print(
                Fore.RED
                + 'Warning: query unique genes file is empty. Check file "ResultDB.tsv" to see if alignment step worked'
                + Style.RESET_ALL
            )
            print(Fore.RED + "Exiting..." + Style.RESET_ALL)
            exit()
        else:
            pass

        # Parse Parsed_ResultDB.tsv file to find unique genes in targetDB
        df_uniqueTarget = ResultDB.drop_duplicates(
            subset=["target"], keep=False
        ).sort_values("query")
        # Error message if no unique genes in the target database are found
        if df_uniqueTarget.empty:
            print(
                Fore.RED
                + 'Warning: target unique genes file is empty. Check file "ResultDB.tsv" to see if alignment step worked'
                + Style.RESET_ALL
            )
            print(Fore.RED + "Exiting..." + Style.RESET_ALL)
            exit()
        else:
            pass

        # Concatenate the 2 dataframes
        df_concat = pd.concat([df_uniqueQuery, df_uniqueTarget]).sort_values("query")

        # Drop query and target columns
        df_concat.drop(columns=["query"], inplace=True)
        df_concat.drop(columns=["target"], inplace=True)

        # Find unique genes in the alignments file
        print("Looking for unique genes...")
        # Only keep one occurrence of those lines that are duplicated. These will be the unique genes
        UniqueGenes = df_concat[df_concat.duplicated(keep="first")]
        # Subset dataframe to only keep column 1, 2, and 5.
        UniqueGenes = UniqueGenes.iloc[:, [0, 1, 4]]
        # Save file with unique genes found in the defined community
        header = ["species_name", "gene_name", "gene_length"]
        UniqueGenes.columns = header
        UniqueGenes.to_csv(
            os.path.join(OUT, r"UniqueGenesList.tsv"),
            index=False,
            sep="\t",
            lineterminator="\n",
        )

    # Worflow if users provide the file 'UniqueGenesList.tsv' as --input
    # Load final genes in a file, if already available from a previous run:
    elif INPUT == "UniqueGenesList.tsv":
        # Load file containing the list of unique genes as a dataframe
        UniqueGenes = pd.read_csv(args.input, sep="\t", lineterminator="\n")

    # Make summary statistics and plot for unique genes detected.
    stats = os.path.join(OUT, r"Stats_pipelineOutput.tsv")
    plot = os.path.join(OUT, r"UniqueGenes.pdf")

    # Check if file 'Stats_pipelineOutput.tsv' exists.
    if os.path.isfile(stats):
        print(
            "File {0} already exists and will be overwritten".format(
                os.path.basename(stats)
            )
        )
        os.remove(stats)
    # Check if file 'UniqueGenes.pdf' exists.
    if os.path.isfile(plot):
        print(
            "File {0} already exists and will be overwritten".format(
                os.path.basename(plot)
            )
        )
        os.remove(plot)
    # Run Rscript to make summary stats and plot
    subprocess.run(["chmod", "755", "./UniqueGenesPlotting.R"]) # Give permissions to run script
    subprocess.call(["Rscript", "./UniqueGenesPlotting.R", OUT, CDS]) # Run R script

    # Remove geneNums file
    os.remove(geneNums)

# Workflow to design taxon-specific primers shared by selected members (i.e. group primers) of the defined microbial community
elif args.primers_type == "group":
    targetList = []
    # Append target bacterial names to a list
    for f in glob(cds_filenames):
        name = os.path.basename(f)
        i = name.index("_cds")
        targetList.append(name[:i])

    ## Creates a dataframe with the number of alignments found for each gene against other organisms/genes

    # Filter dataframe to only keep alignments of genes found in target species
    targetResult = ResultDB[ResultDB["query"].str.contains("|".join(targetList))]
    # Group genes of target species based on the column 'query' (this will show all the columns to the right of query, grouped by query)
    grouped = targetResult.groupby("query")

    # Get number of target species
    numTargets = len(targetList)

    # Create empty list to append dataframes to be concatenated:
    DataIdeal = []
    DataSecond = []
    DataUndesired = []

    # Find which genes could be used to selectively amplify the target species as a group
    for group_name, df_group in grouped:
        notFound = 0  # Variable to check if genes from target species align to unintended targets
        speciesList = (
            []
        )  # Variable to store the names of species to which target genes align to.
        unintendedList = (
            []
        )  # Variable to store the names of unintended species to which target genes align to.
        pidentSum = 0  # Account for %ID. Ideal alignment sums for any one gene would be 100*len(targetList)
        qcovSum = 0  # Account for query coverage. Ideal coverage sum would be 1*len(targetList)
        tcovSum = 0  # Account for target coverage. Ideal coverage sum would be 1*len(targetList)
        # Iterate through rows of each group
        for row_index, row in df_group.iterrows():
            target_col = row["target"]
            pidentity = float(row["pident"])
            qcoverage = float(row["qcov"])
            tcoverage = float(row["tcov"])
            geneLength = int(row["qlen"])
            targetName = target_col[: target_col.index("-")]
            # Define which portion of the gene can be used for primer design.
            sequenceIncluded = [df_group['qstart'].max(), (df_group['qend'].min())] 
            sequenceIncluded_length = sequenceIncluded[1]-sequenceIncluded[0]+1 # This is the length of the sequence we can use to design primers  
            # Check if ANY of the intended target strains are found in each row of the intended_targets column of each group
            if any(substring in target_col for substring in targetList):
                speciesList.append(
                    targetName
                )  # If aligned species is one of the targets, append species name to list
                pidentSum = pidentSum + pidentity
                qcovSum = qcovSum + qcoverage
                tcovSum = tcovSum + tcoverage
            else:
                unintendedList.append(
                    targetName
                )  # If aligned species is NOT one of the targets, append species name this list
                notFound = (
                    notFound + 1
                )  # If not found, increase count by one. A count of 'notFound' higher >= 1 is not desirable and this gene (i.e. group_name) will not be considered.
        # Define more variables for downstream conditional appending
        SpeciesAligned = (
            ",".join(set(speciesList)) if len(speciesList) != 0 else "NA"
        )  # Concatenate list of species found to align to query gene. If none found, print NA
        numTargetsFound = len(set(speciesList))
        numTargetsMissing = numTargets - len(set(speciesList))
        numUnintendedTargets = len(set(unintendedList))
        MissingTargetSpecies = list(
            set(targetList).difference(speciesList)
        )  # Define list of target species that the gene did not align to
        MissingTargets = (
            ",".join(MissingTargetSpecies) if numTargetsMissing != 0 else "NA"
        )  # If non-0 list, concatenate. Otherwise, return NA
        UnintendeSpeciesList = list(
            set(unintendedList).difference(targetList)
        )  # Define list of target species that the gene did not align to
        UnintendedSpecies = (
            ",".join(set(UnintendeSpeciesList))
            if len(UnintendeSpeciesList) != 0
            else "NA"
        )  # Concatenate list of species found to align to query gene

        # Define comments to be attached to each parsed gene:
        ## Comments about number of alignments to intended target genes
        if (len(speciesList) == numTargets) and len(speciesList) == len(
            set(speciesList)
        ):
            NumAlignments = "NUMBER OF ALIGNMENTS:\nThis gene has exactly 1 alignment to each intended target."
        elif len(speciesList) > numTargets and (
            len(speciesList) == len(set(speciesList))
        ):
            NumAlignments = "NUMBER OF ALIGNMENTS:\nThis gene has more than 1 alignment to at least one intended target."
        elif len(speciesList) == numTargets and (
            len(speciesList) > len(set(speciesList))
        ):
            NumAlignments = "NUMBER OF ALIGNMENTS:\nThis gene has more than 1 alignment to at least one intended target and it does not amplify all targets."
        else:
            NumAlignments = "NUMBER OF ALIGNMENTS:\nThis gene does not amplify all the intended targets."

        ## Comments about amplification of unintended species (i.e. non-targets)
        if notFound < 1:
            unintendedComment = "UNINTENDED AMPLIFICATION:\nThis gene only amplifies intended species in the defined community."
        elif notFound >= 1:
            unintendedComment = "UNINTENDED AMPLIFICATION:\nWARNING - This gene amplifies unintended species in the defined community!"

        ## Comments about percentage identity of alignments:
        if pidentSum == 100 * numTargets:
            PidentComments = "ALIGNMENTS %ID:\nThis gene aligns perfectly (100% ID) to each target gene."
        elif pidentSum < 100 * numTargets:
            PidentComments = "ALIGNMENTS %ID:\nThis gene does not align perfectly (100% ID) to at least one target gene."

        ## Comments about alignments query coverage:
        if qcovSum == 1 * numTargets:
            qcovComments = "QUERY COVERAGE:\nThe entire length of this gene (i.e. 100% query coverage) aligns to each target gene."
        else:
            qcovComments = "QUERY COVERAGE:\nOnly a portion of this gene (i.e. <100% query coverage) aligns to at least one target gene."

        ## Comments about alignments target coverage:
        if tcovSum == 1 * numTargets:
            tcovComments = "TARGET COVERAGE:\nThis gene aligns to the entire sequence of each target gene (i.e. 100% target coverage)."
        else:
            tcovComments = "TARGET COVERAGE:\nThis gene does not align to the entire sequence of at least one target gene."

        allCommentsLIST = [
            NumAlignments,
            unintendedComment,
            PidentComments,
            qcovComments,
            tcovComments,
        ]
        allComments = "\n".join(allCommentsLIST)

        # Line to append to each row of dataframe
        LineToAppend = pd.DataFrame(
            [
                {
                    "query": group_name,
                    "gene_length": geneLength,
                    "SequenceIncluded": sequenceIncluded,
                    "SequenceIncluded_length": sequenceIncluded_length,
                    "targets": SpeciesAligned,
                    "num_targets_found": numTargetsFound,
                    "num_targets_missing": numTargetsMissing,
                    "Missing_targets": MissingTargets,
                    "num_unintended_targets": numUnintendedTargets,
                    "unintended_targets": UnintendedSpecies,
                    "Comments": allComments,
                }
            ]
        )

        if (
            (notFound < 1)
            and (len(speciesList) == numTargets)
            and len(speciesList) == len(set(speciesList))
            and (pidentSum == 100 * numTargets)
            and (qcovSum == 1 * numTargets)
            and (tcovSum == 1 * numTargets)
        ):
            DataIdeal.append(LineToAppend)
        elif (notFound >= 1) or (numTargetsFound == 1):  # this situation means that the specific group (i.e. the specific gene aligns to organisms OTHER THAN the ones asked by the user
            DataUndesired.append(LineToAppend)
        else:
            DataSecond.append(LineToAppend)

    # Concatenate list of datatframes
    if len(DataIdeal) > 0:
        DataIdealALL = pd.concat(DataIdeal).sort_values(by=["gene_length"], ascending=False)  
        # Here we don't sort by species_name because since the primers are group-specific, any one gene should amplify all the target species.
        # Thus, the sorting priority here is gene_length
    else:
        print("No ideal group genes identified that amplify all the target organisms.")

    if len(DataSecond) > 0:
        DataSecondALL = pd.concat(DataSecond).sort_values(
                by=["num_targets_found", "SequenceIncluded_length", "num_unintended_targets"],
                ascending=[False, False, True])     # Sort to choose best candidate group genes at the top
    else:
        print('Could not find "second choice" candidate group genes for the target organisms')

    DataUndesiredALL = (
        pd.concat(DataUndesired)
        if len(DataUndesired) > 0
        else print("Could not find genes that target unintended target organisms.")
    )

    # Parse dataframes to separate query into query_species and query_gene:
    partialParsing(DataIdealALL) if len(DataIdeal) > 0 else print(
        "Skipping parsing on optimal genes."
    )
    partialParsing(DataSecondALL) if len(DataSecond) > 0 else print(
        "Skipping parsing on suboptimal genes."
    )
    partialParsing(DataUndesiredALL) if len(DataUndesired) > 0 else print(
        "Skipping parsing on undesired genes."
    )

    if (len(DataIdeal) == 0) and (len(DataSecond) == 0) and (len(DataUndesired) == 0):
        print(
            "We could not categorise genes from the alignment file. Something went wrong. Check target CDSes provided and alignment file"
        )
        print(Fore.RED + "Exiting..." + Style.RESET_ALL)
        exit()
    else:
        pass

    # Save dataframes with gene options to .tsv files
    DataIdealALL.to_csv(
        os.path.join(OUT, r"IdealGroupGenes.tsv"),
        index=False,
        sep="\t",
        lineterminator="\n",
    ) if len(DataIdeal) > 0 else print("Optimal genes not found. Won't ceate file.")
    DataSecondALL.to_csv(
        os.path.join(OUT, r"SecondChoiceGroupGenes.tsv"),
        index=False,
        sep="\t",
        lineterminator="\n",
    ) if len(DataSecond) > 0 else print("Subotimal genes not found. Won't ceate file.")
    DataUndesiredALL.to_csv(
        os.path.join(OUT, r"UndesiredGroupGenes.tsv"),
        index=False,
        sep="\t",
        lineterminator="\n",
    ) if len(DataUndesired) > 0 else print(
        "Undesired genes not found. Won't ceate file."
    )

######################################################## Design Primers ###########################################################################

# Check if folder for individual primer3 files exists. If it doesn't, create it.
if os.path.exists(primer3_files):
    pass
else:
    os.mkdir(primer3_files)

# Design taxon-specific primers for target species defined by the -t flag
count = 1
found = 0
NoUnique = []

# Workflow to design taxon-specific unique primers
if args.primers_type == "unique":
    print(
        "Beginning to design taxon-specific primers for each target organism provided..."
    )

    # Sort unique genes file based on bacterial species and then descending gene length.
    df_sorted = UniqueGenes.sort_values(
        by=["species_name", "gene_length"], ascending=[True, False]
    )

    # filter out genes that are shorter than intended amplicon size
    print(
        f"\nFiltering for genes with length greater or equal to the longest intended amplicon size: {args.product_size_range[1]}\n"
    )
    df_sorted = df_sorted.loc[df_sorted["gene_length"] >= args.product_size_range[1]]

    for f in glob(cds_filenames):
        # Extract species name from CDS filenames
        Index = f.index("_cds")
        speciesName = os.path.basename(f[:Index])
        # Select all the unique genes rows for each target species
        species_rows = df_sorted[df_sorted["species_name"].str.contains(speciesName)]
        # Extract the top N (genes_number flag) rows of each species
        n_rows = species_rows.head(args.genes_number)

        # Check if number of unique genes is less than args.genes_number value
        if len(species_rows) == 0:
            # Skip that species
            print(
                Fore.RED
                + "{0}) Cannot design taxon-specific primers for: {1}. No unique genes found!".format(
                    count, speciesName
                )
                + Style.RESET_ALL
            )
            NoUnique.append(speciesName)
            count += 1
            continue
        elif len(species_rows) < args.genes_number:
            # re-define the args.genes_number
            genes_number = len(species_rows)
            found += 1
            print(
                Fore.YELLOW
                + "{0}) Designing taxon-specific primers for the top {1} unique genes of: {2} ".format(
                    count, genes_number, speciesName
                )
                + Style.RESET_ALL
            )
        else:
            genes_number = args.genes_number
            found += 1
            print(
                f"{count}) Designing taxon-specific primers for the top {genes_number} unique genes of: {speciesName}"
            )

        # Design primers with primer3
        for i in range(genes_number):
            i_length = n_rows.iloc[
                i, 2
            ]  # Return the 3rd column (i.e. the gene_length) for each row.
            i_bacterium = n_rows.iloc[
                i, 0
            ]  # Return the 1st column (i.e. the species_name) for each row.
            i_cds = n_rows.iloc[
                i, 1
            ]  # Return the 3rd column (i.e. the gene_name) for each row.
            record_dictionary = extract_seq_dict(cds_filenames, speciesName)
            i_sequence = extract_seq(i_cds, record_dictionary)
            filename = speciesName + "_" + i_cds
            # Generate N primers (determined by flag: --primers_number) using primer3
            primers = primer3.bindings.designPrimers(
                {
                    "SEQUENCE_ID": filename,
                    "SEQUENCE_TEMPLATE": i_sequence,
                    "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                },
                {
                    "PRIMER_OPT_SIZE": args.optimal_primer_size,
                    "PRIMER_MIN_SIZE": args.min_primer_size,
                    "PRIMER_MAX_SIZE": args.max_primer_size,
                    "PRIMER_OPT_TM": args.optimal_primer_Tm,
                    "PRIMER_MIN_TM": args.min_primer_Tm,
                    "PRIMER_MAX_TM": args.max_primer_Tm,
                    "PRIMER_PAIR_MAX_DIFF_TM": args.max_Tm_diff,
                    "PRIMER_MIN_GC": args.min_primer_gc,
                    "PRIMER_MAX_GC": args.max_primer_gc,
                    "PRIMER_MAX_POLY_X": args.max_poly_x,
                    "PRIMER_PICK_RIGHT_PRIMER": 1,
                    "PRIMER_PICK_LEFT_PRIMER": 1,
                    "PRIMER_PICK_INTERNAL_OLIGO": 0,
                    "PRIMER_PRODUCT_SIZE_RANGE": args.product_size_range,
                    "PRIMER_GC_CLAMP": args.GC_clamp,
                    "PRIMER_NUM_RETURN": args.primers_number,
                },
            )

            # Create a name specific to the bacterium and CDS name
            filenameTSV = filename + ".tsv"
            # Save the dictionary as a .tsv file
            with open(os.path.join(primer3_files, filenameTSV), "w") as f:
                for key in primers.keys():
                    f.write("%s\t %s\n" % (key, primers[key]))
        count += 1

    # Exit with error message if we could not design any taxon-specific primers across the input community
    if len(glob(primer_filenames)) == 0:
        print(
            "\nWe could not design any taxon-specific primers for the input bacterial community"
        )
        print(Fore.RED + "Exiting..." + Style.RESET_ALL)
        os.rmdir(args.outdir)
        exit()
    else:
        print(
            f"\nCompiling a report with the taxon-specific primers of {found}/{count-1} input bacteria..."
        )
        if found != count - 1:
            namesNoUnique = ", ".join(NoUnique)
            print(f"We could not design taxon-specific primers for: {namesNoUnique}.\n")
        elif found == count - 1:
            pass
        # Run every file that finishes on .tsv and store output in a dictionary list
        list_of_df = []
        for f in glob(primer_filenames):
            if "cds" in f:
                Index = f.index("_cds")
            elif "peg" in f:
                Index = f.index("_peg")
            name = os.path.basename(f[:Index])
            list_of_df.append(mainUnique(f, UniqueGenes, name, cds_filenames))
    # Order dataframe by the string length of species and gene sequence.
    df_concatenated = pd.concat(list_of_df).sort_values(
        by=["species", "Full_sequence"], key=lambda x: x.str.len(), ascending=False
    )
    # Store combined primer table as a csv file
    df_concatenated.to_csv(
        os.path.join(OUT, r"UniquePrimerTable.tsv"), index=False, sep="\t"
    )

# Workflow to design taxon-specific group primers
elif args.primers_type == "group":
    print(
        "Beginning to design group-specific primers for the target organisms provided..."
    )

    if len(DataIdeal) > 0: # Design group primers using 'ideal' genes
        print(
            Fore.GREEN
            + "Found genes that amplify every target organism!"
            + Style.RESET_ALL
        )
        time.sleep(1)
        print(
            'Attempting to design primers with these genes (see file "IdealGroupGenes.tsv")'
        )
        df_sorted = DataIdealALL
        # Define which genes to design primers for. The step in the range must be different if df_sorted is from DataIdealALL or DataSecondALL
        step = numTargets


    elif (len(DataIdeal) == 0) and (len(DataSecondALL) > 0): # Design group primers using 'second choice' genes
        print(
            Fore.RED
            + "Could not find genes with perfect homology across all target organisms"
            + Style.RESET_ALL
        )
        time.sleep(2)
        print(
            'Proceeding with "second choice" candidate group genes (see file "SecondChoiceGroupGenes.tsv")'
        )
        df_sorted = DataSecondALL
        # Define which genes to design primers for. The step in the range must be different if df_sorted is from DataIdealALL or DataSecondALL
        step = 1

    elif (len(DataIdeal) == 0) and (len(DataSecondALL) == 0): # Could not identify neither 'ideal' nor 'second cohice' genes. Must exit without designing group primers.
        print(
            'We could not find suitable genes to design group primers for the target organisms requested)'
        )
        print(Fore.RED + "Exiting..." + Style.RESET_ALL)
        exit()


    # filter out genes that are shorter than intended amplicon size
    print(
        f"\nFiltering for genes with length greater or equal to the longest intended amplicon size: {args.product_size_range[1]}\n"
    )
    df_sorted = df_sorted.loc[df_sorted["SequenceIncluded_length"] >= args.product_size_range[1]]

    # Determine how many genes we can design primers for:
    if len(df_sorted) < args.genes_number * step:
        # re-define the args.genes_number
        genes_number = len(df_sorted)
    else:
        genes_number = args.genes_number

    print(df_sorted)
    for i in range(0, genes_number * step, step):
        print("i") 
        print(i) 
        # (start, stop, step). It's important to define step as the same number of target organisms.
        # Due to the nature of this group workflow, we are likely to select for the same genes across different organisms. When sorting by gene length,
        # the same genes (i.e. nucleotide sequences) across the target organisms (genes will be called differently) will be shown consecutively.
        # (same gene across organisms will have the same length and will thus be shown consecutively). Because we ensured only one alignment per target organism must be found,
        # We can confidently use the target of organisms as step to make sure that the same gene but just from a different organism will be used to design group primers.
        speciesName = df_sorted.iloc[i, 0] # Extract species name
        i_cds = df_sorted.iloc[i, 1]  # Return the 2nd column (i.e. the gene name)
        i_sequenceRange = df_sorted.iloc[i, 3]  # Return sequence range to design primers
        i_found=df_sorted.iloc[i, 6] # Return number of targets found
        i_missing=df_sorted.iloc[i, 7] # Return number of targets missing
        # Must now retrieve the CDS filename from the bacterial name:
        ## Define path and name of bacterium to then find which file exists in that path that starts with that bacterial name
        CDS_filename = CDS + "/" + f"{speciesName}*.fna"
        print("CDS_filename")
        print(CDS_filename)
        # Find the exact filename
        for file in glob(CDS_filename):
            print("file")
            print(file)
            CDS_filename = file  # This is the exact filename
            record_dictionary = extract_seq_dict(CDS_filename, speciesName)
            i_sequence = extract_seq(i_cds, record_dictionary)
            print("i_sequence")
            print(i_sequence)
            i_sequenceRangeTemplate = i_sequence[int(i_sequenceRange[0]):int(i_sequenceRange[1])]
            print("i_sequenceRangeTemplate")
            print(i_sequenceRangeTemplate)

            filename = speciesName + "_" + i_cds
            # Generate N primers (determined by flag: --primers_number) using primer3
            primers = primer3.bindings.designPrimers(
                {
                    "SEQUENCE_ID": filename,
                    "SEQUENCE_TEMPLATE": i_sequence,
                    "SEQUENCE_INCLUDED_REGION": [int(i_sequenceRange[0]), int(i_sequenceRange[1]-i_sequenceRange[0])],
                },
                {
                    "PRIMER_OPT_SIZE": args.optimal_primer_size,
                    "PRIMER_MIN_SIZE": args.min_primer_size,
                    "PRIMER_MAX_SIZE": args.max_primer_size,
                    "PRIMER_OPT_TM": args.optimal_primer_Tm,
                    "PRIMER_MIN_TM": args.min_primer_Tm,
                    "PRIMER_MAX_TM": args.max_primer_Tm,
                    "PRIMER_PAIR_MAX_DIFF_TM": args.max_Tm_diff,
                    "PRIMER_MIN_GC": args.min_primer_gc,
                    "PRIMER_MAX_GC": args.max_primer_gc,
                    "PRIMER_MAX_POLY_X": args.max_poly_x,
                    "PRIMER_PICK_RIGHT_PRIMER": 1,
                    "PRIMER_PICK_LEFT_PRIMER": 1,
                    "PRIMER_PICK_INTERNAL_OLIGO": 0,
                    "PRIMER_PRODUCT_SIZE_RANGE": args.product_size_range,
                    "PRIMER_GC_CLAMP": args.GC_clamp,
                    "PRIMER_NUM_RETURN": args.primers_number,
                },
            )

            # Create a name specific to the bacterium and CDS name
            filenameTSV = filename + "-" + str(i_found) + "-" + str(i_missing) + ".tsv"
            # Save the dictionary as a .tsv file
            with open(os.path.join(primer3_files, filenameTSV), "w") as f:
                for key in primers.keys():
                    f.write("%s\t %s\n" % (key, primers[key]))

    # Exit with error message if we could not design any taxon-specific primers across the input community
    if len(glob(primer_filenames)) == 0:
        print(
            "\nWe could not design any group-specific primers for the input targets"
        )
        print(Fore.RED + "Exiting..." + Style.RESET_ALL)
        os.rmdir(args.outdir)
        exit()
    else:
        print("Compiling a report with the group-specific primers designed")

        # Run every file that finishes on .tsv and store output in a dictionary list
        list_of_df = []
        for f in glob(primer_filenames):
            fBase=os.path.basename(f)
            if "cds" in fBase:
                Index = fBase.index("_cds")
                splits = re.split(r'_cds|-|.tsv', fBase)
            elif "peg" in fBase:
                Index = fBase.index("_peg")
                splits = re.split(r'_peg|-|.tsv', fBase)
            name = splits[0]
            i_found = splits[2]
            i_missing = splits[3]
            list_of_df.append(mainGroup(f, df_sorted, name, cds_filenames, i_found,i_missing))
    # Order dataframe by the string length of species and gene sequence.
    df_concatenated = pd.concat(list_of_df).sort_values(
        by=["num_targets_found", "species"], ascending=False
    )
    # Store combined primer table as a .tsv file
    df_concatenated.to_csv(
        os.path.join(OUT, r"GroupPrimerTable.tsv"), index=False, sep="\t"
    )
print(Fore.GREEN + "Done!" + Style.RESET_ALL)