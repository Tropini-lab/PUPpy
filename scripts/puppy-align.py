#!/usr/bin/env python3

# Script: puppy-align
# Part of PUPpy distribution.
# Many-against-many pairwise local sequence alignment of all microbial genes in input microbial community.

# Author = Hans Ghezzi
# Credits = Hans Ghezzi
# Licence = GPL3
# Version = 1.0.0
# Maintainer = Hans Ghezzi
# Email = hans.ghezzi@gmail.com
# Status = Development


# Load libraries
import subprocess
import sys
import pandas as pd
import argparse
from colorama import init
from colorama import Fore, Style
import glob, os
import shutil

######################################################## flags ###########################################################################

parser = argparse.ArgumentParser(
    prog="PROG",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
    PUPpy: A Phylogenetically Unique Primer pipeline in PYthon for high resolution classification of defined microbial communities. Version 1.0
    """,
)

# Setting up flags in python
parser.add_argument(
    "-in",
    "--intended",
    help="Directory with CDS files of intended targets, for which primers should be designed, in the defined microbial community",
    default="",
    required=True,
)

parser.add_argument(
    "-un",
    "--unintended",
    help="Directory with CDS files of unintended targets, for specificity checks, in the defined microbial community",
    default="",
    required=True,
)

parser.add_argument(
    "-o",
    "--outdir",
    help="Relative path to the output folder",
    default="Align_OUT",
)
parser.add_argument(
    "-i",
    "--identity",
    help="Identity thresholds to report sequence alignments by MMseqs2",
    default=0.3,
    type=float,
)

args = parser.parse_args()


######################################################## variables ###########################################################################

# Path to folder with CDS files
cds_intended = os.path.abspath(args.intended)
cds_unintended = os.path.abspath(args.unintended)

# Path to output directory
output = os.path.abspath(args.outdir)
# Check if output folder exists. If it does, exit the program, otherwise create it.
while True:
    if os.path.exists(output):
    # Get user input (yes or no)
        user_input = input("Do you want to OVERWRITE the data in the existing output directory? (yes/no): ").strip().lower()

        # Check the user's response
        if user_input == "yes":
            # Overwrite output directory
            print("Overwriting output directory...")
            directory_to_overwrite = output
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
        os.mkdir(output)
        break #Exit loop if output dir doesn't exist initially

######################################################## Functions ###########################################################################

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


##################################################### RENAME FASTA HEADERS #################################################################

# Execute function 'rename_fasta' to rename FASTA headers in input CDS files.

# List of directories with files to rename
directories = [cds_intended, cds_unintended]

for directory in directories:

    # Get list of files in current directory
    cds_filenames = os.path.join(directory, "*.fna")

    for f in glob.glob(cds_filenames):
        # Extract base name of the file
        name = os.path.basename(f)
        # Find index of "_cds" in filename
        i = name.index("_cds")
        # Function to rename file FASTA headers
        rename_fasta(f, name[:i])

##################################################### Alignments with MMseqs2 #################################################################

# Run MMseqs2.sh script to align all genes in the defined microbial community.
subprocess.run(["chmod", "755", "./MMseqs2.sh"]) # Give permissions to run script

# Define  log file name
log_file = 'align_logfile.txt'
# Specify log file output dir
log_file_path = os.path.join(output, log_file)

# Define the shell command with 'tee' to redirect output
shell_command = f"./MMseqs2.sh {output} {cds_intended} {cds_unintended} {str(args.identity)} | tee {log_file_path}"

# Run the shell command
subprocess.check_call(shell_command, shell=True)

# Print successful completion message
print(Fore.GREEN + "Done!" + Style.RESET_ALL)