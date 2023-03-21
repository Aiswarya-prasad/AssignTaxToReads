#!/env/python3

"""
For this script I start with taxtable made from dada2.
The script assumes that identify_unique_position_set.py is already run and the 
positions dictionary pickled (in database/16S_sequences/unique_position_info.pickle).
This script reads the sequences from a the first column of taxtable and
all the 16S sequences of strains used in the experiment and assigns
the ASV sequence to the closest match using position specific information pickled in 
database/16S_sequences/unique_position_info.pickle

Start with all the strain names in a list. Each time a position of the ASV
is not the same as one of the refs, pop it from the list until all positions are covered.
The output is a csv file where each row is an ASV with the next column being a list of
strains that it matches to (later this can be parsed to make a tax table for easy handling in phyloseq)

Those with multiple copies will be recounted to the mean of all the copies [in a later script].
There are a few cases in which an SDP has multiple copies one if which matches that of another.
In such a case, the count will have to be re-organized based on the other counts of other copies
and that of the single-copy SDP will need to be adjusted. Such operations with be performed 
by a subsequent script and all the samples will be finally merged.

The format of the name of headers in the 16S sequences database files is expected to be:
<strain name>-<# serial number for copy>_<SDP name>.fa
SDP name may contain "_" and "." but not "-"
eg. ESL252-1_firm5_1

Usage:

python3 <path-to-script-dir>/assign_reads_to_strains.py --database_path <path-to-directory> \
                                                   --reads_file <path-to-reads-file> \
                                                   --summary_file_path <path-to-summary-file>
                                                   --outfile_path <path-to-output-dir> \
                                                   --sample <sample-name> \
                                                   --log_file <path-to-log-file> \
                                                   --match_id_cutoff <proportion>
Example:

python3 scripts/assign_reads_to_strains.py --database_path "database/16S_sequences/" \
                                           --reads_file "01_ReadsRenamed/2-12_reads.fastq.gz" \
                                           --summary_file_path "03_assign_reads_to_strain/2-12_summary.txt" \
                                           --outfile_path "03_assign_reads_to_strain/2-12_strain_counts.csv" \
                                           --sample "2-12" \
                                           --log_file  "03_assign_reads_to_strain/2-12_assign_reads_to_strains.log" \
                                           --match_id_cutoff 0.95
"""

import os
import shutil
import sys
import argparse
import gzip
import subprocess
import pickle
import pandas as pd
import math
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import Align
from io import StringIO


def parse_id(id):
    if "_R_" in id:
        name = id.split("_R_")[1]
        name = name.split(":")[0]
    else:
        name = id.split(":")[0]
    return(name)

def SDP_of_id(id):
    id = parse_id(id)
    return("_".join(id.split("_")[1:]))

def strain_name_of_id(id):
    id = parse_id(id)
    return(id.split("-")[0])

# todo add a sanity check to confirm that the other positions match exactly
# AlignIO.write(this_alignment, "temp.aln", "fasta")
def get_strain_matches(aln, positions_dict_given, match_id_cutoff):
    """
    This function reads the ASV sequence aligned at positions
    of interest. At each position it keeps track of the strains that
    the ASV matches to. Finally a list of the strain(s) with maximum
    matches is returned. If the number of matches is < half of the
    positions, then the returned list has the maximum matches strains
    and the string "low".
    since the read was added by mafft, it is always the last
    record. If this ever changes, edit this function accordingly
    aln = this_alignment
    positions_dict_given = positions_dict
    match_id_cutoff = 0.97
    """ 
    read_sequence = aln[-1].seq
    positions_chosen = [x for x in  positions_dict_given[list(positions_dict_given.keys())[0]].keys()]
    names = [x for x in positions_dict_given.keys()]
    matches = {}
    for name in names:
        matches[name] = 0
    for position in positions_chosen:
        base_char = read_sequence[position]
        for name in names:
            if positions_dict_given[name][position] == base_char:
                matches[name] += 1
            else:
                pass
                # list_ref_bases = "".join([positions_dict_given[x][position] for x in positions_dict_given.keys()])
                # print(f"{position}, {base_char} : {list_ref_bases}")
                # print(f"{positions_dict_given[name][position]} and {base_char} in {name}")
    max_matches = max(matches.values())
    strains_matched = [x for x in matches.keys() if matches[x] == max_matches]
    # match_id_cutoff can be 1 because only 16 important positions are chosen
    if max_matches < len(positions_chosen)*float(match_id_cutoff):
        return(["unknown"])
    return(strains_matched)

def progressbar(log_file,current_value,total_value,bar_lengh,progress_char):
    percentage = math.ceil((current_value/total_value)*100)                                                # Percent Completed Calculation 
    progress = math.ceil((bar_lengh * current_value ) / total_value)                                       # Progress Done Calculation 
    loadbar = "Progress: [{:{len}}]{}%".format(progress*progress_char,percentage,len = bar_lengh)    # Progress Bar String
    if log_file == "NA":
        print("\n")
        print(loadbar, end="\r")
    else:
        with open(log_file, "w") as log_fh:
            log_fh.write("\n")
            log_fh.write(loadbar+"\n")

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group("required arguments")
requiredNamed.add_argument("-r", "--reads_file",metavar="output",required=True, help="Path to fastq.gz file with reads for sample", action="store")
requiredNamed.add_argument("-d", "--database_path",metavar="output",required=True, help="Path to directory containing 16 seqeuences", action="store")
requiredNamed.add_argument("-s", "--summary_file_path",metavar="output",required=True, help="Summary file to write summary in", action="store")
requiredNamed.add_argument("-o", "--outfile_path",metavar="output",required=True, help="Path to output file", action="store")
requiredNamed.add_argument("-n", "--sample_name",metavar="output",required=True, help="Name of sample to be mentioned in output", action="store")
requiredNamed.add_argument("-l", "--log_file",metavar="output",required=True, help="Path to log file to write progress into", action="store")
requiredNamed.add_argument("-m", "--match_id_cutoff",metavar="output",required=True, help="Minimum percentage identity (proportion between 0 and 1) to reference positions required", action="store")
args = parser.parse_args()

database_path = args.database_path
input_reads_file = args.reads_file
summary_file_path = args.summary_file_path
outfile_path = args.outfile_path
sample = args.sample_name
log_file_path = args.log_file
match_id_cutoff = args.match_id_cutoff

# database_path = "database/16S_sequences/"
# input_reads_file = "01_ReadsRenamed/2-12_reads.fastq.gz"
# summary_file_path = "temp.summary"
# outfile_path = "temp_strain.csv"
# sample = "2-12"
# log_file_path = "temp.log"
# match_id_cutoff = "0.97"

temp_dir = os.path.join(os.path.dirname(outfile_path), sample+"_temp_files")
# temp_dir = "temp"

try:
    os.makedirs(temp_dir)
except:
    print(f"{temp_dir} exists")

with open(os.path.join(database_path, "unique_position_info.pickle"), "rb") as dict_fh:
    positions_dict = pickle.load(dict_fh)

aligned_16S = os.path.join(temp_dir, "temp_16S_aligned.fasata")
shutil.copyfile(os.path.join(database_path, "16S_aligned.fasta"), aligned_16S)

# # ESL0350 has only 1312 - get complete sequence later
# # Longest is 1562
# # shortest is 1522
# lengths = []
# for record in AlignIO.read(aligned_16S, "fasta"):
#     num_bases = 0
#     for char in record.seq:
#         if char in ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']:
#             num_bases +=1
#     # print(record.id)
#     print(num_bases)
#     lengths.append(num_bases)
# max(lengths)
# lengths.remove(1312)
# min(lengths)
#    
strain_counts_df = pd.DataFrame()
strain_counts_df["strain_16S_copy"] = [x for x in positions_dict.keys()]+ ["unknown"]
with open(summary_file_path, "w") as summary_fh:
    counts_dict = {}
    for name in positions_dict.keys():
        counts_dict[name] = 0
    counts_dict["unknown"] = 0
    shorter = 0
    longer = 0
    within_range = 0
    total_reads = 0
    unknown = 0
    assigned_unambiguous = 0
    assigned_ambiguous = 0
    read_counter = 0
    num_reads = 0
    with gzip.open(input_reads_file, "rt") as reads_fh:
        for record in SeqIO.parse(reads_fh, "fastq"):
            num_reads +=1
    with gzip.open(input_reads_file, "rt") as reads_fh:
        for record in SeqIO.parse(reads_fh, "fastq"):
            # print(f"{record.id}")
            # summary_fh.write(f"{record.id}\n")
            total_reads += 1
            sequence_length = len(record.seq)
            if sequence_length > 1562:
                # print(f"{sequence_length} longer than longest")
                # summary_fh.write(f"{sequence_length} longer than longest\n")
                longer += 1
            if sequence_length < 1522:
                # print(f"{sequence_length} shorter than shortest")
                # summary_fh.write(f"{sequence_length} shorter than shortest\n")
                shorter += 1
            if sequence_length <= 1562 and sequence_length >= 1522:
                # print(f"{sequence_length} within range")
                # summary_fh.write(f"{sequence_length} within range\n")
                within_range += 1
            this_read_file = os.path.join(temp_dir, "temp.fasta")
            success = SeqIO.write([record], this_read_file, "fasta")
            this_id = record.id
            # If the --keeplength option is given, then the alignment length is unchanged.  Insertions at the new sequences are deleted. 
            process = subprocess.Popen(["mafft", "--adjustdirection", "--add", this_read_file, aligned_16S],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            this_alignment = AlignIO.read(StringIO(stdout), "fasta")
            for record in this_alignment:
                record.id = parse_id(record.id)
            matched_names = get_strain_matches(this_alignment, positions_dict, match_id_cutoff)
            matched_strains = [strain_name_of_id(x) for x in matched_names]
            matched_strains_set = set(matched_strains)
            if "unknown" in matched_strains:
                # print(f"{record.id} : unkown")
                unknown += 1
            else:
                if len(matched_strains) == 1:
                    assigned_unambiguous += 1
                else:
                    assigned_ambiguous += 1
            # print(f"{matched_names}")
            # summary_fh.write(f"{matched_names}\n")
            counts_dict.values()
            for name in matched_names:
                counts_dict[name] +=1
            strain_counts_df[sample] = counts_dict.values()
            # print("--------------------------------------------------------------")
            # summary_fh.write("--------------------------------------------------------------\n")
            progressbar(log_file_path,read_counter,num_reads,50,'■')
            read_counter +=1
    summary_fh.write("--------------------------------------------------------------\n")
    summary_fh.write("--------------------------------------------------------------\n")
    summary_fh.write("".join([
        "total_reads =", str(total_reads),"\n",
        "shorter =", str(shorter),"\n",
        "longer =", str(longer), "\n",
        "within_range =", str(within_range), "\n",
        "total_reads =", str(total_reads), "\n",
        "unknown =", str(unknown), "\n",
        "assigned_unambiguous =", str(assigned_unambiguous), "\n",
        "assigned_ambiguous =", str(assigned_ambiguous)
    ]))
strain_counts_df.to_csv(outfile_path)
print("removing temoporary directory")
shutil.rmtree(temp_dir)