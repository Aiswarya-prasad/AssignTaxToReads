#!/env/python3

"""
This script first subsets the 16S sequences of the database to just get the ones defined by 
"pilot strains" if the subset option is set to true and then calculates the set of positions 
that can be used to identify each strain uniquely

For now, I don't get the minimal possible set instead I get positions (including those with
at least < ) that have 3 or 4 different characters depending on which ones yield unique
combinations. I also exclude positions with > a given % gaps. 
The number of unique characters and % of gaps to ask for needs some playing around
to ensure correctness. Eg. I ask for 3 or more and no more than 1 gap

These positions are in a dictionary which is pickled and later used by the script to count
SDPs subsequently
The format of the name of headers in the 16S sequences database files is expected to be:
<strain name>-<# serial number for copy>_<SDP name>.fa
SDP name may contain "_" and "." but not "-"
eg. ESL252-1_firm5_1
Example:
python3 scripts/identify_unique_position_set.py --database_path "database/16S_sequences/" \
                                        --input_fasta "database/16S_sequences/16S_sequences.fasta.all" \
                                        --output_file "database/16S_sequences/unique_position_info.pickle" \
                                        --subset true
"""

import os
import argparse
import gzip
import subprocess
import pickle
from Bio import SeqIO
from Bio import AlignIO
from Bio import Align
from io import StringIO
from collections import Counter
from itertools import combinations

def parse_id(id):
    if "_R_" in id:
        name = id.split("_R_")[1]
        name = name.split(":")[0]
    else:
        name = id.split(":")[0]
    return(name)

def strain_name_of_id(id):
    id = parse_id(id)
    return(id.split("-")[0])

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group("required arguments")
requiredNamed.add_argument("-d", "--database_path",metavar="output",required=True, help="Path to directory containing 16 seqeuences", action="store")
requiredNamed.add_argument("-i", "--input_fasta",metavar="output",required=True, help="Path to directory containing 16 seqeuences", action="store")
requiredNamed.add_argument("-o", "--output_file",metavar="output",required=True, help="Path to output file containing positional information", action="store")
requiredNamed.add_argument("-s", "--subset",metavar="output",required=True, help="Set to true (verbatim) if 16S multi fasta should be subset based on the list defined in the script", action="store")
args = parser.parse_args()

database_path = args.database_path
output_file = args.output_file
input_fasta = args.input_fasta
subset = args.subset

sequences16S_file_all = input_fasta
sequences16S_file = os.path.join(database_path, "16S_sequences.fasta")
aligned_16S = os.path.join(database_path, "16S_aligned.fasta")

pilot_strains = ["ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL0185", "ESL0183", "ESL0184", "ESL0186", "ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0394", "ESL0261", "ESL0262", "ESL0350", "ESL0263"]


if subset == "true":
    print(f"Creating subset {sequences16S_file} from {sequences16S_file_all} containing the strains:")
    print(f"{pilot_strains}")
    selected_records = []
    for record in SeqIO.parse(sequences16S_file_all, "fasta"):
        if strain_name_of_id(record.id) in pilot_strains:
            selected_records.append(record)
    SeqIO.write(selected_records, sequences16S_file, "fasta")
else:
    print(f"Creating {sequences16S_file} with all sequences contained in {sequences16S_file_all}")
    SeqIO.write([record for record in sequences16S_file_all], sequences16S_file, "fasta")

# align 16_sequences from multi fasta
if not os.path.exists(aligned_16S):
    print(f"no aligned 16S sequences found in {database_path}")
    print(f"running alignment in mafft for sequences in {sequences16S_file}")
    process = subprocess.Popen(["mafft", "--auto", "--adjustdirection", sequences16S_file],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = process.communicate()
    AlignIO.write(AlignIO.read(StringIO(stdout), "fasta"), aligned_16S, "fasta")
else:
    print(f"aligned 16S sequences found at {aligned_16S}")

# read alignments
alignments = AlignIO.read(aligned_16S, "fasta")

# get length of alignment
nb_positions = alignments.get_alignment_length()
num_records = len(alignments)
# make a dictionary with a dictionary for each 16S sequences copy containing a set of positions
# that are not all shared and the character in the respective position
# positions_dict {ESLxxxx-1: {1: 'a', ...., 1235: 'g'}, ...., ESLxyyx-4: {1: 'a', ...., 1235: 't'}}


positions_chosen = []
for i in range(nb_positions):
    bad_chars = False
    if i in [493, 1021]:
        # I added this because I knew some positions a priori
        positions_chosen.append(i)
        continue
    for char in alignments[:,i]:
        if char not in ["a", "t", "g", "c", "-"]:
            bad_chars = True
    if bad_chars:
        continue
    # avoid freq filtering
    freq_dict = Counter(alignments[:,i])
    # if more than 70% is - do not include
    # if freq_dict["-"] > 0.7*num_records:
    if freq_dict["-"] > 1:
        continue
    if len(freq_dict.values()) < 3:
        continue
    else:
        positions_chosen.append(i)

positions_dict = {}
for record in alignments:
    name = parse_id(record.id)
    positions_dict[name] = {}
for i in positions_chosen:
    for j, char in enumerate(alignments[:,i]):
                positions_dict[parse_id(alignments[j].id)].update({i: char})

print(f"{len(positions_chosen)} positions make a potentially unique set")

check_positions_dict = {}
for x in positions_dict.keys():
    check_positions_dict[x] = "".join(positions_dict[x].values())
    # "".join(positions_dict[x].values())

# check that these sequences uniquely identify strains
# if it is the same for a pair for strains, check that one of the other copies is ok
# check all pair-wise combinations
name_combinations = combinations(check_positions_dict, 2)
for name, other_name in name_combinations:
    if strain_name_of_id(name) == strain_name_of_id(other_name):
        continue
    if check_positions_dict[name] == check_positions_dict[other_name]:
        print(f"{name} matches {other_name}")

print(f"{len(positions_chosen)} positions chosen")

print("other combinations are unique")
        
# for x in positions_dict.keys():
#     {x: "".join([x for x in positions_dict[x].values()])}

with open(output_file, "wb") as outfile_fh:
    pickle.dump(positions_dict, outfile_fh, protocol=pickle.HIGHEST_PROTOCOL)

print(f"positions dictionary pickled at {output_file}")