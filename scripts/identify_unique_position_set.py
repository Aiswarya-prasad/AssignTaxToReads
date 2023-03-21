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
                                        --list_strains "ESL0825" "ESL0822" "ESL0824" "ESL0827" "ESL0200" "ESL0197" "ESL0295" "ESL0185" "ESL0183" "ESL0184" "ESL0186" "ESL0820" "ESL0170" "ESL0198" "ESL0199" "ESL0819" "ESL0394" "ESL0261" "ESL0262" "ESL0350" "ESL0263"
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
from itertools import chain

def positions_are_unique(alns, test_positions):
    """
    alns = uniq_alignments
    test_positions = positions_chosen
    """
    are_unique = True
    sequences_seen = set()
    sequence_dict = {}
    for aln in alns:
        sequence_name = aln.id
        sequence = "".join([aln[i] for i in test_positions])
        sequence_dict[sequence_name] = sequence
    for name, seq in zip(sequence_dict.keys(), sequence_dict.values()):
        if seq in sequences_seen:
            # print(f"{name} not unique")
            are_unique = False
        sequences_seen.add(seq)
    return(are_unique)

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

def SDP_of_id(id):
    id = parse_id(id)
    return("_".join(id.split("_")[1:]))

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group("required arguments")
requiredNamed.add_argument("-d", "--database_path",metavar="output",required=True, help="Path to directory containing 16 seqeuences", action="store")
requiredNamed.add_argument("-i", "--input_fasta",metavar="output",required=True, help="Path to directory containing 16 seqeuences", action="store")
requiredNamed.add_argument("-o", "--output_file",metavar="output",required=True, help="Path to output file containing positional information", action="store")
requiredNamed.add_argument("-s", "--subset",metavar="output",required=True, help="Set to true (verbatim) if 16S multi fasta should be subset based on the list defined in the script", action="store")
parser.add_argument("-l", "--list_strains",nargs="+",metavar="output",required=True, help="List of genomes names to subset to", action="store")
args = parser.parse_args()

database_path = args.database_path
output_file = args.output_file
input_fasta = args.input_fasta
subset = args.subset
chosen_strains = args.list_strains

# database_path = "database/16S_sequences/"
# output_file = "database/16S_sequences/updated_unique_position_info.pickle"
# input_fasta = "database/16S_sequences/16S_sequences.fasta.all"
# subset = "True"
# chosen_strains = ["ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL0185", "ESL0183", "ESL0184", "ESL0186", "ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0394", "ESL0261", "ESL0262", "ESL0350", "ESL0263"]

sequences16S_file_all = input_fasta
sequences16S_file = os.path.join(database_path, "16S_sequences.fasta")
aligned_16S = os.path.join(database_path, "16S_aligned.fasta")

# chosen_strains = ["ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL0185", "ESL0183", "ESL0184", "ESL0186", "ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0394", "ESL0261", "ESL0262", "ESL0350", "ESL0263"]


if subset == "true":
    print(f"Creating subset {sequences16S_file} from {sequences16S_file_all} containing the strains:")
    print(f"{chosen_strains}")
    selected_records = []
    for record in SeqIO.parse(sequences16S_file_all, "fasta"):
        if strain_name_of_id(record.id) in chosen_strains:
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
alignment_length = alignments.get_alignment_length()
num_records = len(alignments)
# make a dictionary with a dictionary for each 16S sequences copy containing a set of positions
# that are not all shared and the character in the respective position
# positions_dict {ESLxxxx-1: {1: 'a', ...., 1235: 'g'}, ...., ESLxyyx-4: {1: 'a', ...., 1235: 't'}}

uniq_alignments = Align.MultipleSeqAlignment([])
for aln in alignments:
    if any(aln.seq == aln_iter.seq for aln_iter in uniq_alignments):
        pass
    else:
        uniq_alignments.append(aln)

sdp_list = set([SDP_of_id(aln.id) for aln in alignments])
positions_chosen_dict = {}
for sdp in sdp_list:
    # sdp = "firm4_1"
    positions_chosen_dict[sdp] = []
    sdp_uniq_alignments = Align.MultipleSeqAlignment([aln for aln in uniq_alignments if sdp in SDP_of_id(aln.id)])
    complete = False
    invalid_positions = False
    positions_to_consider = []
    for position in range(alignment_length):
        if len(set(sdp_uniq_alignments[:,position])) == 1:
            continue
        else:
            positions_to_consider.append(position)
    for num_positions in range(1, 1+len(positions_to_consider)):
        print(f"testing tuples of length {num_positions}")
        for i, test_positions in enumerate(combinations(positions_to_consider, num_positions)):
            if positions_are_unique(sdp_uniq_alignments, test_positions):
                complete = True
                # adding continue will also pick the earliest possible position
                # not doing so will pick the latest
                # continue
        if complete:
            positions_chosen_dict[sdp] = [pos for pos in test_positions]
            break

positions_chosen = []
for list_x in positions_chosen_dict.values():
    for x in list_x:
        if x not in positions_chosen:
            positions_chosen.append(x)

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

# check that these sequences uniquely identify strains
# if it is the same for a pair for strains, check that one of the other copies is ok
# check all pair-wise combinations
problem_tuples = []
name_combinations = combinations(check_positions_dict, 2)
for name, other_name in name_combinations:
    if strain_name_of_id(name) == strain_name_of_id(other_name):
        continue
    if check_positions_dict[name] == check_positions_dict[other_name]:
        print(f"{name} matches {other_name}")
        if SDP_of_id(name) == SDP_of_id(other_name):
            print("same SDP!")
            # if it is 183 and 262, this is expected
        else:
            if (SDP_of_id(name), SDP_of_id(other_name)) not in problem_tuples:
                problem_tuples.append((SDP_of_id(name), SDP_of_id(other_name)))

print(f"{len(positions_chosen)} positions chosen by sdp-level comparison")

# Since bifidos are too close, also do by phylotype
# This takes long so avoid other phylotypes 
# if sdp-level selection already did a good job
# not done a priori
# problem tuples are identified after the positions
# coming from sdp-level comparison reveal the cross
# sdp pairs that are not uniquely identified
positions_chosen_dict_phylo = {}
for phylo in ["firm5", "firm4", "bifido"]:
# work by pair so it is faster
    # sdps = [x for x in sdp_list if phylo in x]
    for sdps in combinations([x for x in sdp_list if phylo in x], 2):
        if sdps in problem_tuples:
            pass
        else:
            continue
        positions_chosen_dict_phylo[",".join(sdps)] = []
        phylo_uniq_alignments = Align.MultipleSeqAlignment([aln for aln in uniq_alignments if SDP_of_id(aln.id) in sdps])
        # [aln.id for aln in uniq_alignments if SDP_of_id(aln.id) in sdps]
        complete = False
        invalid_positions = False
        positions_to_consider = []
        for position in range(alignment_length):
            chars_in_pos = set(phylo_uniq_alignments[:,position])
            if len(chars_in_pos) == 1 or "-" in chars_in_pos:
                continue
            else:
                positions_to_consider.append(position)
        for num_positions in range(1, 1+len(positions_to_consider)):
            for i, test_positions in enumerate(combinations(positions_to_consider, num_positions)):
                if positions_are_unique(phylo_uniq_alignments, test_positions):
                    complete = True
                    continue
            if complete:
                positions_chosen_dict_phylo[",".join(sdps)] = [pos for pos in test_positions]
                break

for list_x in positions_chosen_dict_phylo.values():
    for x in list_x:
        if x not in positions_chosen:
            positions_chosen.append(x)

positions_dict = {}
for record in alignments:
    name = parse_id(record.id)
    positions_dict[name] = {}
for i in positions_chosen:
    for j, char in enumerate(alignments[:,i]):
                positions_dict[parse_id(alignments[j].id)].update({i: char})

print(f"{len(positions_chosen)} positions make a unique set")

check_positions_dict = {}
for x in positions_dict.keys():
    check_positions_dict[x] = "".join(positions_dict[x].values())

# check that these sequences uniquely identify strains
# if it is the same for a pair for strains, check that one of the other copies is ok
# check all pair-wise combinations
problem_tuples_strains = []
name_combinations = combinations(check_positions_dict, 2)
for name, other_name in name_combinations:
    if strain_name_of_id(name) == strain_name_of_id(other_name):
        continue
    if check_positions_dict[name] == check_positions_dict[other_name]:
        print(f"{name} matches {other_name}")
        if SDP_of_id(name) == SDP_of_id(other_name):
            print("same SDP!")
        else:
            problem_tuples_strains.append((name, other_name))

positions_chosen_dict_strains = {}
for problem_strains in problem_tuples_strains:
    positions_chosen_dict_strains[",".join(problem_strains)] = []
    strain_uniq_alignments = Align.MultipleSeqAlignment([aln for aln in alignments if parse_id(aln.id) in problem_strains])
    # [aln.id for aln in uniq_alignments if SDP_of_id(aln.id) in strains]
    complete = False
    invalid_positions = False
    positions_to_consider = []
    for position in range(alignment_length):
        chars_in_pos = set(strain_uniq_alignments[:,position])
        if len(chars_in_pos) == 1 or "-" in chars_in_pos:
            continue
        else:
            positions_to_consider.append(position)
    for num_positions in range(1, 1+len(positions_to_consider)):
        for i, test_positions in enumerate(combinations(positions_to_consider, num_positions)):
            if positions_are_unique(strain_uniq_alignments, test_positions):
                complete = True
                continue
        if complete:
            positions_chosen_dict_strains[",".join(problem_strains)] = [pos for pos in test_positions]
            break

for list_x in positions_chosen_dict_strains.values():
    for x in list_x:
        if x not in positions_chosen:
            positions_chosen.append(x)

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


print(f"{len(positions_chosen)} positions chosen")

name_combinations = combinations(check_positions_dict, 2)
for name, other_name in name_combinations:
    if strain_name_of_id(name) == strain_name_of_id(other_name):
        continue
    if check_positions_dict[name] == check_positions_dict[other_name]:
        print(f"{name} matches {other_name}")
        if SDP_of_id(name) == SDP_of_id(other_name):
            print("same SDP!")


with open(output_file, "wb") as outfile_fh:
    pickle.dump(positions_dict, outfile_fh, protocol=pickle.HIGHEST_PROTOCOL)

print(f"positions dictionary pickled at {output_file}")