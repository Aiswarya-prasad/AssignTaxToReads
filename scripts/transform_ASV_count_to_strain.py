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

python3 <path-to-script-dir>/transform_ASV_count_to_strain.py --database_dir <path-to-directory> \
                                                   --input_table <path-to-reads> \
                                                   --output_path <path-to-output-directory> \
                                                   --temp_dir <path-to-temp-dir>
Example:

python3 scripts/transform_ASV_count_to_strain.py --database_path "database/16S_sequences/" \
                                      --input_table "04_Dada2ResultTaxonomy/ASVs_pilot_samples.csv" \
                                      --output_path "04_Dada2ResultTaxonomy" \
                                      --temp_dir "04_Dada2ResultTaxonomy/temp"
"""

import os
import argparse
import gzip
import subprocess
import pickle
import pandas as pd
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

def get_relaxed_match(ref_base, read_base):
    """
    match is true when both are one of the 4 valid bases
    if the reference is - then the read being - is important
    otherwise read being - could mean that the read was shorter and gaps were
    added or anything else like that so in that case, we want to consider that
    a match to the base being compared
    i.e. - will always be considered a match so the ref will not be excluded from
    potential matches, however if the ref - 
    To ensure that mismatch is considered only in the case that the reference is
    -
    Use this approach if too many are unassigned
    """
    if read_base == "-":
        return(True)
    else:
        if read_base == ref_base:
            return(True)
        else:
            return(False)

# AlignIO.write(this_alignment, "temp.aln", "fasta")
def get_strain_matches(aln, positions_dict_given):
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
    if max_matches < len(positions_chosen)*0.98:
        return("NA")
        # return(strains_matched.append("low"))
    return(strains_matched)
        
    # todo come back here and get the name of the strain the positions correspond to
    
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group("required arguments")
requiredNamed.add_argument("-d", "--database_path",metavar="output",required=True, help="Path to directory containing 16 seqeuences", action="store")
requiredNamed.add_argument("-i", "--input_table",metavar="output",required=True, help="Path to file ASVs in colnames and counts by sample", action="store")
requiredNamed.add_argument("-o", "--output_path",metavar="output",required=True, help="Path to output directory", action="store")
requiredNamed.add_argument("-t", "--temp_dir",metavar="output",required=True, help="Directory to save temporary output", action="store")
args = parser.parse_args()

database_path = args.database_path
ASV_table = args.input_table
output_path = args.output_path
temp_dir = args.temp_dir

database_path = "database/16S_sequences/"
ASV_table = "04_Dada2ResultTaxonomy/ASVs_pilot_samples.csv"
output_path = "04_Dada2ResultTaxonomy"
temp_dir = "04_Dada2ResultTaxonomy/temp"

with open(os.path.join(database_path, "unique_position_info.pickle"), "rb") as dict_fh:
    positions_dict = pickle.load(dict_fh)

sequences16S_file = os.path.join(database_path, "16S_sequences.fasta")
aligned_16S = os.path.join(database_path, "16S_aligned.fasta")
outfile_path = os.path.join(output_path, "ASV_to_strain_name_counts.csv")
taxtable_path = os.path.join(output_path, "Taxtable_strains.csv")
summary_file_path = os.path.join(output_path, "ASV_to_strain_summary.txt")

strain_names = [x for x in positions_dict.keys()]
# read ASVs from ASV table as seq records their id is the same as their sequence (so it matches the otu table)
ASV_table_df = pd.read_csv(ASV_table)
ASV_table_df.rename(columns={ ASV_table_df.columns[0]: "sample_name" }, inplace = True)
# create dataframe to store counts
strain_16S_df = pd.DataFrame(ASV_table_df["sample_name"])
ASV_records = []
for i, ASV_sequence in enumerate(ASV_table_df.columns[1:]):
    record = SeqIO.SeqRecord(Seq(ASV_sequence, generic_dna), 
                             id="ASV"+str(i))
    ASV_records.append(record)

ASV_table_df.rename(columns={record.seq:record.id for record in ASV_records}, inplace=True)
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

with open(summary_file_path, "w") as summary_fh, open(taxtable_path, "w") as tax_fh:
    tax_fh.write(",".join(["","Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"])+"\n")
    shorter = 0
    longer = 0
    within_range = 0
    total_ASV_sequences = 0
    low_matches = 0
    assigned_unambiguous = 0
    assigned_ambiguous = 0
    # alignment_length = AlignIO.read(aligned_16S, "fasta").get_alignment_length()
    for record in ASV_records:
        print(f"{record.id}")
        summary_fh.write(f"{record.id}\n")
        total_ASV_sequences += 1
        sequence_length = len(record.seq)
        if sequence_length > 1562:
            print(f"{sequence_length} longer than longest")
            summary_fh.write(f"{sequence_length} longer than longest\n")
            longer += 1
        if sequence_length < 1522:
            print(f"{sequence_length} shorter than shortest")
            summary_fh.write(f"{sequence_length} shorter than shortest\n")
            shorter += 1
        if sequence_length <= 1562 and sequence_length >= 1522:
            print(f"{sequence_length} within range")
            summary_fh.write(f"{sequence_length} within range\n")
            within_range += 1 
        # if total_ASV_sequences < 900:
        #     continue
        # if total_ASV_sequences == 1:
        #     break
        this_read_file = os.path.join(temp_dir, "temp.fasta")
        succes = SeqIO.write([record], this_read_file, "fasta")
        this_id = record.id
        # If the --keeplength option is given, then the alignment length is unchanged.  Insertions at the new sequences are deleted. 
        # todo check if adjust direction is truly working and then visualise some alns to get an idea
        process = subprocess.Popen(["mafft", "--adjustdirection", "--add", this_read_file, aligned_16S],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        this_alignment = AlignIO.read(StringIO(stdout), "fasta")
        for record in this_alignment:
            record.id = parse_id(record.id)
        matched_names = get_strain_matches(this_alignment, positions_dict)
        matched_strains = [strain_name_of_id(x) for x in matched_names]
        matched_strains_set = set(matched_strains)
        if len(matched_strains) == 1:
            assigned_unambiguous += 1
        else:
            if "low" in matched_strains:
                low_matches += 1
            else:
                assigned_ambiguous += 1
        print(f"{matched_names}")
        tax_prefix = ",".join(["NA","NA","NA","NA","NA","NA"])
        if any("firm" in x for x in matched_names):
            tax_prefix = ",".join(["Bacteria","Firmicutes","Bacilli","Lactobacillales","Lactobacillaceae","Lactobacillus"])
        if any("bifido" in x for x in matched_names):
            tax_prefix = ",".join(["Bacteria","Actinobacteriota","Actinobacteria","Bifidobacteriales","Bifidobacteriaceae","Bifidobacterium"])
        tax_strain = ";".join([x for x in matched_strains_set])
        seq_formatted = str([x.seq for x in ASV_records if x.id == record.id][0])
        tax_fh.write(",".join([seq_formatted,tax_prefix,SDP_of_id(matched_names[0]),tax_strain])+"\n")
        summary_fh.write(f"{matched_names}\n")
        strains_matched = set([strain_name_of_id(x) for x in matched_names])
        counts_list = list(ASV_table_df[record.id])
        if len(strains_matched) == 1:
            if len(matched_names) == 1:
                # matched uniquely to 1 strain and 1 16S sequence copy
                print(f"easy unique match")
                summary_fh.write(f"easy unique match\n")
                for name in matched_names:
                    if name in strain_16S_df.keys():
                        strain_16S_df[name] = [x + y for x, y in zip(list(strain_16S_df[name]), counts_list)]
                    else:
                        strain_16S_df[name] = counts_list
            else:
                # matched to 1 strain but it has {len(matched_names)} 16S sequence copies
                print(f"easy unique strain multi-copy")
                summary_fh.write(f"easy unique strain multi-copy\n")
                counts_list = [x/len(matched_names) for x in counts_list]
                for name in matched_names:
                    if name in strain_16S_df.keys():
                        strain_16S_df[name] = [x + y for x, y in zip(list(strain_16S_df[name]), counts_list)]
                    else:
                        strain_16S_df[name] = counts_list
        else:
            if set(["ESL0183", "ESL0262"]) == strains_matched:
                # matched to more than 1 strain and because of ESL0183-ESL0262
                print(f"multi-strain multi-copy ESL0262")
                summary_fh.write(f"multi-strain multi-copy ESL0262\n")
                for name in matched_names:
                    if name in strain_16S_df.keys():
                        strain_16S_df[name] = [x + y for x, y in zip(list(strain_16S_df[name]), counts_list)]
                    else:
                        strain_16S_df[name] = counts_list
            else:
                # matched to more than 1 strain and not due to ESL0183-ESL0262
                # add these to "unknown"
                print(f"unknown")
                summary_fh.write(f"unknown\n")
                name = "unknown"
                if name in strain_16S_df.keys():
                    strain_16S_df[name] = [x + y for x, y in zip(list(strain_16S_df[name]), counts_list)]
                else:
                    strain_16S_df[name] = counts_list
        print("--------------------------------------------------------------")
        summary_fh.write("--------------------------------------------------------------\n")

    summary_fh.write("--------------------------------------------------------------\n")
    summary_fh.write("--------------------------------------------------------------\n")
    summary_fh.write("".join([
        "total_ASV_sequences =", str(total_ASV_sequences),"\n",
        "shorter =", str(shorter),"\n",
        "longer =", str(longer), "\n",
        "within_range =", str(within_range), "\n",
        "total_ASV_sequences =", str(total_ASV_sequences), "\n",
        "low_matches =", str(low_matches), "\n",
        "assigned_unambiguous =", str(assigned_unambiguous), "\n",
        "assigned_ambiguous =", str(assigned_ambiguous)
    ]))

strain_16S_df.to_csv(outfile_path)