#!/usr/bin/env python

############################################################
# Argument Options
# Header needs to have phylum -> species all separated by '_'

import argparse
parser = argparse.ArgumentParser("Reads and writes each entry as a single file.")
parser.add_argument("-i1",
                    action = "store", 
                    dest = "input1", 
                    metavar = "input1",
                    help = "[REQUIRED] Master files", 
                    required = True)
parser.add_argument("-i2",
                    action = "store", 
                    dest = "input2", 
                    metavar = "input2",
                    help = "[REQUIRED] Tables to add", 
                    required = True)
options = parser.parse_args()

############################################################

from Bio import SeqIO, AlignIO
import pandas as pd


outfile_table    = open("ASVs_counts_combined.txt", "w")
outfile_taxonomy = open("ASVs_taxonomy_combined.txt", "w")
outfile_fasta    = open("ASVs_combined.fasta", "w")

# 1
infile_fasta_1 = open(options.input1.split(",")[0], "r")
infile_table_1 = open(options.input1.split(",")[1], "r")
infile_taxonomy_1 = open(options.input1.split(",")[2], "r")

UID2seq_1 = {}
for record in SeqIO.parse(infile_fasta_1, "fasta"):
    UID2seq_1[str(record.description)] = str(record.seq)
seq2UID_1 = {v: k for k, v in UID2seq_1.items()}

table_1 = pd.read_csv(infile_table_1, sep = "\t")
table_1 = table_1.rename(index = UID2seq_1)

taxonomy_1 = pd.read_csv(infile_taxonomy_1, sep = "\t", names = ["id", "taxonomy", "score"])
taxonomy_1.set_index("id", inplace = True)
taxonomy_1 = taxonomy_1.rename(index = UID2seq_1)


# 2
infile_fasta_2 = open(options.input2.split(",")[0], "r")
infile_table_2 = open(options.input2.split(",")[1], "r")
infile_taxonomy_2 = open(options.input2.split(",")[2], "r")

UID2seq_2 = {}
for record in SeqIO.parse(infile_fasta_2, "fasta"):
    UID2seq_2[str(record.description)] = str(record.seq)
seq2UID_2 = {v: k for k, v in UID2seq_2.items()}

table_2 = pd.read_csv(infile_table_2, sep = "\t")
table_2 = table_2.rename(index = UID2seq_2)

taxonomy_2 = pd.read_csv(infile_taxonomy_2, sep = "\t", names = ["id", "taxonomy", "score"])
taxonomy_2.set_index("id", inplace = True)
taxonomy_2 = taxonomy_2.rename(index = UID2seq_2)


# Join two tables
table_joined = table_1.join(table_2, how = "outer")


# Dictionary with sequence:UID (new)
seqsBeingAdded = set(UID2seq_2.values()) - set(UID2seq_1.values()) # New sequences in table_2
startingUIDIndex = max([int(i.split("ASV_")[1]) for i in seq2UID_1.values()]) + 1
seq2newUID = {}
for seq in seqsBeingAdded:
    newUID = "ASV_%010d" % startingUIDIndex
    startingUIDIndex += 1
    seq2newUID[seq] = newUID

# Final seq2UID
seq2UID_updated = {**seq2UID_1, **seq2newUID}

# Rename Table
table_joined = table_joined.rename(index = seq2UID_updated)
table_joined = table_joined.fillna(0)
table_joined.sort_index(inplace=True)


# Deal with taxonomy
taxonomyBeingAdded = taxonomy_2[taxonomy_2.index.isin(seqsBeingAdded)]
taxonomy_joined = pd.concat([taxonomy_1, taxonomyBeingAdded])
taxonomy_joined = taxonomy_joined.rename(index = seq2UID_updated)
taxonomy_joined.sort_index(inplace=True)


#############
## Outputs ##
#############

# Output ASVs_counts_combined.txt
table_joined.to_csv("ASVs_counts_combined.txt", index = True, sep = "\t")

# Output ASVs_counts_combined.txt
taxonomy_joined.to_csv("ASVs_taxonomy_combined.txt", index = True, sep = "\t", header = False)

# Output ASV_combined.fasta
sorted(seq2UID_updated.items(), key=lambda x: x[1])
for key, value in seq2UID_updated.items():
    outfile_fasta.write(">" + value + "\n")
    outfile_fasta.write(key + "\n")

exit(0)



table_joined = table_joined.rename(index = seq2newUID)

taxonomy_joined = taxonomy_joined.rename(index = seq2newUID)
print(taxonomy_joined)


exit(0)
#taxonomy_joined = pd.merge(taxonomy_1, taxonomy_2, how = "left")

#.fillna(0)

print(table_joined)
print(table_joined2)


print(taxonomy_1)

print(taxonomy_joined)
exit(0)

taxonomy_joined = taxonomy_joined.rename(index = seq2UID_1)


print(difference)

# print(taxonomy_1)
# print(taxonomy_2)
# print(taxonomy_joined)
exit(0)

print(taxonomy_joined)


print(newASVstartingIndex)



exit(0)



print(seq2UID_1.values())





exit(0)

print(UID2seq_1.values())
print(UID2seq_2.values())


print("")
print(difference)
exit(0)



print(table_1)
print(table_2)


print(result)

result1 = result.rename(index = seq2UID_1)

print(result1)

exit(0)


#outfile = open(options.outfile, "w")






infile.close()
outfile.close()
