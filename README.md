[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


# HONEYPI
### :sunflower: -> :honeybee: -> :honey_pot: -> :+1:

### An automated pipeline to process ITS2 sequences from honey and plants sequences from National Honey Monitoring Scheme ###

** Current version: 1.0 (updated 2020-03-31) **


## A. Installation

### HONEYPI

In your home directory, copy and paste the following (line by line):

```
cd ~
conda create -n honeypi_env -y python=3.6 progressbar2 requests rdptools itsx vsearch trim-galore bbmap seqkit pandas biopython -c bioconda -c conda-forge -c anaconda
source activate honeypi_env
git clone https://github.com/hsgweon/honeypi.git
pip install ./honeypi
conda deactivate
```

### R DADA2 package

Within R, check to see if **dada2** is installed:

```
R
library(dada2)
```

If not, then install them with:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("dada2", ask = FALSE)
```
If it asks "Would you like to use a personal library instead?", say yes. It will ask for a location, and again just say yes.
When all done, type "quit()" and say "n".

Check again to see if the packages are all in place.

*That's it!*

### Test HONEYPI before you run it on your data

It is highly recommended that you test HONEYPI (see 4. Testing HONEYPI) before you run it on your dataset.

## 2. Running HONEYPI

Since we just created a sandbox "honeypi_env" in which all tools (except R and its packages) are installed, we need to be running honeypi and processing data within the environment:

```
source activate honeypi_env
```

Then, go to a directory where with your rawdata directory is located, and create a readpairslist file. This file is needed to ensure all files and sample names are correctly labelled. It does some internal checks to make sure there are no human errors with samples names etc. 

N.B. **rawdata** is your directory with paired-end sequences*

```
cd honeypi/testdata
honeypi_createreadpairslist -i rawdata -o readpairslist.txt
honeypi -i rawdata -o honeypi_output --amplicontype ITS2 -l readpairslist.txt
```

Done... simple... isn't it? Ah, one more thing - after you finished with HONEYPI, don't forget to get out of the sandbox by:

```
conda deactivate
```


## 3. To uninstall HONEYPI completely:

```
source activate honeypi_env
pip uninstall honeypi
conda deactivate
conda env remove --name honeypi -y
```


## 4. Testing HONEYPI

```
source activate honeypi_env

cd ~/honeypi/testdata
honeypi_createreadpairslist -i rawdata -o readpairslist.txt
honeypi -i rawdata -o honeypi_output --amplicontype ITS2 -l readpairslist.txt
```

**Look inside the output directory ("honeypi_output"), and find:**

1. ASVs_counts.txt
2. ASVs.fasta
3. ASVs_taxonomy.txt
4. summary.log


**(Misc) Check to see if your files have:**


*ASVs_counts.txt*

```
H12BB43 H12DSPINHK9
ASV_0000000001  0       2853
ASV_0000000002  0       1836
ASV_0000000003  0       889
ASV_0000000004  0       345
ASV_0000000005  0       340
ASV_0000000006  50      200
ASV_0000000007  0       218
ASV_0000000008  0       172
ASV_0000000009  100     0
ASV_0000000010  0       69
ASV_0000000011  0       62
ASV_0000000012  0       34
ASV_0000000013  0       31
ASV_0000000014  0       14
ASV_0000000015  0       10
ASV_0000000016  0       9
ASV_0000000017  0       5
ASV_0000000018  0       3
ASV_0000000019  0       3
ASV_0000000020  0       3
ASV_0000000021  0       2
```

*ASVs.fasta

```
>ASV_0000000001
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCTTCTGGCCGAGGGCACGTCTGCCTGGGTGTCACAAATCGTCGTCCCCCCCATCCTCTCGAGGATATGGGACGGAAGCTGGTCTCCCGTGTGTTACCGCACGCGGTTGGCCAAAATCCGAGCTAAGGGCGCCAGGAGCGTCTCGACATGCGGTGGTGAATTCAAGCCTCGTAATATCGTCGGTCGTTCCGGTCCAAAAGCTCTCGATGACCCAAAGTCCTCAACGCGACCCCAGGTCAGGCGGGATCACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000002
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGATGCCATTAGGCTGAAGGCACGTCTGCCTGGGCGTCACGCATGGCGTCGCCCACTCACCCCGTGCCTCTGTGGGCGGAAGGTGTGTGAGCGGATATTGGCCCCCCGTTCACGTTCGTGCTCGGTCGGTCTAAAAGGAAAGTCCCCAACGACGGACATCACGGCGAGTGGTGGTTGCCAGACCGTCCCGACGCGTCGTGCATGCTGTTCTTTGTCGTTGGCCGGCTCATCGACCCCCGAGTACCGTCAGGTACTCGGTACCTCGACTGCGACCCCAGGTCAGGCGGGATCACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000003
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCTGAAGCCATCAGGTTGAAGGCACGTCTGCCTGGGCGTCACGCATTGCGTCGCTCACTCACCCCGTGCATCATTGGGCGGGCAAGTGTGTGGGCGGATATTGGCCCCCCGTTCACATTTGTGCTCGGTCGGCCTAAAAAGAAGGTCCTTGATGACGGACATCACAACAAGTGGTGGTTGCTAAACCGTCGCGCCATGTTGTGCATTATACTCCGTCGTCGGTTGCCTCATTGACCCTTAAGTGCCATTGAACTTGGTACCTCAACTGCGACCCCAGGTCAGGCGGGATTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000004
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAAGCCATTAGGCCGAGGGCACGCCTGCCTGGGCGTCACACGTCGTTGCCCCCCCCCAAACCCCTCGGGAGTTGGGCGGGACGGATGATGGCCTCCCGTGTGCTCTGTCATGCGGTTGGCATAAAAACAAGTCCTCGGCGACTAACGCCACGACAATCGGTGGTTGTCAAACCTCTGTTGCCTATCGTGTGCGCGTGTCGAGCGAGGGCTCAACAAACCATGTTGCATCGATTCGTCGATGCTTTCAACGCGACCCCAGGTCAGGCGGGGTTACCCGCTGAATTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000005
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCATCCGGTTGAGGGCACGCCTGCCTGGGCGTCACGCTTCGCATCGCCCCCCACCATACATACCCAACGGGTACTAATGGTGTTTGGGGCGGAGATTGGCCTCCCGCACCTCTGATGCGGTTGGCCTAAAAATGAGTCCCCTTCAGCGGACACACGACTAGTGGTGGTTGAACAGACCCTCGTCCTTATCGTGTGTCGTGAGCTGCAAGGGAAACCCTCACCAAAGACCCTATTGCATTGTTTTTTGGACAATGCTTCGACCGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000006
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGATGCCATTAGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGAAGCCTCTCGCCAATTTCCTATATTGATAGGGGTATTGTGCAGGGCGAATGTTGGCCTCCCGTGAGCTTTATTGCCTCATGGTTGGTTGAAAATCGAGACCTTGGTAGGGTGTGCCATGATAGGTGGTGGCTGTGTTACGCACGAGACCAAGTAAGTCATGTGCTGCTCTATTGAATTTAGGCCTCTTTTACCCACATGCGTTTCGAAACGCTCGTGATGAGACCTCAGGTCAGGCGGGGCTACCCGCTGAATTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000007
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCTGATGCCATTAGGTTGAGGGCACGTCTGCCTGGGCGTCACATATCGAAGCCTCTTGCCAATTTCCTATTGATTGGTATTGTGCAAGATGATGTTGGCCTCCCGTGAGCACCATCGCCTCATGGTTGGTTGAAAATCGAGACCTTGGTAGAGTGTGCCATGATAAATGGTGCATGTGTTAAGCACGAGACCAAACAATCATGTGCTGCTCTATTGAATTTAGCCTCTTTTACCCACATGCGTGTCTAAACGCTCGTGATGAGACCTCAGGTCAGGCGGGGCTACCCGCTGAATTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000008
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGATGCCATTAGGTTGAGGGCACGTCTGCCTGGGTGTCACATATCGAAGCCTCCTTGCCAATTTCCCTGATTATTGTGCAGGGTGGATGTTGGCCTCCCGTGAGCTCTTTCGTCTCATGGTTGGTTGAAAATTGAGACCTTGGTAGGGTGTGCCATGATAGATGGTGGTTGTGTGACCCACGAGACCAATCATGCGCTGCTCTATTGAATTTGGCCTCCTTTACCCATATGCGTTTCCAAACGCTCGTGATGAGACCTCAGGTCAGGCGGGGCTACCCGCTGAATTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000009
ATGCGATACTTGGTGTGAATTGCAGAATTCCGTGAATCATCGAATCTTTGAACGCACATTGCGCCCCTTGGTATTCCAGGGGGCATGCCTGTTTGAGCGTCATTTCCCTCTCAAACGCTTGCGTTTGGTAGTGAGCGATACTCTTTTTGTGTGTATCTCTGAGGAGTTTGCTTGAAAGTGGGAGGCCATAGGCGGAGCCTAGCTTGAGCGTGTGGTGGAGGAACTGTGCCGAGAGGTGCAGGGCCGCGCTGCAACGCCTGGCCACGAAAACGAAGTCGTATTAGGTTTTACCGACTCGGCGAAGGAAGTAGTGGACGGGGGGAAAAGAGCGGAGCTCTCTTTTTTGTTTTGTTTGTTGATGATACGACGAGCAAGAGCAGCAGAGCCTGGCTTGAGAGAATTCACAAAGTTTGACCTCAAATCAGGTAGGATTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000010
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGACGCCTTCGGGCTGAGGGCACGCCTGCCTGGGCGTCACGCATCGCGTCGTCCCCTCCCATTCCCTCACGGGTTTGGTTATGGGACGGATAATGGCTTCCCGTTAGCTCGGTTAGCCCAAAAAGGATCCCTCATCGACGGATGTCACAACCAGTGGTGGTTGAAAGATCATTGGTGCTGTTGTGCTTCACCCTGTCGCTTGCTAGGGCATCGTCATAAACTAACGGCGTGTAATGCGCCTTCGACCGCGACCCCAGGTCAGACGGGACTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000011
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCTTCTGGCCGAGGGCACGTCTGCCTGGGTGTCACAAATCGTCGTCCCCCCATCCTCTCGAGGATATGGGACGGAAGCTGGTCTCCCGTGTGTTACCGCACGCGGTTGGCCAAAATCCGAGCTAAGGACGTTTTGGAGCGTCTCGACATGCGGTGGTGAATTGTAACCTCGTCATATTGTCGGTCGTTCCGGTTCAAAAGCTCTTGATGACCCAAAGTCCTCAACGCGACCCCAGGTCAGGCGGGATCACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000012
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCATCCGGCCGAGGGCACGCCTGCCTGGGCGTCACGCATCGCGTCGCCCCCACCAAATTTCCAAATCTGGTTGGGGGCGGAGATTGGCCTCCCGTACCTGTTGTGGTTGGCCTAAAAAGGAGTCCCCTTCGGTGGACACACGACTAGTGGTGGTTGAACAGACCCTCGTCTTTATTGTGTGTCATGAGCTGCTAGGGAGCCCTCATCAAAGACCCTTTGTATCGTTTTCGGACGGTGCTTCGACCGCGACCCCAGGTCAGACGGGACTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000013
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCATTCGGCCGAGGGCACGTCTGCCTGGGCGTCACGCATTGCGTCGCCCCAGACTACGCCTCCCCAACGGGGATGCGTTCGACTGGGGCGGAGAATGGTCTCCCGTGTCGTCGGCGTGGTTGGCCTAAAAAGGAGTCCCCTTCGGCGGACGCACGGCTAGTGGTGGTTGTTAAGGCCTTCGTATCGAGCCGTGTGTCGTTAGCCGCAAGGGAAGCACTCTTTAAAGACCCCAATGTGTCGTCTCGTGACGACGCTTCGACCGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000014
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCTAAGCCATTTGGCCGAGGGCACGCCTGCCTGGGCGTCAATCATCTATTCGTCACCCCAACCTCTGCTCCCCATAAAGGAGCTCGGGTCCTGGTTACGGAAGTTGGCCTCCCGTGGTCTCGAAGCGCGGCTGGCCTAAAATTGAGCATCGGGTTGGTGATCTCCGAGGCACGCGGTGGTTGTTCATTCTTACCTCGTGATGTTGCCCCGGGGCATCTTCCACAAGAAGCTCCACGACCCTAGATACATATCGATGCGACCCCAGGTCAGGCGGGGCCACCCGCTGAATTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000015
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCAAGTCTTTGAACGCAAGTTGCGCCCGAAGCCATTAGGCCGAGGGCACGTCTGCCTGGGCGTCACACGTCGTTGCCCCCCCCCAACCCCCTCGGGAGTTGGATGGGACGGATGATGGCCTCCCGTGTGCTCAGTCACGCGGTTGGCATAAATACCAAGTCCTCGGCGACCAACGCCACGACAATCGGTGGTTGTCAAACCTCGGTTTCCTGTCGTGCGCGCGTGTTGATCGAGTGCTTTCTTAAACAATGCGTGTCGATCCGTCGATGCTTACAACGCGACCCCAGGTCAGGCGGGGTTACCCGCTGAATTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000016
ATGCGATACTTGGTGTGAATTGCATAATCCCGAAGTGACGGCACGACCGAACAAAGCCCGAGCGGTAGCGGCGGAGACGTCGTGCCCTCGGAAACGCGACCCCAGGTCAGGCGGGGCCACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000017
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAAGCCATTAGGCTGAGGGCACGTCTGCCTGGGCGTCACGCATGACGTCGCCCCCCAACCTCGCTCTCACTCGTGGGAGTTGTTGCGGAGGGGCGGATACTGGCCTCCCGTGCCTCATCGTATGGTTGGCCCAAATGTGAGTCCTTGGCGACGGACGTCACGACAAGTGGTGGTTGTAAAAAGCCCTCTTCTCCTGTCGTGCGGTGGCGCGTCGCCAGCAAGAACTCTCGTGACCCTGTTGTGCCGTTGTCAACGCGCACTCCGACTGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000018
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCTTCTGGCCGAGGGCACGTCTGCCTGGGTGTCACAAAAGCTCTCGATGACCCAAAGTCCTCAAAGCGACCCCAGGTCAGGCGGGATCACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000019
AGGCGATACTTGGTGTGAATGACCGTCACCTTAGTTAGCTCAACGACCCTTACACGCCACAAACTTTGTGCGCTTCGATTGTGACCCCATGTCAGGCGGGATTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000020
ATGCGATACTTGGTGTGAATTGAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCTGAAGCCATCAGGTTGAAGGCACGTCTGCCTGGGTGTCACAAATCGTCGTCCCCCCCATCCTCTCGAGGATATGGGACGGAAGCTGGTCTCCCGTGTGTTACCGCACGCGGTTGGCCAAAATCCGAGCTAAGGATGCCAGGAGCGTCTTGACATGCGGTGGTGAATTCAATCTCCTCGTCATATCGTCGGTCGTTCCGGTCCAAAAGCTCTCGATGACCCAAAGTCCTCAACGTGACCCCAGGTCAGGCGGGATCACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
>ASV_0000000021
ATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAAGCCATTAGGTCGAGGGCACGCCTGCCTGGGCGTCACGCATCGCGTCTCCCCCCAACCACCCTGCGTGGATTGGGAGGAGGATGATGGCCTCCCATGCCTCACCGGGCGTGGATGGCCTAAATAAGGAGCCCCCGGTTACGAAGTGCCGCGGCGATTGGTGGAATACAAGGCCTAGCCTAGGACGAAATCGAAGTCGCGCACATCGTAGCTCTTGAGGACTCGCAGGACCCTAACTTGTTTGCCCCTAGGGGCGGCAAAACCGTTGCGACCCCAGGTCAGGCGGGGCTACCCGCTGAGTTTAAGCATATCAATAAGCGGAGGA
```

*ASVs_taxonomy.txt*

```
ASV_0000000001  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Brassicales; f__Brassicaceae; g__Brassica; s__Brassica_nigra    0.79
ASV_0000000002  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Ericales; f__Ericaceae; g__Calluna; s__Calluna_vulgaris 1.0
ASV_0000000003  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Ericales; f__Ericaceae; g__Erica; s__Erica_cinerea      1.0
ASV_0000000004  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Rosaceae; g__Rubus  1.0
ASV_0000000005  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Asterales; f__Asteraceae; g__Scorzoneroides; s__Scorzoneroides_autumnalis       1.0
ASV_0000000006  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Fabales; f__Fabaceae; g__Trifolium; s__Trifolium_repens 1.0
ASV_0000000007  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Fabales; f__Fabaceae; g__Trifolium; s__Trifolium_pratense       1.0
ASV_0000000008  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Fabales; f__Fabaceae; g__Vicia; s__Vicia_faba   1.0
ASV_0000000009  k__Fungi; p__Ascomycota; c__Saccharomycetes; o__Saccharomycetales; f__Saccharomycetaceae; g__Zygosaccharomyces; s__Zygosaccharomyces_mellis     1.0
ASV_0000000010  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Lamiales; f__Plantaginaceae; g__Plantago; s__Plantago_lanceolata        1.0
ASV_0000000011  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Brassicales; f__Brassicaceae; g__Sinapis; s__Sinapis_alba       1.0
ASV_0000000012  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Asterales; f__Asteraceae; g__Hypochaeris; s__Hypochaeris_radicata       0.98
ASV_0000000013  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Asterales; f__Asteraceae; g__Cirsium; s__Cirsium_arvense        1.0
ASV_0000000014  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Myrtales; f__Onagraceae; g__Chamaenerion; s__Chamaenerion_angustifolium 1.0
ASV_0000000015  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Rosaceae; g__Rosa   0.81
ASV_0000000016  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Asparagales; f__Asphodelaceae; g__Phormium; s__Phormium_cookianum       0.5
ASV_0000000017  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Apiales; f__Araliaceae; g__Hedera; s__Hedera_helix      0.98
ASV_0000000018  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Ranunculales; f__Menispermaceae; g__Tinospora; s__Tinospora_malabarica  0.5
ASV_0000000019  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Apiales; f__Apiaceae; g__Conopodium; s__Conopodium_majus        0.74
ASV_0000000020  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Brassicales; f__Brassicaceae; g__Raphanus; s__Raphanus_sativus  0.7
ASV_0000000021  k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Caryophyllales; f__Chenopodiaceae; g__Atriplex; s__Atriplex_patula      0.77
```

*summary.log*

```
Number of reads: 10441
Number of joined reads: 9379
Number of ASVs (after dada2): 21
Number of ASVs (after ITSx): 18
```


## 4. Combining two sets of honeypi results

You need two sets of results.

First honeypi results:

```
ASVs_1.fasta
ASVs_counts_1.txt
ASVs_taxonomy_1.txt
```

Second honeypi results to add to the first:

```
ASVs_2.fasta
ASVs_counts_2.txt
ASVs_taxonomy_2.txt
```

Then execute the following commands - watch out for the order of files. 

```
source activate honeypi_env
honeypi_joinTwoResults.py -i1 ASVs_1.fasta,ASVs_counts_1.txt,ASVs_taxonomy_1.txt -i2 ASVs_2.fasta,ASVs_counts_2.txt,ASVs_taxonomy_2.txt
```
