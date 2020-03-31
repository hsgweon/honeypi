#!/usr/bin/env python

import argparse, sys, os, argparse, shutil, subprocess, gzip, bz2

from honeypi.logger import *
from honeypi.runcmd import *

__version__    = 0.1

__author__     = "Hyun Soon Gweon"
__copyright__  = "Copyright 2015"
__credits__    = ["Hyun Soon Gweon"]
__license__    = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__      = "h.s.gweon@reading.ac.uk"

PEAR                       = "pear"
VSEARCH                    = "vsearch"
FASTQJOIN                  = "fastq-join"
FASTX_FASTQ_QUALITY_FILTER = "fastq_quality_filter"
FASTX_FASTQ_TO_FASTA       = "fastq_to_fasta"


def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def getFileLineCount(filename, extensionType = "uncompressed"):

    import gzip, bz2
    
    if extensionType == "gz":
        f = gzip.open(filename, "r")
        return sum(bl.count(b"\n") for bl in blocks(f))
    elif extensionType == "bz2":
        f = bz2.open(filename, "r")
        return sum(bl.count(b"\n") for bl in blocks(f))
    else:
        f = open(filename, "r")
        return sum(bl.count("\n") for bl in blocks(f))


def get_md5(filename):
    
    import hashlib

    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def downloadDB(
    url,
    md5,
    output_dir,
    logging_file, 
    summary_file, 
    verbose):

    import progressbar
    import requests
    import tarfile

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # else:
    #     shutil.rmtree(output_dir)
    #     os.mkdir(output_dir)

    filename     = url.split("/")[-1]

    # Check if unpacked file exists
    DOWNLOADDB = True
    # print(output_dir + "/" + filename)
    if os.path.isfile(output_dir + "/" + filename):
        if get_md5(output_dir + "/" + filename) != md5:
            logger("... File exits, but seems to be different, so re-downloading...", logging_file, display = True)
            DOWNLOADDB = True
        else:
            logger("... File exits. Skipping downloading and preparing to unpack the existing file...", logging_file, display = True)
            DOWNLOADDB = False

    if DOWNLOADDB:

        request = requests.get(url, stream=True)

        file         = open(output_dir + "/" + filename, 'wb')
        file_size    = int(request.headers['Content-Length'])
        file_size_mb = round(file_size / 1024 / 1024,2)

        chunk_sz = 512

        widgets = [progressbar.Bar(marker="#",left="[",right="]"),
                   progressbar.Percentage()," | ",
                   progressbar.FileTransferSpeed()," | ",
                   progressbar.SimpleProgress()," | ",
                   progressbar.ETA()]

        bar = progressbar.ProgressBar(widgets=widgets, maxval=file_size).start()

        i = 0

        for chunk in request.iter_content(chunk_size=chunk_sz):
            file.write(chunk)
            i += len(chunk)
            bar.update(i)

        bar.finish()
        print('File size: {0} MB'.format(file_size_mb))

        file.close()

        if get_md5(output_dir + "/" + filename) != md5:
            print("Downloaded data is corrupt. Get in touch with Soon!. Exiting...")
            exit(1)

    # Unpack
    logger("... Unpacking", logging_file, display = True)
    tar = tarfile.open(output_dir + "/" + filename)
    tar.extractall(path=output_dir)
    tar.close()


def count_sequences(
    input_dir, 
    filenames_list, 
    logging_file, 
    summary_file, 
    verbose):

    filenameextensions = []
    for filename in filenames_list:
        filenameextensions.append(filename.split(".")[-1].rstrip())
    if len(set(filenameextensions)) > 1:
        logger("Error: More than two types of extensions", logging_file)
        exit(1)
    extensionType = next(iter(filenameextensions))

    # Count
    numberofsequences = 0
    for filename in filenames_list:
        numberofsequences += int(getFileLineCount(input_dir  + "/" + filename, extensionType) / 4)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences in the rawdata!", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of reads: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of reads: " + str(numberofsequences) + "\n")



# Join paired-end reads
def run_trimgalore(
    input_dir,
    output_dir,
    fastqs_f,
    fastqs_r,
    sampleids_list,
    logging_file,
    summary_file,
    verbose):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    for i in range(len(sampleids_list)):

        cmd = " ".join([
            "trim_galore --phred33 --paired --illumina -q 25 --retain_unpaired --length 100 --max_n 1",
            input_dir + "/" + fastqs_f[i],
            input_dir + "/" + fastqs_r[i],
            "-o", output_dir
            ])

        # If empty, then create empty outputs
        if os.stat(input_dir + "/" + fastqs_f[i]).st_size == 0 or os.stat(input_dir + "/" + fastqs_r[i]).st_size == 0:
            open(output_dir + "/" + sampleids_list[i] + ".fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".discarded.fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".unassembled.forward.fastq", 'a').close()
            open(output_dir + "/" + sampleids_list[i] + ".unassembled.reverse.fastq", 'a').close()
            continue
        
        run_cmd(cmd, logging_file, verbose)

    # Count
    numberofsequences = 0
    for i in range(len(sampleids_list)):
        filename = output_dir  + "/" + fastqs_f[i].split(".fastq.gz")[0]+ "_val_1.fq.gz"
        numberofsequences += int(getFileLineCount(filename, "gz") / 4)
    
    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences after joining. Something is not right.", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of joined reads: " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of joined reads: " + str(numberofsequences) + "\n")


def run_dada2(
    input_dir,
    output_dir,
    logging_file,
    summary_file,
    verbose):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    cmd = " ".join(["honeypi_dada2", input_dir, output_dir])

    run_cmd(cmd, logging_file, verbose)


    # Count
    numberofsequences = 0
    filename = output_dir  + "/ASVs.fasta"
    numberofsequences += int(getFileLineCount(filename) / 2)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences after running dada2...", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of ASVs (after dada2): " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of ASVs (after dada2): " + str(numberofsequences) + "\n")


def run_ITSx(
    input_dir,
    output_dir,
    its_region,
    logging_file,
    summary_file,
    verbose):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    cmd = " ".join([
            "ITSx",
            "-i", input_dir + "/ASVs.fasta", 
            "-o", output_dir + "/ASVs_ITSx",
            "--preserve T"])
    run_cmd(cmd, logging_file, verbose)


    if its_region == "ITS2":
        filename = output_dir  + "/ASVs_ITSx.ITS2.fasta"


    # Remove short sequences (below 100bp)
    cmd = " ".join([
            "vsearch",
            "--fastx_filter", filename, 
            "--fastq_minlen 100",
            "--fasta_width 0",
            "--fastaout", output_dir + "/ASVs_ITSxed.fasta"])
    run_cmd(cmd, logging_file, verbose)


    # Count
    numberofsequences = 0
    numberofsequences += int(getFileLineCount(output_dir + "/ASVs_ITSxed.fasta") / 2)

    if numberofsequences == 0: 
        logger("ERROR: You have 0 sequences after running ITSx...", logging_file, display = True)
        exit(1)
    else:
        logger(BLUE + "... number of ASVs (after ITSx): " + str(numberofsequences) + ENDC, logging_file, display = True)
        summary_file.write("Number of ASVs (after ITSx): " + str(numberofsequences) + "\n")


def consolidateSequences(
    output_dir,
    logging_file,
    summary_file,
    verbose):

    # filterbyname.sh in=honeypi_output/dada2/ASVs.fasta out=foo.fasta names=honeypi_output/ITSx/ITSxed.list include=t

    cmd = "grep '>' " + output_dir + "/ITSx/ASVs_ITSxed.fasta | cut -c 2- > " + output_dir + "/ITSx/ITSxed.list"
    run_cmd(cmd, logging_file, verbose)

    cmd = "filterbyname.sh in=" + output_dir + "/dada2/ASVs.fasta out=" + output_dir + "/ITSx/ASVs_ITSxed_NOT.fasta names=" + output_dir + "/ITSx/ITSxed.list include=f overwrite=t -Xmx500g"
    run_cmd(cmd, logging_file, verbose)

    cmd = "cat " + output_dir + "/ITSx/ASVs_ITSxed.fasta " + output_dir + "/ITSx/ASVs_ITSxed_NOT.fasta | seqkit sort --quiet > " + output_dir + "/ITSx/ASVs.fasta"
    run_cmd(cmd, logging_file, verbose)


def run_RDPClassifier(
    input_fasta,
    output_dir,
    rdpclassifier_properties,
    logging_file,
    summary_file,
    verbose):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    cmd = " ".join(["classifier",
                    "-Xmx16G",
                    "classify", 
                    "-t", rdpclassifier_properties,
                    "-o", output_dir + "/assigned_taxonomy_prelim.txt", 
                    input_fasta])
    run_cmd(cmd, logging_file, verbose)

    cmd = " ".join(["honeypi_reformatAssignedTaxonomy", 
                    "-i", output_dir + "/assigned_taxonomy_prelim.txt", 
                    "-o", output_dir + "/assigned_taxonomy.txt", 
                    "-c", "0.5"])
    run_cmd(cmd, logging_file, verbose)


def filterASVtable(
    input_table,
    input_fasta,
    output_dir,
    logging_file,
    summary_file,
    verbose):

    cmd = " ".join(["honeypi_filterASVtable",
                    "--otu", input_table,
                    "--fas", input_fasta,
                    "--out", output_dir + "/ASVs_counts.txt"])
    run_cmd(cmd, logging_file, verbose)






