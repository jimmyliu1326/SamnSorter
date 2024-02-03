#!/usr/bin/env Rscript

# load pkgs
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(digest))

# get version
version <- tryCatch(
    version <- readLines("Version"),
    error = function(e) {
        message("Cannot read tool version")
        q(status = 1)
    }
)

# cli argument parser
parser <- ArgumentParser()
parser$add_argument("file", nargs='+')
parser$add_argument("-o", "--outdir", type="character", default="./samnsorter",
    help="Directory path to where output files will be written to [%(default)s]")
parser$add_argument("--tmpdir", type="character", default="",
    help="Directory path to where temporary files will be stored [%(default)s]")
parser$add_argument("--threshold", type="integer", default=1538,
    help="Clustering distance threshold for cluster assignments [%(default)s]")
parser$add_argument("-t", metavar="THREADS", dest="threads", type="double", default=4,
    help="Number of threads to use [%(default)s]")
parser$add_argument("-v", "--version", action="version", version=paste0("SamnSorter v", version),
    help="Print version information")
args <- parser$parse_args()

# input validation
check_args <- function(args=NULL) {
    cat("Validating inputs...\n")
    # check if input files exist
    purrr::walk(args$file, function(f) {
        if (!file.exists(f)) {
            message(f, " does not exist, please verify input file path(s)")
            q(status = 1)
        }
    })
    # set up tmpdir
    hash <- digest(Sys.time(), algo='crc32')
    if (nchar(args$tmpdir) == 0) {
        tmpdir <<- file.path(args$outdir, hash) 
    } else {
        if (!file.exists(args$tmpdir)) { 
            message(args$tmpdir, "does not exist, please verify the directory path")
            q(status = 1)
        }
        tmpdir <<- file.path(args$tmpdir, hash)
    }
    cat("Created temporary file directory at:", tmpdir, "\n")
    system(paste('mkdir -p', tmpdir))
    # archive logs in outdir
    logdir <<- file.path(args$outdir, "logs")
    cat("Analysis logs are saved to:", logdir, "\n")
    system(paste('mkdir -p', logdir))

}

# set up input directory
input_setup <- function (
    file = NULL, # a vector of paths to genomes
    tmpdir = NULL # tmp directory path where to run analysis
) { 
    cat("Setting up analysis input directory...\n")
    system(paste('mkdir -p', file.path(tmpdir, "input")))
    walk(file, function(f) {
        realpath <- file.path(normalizePath(f))
        cmd <- paste("ln -s", realpath, file.path(tmpdir, "input", basename(realpath)))
        cat(cmd, "\n")
        system(cmd)
    })
}

# cgMLST with ChewBBACA
cgmlst <- function(
    tmpdir = NULL,
    logdir = NULL,
    outdir = NULL,
    threads = 4
) {
    cat("Running cgMLST...\n")
    cmd <- paste("chewBBACA.py AlleleCall",
    "-i", file.path(tmpdir, "input"),
    "-g $REF_PATH",
    "--ptf $TRAINING_PATH",
    "-o", file.path(tmpdir, "cgmlst"),
    "--hash-profiles sha1",
    "--cpu", threads,
    "--no-inferred",
    ">", file.path(logdir, "chewBBACA.log"), "2>&1") # save log to file
    cat(cmd, "\n")
    system(cmd)
    
    # publish hashed and unhashed profiles in outdir
    cat("Publishing cgMLST results...\n")
    cmd <- paste("cp -r", file.path(tmpdir, "cgmlst"),
                 outdir)
    cat(cmd, "\n")
    system(cmd)
}

# calculate query distance in respect to references
cgmlst_dist_query <- function(
    tmpdir = NULL,
    logdir = NULL,
    outdir = NULL
) {
    cat("Calculating allelic distance between query and reference...\n")
    cmd <- paste("$CGMLST_DISTS -H", 
                 "$REF_ALLELES", # reference profiles
                 file.path(tmpdir, "cgmlst", "results_alleles_hashed.tsv"), # query profiles
                 ">", file.path(outdir, "cgmlst_dist.tsv")) # output distance matrix
    cat(cmd, "\n")
    system(cmd)
}

# main workflow
tic()
cat("This is SamnSorter", paste0("v", version, "\n"))
check_args(args)
input_setup(args$file, tmpdir)
tmpdir <- file.path(args$outdir, "6cbdab8f")
#cgmlst(tmpdir, logdir, args$outdir, args$threads)
cgmlst_dist_query(tmpdir, logdir, args$outdir)
cat("Workflow completed successfully.\n")
toc(log = T)

