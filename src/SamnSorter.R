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

# helper func for running shell cmds
shell <- function(command, echo = F) {
    if (echo) {
        system(command)
        cat(command, "\n")
    } else {
        system(command)
    }
}

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
            message(args$tmpdir, " does not exist, please verify the directory path")
            q(status = 1)
        }
        tmpdir <<- file.path(args$tmpdir, hash)
    }
    cat("Created temporary file directory at:", tmpdir, "\n")
    shell(paste('mkdir -p', tmpdir))
    # archive logs in outdir
    logdir <<- file.path(args$outdir, "logs")
    cat("Analysis logs are saved to:", logdir, "\n")
    shell(paste('mkdir -p', logdir))

}

# set up input directory
input_setup <- function (
    file = NULL, # a vector of paths to genomes
    tmpdir = NULL # tmp directory path where to run analysis
) { 
    cat("Setting up analysis input directory...\n")
    shell(paste('mkdir -p', file.path(tmpdir, "input")))
    walk(file, function(f) {
        realpath <- file.path(normalizePath(f))
        cmd <- paste("ln -s", realpath, file.path(tmpdir, "input", basename(realpath)))
        shell(cmd, echo = T)
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
    shell(cmd, echo = T)
    
    # publish hashed and unhashed profiles in outdir
    cat("Publishing cgMLST results...\n")
    cmd <- paste("cp -r", file.path(tmpdir, "cgmlst"),
                 outdir)
    shell(cmd, echo = T)
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
                 ">", file.path(tmpdir, "cgmlst_dist.tsv")) # output distance matrix
    shell(cmd, echo = T)
    # reformat distance matrix
    dist <- fread(file.path(tmpdir, "cgmlst_dist.tsv"))
    dist_out <<- dist %>% column_to_rownames("cgmlst-dists")
    rownames(dist_out) <- str_replace_all(rownames(dist_out), "\\.1|\\.2", "")
    colnames(dist_out) <- str_replace_all(colnames(dist_out), "\\.1|\\.2", "")
    write.table(dist_out, file.path(outdir, "cgmlst_dist.tsv"),
                sep = "\t", quote = F, row.names = T, col.names = NA)

}

# best hit search
bh_search <- function(
    tmpdir = NULL,
    outdir = NULL) {
        cat("Identifying best hit...\n")
        bh <- apply(dist_out, 1, function(x) { names(which(x == min(x)))[1] })
        bh_res <- data.frame("id" = names(bh),
                             "best_hit" = bh)
        write.table(bh_res, file.path(tmpdir, "best_hit.tsv"),
                    sep = "\t", row.names = F, quote = F)

}

# main search 
search <- function(
    tmpdir = NULL,
    outdir = NULL) {

        # run best hit
        bh_search(tmpdir, outdir)
}

# main workflow
tic()
cat("This is SamnSorter", paste0("v", version, "\n"))
check_args(args)
tmpdir <- file.path(args$outdir, "6cbdab8f")
#tmpdir <- file.path(args$outdir, "demo", "b5329061")
#input_setup(args$file, tmpdir)
#cgmlst(tmpdir, logdir, args$outdir, args$threads)
cgmlst_dist_query(tmpdir, logdir, args$outdir)
search(tmpdir, args$outdir)
cat("Workflow completed successfully.\n")
toc(log = T)

