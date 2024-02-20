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
parser$add_argument("--min_lr", type="double", default=0.25,
    help="Minimum likelihood weighted ratio for phylogenetic placements [%(default)s]")
parser$add_argument("-t", metavar="THREADS", dest="threads", type="integer", default=4,
    help="Number of threads to use [%(default)s]")
parser$add_argument("-v", "--version", action="version", version=paste0("SamnSorter v", version),
    help="Print version information")
args <- parser$parse_args()

# helper func for running sh cmds
sh <- function(command, echo = F) {
    if (echo) {
        cat(command, "\n")
        ec <- system(command)
    } else {
        ec <- system(command)
    }
    if (ec != 0) {
        cat("\nNonzero exit code detected: check the log for error messages.\n")
        q(status = ec)
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
    sh(paste('mkdir -p', tmpdir))
    # archive logs in outdir
    logdir <<- file.path(args$outdir, "logs")
    cat("Analysis logs are saved to:", logdir, "\n")
    sh(paste('mkdir -p', logdir))

}

# set up input directory
input_setup <- function (
    file = NULL, # a vector of paths to genomes
    tmpdir = NULL # tmp directory path where to run analysis
) { 
    cat("Setting up analysis input directory...\n")
    sh(paste('mkdir -p', file.path(tmpdir, "input")))
    walk(file, function(f) {
        realpath <- file.path(normalizePath(f))
        cmd <- paste("ln -s", realpath, file.path(tmpdir, "input", basename(realpath)))
        sh(cmd, echo = T)
    })
}

# cgMLST with ChewBBACA
cgmlst <- function(
    tmpdir = NULL, # tmp directory path where to run analysis
    logdir = NULL, # log directory path for storing stdout and stderr
    outdir = NULL, # output dir
    threads = 4    # number of threads
) {
    cat("Running cgMLST...\n")
    cmd <- paste(
    "chewBBACA.py",
    "AlleleCall",
    "-i", file.path(tmpdir, "input"),
    "-g $REF_PATH",
    "--ptf $TRAINING_PATH",
    "-o", file.path(tmpdir, "cgmlst"),
    "--hash-profiles sha1",
    "--cpu", threads,
    "--no-inferred",
    ">", file.path(logdir, "chewBBACA.log"), "2>&1") # save log to file
    sh(cmd, echo = T)
    
    # publish hashed and unhashed profiles in outdir
    cat("Publishing cgMLST results...\n")
    cmd <- paste("cp -r", file.path(tmpdir, "cgmlst"),
                 outdir)
    sh(cmd, echo = T)
}

# calculate query distance in respect to references
cgmlst_dist_query <- function(
    tmpdir = NULL, # tmp directory path where to run analysis
    logdir = NULL, # log directory path for storing stdout and stderr
    outdir = NULL # output dir
) {
    cat("Calculating allelic distance between query and reference...\n")
    cmd <- paste("$CGMLST_DISTS", "-H", 
                 "$REF_ALLELES", # reference profiles
                 file.path(tmpdir, "cgmlst", "results_alleles_hashed.tsv"), # query profiles
                 ">", file.path(tmpdir, "cgmlst_dist.tsv")) # output distance matrix
    sh(cmd, echo = T)
    # reformat distance matrix
    dist <- fread(file.path(tmpdir, "cgmlst_dist.tsv"), sep = "\t")
    dist_out <<- dist %>% column_to_rownames("cgmlst-dists")
    #rownames(dist_out) <- str_replace_all(rownames(dist_out), "\\.1|\\.2", "")
    #colnames(dist_out) <- str_replace_all(colnames(dist_out), "\\.1|\\.2", "")
    write.table(dist_out, file.path(outdir, "cgmlst_dist.tsv"),
                sep = "\t", quote = F, row.names = T, col.names = NA)

}

# best hit search
bh_search <- function(
    tmpdir = NULL # tmp directory path where to run analysis
) {
    cat("Identifying best hit...\n")
    bh <- apply(dist_out, 1, function(x) { names(which(x == min(x)))[1] })
    # get cluster ID of the bh
    clust_path <- Sys.getenv("REF_CLUSTERS")  
    clust <- fread(clust_path, sep = "\t")
    bh_clust <- clust$clust[which(clust$id == bh)]
    bh_res <- data.frame("id" = names(bh),
                          "best_hit" = bh_clust)
    write.table(bh_res, file.path(tmpdir, "best_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
    return(bh_res)
}

# phylogenetic placement using APPLES
pp_search <- function(
    tmpdir = NULL, # tmp directory path where to run analysis
    logdir = NULL, # log directory path for storing stdout and stderr
    outdir = NULL, # output dir
    threads = 4,   # number of threads
    min_lr = 0.25  # minimum placement likelihood ratio

) {
    cat("Running phylogenetic placement...\n")
    cmd <- paste(
        "run_apples.py", 
        "-t", "$REF_NWK", # backbone tree
        "-d", file.path(outdir, "cgmlst_dist.tsv"), # query distance matrix
        "-o", file.path(tmpdir, "query.jplace"), # jplace output
        "-T", threads,
        ">", file.path(logdir, "APPLES.log"), "2>&1" 
    )
    sh(cmd, echo = T)
    cmd <- paste(
        "gappa", "examine", "assign",
        "--jplace-path", file.path(tmpdir, "query.jplace"), # input jplace
        "--taxon-file", "$REF_TAXONOMY", # input reference taxonomy
        "--allow-file-overwriting",
        "--per-query-results",
        "--out-dir", file.path(tmpdir), # output dir
        "--threads", threads,
        ">", file.path(logdir, "gappa.log"), "2>&1" # write log
    )
    #sh(cmd, echo = T)
    # analyze pp results
    invalid_tax <- c("DISTANT", "enterica", "bongori", "bongori;outgroup")
    asgmnts <- fread(file.path(tmpdir, "per_query.tsv"), sep = "\t")
    asgmnts.val <- filter(asgmnts, LWR >= min_lr, !(taxopath %in% invalid_tax))
    pp_res <- asgmnts.val %>%
        rename("id" = "name",
               "pp_clust" = "taxopath") %>%
        mutate(pp_clust = str_replace_all(pp_clust, ".*;", "")) %>%
        select(id, pp_clust)
    write.table(pp_res, file.path(tmpdir, "pp_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
    return(pp_res)
}

# main search 
search <- function(
    tmpdir = NULL, # tmp directory path where to run analysis
    logdir = NULL, # log directory path for storing stdout and stderr
    args = NULL    # arg parameters
) {
        # run best hit
        bh_res <- bh_search(tmpdir)
        # run pp
        pp_res <- pp_search(tmpdir, logdir, args$outdir, args$threads, args$min_lr)
        # merge results
        search_res <- full_join(pp_res, bh_res, by = "id") %>%
            mutate(pp_clust = if_else(is.na(pp_clust), "NOVEL", pp_clust))
        # final prediction by majority voting
        final_preds <- apply(search_res[,2:ncol(search_res)], 1,
                             function(x) {
                                if (length(x) == length(unique(x))) {
                                    return("NOVEL")
                                } else {
                                    return(as.character(x[which.max(table(x))]))
                                }
                             })
        search_res <- cbind(search_res, final_clust = final_preds)
        cat("Writing cluster assignment results...\n")
        write.table(search_res, file.path(tmpdir, "samnsorter_res.tsv"),
                sep = "\t", row.names = F, quote = F)
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
search(tmpdir, logdir, args)
cat("Workflow completed successfully.\n")
toc(log = T)

