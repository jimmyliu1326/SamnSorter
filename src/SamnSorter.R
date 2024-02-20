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

# workflow metadata class
setClass("Workflow metadata", 
         slots=list(
            tmpdir="character", 
            logdir="character", 
            outdir="character",
            args="list"
        )
)

# S4 generics
invisible(setGeneric("input_setup", function(meta, args) standardGeneric("input_setup")))
invisible(setGeneric("cgmlst", function(meta, args) standardGeneric("cgmlst")))
invisible(setGeneric("cgmlst_dist_query", function(meta, args) standardGeneric("cgmlst_dist_query")))
invisible(setGeneric("main_search", function(meta, args) standardGeneric("main_search")))
invisible(setGeneric("pp_search", function(meta, args) standardGeneric("pp_search")))
invisible(setGeneric("bh_search", function(meta, args) standardGeneric("bh_search")))

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
        wf_meta@tmpdir <- file.path(args$outdir, hash) 
    } else {
        if (!file.exists(args$tmpdir)) { 
            message(args$tmpdir, " does not exist, please verify the directory path")
            q(status = 1)
        }
        wf_meta@tmpdir <- file.path(args$tmpdir, hash)
    }
    cat("Created temporary file directory at:", wf_meta@tmpdir, "\n")
    sh(paste('mkdir -p', wf_meta@tmpdir))
    # archive logs in outdir
    wf_meta@logdir <- file.path(args$outdir, "logs")
    cat("Analysis logs are saved to:", wf_meta@logdir, "\n")
    sh(paste('mkdir -p', wf_meta@logdir))
    # set outdir in wf meta
    wf_meta@outdir <- args$outdir
    return(wf_meta)
}

# set up input directory
setMethod("input_setup",
signature="Workflow metadata",
function(
    meta,
    args
) { 
    cat("Setting up analysis input directory...\n")
    sh(paste('mkdir -p', file.path(meta@tmpdir, "input")))
    walk(args$file, function(f) {
        realpath <- file.path(normalizePath(f))
        cmd <- paste("ln -s", realpath, file.path(tmpdir, "input", basename(realpath)))
        sh(cmd, echo = T)
    })
})

# cgMLST with ChewBBACA
setMethod("cgmlst",
signature="Workflow metadata",
function(
    meta,
    args
) {
    cat("Running cgMLST...\n")
    cmd <- paste(
    "chewBBACA.py",
    "AlleleCall",
    "-i", file.path(meta@tmpdir, "input"),
    "-g $REF_PATH",
    "--ptf $TRAINING_PATH",
    "-o", file.path(meta@tmpdir, "cgmlst"),
    "--hash-profiles sha1",
    "--cpu", args$threads,
    "--no-inferred",
    ">", file.path(meta@logdir, "chewBBACA.log"), "2>&1") # save log to file
    sh(cmd, echo = T)
    
    # publish hashed and unhashed profiles in outdir
    cat("Publishing cgMLST results...\n")
    cmd <- paste("cp -r", 
                 file.path(meta@tmpdir, "cgmlst"),
                 meta@outdir)
    sh(cmd, echo = T)
})

# calculate query distance in respect to references
setMethod("cgmlst_dist_query",
signature="Workflow metadata",
function(
    meta,
    args
) {
    cat("Calculating allelic distance between query and reference...\n")
    cmd <- paste("$CGMLST_DISTS", "-H", 
                 "$REF_ALLELES", # reference profiles
                 file.path(meta@tmpdir, "cgmlst", "results_alleles_hashed.tsv"), # query profiles
                 ">", file.path(meta@tmpdir, "cgmlst_dist.tsv")) # output distance matrix
    sh(cmd, echo = T)
    # reformat distance matrix
    dist <- fread(file.path(meta@tmpdir, "cgmlst_dist.tsv"), sep = "\t")
    dist_out <<- dist %>% column_to_rownames("cgmlst-dists")
    #rownames(dist_out) <- str_replace_all(rownames(dist_out), "\\.1|\\.2", "")
    #colnames(dist_out) <- str_replace_all(colnames(dist_out), "\\.1|\\.2", "")
    cat("Saving distance matrix to:", file.path(meta@outdir, "cgmlst_dist.tsv"), "\n")
    write.table(dist_out, file.path(meta@outdir, "cgmlst_dist.tsv"),
                sep = "\t", quote = F, row.names = T, col.names = NA)
})

# best hit search
setMethod("bh_search",
signature="Workflow metadata",
function(
    meta,
    args
) {
    cat("Identifying best hit...\n")
    bh <- apply(dist_out, 1, function(x) { names(which.min(x))[1] })
    # get cluster ID of the bh
    clust_path <- Sys.getenv("REF_CLUSTERS")  
    clust <- fread(clust_path, sep = "\t")
    bh_clust <- clust$clust[which(clust$id == bh)]
    bh_res <- data.frame("id" = names(bh),
                          "best_hit" = bh_clust)
    write.table(bh_res, file.path(meta@tmpdir, "best_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
    return(bh_res)
})

# phylogenetic placement using APPLES
setMethod("pp_search",
signature="Workflow metadata",
function(
    meta,
    args
) {
    cat("Running phylogenetic placement...\n")
    cmd <- paste(
        "run_apples.py", 
        "-t", "$REF_NWK", # backbone tree
        "-d", file.path(meta@outdir, "cgmlst_dist.tsv"), # query distance matrix
        "-o", file.path(meta@tmpdir, "query.jplace"), # jplace output
        "-T", args$threads,
        ">", file.path(meta@logdir, "APPLES.log"), "2>&1"
    )
    sh(cmd, echo = T)
    cmd <- paste(
        "gappa", "examine", "assign",
        "--jplace-path", file.path(meta@tmpdir, "query.jplace"), # input jplace
        "--taxon-file", "$REF_TAXONOMY", # input reference taxonomy
        "--allow-file-overwriting",
        "--per-query-results",
        "--distant-label",
        "--out-dir", file.path(meta@tmpdir), # output dir
        "--threads", args$threads,
        ">", file.path(meta@logdir, "gappa.log"), "2>&1" # write log
    )
    sh(cmd, echo = T)
    # analyze pp results
    invalid_tax <- c("DISTANT", "enterica", "bongori", "bongori;outgroup")
    asgmnts <- fread(file.path(wf_meta@tmpdir, "per_query.tsv"), sep = "\t")
    asgmnts.val <- filter(asgmnts, LWR >= args$min_lr, !(taxopath %in% invalid_tax))
    pp_res <- asgmnts.val %>%
        rename("id" = "name",
               "pp_clust" = "taxopath") %>%
        mutate(pp_clust = str_replace_all(pp_clust, ".*;", "")) %>%
        select(id, pp_clust)
    write.table(pp_res, file.path(meta@tmpdir, "pp_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
    return(pp_res)
})

# main search 
setMethod("main_search",
signature="Workflow metadata",
function(
    meta,
    args
) {
        # run best hit
        bh_res <- bh_search(meta, args)
        # run pp
        pp_res <- pp_search(meta, args)
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
        write.table(search_res, file.path(meta@tmpdir, "samnsorter_res.tsv"),
                sep = "\t", row.names = F, quote = F)
})

# main workflow
tic()
cat("This is SamnSorter", paste0("v", version, "\n"))
wf_meta <- new("Workflow metadata")
wf_meta <- check_args(args)
tmpdir <- file.path(args$outdir, "6cbdab8f")
wf_meta@tmpdir <- tmpdir
#tmpdir <- file.path(args$outdir, "demo", "b5329061")
#input_setup(wf_meta, args)
cgmlst(wf_meta, args)
cgmlst_dist_query(wf_meta, args)
main_search(wf_meta, args)
cat("Workflow completed successfully.\n")
toc(log = T)

