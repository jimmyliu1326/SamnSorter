#!/usr/bin/env Rscript

# load pkgs
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(uwot))
suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(kknn))
suppressPackageStartupMessages(library(parsnip))

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
            args="list",
            query_dist="list"
        )
)

# S4 generics
invisible(setGeneric("input_setup", function(meta, args) standardGeneric("input_setup")))
invisible(setGeneric("cgmlst", function(meta, args) standardGeneric("cgmlst")))
invisible(setGeneric("cgmlst_dist_query", function(meta, cgmlst_path) standardGeneric("cgmlst_dist_query")))
invisible(setGeneric("main_search", function(meta, args) standardGeneric("main_search")))

# cli argument parser
parser <- ArgumentParser()
parser$add_argument("file", nargs='+')
parser$add_argument("-o", "--outdir", type="character", default="./samnsorter",
    help="Directory path to where output files will be written to [%(default)s]")
parser$add_argument("--tmpdir", type="character", default="",
    help="Directory path to where temporary files will be stored [%(default)s]")
parser$add_argument("--min_lr", type="double", default=0.25,
    help="Minimum likelihood weighted ratio for phylogenetic placements [%(default)s]")
parser$add_argument("--min_lof", type="double", default=4.0,
    help="Minimum local outlier factor score [%(default)s]")
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
check_args <- function(args) {
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
        cmd <- paste("ln -s", realpath, file.path(meta@tmpdir, "input", basename(realpath)))
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
    cgmlst_profile <- file.path(meta@outdir, "cgmlst", "results_alleles_hashed.tsv")
    return(cgmlst_profile)
})

# calculate query distance in respect to references
setMethod("cgmlst_dist_query",
signature="Workflow metadata",
function(
    meta,
    cgmlst_path
) {
    cat("Calculating allelic distance between query and reference...\n")
    cmd <- paste("$CGMLST_DISTS", "-H", 
                 "$REF_ALLELES", # reference profiles
                 cgmlst_path, # query profiles
                 ">", file.path(meta@tmpdir, "cgmlst_dist.tsv"), # output distance matrix
                 "2>", file.path(meta@logdir, "cgmlst-dists.log"))
    sh(cmd, echo = T)
    # reformat distance matrix
    dist <- fread(file.path(meta@tmpdir, "cgmlst_dist.tsv"), sep = "\t")
    dist_out <- dist %>% column_to_rownames("cgmlst-dists")
    dist_out.path <- file.path(meta@outdir, "cgmlst_dist.tsv")
    #rownames(dist_out) <- str_replace_all(rownames(dist_out), "\\.1|\\.2", "")
    #colnames(dist_out) <- str_replace_all(colnames(dist_out), "\\.1|\\.2", "")
    cat("Saving distance matrix to:", dist_out.path, "\n")
    write.table(dist_out, dist_out.path, sep = "\t", 
                quote = F, row.names = T, col.names = NA)
    return(list("mat" = dist_out, "path" = dist_out.path))
})

# best hit search
bh_search <- function(
    dist # query-to-ref distance data frame
) {
    cat("Identifying best hit...\n")
    bh <- apply(dist[,1:ncol(dist)-1], # drop last column to ignore outgroup
                1, function(x) { 
        idx <- which.min(x)[1]
        list("dist" = x[idx],
             "ref_id" = names(x)[idx]
        )
    })
    bh_dist <- map_dbl(bh, ~return(.$dist))
    bh_id <- map_chr(bh, ~return(.$ref_id))
    # get cluster ID of bh
    clust_path <- Sys.getenv("REF_CLUSTERS")
    clust <- fread(clust_path, sep = "\t")
    bh_clust <- clust$clust[match(bh_id, clust$id)]
    bh_res <- data.frame("id" = rownames(dist),
                          "best_hit" = bh_clust,
                          "best_hit_dist" = bh_dist)
    return(bh_res)
}

# get nearest neighbour for asymmetrical distance matrix
query_nn <- function(X, k, include_self = TRUE) {
  X <- as.matrix(X)
  if (nrow(X) == 1) {
    nn_idx <- t(as.matrix(t(apply(X, 1, order))[, 1:k]))
  } else {
    nn_idx <- t(apply(X, 1, order))[, 1:k]
  }
  nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
  for (i in seq_len(nrow(nn_idx))) {
    nn_dist[i, ] <- X[i, nn_idx[i, ]]
  }
  attr(nn_idx, "dimnames") <- NULL
  attr(nn_dist, "dimnames") <- NULL
  list(idx = nn_idx, dist = nn_dist)
}

# kkNN prediction on UMAP embeddings
knn_search <- function(
    dist,       # query-to-ref distance
    umap_model, # path to umap model
    knn_model,  # path to knn model
    min_lof     # outlier score threshold
) {
    set.seed(123)
    # compute nearest neighbour for new data
    query.nn <- query_nn(dist[,1:ncol(dist)-1], k = 15)
    # transform new data into embedding
    query.embed <- umap_transform(X = dist[,1:ncol(dist)-1],
                                  model = umap_model,
                                  nn_method = query.nn)
    # compute local outlier factor (lof) for query embedding
    query.lof <- map_dbl(seq_along(rownames(dist)), function(x) {
        lof(rbind(umap_model$embedding, query.embed[x,]))[nrow(umap_model$embedding)+1]
    })
    # to detect query outliers
    # kNN classification on transformed query
    query.clust <- predict(knn_model, new_data = as.data.frame(query.embed))
    # ignore predictions with lof >= min_lof
    query.clust[which(query.lof >= min_lof)] <- "NOVEL"
    # generate output table
    knn_res <- data.frame('id' = rownames(dist),
                          'knn_clust' = query.clust$.pred_class,
                          'lof_score' = query.lof)
    return(knn_res)
}

# phylogenetic placement using APPLES
pp_search <- function(
    dist_file, # path to query-ref distance matrix
    tmpdir,    # path to tmp dir
    logdir,    # path to log dir
    threads,   # number of threads
    min_lr     # minimum likelihood ratio
) {

    cat("Running phylogenetic placement...\n")
    cmd <- paste(
        "run_apples.py", 
        "-t", "$REF_NWK", # backbone tree
        "-d", dist_file, # query distance matrix
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
        "--distant-label",
        "--out-dir", tmpdir, # output dir
        "--threads", threads,
        ">", file.path(logdir, "gappa.log"), "2>&1" # write log
    )
    sh(cmd, echo = T)
    # analyze pp results
    invalid_tax <- c("DISTANT", "enterica", "bongori", "bongori;outgroup")
    asgmnts <- fread(file.path(tmpdir, "per_query.tsv"), sep = "\t")
    pp_res <- asgmnts %>%
        group_by(name) %>%
        split(f = as.factor(.$name)) %>%
        map_dfr(function(x) {
            # get sample id
            id <- x$name[1]
            # get best LWR
            lwr <- max(x$LWR)[1]
            # get assignment
            val <- str_replace_all(pull(filter(x, LWR >= min_lr, !(taxopath %in% invalid_tax)), taxopath),
                                   ".*;", "")
            if (length(val) == 0) val <- "NOVEL"
            # return
            data.frame("id" = id,
                       "pp_clust" = val,
                       "LWR" = lwr)
        })
    return(pp_res)
}

# main search 
setMethod("main_search",
signature="Workflow metadata",
function(
    meta,
    args
) {
        # run best hit
        bh_res <- bh_search(meta@query_dist$mat)
        write.table(bh_res, file.path(meta@tmpdir, "best_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
        # run pp
        pp_res <- pp_search(meta@query_dist$path, 
                            tmpdir = meta@tmpdir,
                            logdir = meta@logdir,
                            threads = args$threads,
                            min_lr = args$min_lr
        )
        write.table(pp_res, file.path(meta@tmpdir, "pp_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
        # run knn classification
        cat("Cluster prediction using kNN on UMAP embedding...\n")
        model_dir <- Sys.getenv("MODEL_DIR")
        # load umap model
        umap_model <- readRDS(file.path(model_dir, "umap_k15_a1_b1_dim500_8154.RDS"))
        # load knn model
        knn_model <- readRDS(file.path(model_dir,"kknn_k5_k15_a1_b1_dim500_8154.RDS"))
        knn_res <- map_dfr(seq_along(rownames(meta@query_dist$mat)), 
                           ~knn_search(
                                dist=meta@query_dist$mat[.,], # query-to-ref distance
                                umap_model,
                                knn_model,
                                args$min_lof # outlier score threshold
                           )
        )
        write.table(knn_res, file.path(meta@tmpdir, "knn_hit.tsv"),
                sep = "\t", row.names = F, quote = F)
        # merge results
        search_res <- full_join(pp_res, bh_res, by = "id") %>%
            full_join(knn_res, by = "id") %>%
            mutate(pp_clust = if_else(is.na(pp_clust), "NOVEL", pp_clust),
                   knn_clust = if_else(is.na(knn_clust), "NOVEL", knn_clust)
            )
        # final prediction by majority voting
        final_preds <- apply(search_res[,c("best_hit", "pp_clust", "knn_clust")], 1,
                             function(x) {
                                if (length(x) == length(unique(x))) {
                                    return("NOVEL")
                                } else {
                                    return(as.character(x[which.max(table(x))]))
                                }
                             })
        search_res <- cbind(search_res, final_clust = final_preds) %>%
            select(id, best_hit, pp_clust, knn_clust, final_clust, everything())
        cat("Writing cluster assignment results...\n")
        write.table(search_res, file.path(meta@tmpdir, "samnsorter_res.tsv"),
                sep = "\t", row.names = F, quote = F)
})

# main workflow
tic()
cat("This is SamnSorter", paste0("v", version, "\n"))
wf_meta <- new("Workflow metadata")
wf_meta <- check_args(args)
tmpdir <- file.path(wf_meta@outdir, "6cbdab8f")
wf_meta@tmpdir <- tmpdir
#input_setup(wf_meta, args)
#hashed_profiles <- cgmlst(wf_meta, args)
hashed_profiles <- file.path(wf_meta@tmpdir, "cgmlst", "results_alleles_hashed.tsv")
wf_meta@query_dist <- cgmlst_dist_query(wf_meta, hashed_profiles)
main_search(wf_meta, args)
cat("Workflow completed successfully.\n")
toc(log = T)