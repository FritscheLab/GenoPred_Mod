#!/usr/bin/Rscript
# This script was added by Lars Fritsche inpired by the original code by Oliver Pain.
# Modified to generate raw, unscaled PGS for all individuals regardless of ancestry
start.time <- Sys.time()
library("optparse")

option_list <- list(
    make_option("--config",
        action = "store", default = NULL, type = "character",
        help = "Pipeline configuration file [required]"
    ),
    make_option("--name",
        action = "store", default = NULL, type = "character",
        help = "Name of target sample [required]"
    ),
    make_option("--score",
        action = "store", default = NULL, type = "character",
        help = "Score to be used [required]"
    ),
    make_option("--plink2",
        action = "store", default = "plink2", type = "character",
        help = "Path PLINK v2 software binary [optional]"
    ),
    make_option("--n_cores",
        action = "store", default = 1, type = "numeric",
        help = "Number of cores to use [optional]"
    ),
    make_option("--test",
        action = "store", default = NA, type = "character",
        help = "Specify number of SNPs to include [optional]"
    ),
    make_option("--memory",
        action = "store", default = 5000, type = "numeric",
        help = "Memory limit [optional]"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
library(bigstatsr)
source("../functions/misc.R")
source_all("../functions")
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

# Check required inputs
if (is.null(opt$config)) {
    stop("--config must be specified.\n")
}
if (is.null(opt$name)) {
    stop("--name must be specified.\n")
}

# Read in outdir
outdir <- read_param(config = opt$config, param = "outdir", return_obj = F)

# Create output directory for raw PGS
opt$output <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry")
system(paste0("mkdir -p ", opt$output))

# Create logs directory
logs_dir <- paste0(outdir, "/", opt$name, "/logs")
system(paste0("mkdir -p ", logs_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(logs_dir, "/target_pgs_raw_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log")
log_header(log_file = log_file, opt = opt, script = "target_scoring_pipeline_raw.R", start.time = start.time)

# If testing, change CHROMS to chr value
if (!is.na(opt$test) && opt$test == "NA") {
    opt$test <- NA
}
if (!is.na(opt$test)) {
    CHROMS <- as.numeric(gsub("chr", "", opt$test))
}

# Identify score files to be combined
score_files <- list_score_files(opt$config)

if (is.null(score_files) || nrow(score_files) == 0) {
    log_add(log_file = log_file, message = paste0("No score files found for raw PGS generation."))
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    sink(file = log_file, append = T)
    cat("Analysis finished at", as.character(end.time), "\n")
    cat("Analysis duration was", as.character(round(time.taken, 2)), attr(time.taken, "units"), "\n")
    sink()
    quit(save = "no", status = 0)
}

# Check whether score files or target genetic data are newer than target pgs
ancestry_reporter_file <- paste0(outdir, "/reference/target_checks/", opt$name, "/ancestry_reporter.done")
ancestry_reporter_file_time <- file.info(ancestry_reporter_file)$mtime
score_files_to_do <- data.table()
for (i in 1:nrow(score_files)) {
    pgs_i <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i], "/", opt$name, "-", score_files$name[i], "-AllAncestry-raw.profiles")
    score_i <- paste0(outdir, "/reference/pgs_score_files/", score_files$method[i], "/", score_files$name[i], "/ref-", score_files$name[i], ".score.gz")
    if (!file.exists(pgs_i)) {
        score_files_to_do <- rbind(score_files_to_do, score_files[i, ])
    } else {
        score_i_time <- file.info(score_i)$mtime
        pgs_i_time <- file.info(pgs_i)$mtime
        if (score_i_time > pgs_i_time | ancestry_reporter_file_time > pgs_i_time) {
            score_files_to_do <- rbind(score_files_to_do, score_files[i, ])
            system(paste0("rm ", pgs_i))
        }
    }
}
log_add(log_file = log_file, message = paste0("After checking timestamps, ", nrow(score_files_to_do), "/", nrow(score_files), " score files will be used for raw target scoring."))
score_files <- score_files_to_do

if (nrow(score_files) == 0) {
    log_add(log_file = log_file, message = paste0("No score files to be processed for raw PGS generation."))
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    sink(file = log_file, append = T)
    cat("Analysis finished at", as.character(end.time), "\n")
    cat("Analysis duration was", as.character(round(time.taken, 2)), attr(time.taken, "units"), "\n")
    sink()
    quit(save = "no", status = 0)
}

# Subset score files
if (!is.null(opt$score)) {
    if (all(score_files$name != opt$score)) {
        stop("Requested score files not present in gwas_list or score_list")
    }
    score_files <- score_files[score_files$name == opt$score, ]
}

# Read in target_list
target_list <- read_param(config = opt$config, param = "target_list", return_obj = T)

# Set params for plink_score - use ALL individuals (no population filtering)
opt$target_plink_chr <- paste0(outdir, "/", opt$name, "/geno/", opt$name, ".ref.chr")
opt$target_keep <- NULL # No population filtering - include ALL individuals

refdir <- read_param(config = opt$config, param = "refdir", return_obj = F)

# Read in reference SNP data
ref <- read_pvar(paste0(refdir, "/ref.chr"), chr = CHROMS)[, c("CHR", "SNP", "A1", "A2"), with = F]

log_add(log_file = log_file, message = paste0("Generating raw, unscaled PGS for ALL individuals (", nrow(score_files), " score files)"))

# We will process score files and perform target scoring for one chromosome for efficiency
for (chr_i in CHROMS) {
    log_add(log_file = log_file, message = "########################")
    log_add(log_file = log_file, message = paste0("Processing chromosome ", chr_i, ":"))

    #####
    # Combine score files
    #####
    # Create row number index to subset score files by chromosome
    row_index <- format(which(ref$CHR == chr_i) + 1, scientific = FALSE)
    write.table(row_index, paste0(tmp_dir, "/row_index.txt"), row.names = F, quote = F, col.names = F)

    # Create file containing SNP, A1, and A2 information for each chromosome
    fwrite(ref[ref$CHR == chr_i, c("SNP", "A1", "A2"), with = F], paste0(tmp_dir, "/map.txt"), row.names = F, quote = F, sep = " ")

    # Extract process score files for each name (gwas/score) in parallel
    foreach(i = 1:nrow(score_files), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
        system(paste0(
            "zcat ", outdir, "/reference/pgs_score_files/", score_files$method[i], "/", score_files$name[i], "/ref-", score_files$name[i], ".score.gz | ",
            "awk 'NR==FNR {rows[$1]; next} FNR==1 || FNR in rows' ", paste0(tmp_dir, "/row_index.txt"), " - | ", # Corrected to retain the header and process indexed rows
            "cut -d' ' --complement -f1-3 | ", # Keep relevant columns, remove first 3
            "sed '1 s/SCORE_/", paste0("score_file_", i, "."), "/g' > ", # Replace SCORE in the header
            tmp_dir, "/tmp_score.", paste0(score_files$method[i], ".", score_files$name[i]), ".txt"
        ))
    }

    # Paste files together in batches
    # Set number of batches according to the number of score files to combine
    num_batches <- max(c(1, min(c(opt$n_cores, floor(nrow(score_files) / 2)))))
    tmp_score_files <- paste0(tmp_dir, "/tmp_score.", score_files$method, ".", score_files$name, ".txt")
    set.seed(1)
    batches <- split(sample(tmp_score_files), rep(1:num_batches, length.out = length(tmp_score_files)))
    log_add(log_file = log_file, message = paste0("Aggregating score files in ", num_batches, " batches."))
    foreach(i = 1:length(batches), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
        system(paste0("paste -d ' ' ", paste(batches[[i]], collapse = " "), " > ", tmp_dir, "/tmp_batch_", i))
        system(paste0("rm ", paste(batches[[i]], collapse = " ")))
    }

    # Paste batches together
    log_add(log_file = log_file, message = paste0("Aggregating batched score files."))
    tmp_batch_files <- paste0(tmp_dir, "/tmp_batch_", 1:length(batches))
    system(paste0("paste -d ' ' ", tmp_dir, "/map.txt ", paste(tmp_batch_files, collapse = " "), " > ", tmp_dir, "/all_score.txt"))
    system(paste0("rm ", paste(tmp_batch_files, collapse = " ")))

    # Perform polygenic risk scoring
    scores_i <-
        plink_score(
            pfile = opt$target_plink_chr,
            chr = chr_i,
            plink2 = opt$plink2,
            score = paste0(tmp_dir, "/all_score.txt"),
            keep = opt$target_keep, # This is NULL - includes ALL individuals
            threads = opt$n_cores,
            fbm = T
        )

    # Sum scores across chromosomes
    if (chr_i == CHROMS[1]) {
        scores_ids <- scores_i$ids
        cols <- scores_i$cols

        # Initialize a FBM (backed on disk) for running PGS sum
        file.remove(paste0(tmp_dir, "/PGS_raw_fbm.bk"))
        scores <- FBM(
            nrow = nrow(scores_ids),
            ncol = length(cols),
            backingfile = paste0(tmp_dir, "/PGS_raw_fbm"),
            init = 0
        )
    }

    # In-place addition: for each score column
    for (j in cols) {
        scores[, which(cols == j)] <- scores[, which(cols == j)] + scores_i$scores[, which(scores_i$cols == j)]
    }

    file.remove(
        scores_i$scores$backingfile,
        scores_i$scores$rds
    )
    rm(scores_i)
    gc()
}

# Combine score with IDs
scores <- as.data.table(matrix(scores[, ], ncol = length(cols)))
setnames(scores, cols)
scores <- cbind(scores_ids, scores)

###
# Save the allscore file
###
log_add(log_file = log_file, message = "Saving allscore file for raw PGS...")
dir.create(paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/"), recursive = T)
fwrite(scores, paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", opt$name, "-AllAncestry-raw.profiles"), sep = " ", na = "NA", quote = F)


###
# Save RAW, UNSCALED polygenic scores - NO SCALING APPLIED
###

log_add(log_file = log_file, message = "Saving raw, unscaled polygenic scores for ALL individuals (no ancestry filtering, no scaling).")

for (i in 1:nrow(score_files)) {
    scores_i <- scores[, c("FID", "IID", names(scores)[grepl(paste0("^score_file_", i, "\\."), names(scores))]), with = F]
    names(scores_i) <- gsub(paste0("^score_file_", i, "\\."), paste0(score_files$name[i], "_"), names(scores_i))
    dir.create(paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i]), recursive = T)
    fwrite(scores_i, paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i], "/", opt$name, "-", score_files$name[i], "-AllAncestry-raw.profiles"), sep = " ", na = "NA", quote = F)
}

log_add(log_file = log_file, message = paste0("Saved raw, unscaled polygenic scores for ALL individuals to: ", outdir, "/", opt$name, "/pgs_raw/AllAncestry/"))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat("Analysis finished at", as.character(end.time), "\n")
cat("Analysis duration was", as.character(round(time.taken, 2)), attr(time.taken, "units"), "\n")
sink()
