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
raw_pgs_dir <- read_param(config = opt$config, param = "raw_pgs_dir", return_obj = F)
if (!is.null(raw_pgs_dir) && length(raw_pgs_dir) > 1) {
    raw_pgs_dir <- raw_pgs_dir[1]
}

# Create output directory for raw PGS
opt$output <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry")
system(paste0("mkdir -p ", opt$output))

# Store the combined plink --score weight file alongside the raw profiles
weights_dir <- opt$output
dir.create(weights_dir, recursive = T, showWarnings = F)
weight_file <- paste0(weights_dir, "/", opt$name, "-AllAncestry-raw.weights.score")

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

# Sort score files to ensure consistent processing order
if (!is.null(score_files)) {
    score_files <- score_files[order(score_files$method, score_files$name), ]
}
score_files_all <- copy(score_files)

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
    pgs_i <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i], "/", opt$name, "-", score_files$method[i], "-", score_files$name[i], "-AllAncestry-raw.profiles")
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
if (!file.exists(weight_file)) {
    log_add(log_file = log_file, message = "Combined weight file missing; regenerating raw scores to rebuild weight file.")
    score_files_to_do <- copy(score_files)
}
log_add(log_file = log_file, message = paste0("After checking timestamps, ", nrow(score_files_to_do), "/", nrow(score_files), " score files will be used for raw target scoring."))
# Sort score files to ensure consistent processing order
score_files <- score_files_to_do[order(score_files_to_do$method, score_files_to_do$name), ]

if (nrow(score_files) == 0) {
    log_add(log_file = log_file, message = paste0("No score files to be processed for raw PGS generation. Ensuring method-level outputs exist."))
    if (file.exists(weight_file)) {
        weights_header <- fread(weight_file, nrows = 0)
        for (i in 1:nrow(score_files_all)) {
            cols_i <- names(weights_header)[grepl(paste0("^score_file_", i, "\\."), names(weights_header))]
            if (length(cols_i) == 0) {
                next
            }
            weights_subset <- fread(weight_file, select = c("SNP", "A1", "A2", cols_i))
            weights_out_dir <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files_all$method[i], "/", score_files_all$name[i])
            dir.create(weights_out_dir, recursive = TRUE, showWarnings = FALSE)
            weights_out_new <- paste0(weights_out_dir, "/", opt$name, "-", score_files_all$method[i], "-", score_files_all$name[i], "-AllAncestry-raw.weights.txt")
            setnames(weights_subset, old = cols_i, new = gsub(paste0("^score_file_", i, "\\."), paste0(score_files_all$method[i], "_", score_files_all$name[i], "_"), cols_i))
            fwrite(weights_subset, weights_out_new, sep = " ", na = "NA", quote = FALSE)
            log_add(log_file = log_file, message = paste0("Saved weight file for ", score_files_all$method[i], "/", score_files_all$name[i], ": ", weights_out_new))
        }
    }
    # Ensure method name appears in per-method profile filenames
    for (i in 1:nrow(score_files_all)) {
        profiles_dir_i <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files_all$method[i], "/", score_files_all$name[i])
        method_profile <- paste0(profiles_dir_i, "/", opt$name, "-", score_files_all$method[i], "-", score_files_all$name[i], "-AllAncestry-raw.profiles")
        if (!file.exists(method_profile)) {
            dir.create(profiles_dir_i, recursive = TRUE, showWarnings = FALSE)
        }
    }
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

# Ensure we start with a fresh combined weight file when running
if (file.exists(weight_file)) {
    file.remove(weight_file)
}

# We will process score files and perform target scoring for one chromosome for efficiency
for (chr_i in CHROMS) {
    log_add(log_file = log_file, message = "########################")
    log_add(log_file = log_file, message = paste0("Processing chromosome ", chr_i, ":"))
    tmp_all_score <- paste0(tmp_dir, "/all_score_chr", chr_i, ".txt")
    if (file.exists(tmp_all_score)) {
        file.remove(tmp_all_score)
    }

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
    # Use deterministic ordering instead of random sampling to ensure consistent column order
    batches <- split(tmp_score_files, rep(1:num_batches, length.out = length(tmp_score_files)))
    log_add(log_file = log_file, message = paste0("Aggregating score files in ", num_batches, " batches."))
    foreach(i = 1:length(batches), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
        system(paste0("paste -d ' ' ", paste(batches[[i]], collapse = " "), " > ", tmp_dir, "/tmp_batch_", i))
        system(paste0("rm ", paste(batches[[i]], collapse = " ")))
    }

    # Paste batches together
    log_add(log_file = log_file, message = paste0("Aggregating batched score files."))
    tmp_batch_files <- paste0(tmp_dir, "/tmp_batch_", 1:length(batches))
    system(paste0("paste -d ' ' ", tmp_dir, "/map.txt ", paste(tmp_batch_files, collapse = " "), " > ", tmp_all_score))
    system(paste0("rm ", paste(tmp_batch_files, collapse = " ")))
    weight_stat <- file.info(tmp_all_score)
    if (!file.exists(tmp_all_score) || is.na(weight_stat$size) || weight_stat$size == 0) {
        stop(paste0("Weight file missing or empty: ", tmp_all_score))
    }
    # Append to run-level weight file (header only once)
    if (!file.exists(weight_file)) {
        file.copy(tmp_all_score, weight_file, overwrite = TRUE)
    } else {
        system(paste0("tail -n +2 ", shQuote(tmp_all_score), " >> ", shQuote(weight_file)))
    }
    log_add(log_file = log_file, message = paste0("Saved plink --score weight file for chr ", chr_i, ": ", tmp_all_score))

    # Perform polygenic risk scoring
    scores_i <-
        plink_score(
            pfile = opt$target_plink_chr,
            chr = chr_i,
            plink2 = opt$plink2,
            score = tmp_all_score,
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
combined_profiles_file <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", opt$name, "-AllAncestry-raw.profiles")
fwrite(scores, combined_profiles_file, sep = " ", na = "NA", quote = F)


###
# Save RAW, UNSCALED polygenic scores - NO SCALING APPLIED
###

log_add(log_file = log_file, message = "Saving raw, unscaled polygenic scores for ALL individuals (no ancestry filtering, no scaling).")

for (i in 1:nrow(score_files)) {
    scores_i <- scores[, c("FID", "IID", names(scores)[grepl(paste0("^score_file_", i, "\\."), names(scores))]), with = F]
    names(scores_i) <- gsub(paste0("^score_file_", i, "\\."), paste0(score_files$method[i], "_", score_files$name[i], "_"), names(scores_i))
    profiles_dir_i <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i])
    dir.create(profiles_dir_i, recursive = T)
    profile_file_new <- paste0(profiles_dir_i, "/", opt$name, "-", score_files$method[i], "-", score_files$name[i], "-AllAncestry-raw.profiles")
    fwrite(scores_i, profile_file_new, sep = " ", na = "NA", quote = F)
}

log_add(log_file = log_file, message = paste0("Saved raw, unscaled polygenic scores for ALL individuals to: ", outdir, "/", opt$name, "/pgs_raw/AllAncestry/"))

weight_stat <- file.info(weight_file)
if (!file.exists(weight_file) || is.na(weight_stat$size) || weight_stat$size == 0) {
    stop(paste0("Missing weight file after scoring: ", weight_file))
}
# Sanity check: ensure combined weight file has rows for all reference SNPs across chromosomes
weight_rows <- nrow(fread(weight_file, select = "SNP"))
expected_rows <- nrow(ref)
if (weight_rows < expected_rows) {
    stop(paste0("Combined weight file is incomplete: found ", weight_rows, " rows, expected ", expected_rows, "."))
}
log_add(log_file = log_file, message = paste0("All plink --score weights written to ", weight_file))

# Save method-level weight files next to their profiles for transparency
weights_header <- fread(weight_file, nrows = 0)
for (i in 1:nrow(score_files)) {
    cols_i <- names(weights_header)[grepl(paste0("^score_file_", i, "\\."), names(weights_header))]
    if (length(cols_i) == 0) {
        next
    }
    weights_subset <- fread(weight_file, select = c("SNP", "A1", "A2", cols_i))
    weights_out_dir <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i])
    dir.create(weights_out_dir, recursive = TRUE, showWarnings = FALSE)
    weights_out_new <- paste0(weights_out_dir, "/", opt$name, "-", score_files$method[i], "-", score_files$name[i], "-AllAncestry-raw.weights.txt")
    setnames(weights_subset, old = cols_i, new = gsub(paste0("^score_file_", i, "\\."), paste0(score_files$method[i], "_", score_files$name[i], "_"), cols_i))
    fwrite(weights_subset, weights_out_new, sep = " ", na = "NA", quote = FALSE)
    log_add(log_file = log_file, message = paste0("Saved weight file for ", score_files$method[i], "/", score_files$name[i], ": ", weights_out_new))
}

# Optionally copy raw outputs to an external raw_pgs directory specified in the config
if (!is.null(raw_pgs_dir) && !is.na(raw_pgs_dir)) {
    copy_base <- paste0(raw_pgs_dir, "/", opt$name, "/pgs_raw/AllAncestry")
    dir.create(copy_base, recursive = TRUE, showWarnings = FALSE)

    # Copy combined outputs
    file.copy(combined_profiles_file, paste0(copy_base, "/", basename(combined_profiles_file)), overwrite = TRUE)
    file.copy(weight_file, paste0(copy_base, "/", basename(weight_file)), overwrite = TRUE)

    # Copy per-method outputs
    for (i in 1:nrow(score_files)) {
        src_dir <- paste0(outdir, "/", opt$name, "/pgs_raw/AllAncestry/", score_files$method[i], "/", score_files$name[i])
        dest_dir <- paste0(copy_base, "/", score_files$method[i], "/", score_files$name[i])
        dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

        src_profile <- paste0(src_dir, "/", opt$name, "-", score_files$method[i], "-", score_files$name[i], "-AllAncestry-raw.profiles")
        src_weight <- paste0(src_dir, "/", opt$name, "-", score_files$method[i], "-", score_files$name[i], "-AllAncestry-raw.weights.txt")

        if (file.exists(src_profile)) {
            file.copy(src_profile, paste0(dest_dir, "/", basename(src_profile)), overwrite = TRUE)
        }
        if (file.exists(src_weight)) {
            file.copy(src_weight, paste0(dest_dir, "/", basename(src_weight)), overwrite = TRUE)
        }
    }
    log_add(log_file = log_file, message = paste0("Copied raw PGS outputs to ", copy_base))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat("Analysis finished at", as.character(end.time), "\n")
cat("Analysis duration was", as.character(round(time.taken, 2)), attr(time.taken, "units"), "\n")
sink()
