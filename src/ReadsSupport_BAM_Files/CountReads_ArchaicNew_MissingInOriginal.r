#!/usr/bin/env Rscript
# ======================================================================
# README — Counts at “new vs Original” sites in Archaic New ROIs
# ======================================================================
# Purpose
#   For each individual and each imputed coverage label, find sites in
#   Archaic New regions that are “new” relative to the Original VCF
#   (absent in Original OR present but Original GT is missing), then
#   pile up BAM reads at those sites to quantify A/C/G/T/N and derive
#   ref/alt counts and proportions.
#
# Inputs (set below in Config)
#   - VCFs (bgzip + tabix): one per coverage label + one "Original".
#   - BAMs (+ .bai): one per individual (truth reads for pileup).
#   - MergeFragmentsArchaic_*.csv: introgression map with segments
#     annotated by Individual, Coverage, Filter, segment_type, chr/start/end.
#
# Outputs
#   - OutCounts_perIndividual_CorrectOriginal/
#       merged_counts_gt_<IND>_AllCoverages.tsv.gz
#     Columns include:
#       Individual, Sample, label, cov, chr, pos, ref, alt, GT,
#       A, C, G, T, N, dp_total, ref_count, alt_count, prop_ref, prop_alt
#     Note: dp_total = A+C+G+T (N excluded by design).
#
# Requirements
#   - Same genome build across VCFs, BAMs, and ROI coordinates.
#   - VCFs are bgzipped (.vcf.gz) and tabix-indexed (.tbi).
#   - BAMs are indexed (.bai). Indexing is auto-created if missing.
#   - R packages: VariantAnnotation, GenomicRanges, GenomeInfoDb, IRanges,
#                 Rsamtools, data.table, dplyr, BiocParallel.
#
# Quick usage
#   $ Rscript this_script.R
#   (Optionally edit Config to switch GP99 vs AllSNPs, paths, and workers.)
#
# What counts as “new vs Original” (per coverage label, within that label’s ROI)
#   (a) Site exists in Coverage VCF but NOT in Original VCF by chr:pos, OR
#   (b) Site exists in BOTH, but Original GT is missing (NA, "./.", or ".|.").
#
# High-level flow
#   1) Read introgression map; identify Individuals with BAMs.
#   2) For each Individual (parallel over individuals):
#      a. Check BAM exists & index if needed.
#      b. Build ROIs: Kept ∩ (Archaic New) ∩ labels_target; reduce overlaps.
#      c. Open all VCFs once; resolve correct sample names (aliases allowed).
#      d. For each coverage label:
#           - Fetch fast SNP keys (chr,pos,ref,alt) in ROI via tabix.
#           - Compute “new” positions by (a) absence in Original OR (b) GT missing.
#           - Attach coverage ref/alt; tag with cov/label.
#      e. Stack per-label results → new_all; fetch GTs per label at those sites.
#      f. Union all sites once; run Rsamtools::pileup on BAM (parallel by chunks).
#      g. Join counts back; compute ref_count/alt_count and prop_ref/prop_alt.
#      h. Write merged TSV.GZ to output directory.
#
# Parallelism knobs
#   - workers_indiv   : outer parallelism (how many individuals at once).
#   - workers_pileup  : inner parallelism per individual (pileup chunks).
#   - target_tasks_factor × workers_pileup ≈ number of chunks per individual.
#   - chunk_min       : minimum sites per chunk (reduces overhead on tiny chunks).
#   Be mindful of storage I/O: effective concurrency ≈ workers_indiv × workers_pileup.
#
# Quality thresholds
#   - min_mapq, min_baseq are passed to PileupParam. Currently 0 (permissive).
#     Typical stricter values: min_mapq=20–30, min_baseq=20.
#
# Contig naming
#   - The script auto-normalizes “chr” prefixes for both VCF and BAM access.
#     Output chr column is normalized WITHOUT “chr”.
#
# Notes & gotchas
#   - Only single-base A/C/G/T SNPs are considered (indels/MNPs dropped).
#   - If a base never appears at a site, its column is added with 0.
#   - Errors for an individual are caught and logged; other individuals continue.
#   - Consider raising thresholds or additional filters if you observe excess noise.
# ======================================================================

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(data.table)
  library(Rsamtools)
  library(dplyr)
  library(BiocParallel)
})

options(stringsAsFactors = FALSE)

setwd("/st1hdd/Huerta-Sanchez/Marco/Analyses/ImputationTest/1.TestImputation/DownSamples/Selection/DenisovaUstIshim/")

# ── Config ─────────────────────────────────────────────────────────────────────
labels_target <- c("2X","1X","0.5X","0.0625X","ImputedOC")

# Choose ONE block (comment the other). Here: GP99; swap to AllSNPs if needed
#vcf_paths <- c(
#  "2X"        = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed2X_GP99_Phased.vcf.gz",
#  "1X"        = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed_1X_GP99_Phased.vcf.gz",
#  "0.0625X"   = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed_0625X_GP99_Phased.vcf.gz",
#  "0.5X"      = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed_05X_GP99_Phased.vcf.gz",
#  "ImputedOC" = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/ImputedHC_GP99_phased.vcf.gz",
#  "Original"  = "/st1hdd/Huerta-Sanchez/Marco/Analyses/ImputationTest/1.TestImputation/DownSamples/Dataset/HighCoverage/highCoverage_Anno_NoInfo.vcf.gz"
#)

# # AllSNPs variant:
vcf_paths <- c(
   "2X"        = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed2_Raw_Phased.vcf.gz",
   "1X"        = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed1_Raw_Phased.vcf.gz",
   "0.0625X"   = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed06_Raw_Phased.vcf.gz",
   "0.5X"      = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/Imputed05_Raw_Phased.vcf.gz",
   "ImputedOC" = "/st1hdd/Huerta-Sanchez/Downsampling/Phased/ImpHC_Raw_Phased.vcf.gz",
   "Original"  = "/st1hdd/Huerta-Sanchez/Marco/Analyses/ImputationTest/1.TestImputation/DownSamples/Dataset/HighCoverage/highCoverage_Anno_NoInfo.vcf.gz"
 )

# BAMs
bam_dir <- "/st1hdd/Huerta-Sanchez/Downsampling/Original/BAM"
bam_map <- c(
  "2H10"      = "2H10.collapsed.sorted.grouped.indel.mq25.clip2bp.34bp.duprm.bam",
  "2H11"      = "2H11.collapsed.sorted.grouped.indel.mq25.clip2bp.34bp.duprm.bam",
  "WC1"       = "WC1.all_SG_join.Mkdup.len.realg.reheader.mq25.clip2bp.34bp.bam",
  "BA64"      = "BA64-A-EX2-LIB1.trimmed25bp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam",
  "USR1"      = "USR1.reheader.clip2bp.mq37.30bp.mq25.clip2bp.34bp.bam",
  "KK1"       = "KK1.trimmed25bp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam",
  "sf12"      = "sf12_90perc_libr_160201.merge.hs37d5.fa.grouped.duprm.indel.35bp.mq25.clip2bp.35bp.bam",
  "NE1"       = "NE1.trimmed25bp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam",
  "Stuttgart" = "LBK.hg19_1000g.bam.reheader.rmdup.indel.mq25.clip2bp.34bp.bam",
  "Loschbour" = "Loschbour.hg19_1000g.bam.reheader.rmdup.indel.mq25.clip2bp.34bp.bam",
  "atp016"    = "atp016.hs37d5.fa.cons.90perc.grouped.duprm.indel.35bp.mq25.clip2bp.35bp.bam",
  "JP14"      = "JP14-A2-EX3.trimmed25bp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam",
  "PB675"     = "PB675-A1-EX3.trimmed25bp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam",
  "SRA62"     = "SRA62.trimmed25bp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam",
  "Kolyma1"   = "Kolyma_River.sort.rmdup.realign.md.reheader.clip2bp.mq37.34bp.bam",
  "Yana1"     = "Yana_old.sort.rmdup.realign.md.reheader.clip2bp.mq37.34bp.bam",
  "BOT2016"   = "BOT2016.realigned.calmd.readsadded.reheader.clip2bp.mq37.clip2bp.34bp.bam",
  "Yamnaya"   = "Yamnaya.realigned.calmd.readsadded.reheader.clip2bp.mq37.clip2bp.34bp.bam",
  "Ust_Ishim" = "Ust_Ishim.hg19_1000g.all.reheader.mq25.clip2bp.34bp.bam",
  "Mota"      = "Mota.trimmed25bp25q.fish.sorted.grouped.duprm.indel.mq20.clip2.34bp.bam"
)

# Aliases for sample name resolution in VCFs
sample_aliases <- list(
  Stuttgart = c("LBK"),
  Kolyma1   = c("Kolyma_River","Kolima1"),
  BA64      = c("BA64-A-EX2-LIB1"),
  JP14      = c("JP14-A2-EX3")
)

# Parallel settings
workers_indiv  <- 6L     # outer parallelism across individuals (start modest on HDD)
workers_pileup <- 10L     # inner pileup workers per individual
target_tasks_factor <- 6L  # make ~4×workers_pileup chunks
chunk_min      <- 2000L  # minimum positions per chunk

min_mapq <- 0L
min_baseq <- 0L

out_dir_base <- file.path(getwd(), "OutCounts_perIndividual_CorrectOriginal")
dir.create(out_dir_base, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ── IO helpers ─────────────────────────────────────────────────────────────────
open_tabix <- function(paths) {
  tfs <- lapply(paths, function(p) { tf <- TabixFile(p); open(tf); tf })
  names(tfs) <- names(paths); tfs
}
close_tabix <- function(tfs) invisible(lapply(tfs, close))

resolve_sample_name <- function(vcf_file, individual) {
  hdr <- VariantAnnotation::scanVcfHeader(vcf_file)
  smps <- VariantAnnotation::samples(hdr)
  cand <- c(individual, unlist(sample_aliases[[individual]] %||% character()))
  # exact
  for (x in cand) if (x %in% smps) return(x)
  # case-insensitive exact
  lsm <- tolower(smps)
  for (x in tolower(cand)) {
    hit <- which(lsm == x); if (length(hit)) return(smps[hit[1]])
  }
  # startsWith (case-insensitive)
  for (x in tolower(cand)) {
    hit <- which(startsWith(lsm, x)); if (length(hit)) return(smps[hit[1]])
  }
  if (length(smps) == 1L) return(smps[1])
  stop(sprintf("Cannot resolve sample for %s in %s", individual, basename(vcf_file)))
}

# Fast key extraction via Tabix (biallelic SNPs only)
vcf_keys_fast <- function(tf, gr) {
  if (length(gr) == 0L) return(data.table(chr=character(), pos=integer(), ref=character(), alt=character()))
  gr <- reduce(gr)
  lines <- unlist(scanTabix(tf, param = gr), use.names = FALSE)
  if (!length(lines)) return(data.table(chr=character(), pos=integer(), ref=character(), alt=character()))
  parts <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
  dt <- data.table(chr = parts[[1]], pos = as.integer(parts[[2]]), ref = parts[[4]], alt = parts[[5]])
  dt <- dt[nchar(ref) == 1L]
  dt[, alt := strsplit(alt, ",", fixed = TRUE)]
  dt <- dt[, .(alt = unlist(alt, use.names = FALSE)), by = .(chr, pos, ref)]
  dt <- dt[nchar(alt) == 1L &
           ref %chin% c("A","C","G","T") &
           alt %chin% c("A","C","G","T")]
  unique(dt)
}

# GT only for selected positions (fast)
gt_for_positions <- function(vcf_file, gr_positions, sample_name) {
  if (length(gr_positions) == 0L) return(data.table(chr=character(), pos=integer(), GT=character()))
  hdr <- scanVcfHeader(vcf_file)
  si  <- seqinfo(hdr)
  has_chr <- any(grepl("^chr", seqlevels(si)))
  gr2 <- GRanges(
    seqnames = if (has_chr) paste0("chr", sub("^chr","", as.character(seqnames(gr_positions))))
               else          sub("^chr","", as.character(seqnames(gr_positions))),
    ranges = ranges(gr_positions)
  )
  param <- ScanVcfParam(which = gr2, samples = as.character(sample_name), info = character(), geno = "GT")
  v <- readVcf(vcf_file, genome = si, param = param)
  if (length(v) == 0L) return(data.table(chr=character(), pos=integer(), GT=character()))
  rr <- rowRanges(v)
  gt <- geno(v)$GT
  data.table(chr = sub("^chr","", as.character(seqnames(rr))), pos = as.integer(start(rr)), GT = as.character(gt[,1]))
}

# ── Pileup helpers ─────────────────────────────────────────────────────────────
make_gr_for_bam <- function(dt_positions, bam_has_chr) {
  sn <- dt_positions$chr
  sn <- if (bam_has_chr) ifelse(grepl("^chr", sn), sn, paste0("chr", sn)) else sub("^chr", "", sn)
  GRanges(seqnames = sn, ranges = IRanges(start = dt_positions$pos, width = 1L))
}

pileup_counts <- function(bam_path, sites_dt, workers = workers_pileup) {
  if (nrow(sites_dt) == 0L) return(data.table(chr=character(), pos=integer(), A=integer(), C=integer(), G=integer(), T=integer(), N=integer(), dp_total=integer()))
  hdr <- scanBamHeader(bam_path)[[1]]$targets
  contigs <- names(hdr)
  bam_has_chr <- any(grepl("^chr", contigs))
  valid_chr <- sub("^chr", "", contigs)
  sites_dt <- sites_dt[chr %in% valid_chr]

  # Chunking: ~4×workers tasks
  n <- nrow(sites_dt)
  target_tasks <- max(target_tasks_factor * workers, workers)
  chunk_size  <- max(chunk_min, ceiling(n / target_tasks))
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  chunks <- lapply(idx, function(ii) sites_dt[ii])

  scan_flag  <- scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE)
  pile_param <- PileupParam(
    distinguish_strands     = FALSE,
    distinguish_nucleotides = TRUE,
    include_deletions       = FALSE,
    include_insertions      = FALSE,
    min_mapq                = min_mapq,
    min_base_quality        = min_baseq
  )

  BPP <- if (.Platform$OS.type == "unix") MulticoreParam(workers = workers, progressbar = TRUE) else SnowParam(workers = workers, type = "SOCK", progressbar = TRUE)

  pu_one <- function(dt_pos) {
    if (nrow(dt_pos) == 0L) return(data.table(chr=character(), pos=integer(), nucleotide=character(), count=integer()))
    gr <- make_gr_for_bam(dt_pos, bam_has_chr)
    sbp <- ScanBamParam(which = gr, flag = scan_flag)
    pu  <- Rsamtools::pileup(file = bam_path, scanBamParam = sbp, pileupParam = pile_param)
    if (nrow(pu) == 0L) return(data.table(chr=character(), pos=integer(), nucleotide=character(), count=integer()))
    dt <- as.data.table(pu)
    dt[, seqnames := sub("^chr","", seqnames)]
    setnames(dt, c("seqnames","pos","nucleotide","count"), c("chr","pos","nucleotide","count"))
    dt[, chr := as.character(chr)][, pos := as.integer(pos)]
    dt
  }

  parts <- bplapply(chunks, pu_one, BPPARAM = BPP)
  long  <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  wide  <- dcast(long, chr + pos ~ nucleotide, value.var = "count", fun.aggregate = sum)
  for (b in c("A","C","G","T","N")) if (!b %in% names(wide)) wide[, (b) := 0L]
  wide[, dp_total := A + C + G + T]
  setorder(wide, chr, pos)
  wide[]
}

# ── Data ───────────────────────────────────────────────────────────────────────
MergeFragmentsArchaic <- read.csv("MergeFragmentsArchaic_SimilarityDistance_Sharing_Filtering_Masking_Sep04.csv", header = TRUE) |> as.data.frame()

all_inds <- sort(unique(MergeFragmentsArchaic$Individual))
inds_with_bam <- intersect(all_inds, names(bam_map))

message("Individuals (with BAM): ", paste(inds_with_bam, collapse = ", "))
# Consider these patterns as GT-missing in the Original sample
is_gt_missing <- function(GT) is.na(GT) || GT %chin% c("./.", ".|.")

# ── Per-individual worker ──────────────────────────────────────────────────────
process_individual <- function(individual) {
  message("\n→ ", individual)

  # BAM
  bam_path <- file.path(bam_dir, bam_map[[individual]])
  if (!file.exists(bam_path)) { warning("Missing BAM for ", individual); return(invisible(NULL)) }
  if (!file.exists(paste0(bam_path, ".bai"))) indexBam(bam_path)

  # ROIs from map (Kept, Archaic New, by coverage)
  roi_df <- MergeFragmentsArchaic |>
    dplyr::filter(Filter == "Kept",
                  Individual == !!individual,
                  segment_type == "Archaic New",
                  Coverage %in% labels_target) |>
    transmute(Coverage,
              chr   = sub("^chr","", as.character(chr)),
              start = as.integer(start),
              end   = as.integer(end)) |>
    distinct()

  if (!nrow(roi_df)) { message("   No ROIs."); return(invisible(NULL)) }

  roi_by_cov <- lapply(split(roi_df, roi_df$Coverage), function(df) reduce(GRanges(df$chr, IRanges(df$start, df$end))))
  # Ensure all labels exist in the list
  for (lbl in labels_target) if (is.null(roi_by_cov[[lbl]])) roi_by_cov[[lbl]] <- GRanges()

  # Open tabix once per individual
  tfs <- open_tabix(vcf_paths); on.exit(close_tabix(tfs), add = TRUE)

  # Resolve sample names per VCF
  samp_cov <- vapply(labels_target, function(lbl) resolve_sample_name(vcf_paths[[lbl]], individual), character(1))
  samp_ori <- resolve_sample_name(vcf_paths[["Original"]], individual)
  for (lbl in labels_target) message(sprintf("   [%s] sample: %s | Original: %s", lbl, samp_cov[[lbl]], samp_ori))

	# For each label: positions in label minus Original within SAME coverage ROI (by chr:pos)
	# For each label: coverage positions are "new" if:
	# (a) Position is absent in Original, OR
	# (b) Position exists in Original but Original GT is missing (./. or .|.)
	new_list <- lapply(labels_target, function(lbl) {
		grL <- roi_by_cov[[lbl]]

		# Pull keys (fast, tabix): bi-allelic A/C/G/T only
		cov_keys <- vcf_keys_fast(tfs[[lbl]],        grL)
		ori_keys <- vcf_keys_fast(tfs[["Original"]], grL)

		# No coverage records in ROI -> nothing to do
		if (!nrow(cov_keys)) {
			return(data.table(cov = lbl, label = lbl, chr = character(), pos = integer())[0])
		}

		# Normalize to plain "1..22,X,Y" (no "chr" prefix) for set ops
		cov_pos <- unique(data.table(chr = sub("^chr","", cov_keys$chr), pos = cov_keys$pos))
		ori_pos <- unique(data.table(chr = sub("^chr","", ori_keys$chr), pos = ori_keys$pos))

		# (a) Absent in Original by position (cheap, no GT needed)
		newpos_absent <- fsetdiff(cov_pos, ori_pos)

		# (b) Present in both, but Original GT is missing
		#     -> check GT only for the intersection (targeted, not the whole ROI)
		shared <- fintersect(cov_pos, ori_pos)
		if (nrow(shared)) {
			gr_shared <- GRanges(shared$chr, IRanges(shared$pos, shared$pos))
			ori_gt_dt <- gt_for_positions(vcf_paths[["Original"]], gr_shared, samp_ori)
			# Keep only GT-missing in Original
			ori_missing_pos <- ori_gt_dt[is.na(GT) | GT %chin% c("./.", ".|."), .(chr, pos)]
		} else {
			ori_missing_pos <- shared[0]
		}

		# Final "new" positions = absent OR GT-missing in Original
		newpos <- unique(rbind(newpos_absent, ori_missing_pos, use.names = TRUE, fill = TRUE))
		if (!nrow(newpos)) {
			message("   ", lbl, ": 0 new sites (absent+GT-missing)")
			return(data.table(cov = lbl, label = lbl, chr = character(), pos = integer())[0])
		}

		# Attach ref/alt from the coverage keys (what you pile up against)
		refalt <- unique(data.table(
			chr = sub("^chr","", cov_keys$chr),
			pos = cov_keys$pos,
			ref = cov_keys$ref,
			alt = cov_keys$alt
		))
		ans <- merge(newpos, refalt, by = c("chr","pos"), all.x = TRUE)
		ans[, `:=`(cov = lbl, label = lbl)]
		ans[]
	})
	names(new_list) <- labels_target


  new_all <- rbindlist(new_list, use.names = TRUE, fill = TRUE)
  if (!nrow(new_all)) { message("   No new positions vs Original."); return(invisible(NULL)) }

  # GT per coverage at those positions
  gt_list <- lapply(labels_target, function(lbl) {
    pos <- new_all[label == lbl, .(chr, pos)]
    pos <- unique(pos)
    if (!nrow(pos)) return(data.table(chr=character(), pos=integer(), GT=character())[0])
    gr <- GRanges(pos$chr, IRanges(pos$pos, pos$pos))
    gt <- gt_for_positions(vcf_paths[[lbl]], gr, samp_cov[[lbl]])
    gt[, label := lbl]
    gt[]
  })
  gt_all <- rbindlist(gt_list, use.names = TRUE, fill = TRUE)

  new_all_gt <- merge(new_all, gt_all, by = c("chr","pos","label"), all.x = TRUE)
  new_all_gt[, Individual := individual]
  # Carry sample used for that label
  new_all_gt[, Sample := samp_cov[label]]

  # Union positions for pileup
  union_sites <- unique(new_all_gt[, .(chr, pos)])
  setorder(union_sites, chr, pos)

  # Pileup
  base_counts <- pileup_counts(bam_path, union_sites, workers = workers_pileup)

  # Join counts and compute ref/alt proportions
  x <- merge(new_all_gt, base_counts, by = c("chr","pos"), all.x = TRUE)
  for (b in c("A","C","G","T","N","dp_total")) x[is.na(get(b)), (b) := 0L]

  mat <- as.matrix(x[, .(A, C, G, T)])
  idx_ref <- match(x$ref, c("A","C","G","T"))
  idx_alt <- match(x$alt, c("A","C","G","T"))
  x[, `:=`(
    ref_count = ifelse(is.na(idx_ref), 0L, mat[cbind(seq_len(.N), idx_ref)]),
    alt_count = ifelse(is.na(idx_alt), 0L, mat[cbind(seq_len(.N), idx_alt)])
  )]
  x[, `:=`(
    prop_ref = fifelse(dp_total > 0, ref_count / dp_total, NA_real_),
    prop_alt = fifelse(dp_total > 0, alt_count / dp_total, NA_real_)
  )]

  setcolorder(x, c("Individual","Sample","label","cov","chr","pos","ref","alt","GT",
                   "A","C","G","T","N","dp_total","ref_count","alt_count","prop_ref","prop_alt"))
  setorder(x, label, chr, pos)

  # Write
  of <- file.path(out_dir_base, sprintf("merged_counts_gt_%s_AllCoverages.tsv.gz", individual))
  fwrite(x, of, sep = "\t")
  message("   wrote: ", of)

  invisible(TRUE)
}

# ── Run (outer parallel over individuals) ──────────────────────────────────────
BPOUT <- if (.Platform$OS.type == "unix") MulticoreParam(workers = workers_indiv, progressbar = TRUE) else SnowParam(workers = workers_indiv, type = "SOCK", progressbar = TRUE)

BiocParallel::bplapply(inds_with_bam, function(ind) {
  tryCatch(process_individual(ind),
           error = function(e) { message("Error in ", ind, ": ", conditionMessage(e)) })
}, BPPARAM = BPOUT)

message("\nDone.")
