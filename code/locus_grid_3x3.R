#!/usr/bin/env Rscript
# locus_grid_3x3.R
# 3x3 grid (rows=baseline gene, cols=conditioning gene)
# each cell: baseline / cond on top SNP / cond on expression
#
# Usage:
#   Rscript locus_grid_3x3.R <LDGZ> <SNPLIST> <ERAP2_LEAD> <ERAP1_TOP> <LNPEP_TOP> <OUTPNG> [--no_vline_expr] [--width 18] [--height 26] [--dpi 300]
#
# Assumed assoc filenames in cwd:
#   ERAP2.base_pm500k.assoc.linear
#   ERAP1.base.assoc.linear
#   LNPEP.base.assoc.linear
#   <OUTCOME>_cond_<SNP>.assoc.linear
#   <OUTCOME>_cond_on_<GENE>expr.assoc.linear   (e.g. ERAP2_cond_on_ERAP1expr.assoc.linear)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

parse_args <- function(args) {
  if (length(args) < 6) {
    stop("Need at least 6 args: LDGZ SNPLIST ERAP2_LEAD ERAP1_TOP LNPEP_TOP OUTPNG", call.=FALSE)
  }
  opt <- list(
    ld_gz = args[1],
    snplist = args[2],
    erap2 = args[3],
    erap1 = args[4],
    lnpep = args[5],
    outpng = args[6],
    no_vline_expr = FALSE,
    width = 18,
    height = 26,
    dpi = 300
  )
  if (length(args) > 6) {
    i <- 7
    while (i <= length(args)) {
      a <- args[i]
      if (a == "--no_vline_expr") {
        opt$no_vline_expr <- TRUE
        i <- i + 1
      } else if (a == "--width") {
        opt$width <- as.numeric(args[i+1]); i <- i + 2
      } else if (a == "--height") {
        opt$height <- as.numeric(args[i+1]); i <- i + 2
      } else if (a == "--dpi") {
        opt$dpi <- as.integer(args[i+1]); i <- i + 2
      } else {
        stop(paste0("Unknown option: ", a), call.=FALSE)
      }
    }
  }
  opt
}

read_snplist <- function(path) {
  snps <- scan(path, what=character(), quiet=TRUE)
  snps <- snps[!is.na(snps) & snps != ""]
  if (length(snps) >= 1 && tolower(snps[1]) %in% c("snp","rsid","id")) snps <- snps[-1]
  if (length(snps) < 2) stop("SNPLIST seems empty.", call.=FALSE)
  snps
}

read_ld_matrix <- function(ld_gz) {
  ld <- fread(ld_gz, header=FALSE, data.table=FALSE, showProgress=FALSE)
  as.matrix(ld)
}

read_plink_linear_add <- function(path) {
  dt <- fread(path, data.table=FALSE, showProgress=FALSE)
  if (!all(c("SNP","BP","TEST","P") %in% colnames(dt))) {
    stop(paste0("Unexpected columns in: ", path), call.=FALSE)
  }
  dt <- dt[dt$TEST == "ADD", , drop=FALSE]
  if (!("BETA" %in% colnames(dt))) dt$BETA <- NA_real_
  dt$BP <- as.numeric(dt$BP)
  dt$P <- suppressWarnings(as.numeric(dt$P))
  dt$BETA <- suppressWarnings(as.numeric(dt$BETA))
  dt
}

align_to_snplist <- function(dt, snps) {
  idx <- match(snps, dt$SNP)
  out <- dt[idx, , drop=FALSE]
  out$SNP <- snps
  out
}

safe_assoc <- function(path, base_aligned) {
  if (!file.exists(path)) {
    x <- base_aligned
    x$P <- NA_real_
    x$BETA <- NA_real_
    return(x)
  }
  dt <- read_plink_linear_add(path)
  align_to_snplist(dt, base_aligned$SNP)
}

make_cell_plot <- function(outcome, covar, base_aligned, assoc_cond_snp, assoc_cond_expr,
                           r2_vec, ref_snp, ref_bp, no_vline_expr, y_max_outcome) {

  base <- base_aligned
  base$STATE_ID <- "baseline"
  cond1 <- assoc_cond_snp
  cond1$STATE_ID <- "cond_snp"
  cond2 <- assoc_cond_expr
  cond2$STATE_ID <- "cond_expr"

  df <- rbind(base, cond1, cond2)
  df$MB <- df$BP / 1e6
  df$LOGP <- suppressWarnings(-log10(df$P))
  df$r2 <- rep(r2_vec, times=3)

  # panel labels
  state_levels <- c("baseline","cond_snp","cond_expr")
  state_labels <- c(
    "baseline",
    sprintf("cond %s top\n(%s)", covar, ref_snp),
    sprintf("cond on %s expression", covar)
  )
  df$STATE_ID <- factor(df$STATE_ID, levels=state_levels, labels=state_labels)

  # vline per panel
  v_states <- c("baseline", sprintf("cond %s top\n(%s)", covar, ref_snp))
  if (!no_vline_expr) v_states <- c(v_states, sprintf("cond on %s expression", covar))
  vline_df <- data.frame(
    STATE_ID=factor(v_states, levels=levels(df$STATE_ID)),
    MB=ref_bp/1e6
  )

  # label only in baseline panel
  label_df <- data.frame(
    STATE_ID=factor("baseline", levels=levels(df$STATE_ID)),
    MB=ref_bp/1e6,
    y=y_max_outcome * 0.98,
    lab=ref_snp
  )

  # NA annotation for self-expression covariate (diagonal)
  # -> show "NA" text in the expression facet even if all points removed
  expr_label <- sprintf("cond on %s expression", covar)
  xmid <- mean(range(df$MB, na.rm=TRUE))
  if (!is.finite(xmid) || is.na(xmid)) xmid <- ref_bp/1e6
  na_df <- data.frame(
    STATE_ID=factor(expr_label, levels=levels(df$STATE_ID)),
    MB=xmid,
    y=y_max_outcome * 0.55,
    lab=if (outcome == covar) "NA\n(self expression covariate)" else NA_character_
  )
  na_df <- na_df[!is.na(na_df$lab), , drop=FALSE]

  p <- ggplot(df, aes(x=MB, y=LOGP, color=r2)) +
    geom_hline(yintercept=-log10(1e-6), linetype="dashed", linewidth=0.3, alpha=0.6) +
    geom_hline(yintercept=-log10(1e-4), linetype="dashed", linewidth=0.3, alpha=0.4) +
    geom_vline(data=vline_df, aes(xintercept=MB), inherit.aes=FALSE,
               linetype="dashed", linewidth=0.4, alpha=0.9) +
    geom_point(size=0.75, alpha=0.95, na.rm=TRUE) +
    geom_text(data=label_df, aes(x=MB, y=y, label=lab), inherit.aes=FALSE,
              size=2.6, vjust=0, fontface="bold") +
    geom_text(data=na_df, aes(x=MB, y=y, label=lab), inherit.aes=FALSE,
              size=3.0, fontface="bold", lineheight=0.95) +
    facet_wrap(~STATE_ID, ncol=1, scales="fixed") +
    scale_color_gradientn(
      colours=c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
      limits=c(0,1), oob=squish, name=expression(r^2)
    ) +
    coord_cartesian(ylim=c(0, y_max_outcome)) +
    labs(
      title=sprintf("%s baseline | cond %s", outcome, covar),
      x="Position (Mb, chr5)",
      y=expression(-log[10]~P)
    ) +
    theme_bw(base_size=9) +
    theme(
      plot.title=element_text(size=10, face="bold", hjust=0.5),
      strip.background=element_rect(fill="grey95"),
      strip.text=element_text(size=8, face="bold"),
      legend.key.height=unit(0.45, "cm"),
      legend.key.width=unit(0.45, "cm"),
      panel.grid.minor=element_blank()
    )

  p
}

args <- commandArgs(trailingOnly=TRUE)
opt <- parse_args(args)

genes <- c("ERAP2","ERAP1","LNPEP")
ref_snp_by_gene <- c(ERAP2=opt$erap2, ERAP1=opt$erap1, LNPEP=opt$lnpep)

base_file <- c(
  ERAP2="ERAP2.base_pm500k.assoc.linear",
  ERAP1="ERAP1.base.assoc.linear",
  LNPEP="LNPEP.base.assoc.linear"
)

snps <- read_snplist(opt$snplist)
ld <- read_ld_matrix(opt$ld_gz)

if (nrow(ld) != length(snps) || ncol(ld) != length(snps)) {
  stop(sprintf("LD matrix dim (%d x %d) != length(SNPLIST) (%d).", nrow(ld), ncol(ld), length(snps)), call.=FALSE)
}

# r^2 vectors per column gene (by column ref SNP)
r2_by_col <- list()
ref_bp_by_col <- c()
for (g in genes) {
  rs <- ref_snp_by_gene[[g]]
  idx <- match(rs, snps)
  if (is.na(idx)) stop(sprintf("Ref SNP %s (%s) not found in SNPLIST.", rs, g), call.=FALSE)
  r2_by_col[[g]] <- pmin(1, pmax(0, (ld[, idx])^2))
  ref_bp_by_col[[g]] <- NA_real_
}

# read & align baselines first (also gives BP)
base_aligned <- list()
for (out in genes) {
  if (!file.exists(base_file[[out]])) stop(sprintf("Missing baseline file: %s", base_file[[out]]), call.=FALSE)
  dtb <- read_plink_linear_add(base_file[[out]])
  base_aligned[[out]] <- align_to_snplist(dtb, snps)
}

# fill ref BP using any baseline that contains it
for (g in genes) {
  rs <- ref_snp_by_gene[[g]]
  bp <- NA_real_
  for (out in genes) {
    hit <- base_aligned[[out]][base_aligned[[out]]$SNP == rs, "BP"]
    if (length(hit) == 1 && !is.na(hit)) { bp <- hit; break }
  }
  if (is.na(bp)) stop(sprintf("Could not find BP for ref SNP %s in any baseline file.", rs), call.=FALSE)
  ref_bp_by_col[[g]] <- bp
}

# y max per outcome gene
ymax_out <- c()
for (out in genes) {
  yvals <- c(-log10(base_aligned[[out]]$P))
  for (cov in genes) {
    cond_snp_path <- sprintf("%s_cond_%s.assoc.linear", out, ref_snp_by_gene[[cov]])
    cond_expr_path <- sprintf("%s_cond_on_%sexpr.assoc.linear", out, cov)
    d1 <- safe_assoc(cond_snp_path, base_aligned[[out]])
    d2 <- safe_assoc(cond_expr_path, base_aligned[[out]])
    yvals <- c(yvals, -log10(d1$P), -log10(d2$P))
  }
  ymax <- suppressWarnings(max(yvals, na.rm=TRUE))
  if (!is.finite(ymax) || is.na(ymax)) ymax <- 1
  ymax_out[[out]] <- ymax * 1.05
}

# build 9 cell plots
cell <- list()
for (out in genes) {
  for (cov in genes) {
    cond_snp_path <- sprintf("%s_cond_%s.assoc.linear", out, ref_snp_by_gene[[cov]])
    cond_expr_path <- sprintf("%s_cond_on_%sexpr.assoc.linear", out, cov)

    d1 <- safe_assoc(cond_snp_path, base_aligned[[out]])
    d2 <- safe_assoc(cond_expr_path, base_aligned[[out]])

    p <- make_cell_plot(
      outcome=out, covar=cov,
      base_aligned=base_aligned[[out]],
      assoc_cond_snp=d1,
      assoc_cond_expr=d2,
      r2_vec=r2_by_col[[cov]],
      ref_snp=ref_snp_by_gene[[cov]],
      ref_bp=ref_bp_by_col[[cov]],
      no_vline_expr=opt$no_vline_expr,
      y_max_outcome=ymax_out[[out]]
    )

    if (cov != "ERAP2") p <- p + theme(axis.title.y=element_blank())
    if (out != "LNPEP") p <- p + theme(axis.title.x=element_blank())

    cell[[paste(out, cov, sep="_")]] <- p
  }
}

pgrid <- (cell[["ERAP2_ERAP2"]] | cell[["ERAP2_ERAP1"]] | cell[["ERAP2_LNPEP"]]) /
         (cell[["ERAP1_ERAP2"]] | cell[["ERAP1_ERAP1"]] | cell[["ERAP1_LNPEP"]]) /
         (cell[["LNPEP_ERAP2"]] | cell[["LNPEP_ERAP1"]] | cell[["LNPEP_LNPEP"]])

topline <- sprintf("Colours: r^2 to column ref SNP | ERAP2=%s | ERAP1=%s | LNPEP=%s",
                   opt$erap2, opt$erap1, opt$lnpep)

pgrid <- pgrid +
  plot_layout(guides="collect") &
  theme(legend.position="right")

pgrid <- pgrid + plot_annotation(
  title=topline,
  theme=theme(plot.title=element_text(size=10, face="bold", hjust=0.5))
)

ggsave(opt$outpng, pgrid, width=opt$width, height=opt$height, dpi=opt$dpi, bg="white")
