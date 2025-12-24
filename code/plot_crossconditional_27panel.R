#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  library(ggplot2)
}))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9) {
  cat("usage: plot_crossconditional_27panel.R <runs.tsv> <signals.tsv> <ld_erap2.ld.gz> <ld_erap1.ld.gz> <ld_lnpep.ld.gz> <winkb> <out_pdf> <out_png> <out_table>\n", file=stderr())
  quit(status=2)
}

runs_path    <- args[1]
signals_path <- args[2]
ld_erap2     <- args[3]
ld_erap1     <- args[4]
ld_lnpep     <- args[5]
winkb        <- as.numeric(args[6])
out_pdf      <- args[7]
out_png      <- args[8]
out_table    <- args[9]

read_tsv <- function(p){
  read.table(p, header=TRUE, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
}

read_assoc <- function(p){
  x <- read.table(p, header=TRUE, stringsAsFactors=FALSE)
  x <- x[x$TEST=="ADD", c("CHR","SNP","BP","BETA","STAT","P")]
  x$P <- suppressWarnings(as.numeric(x$P))
  x$neglog10p <- ifelse(is.na(x$P) | x$P<=0, NA, -log10(x$P))
  x
}

read_ld <- function(p){
  con <- gzfile(p, "rt")
  ld <- read.table(con, header=TRUE, stringsAsFactors=FALSE)
  close(con)
  ld <- ld[, c("SNP_B","R2")]
  colnames(ld) <- c("SNP","r2")
  aggregate(r2 ~ SNP, data=ld, FUN=max)
}

runs <- read_tsv(runs_path)
signals <- read_tsv(signals_path)

lead_of <- function(g){
  s <- signals[signals$gene==g & signals$signal_id=="signal1", ]
  if (nrow(s)==0) stop(paste("missing signal1 for", g))
  s$lead[1]
}
lead_ERAP2 <- lead_of("ERAP2")
lead_ERAP1 <- lead_of("ERAP1")
lead_LNPEP <- lead_of("LNPEP")

ld_map <- list(
  ERAP2 = read_ld(ld_erap2),
  ERAP1 = read_ld(ld_erap1),
  LNPEP = read_ld(ld_lnpep)
)

assoc_path <- function(outcome, run_id){
  r <- runs[runs$outcome==outcome & runs$run_id==run_id, ]
  if (nrow(r)==0) return(NA)
  r$assoc_path[1]
}

panel_data_for_outcome <- function(outcome){
  p_base <- assoc_path(outcome, "baseline")
  if (is.na(p_base)) stop(paste("baseline not found for", outcome))
  base <- read_assoc(p_base)
  base$state <- "baseline"

  cols <- data.frame(
    cond_gene=c("ERAP2","ERAP1","LNPEP"),
    ref_snp=c(lead_ERAP2, lead_ERAP1, lead_LNPEP),
    stringsAsFactors=FALSE
  )

  out_all <- NULL
  for (i in 1:nrow(cols)){
    cg <- cols$cond_gene[i]
    ref <- cols$ref_snp[i]

    run_snp <- paste0("cond_snp_", ref)
    p_snp <- assoc_path(outcome, run_snp)
    if (!is.na(p_snp)) {
      snp <- read_assoc(p_snp)
      snp$state <- paste0("cond on ", cg, " lead")
    } else {
      snp <- base; snp$neglog10p <- NA
      snp$state <- paste0("cond on ", cg, " lead")
    }

    run_expr <- paste0("cov_expr_", cg)
    p_expr <- assoc_path(outcome, run_expr)
    if (outcome == cg) {
      expr <- base; expr$neglog10p <- NA
      expr$state <- "cond on self expression (NA)"
    } else if (!is.na(p_expr)) {
      expr <- read_assoc(p_expr)
      expr$state <- paste0("cond on ", cg, " expression")
    } else {
      expr <- base; expr$neglog10p <- NA
      expr$state <- paste0("cond on ", cg, " expression")
    }

    b <- base; b$state <- "baseline"
    dd <- rbind(b, snp, expr)
    dd$cond_gene <- cg
    dd$ref_snp <- ref

    dd <- merge(dd, ld_map[[cg]], by="SNP", all.x=TRUE)
    dd$r2[is.na(dd$r2)] <- 0
    out_all <- rbind(out_all, dd)
  }

  ord <- c("baseline",
           paste0("cond on ", c("ERAP2","ERAP1","LNPEP"), " lead"),
           paste0("cond on ", c("ERAP2","ERAP1","LNPEP"), " expression"),
           "cond on self expression (NA)")
  out_all$state <- factor(as.character(out_all$state), levels=ord)
  out_all
}

d1 <- panel_data_for_outcome("ERAP2"); d1$outcome <- "ERAP2"
d2 <- panel_data_for_outcome("ERAP1"); d2$outcome <- "ERAP1"
d3 <- panel_data_for_outcome("LNPEP"); d3$outcome <- "LNPEP"
D <- rbind(d1,d2,d3)

get_ref_bp <- function(ref){
  rr <- D[D$SNP==ref, c("cond_gene","BP")]
  rr <- unique(rr)
  if (nrow(rr)==0) return(data.frame(cond_gene=NA, BP=NA))
  rr
}
vlines <- rbind(
  transform(get_ref_bp(lead_ERAP2), cond_gene="ERAP2"),
  transform(get_ref_bp(lead_ERAP1), cond_gene="ERAP1"),
  transform(get_ref_bp(lead_LNPEP), cond_gene="LNPEP")
)
vlines <- unique(vlines[!is.na(vlines$BP), c("cond_gene","BP")])

sigline <- -log10(5e-8)

p <- ggplot(D, aes(x=BP/1e6, y=neglog10p, color=r2)) +
  geom_point(size=0.55, alpha=0.9, na.rm=TRUE) +
  geom_hline(yintercept=sigline, linetype="dashed", linewidth=0.3) +
  geom_vline(data=vlines, aes(xintercept=BP/1e6), linetype="dotted", linewidth=0.35, inherit.aes=FALSE) +
  facet_grid(outcome + state ~ cond_gene, scales="free_y") +
  labs(x="Position (Mb, chr5)", y=expression(-log[10](P)), color=expression(r^2)) +
  theme_bw(base_size=9) +
  theme(
    legend.position="right",
    strip.text.y = element_text(size=7),
    strip.text.x = element_text(size=8),
    panel.grid.minor = element_blank()
  )

ggsave(out_pdf, p, width=9.5, height=12.5, units="in", dpi=300)
ggsave(out_png, p, width=9.5, height=12.5, units="in", dpi=300)

att_path <- sub("runs.tsv$", "attenuation_summary.tsv", runs_path)
att <- read_tsv(att_path)

att_s1 <- att[
  (att$outcome=="ERAP2" & att$signal_id=="signal1") |
  (att$outcome=="ERAP1" & att$signal_id=="signal1") |
  (att$outcome=="LNPEP" & att$signal_id=="signal1"), ]

keep_run <- function(x){ x=="baseline" || grepl("^cond_snp_", x) || grepl("^cov_expr_", x) }
att_s1 <- att_s1[ keep_run(att_s1$run_id), ]

att_s1$conditioning_gene <- NA
att_s1$conditioning_type <- NA
att_s1$conditioning_item <- NA

for (i in 1:nrow(att_s1)){
  rid <- att_s1$run_id[i]
  if (rid=="baseline"){
    att_s1$conditioning_gene[i] <- "baseline"
    att_s1$conditioning_type[i] <- "baseline"
    att_s1$conditioning_item[i] <- "baseline"
  } else if (grepl("^cond_snp_", rid)){
    snp <- sub("^cond_snp_", "", rid)
    g <- if (snp==lead_ERAP2) "ERAP2" else if (snp==lead_ERAP1) "ERAP1" else if (snp==lead_LNPEP) "LNPEP" else "other"
    att_s1$conditioning_gene[i] <- g
    att_s1$conditioning_type[i] <- "SNP"
    att_s1$conditioning_item[i] <- snp
  } else if (grepl("^cov_expr_", rid)){
    g <- sub("^cov_expr_", "", rid)
    att_s1$conditioning_gene[i] <- g
    att_s1$conditioning_type[i] <- "expression"
    att_s1$conditioning_item[i] <- paste0("expr_", g)
  }
}

tab <- att_s1[att_s1$run_id!="baseline", c(
  "outcome","lead_snp","conditioning_gene","conditioning_type","conditioning_item",
  "beta_base","p_base","neglog10p_base",
  "beta_cond","p_cond","neglog10p_cond",
  "delta_beta_pct","delta_neglog10p_pct"
)]
write.table(tab, out_table, sep="\t", quote=FALSE, row.names=FALSE)
