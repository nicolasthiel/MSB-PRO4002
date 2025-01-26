library(DESeq2)
library(clusterProfiler)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(edgeR)
library(tidyr)
library(enrichplot)
library(gt)
library(gtsummary)


# set theme
theme_set(theme_minimal(base_size = 15))

# set etiology
etiology <- "HCM"

# set output path
path_out <- paste("out/", etiology, "/", sep = "")

# -----------------------------------------------------------------------------#
# pre-processing
# -----------------------------------------------------------------------------#

# read in the tabular file
data_expr <- read.delim(paste("data/", etiology, "/counts.tsv", sep = ""), sep = "\t", row.names = 1)

data_sample_info <- read.delim("data/sampleData_new.txt", row.names = 1)
data_sample_info2 <- read.delim("data/MAGNET_SampleData_18112022.csv", sep = ",", row.names = 1)

data_sample_info$RIN <- data_sample_info2$RIN

data_sample_info <- data_sample_info[rownames(data_sample_info) %in% colnames(data_expr), ]

# filter out samples with missing afib
data_sample_info <- data_sample_info[!is.na(data_sample_info$afib), ]
data_sample_info <- data_sample_info[!is.na(data_sample_info$RIN), ]


data_sample_info <- data_sample_info[rownames(data_sample_info) %in% colnames(data_expr), ]
data_counts <- data_expr[, rownames(data_sample_info)]

# convert variables to factors
factor_vars <- c("Library.Pool", "gender", "Hypertension", "VTVF", "race", "afib", "Diabetes", "disease_race")
data_sample_info[factor_vars] <- lapply(data_sample_info[factor_vars], as.factor)


# -----------------------------------------------------------------------------#
# Data summary
# -----------------------------------------------------------------------------#

# create table one comparing across the 4 etiologies
table_two <- data_sample_info %>%
    select(
        # demographic
        age,
        gender,
        race,
        weight,
        height,
        # clinical
        afib,
        LVEF,
        Diabetes,
        Hypertension
    ) %>%
    tbl_summary(
        by = afib,
        missing = "no",
        statistic = list(
            all_continuous() ~ "{mean} ± {sd}"
        )
    ) %>%
    add_p() %>%
    add_overall() %>%
    modify_header(label = "**Variable**") %>%
    modify_caption("**Table 3. HCM Patients - afib vs no afib**") %>%
    bold_labels()

# save as HTML
as_gt(table_two) %>%
    gtsave(paste(path_out, "table_two.html", sep = ""))


# -----------------------------------------------------------------------------#
# Exploratory data analysis
# -----------------------------------------------------------------------------#

# Calculate log CPM (Counts Per Million)
data_cpm <- cpm(data_counts, log = FALSE)
data_cpm <- log2(data_cpm + 1)
data_cpm <- as.data.frame(data_cpm)
data_cpm_gathered <- gather(data_cpm, key = "SampleID", value = "CPM")
data_cpm_gathered$afib <- data_sample_info$afib[match(data_cpm_gathered$SampleID, colnames(data_cpm))]


# -----------------------------------------------------------------------------#
# DESeq2 analysis
# -----------------------------------------------------------------------------#

# making sure the colnames of data_counts match the rownames of data_sample_info
all(colnames(data_counts) == rownames(data_sample_info))

# scale and center age
data_sample_info$age <- scale(data_sample_info$age, center = TRUE, scale = TRUE)

# identify covariates to afib
covariates <- colnames(data_sample_info)[colnames(data_sample_info) != "afib"]
p_values <- sapply(covariates, function(covariate) {
    if (is.numeric(data_sample_info[[covariate]])) {
        # Perform t-test for continuous variables
        t_test <- t.test(data_sample_info[[covariate]] ~ data_sample_info$afib)
        return(t_test$p.value)
    } else {
        # Perform chi-squared test for categorical variables
        chi_sq_test <- chisq.test(table(data_sample_info[[covariate]], data_sample_info$afib))
        return(chi_sq_test$p.value)
    }
})
significant_covariates <- names(p_values)[p_values < 0.05]
print(significant_covariates)

# construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
    countData = data_counts,
    colData = data_sample_info,
    design = ~ Library.Pool + RIN + age + race + gender + afib
)

# set the factor level
dds$afib <- relevel(dds$afib, ref = "No")

# run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, pAdjustMethod = "none")
summary(res)

df_deg <- as.data.frame(res)
write.csv(df_deg, paste(path_out, "DEG_results.csv", sep = ""), row.names = TRUE)


# -----------------------------------------------------------------------------#
# annotation of results
# -----------------------------------------------------------------------------#

df_deg$diffexpressed <- "NO"
df_deg$diffexpressed[df_deg$pvalue < 0.05 & df_deg$log2FoldChange > 0] <- "UP"
df_deg$diffexpressed[df_deg$pvalue < 0.05 & df_deg$log2FoldChange < -0] <- "DOWN"

# get gene annotations
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
    filters = "entrezgene_id",
    values = rownames(data_counts),
    mart = ensembl
)
gene_annotations <- gene_annotations[!duplicated(gene_annotations$entrezgene_id), ]

df_deg <- merge(df_deg, gene_annotations, by.x = "row.names", by.y = "entrezgene_id", all.x = TRUE)
rownames(df_deg) <- df_deg$Row.names
df_deg$Row.names <- NULL

top20_degs <- head(df_deg[order(df_deg$pvalue), "hgnc_symbol"], 20)
df_deg$delabel <- ifelse(df_deg$hgnc_symbol %in% top20_degs, df_deg$hgnc_symbol, NA)


# -----------------------------------------------------------------------------#
# Plots
# -----------------------------------------------------------------------------#

# volcano plot
plot_volcano <- ggplot(data = df_deg, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = 0, col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
    geom_point(alpha = 0.5) +
    scale_color_manual(
        values = c("#00AFBB", "grey", "#d04b39"),
        labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    labs(
        color = "diffexpressed", "Volcano plot",
        x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value")
    ) +
    ggtitle(paste("Differentially Expressed Genes in ", etiology, " Patients with afib vs without afib")) +
    geom_text_repel(aes(label = ifelse(is.na(delabel), "", delabel)), max.overlaps = Inf)
ggsave(path = path_out, "plot_volcano.jpg", plot = plot_volcano, width = 10, height = 8, dpi = 500)

jpeg(file = paste(path_out, "plot_MA.jpg", sep = ""), width = 800, height = 600, quality = 100)
plot_MA <- DESeq2::plotMA(dds, ylim = c(-1, 1), colSig = "red")
dev.off()


# -----------------------------------------------------------------------------#
# Enrichment analysis
# -----------------------------------------------------------------------------#

logFC <- df_deg$log2FoldChange
names(logFC) <- rownames(df_deg)
logFC <- sort(logFC, decreasing = TRUE)

### GO gene set enrichment analysis ###
gsea_go <- gseGO(
    geneList = logFC,
    OrgDb = "org.Hs.eg.db"
)
df_gsea_go <- as.data.frame(gsea_go)
write.csv(df_gsea_go, paste(path_out, "gsea_go_results.csv", sep = ""), row.names = FALSE)

# dot plot
plot_dot_gsea_go <- enrichplot::dotplot(
    gsea_go,
    showCategory = 20,
    title = paste("GO terms enriched in ", etiology, " ordered by gene ratio", sep = "")
)
ggsave(paste(path_out, "plot_dot_gsea_go.jpg", sep = ""), plot = plot_dot_gsea_go, width = 10, height = 14, dpi = 500)

# tree plot
sim_matrix <- pairwise_termsim(gsea_go)
plot_tree_gsea_go <- enrichplot::treeplot(
    sim_matrix,
    title = paste("GO terms enriched in ", etiology, " ordered by gene ratio", sep = "")
)
ggsave(paste(path_out, "plot_tree_gsea_go.jpg", sep = ""), plot = plot_tree_gsea_go, width = 23, height = 10, dpi = 500)


### WP gene set enrichment analysis ###
gsea_wp <- gseWP(
    geneList = logFC,
    organism = "Homo sapiens"
)
df_gsea_wp <- as.data.frame(gsea_wp)
write.csv(df_gsea_wp, paste(path_out, "gsea_wp_results.csv", sep = ""), row.names = FALSE)

# dot plot
plot_dot_gsea_wp <- enrichplot::dotplot(
    gsea_wp,
    showCategory = 20,
    title = paste("WP pathways enriched in ", etiology, " ordered by gene ratio", sep = "")
)
ggsave(paste(path_out, "plot_dot_gsea_wp.jpg", sep = ""), plot = plot_dot_gsea_wp, width = 10, height = 10, dpi = 500)

# cnet plot
foldChange <- df_deg$log2FoldChange[df_deg$diffexpressed != "NO"]
names(foldChange) <- df_deg$hgnc_symbol[df_deg$diffexpressed != "NO"]
gsea_wp <- setReadable(gsea_wp, "org.Hs.eg.db", "ENTREZID")
plot_cnet_gsea_wp <- enrichplot::cnetplot(
    gsea_wp,
    foldChange = foldChange,
    layout = "kk",
    colorEdge = TRUE,
    title = paste("Top 5 WP pathways enriched in ", etiology, sep = "")
)
ggsave(paste(path_out, "plot_cnet_gsea_wp.jpg", sep = ""), plot = plot_cnet_gsea_wp, width = 18, height = 15, dpi = 500)

# tree plot
sim_matrix <- pairwise_termsim(gsea_wp)
plot_tree_gsea_wp <- enrichplot::treeplot(sim_matrix)
ggsave(paste(path_out, "plot_tree_gsea_wp.jpg", sep = ""), plot = plot_tree_gsea_wp, width = 18, height = 10, dpi = 500)


### KEGG gene set enrichment analysis ###
gsea_kegg <- gseKEGG(
    geneList = logFC,
    organism = "hsa"
)
df_gsea_kegg <- as.data.frame(gsea_kegg)
write.csv(df_gsea_kegg, paste(path_out, "gsea_kegg_results.csv", sep = ""), row.names = FALSE)

# dot plot
plot_dot_gsea_kegg <- enrichplot::dotplot(gsea_kegg, showCategory = 20)
ggsave(paste(path_out, "plot_dot_gsea_kegg.jpg", sep = ""), plot = plot_dot_gsea_kegg, width = 10, height = 10, dpi = 500)

# cnet plot
foldChange <- df_deg$log2FoldChange
names(foldChange) <- df_deg$hgnc_symbol
gsea_kegg <- setReadable(gsea_kegg, "org.Hs.eg.db", "ENTREZID")
plot_cnet_gsea_kegg <- enrichplot::cnetplot(gsea_kegg, foldChange = foldChange, layout = "kk", colorEdge = TRUE)
ggsave(paste(path_out, "plot_cnet_gsea_kegg.jpg", sep = ""), plot = plot_cnet_gsea_kegg, width = 18, height = 15, dpi = 500)
