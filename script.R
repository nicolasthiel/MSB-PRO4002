#-----------------------------------------------------------------------------#
# load and pre-process data
#-----------------------------------------------------------------------------#

# load data
data_cpm <- read.delim("data/MAGNET_GeneExpressionData_CPM_19112020.txt", sep = "\t", row.names = 1)
data_exon_lengths <- read.delim("data/MAGNET_exonLengths.txt", as.is = TRUE, row.names = 1)
data_sample_info <- read.csv("data/MAGNET_SampleData_18112022.csv", row.names = 1)

# convert binary columns to logical
binary_columns <- c("afib", "VTVF", "Diabetes", "Hypertension")
for (col in binary_columns) {
    data_sample_info[[col]] <- data_sample_info[[col]] == "Yes"
}

# filter for samples with DCM or HCM
filtered_sample_info <- data_sample_info[data_sample_info$etiology %in% c("DCM", "HCM"), ]
filtered_cpm <- data_cpm[, rownames(filtered_sample_info)]
