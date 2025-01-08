required_packages <- c("BiocManager")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

data_cpm <- read.delim("data/MAGNET_GeneExpressionData_CPM_19112020.txt", sep = "\t", row.names = 1)
data_exon_lengths <- read.delim("data/MAGNET_exonLengths.txt", as.is = TRUE, row.names = 1)
data_sample_info <- read.csv("data/MAGNET_SampleData_18112022.csv", row.names = 1)

#Hello
