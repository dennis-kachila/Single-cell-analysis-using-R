# This script processes single-cell RNA sequencing data using the Seurat package.
# 
# Libraries:
# - Seurat: A toolkit for single-cell genomics analysis.
#
# Data:
# - Female1: Reads the filtered feature-barcode matrix for the Female1 sample 
#   from the specified directory.
# - Male9: Reads the filtered feature-barcode matrix for the Male9 sample 
#   from the specified directory.
#
# Filepaths:
# - The script specifies the file paths for the input data directories 
#   containing the filtered feature-barcode matrices for each sample.
library(Seurat)

Female1 <- Read10X(data.dir = "C:/Users/HomePC/AppData/Roaming/SPB_Data/Single Cell/Assignment/F1/filtered_feature_bc_matrix/")
Male9 <-   Read10X(data.dir = "C:/Users/HomePC/AppData/Roaming/SPB_Data/Single Cell/Assignment/M9/filtered_feature_bc_matrix/")



# This script performs the following operations:
# 1. Checks the class of the objects `Female1` and `Male9`.
# 2. Extracts specific gene expression data for genes "Gapdh", "Jag1", "Notch2", and "Wnt7b" 
#    from the first 50 cells of the `Female1` dataset.
# 3. Creates Seurat objects for `Female1` and `Male9` datasets with a minimum of 3 cells 
#    and 300 features, assigning project names "F1" and "M9" respectively.
# 4. Displays the first few rows of the newly created Seurat objects for `Female1` and `Male9`.
# 5. Merges the `Female1` and `Male9` Seurat objects into a single object named `Early_Merge`, 
#    adding cell identifiers "F1" and "M9" to distinguish between the two datasets.
class(Female1)

class(Male9)

Female1[c("Gapdh", "Jag1", "Notch2", "Wnt7b"), 1:50]


Female1 <- CreateSeuratObject(Female1, min.cells = 3, min.features = 300, project = "F1")
Male9 <- CreateSeuratObject(Male9, min.cells = 3, min.features = 300, project = "M9")

head(Female1)

head(Male9)

Early_Merge <- merge(x = Female1, y = Male9, add.cell.ids = c("F1","M9"))

