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

# 1. Read in filtered feature bc matrix files
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

# 2. Identify the type of file, and look inside
class(Female1)
class(Male9)
Female1[c("Gapdh", "Jag1", "Notch2", "Wnt7b"), 1:50]


# 3. Create Seurat Objects from each matrix file
Female1 <- CreateSeuratObject(Female1, min.cells = 3, min.features = 300, project = "F1")
Male9 <- CreateSeuratObject(Male9, min.cells = 3, min.features = 300, project = "M9")
head(Female1)
head(Male9)


# Section 2. Merge separate files

Early_Merge <- merge(x = Female1, y = Male9, add.cell.ids = c("F1","M9"))
tail(Early_Merge@meta.data)

replicate <- substring(rownames(Early_Merge@meta.data), 1,2)
names(replicate) <- rownames(Early_Merge@meta.data)
Early_Merge <- AddMetaData(Early_Merge, metadata = replicate, col.name = "replicate") 
head(Early_Merge@meta.data)


sex <- substring(rownames(Early_Merge@meta.data), 1,1)
names(sex) <- rownames(Early_Merge@meta.data)
Early_Merge <- AddMetaData(Early_Merge, metadata = sex, col.name = "sex") 
head(Early_Merge@meta.data)
tail(Early_Merge@meta.data)



# Section 3. Quality Control Metrics to remove low quality cells
Early_Merge [["percent.mt"]] <- PercentageFeatureSet(Early_Merge , pattern = "^mt-")
#visualize spread of quality metrix to determine cutoffs
VlnPlot(Early_Merge, features = c("percent.mt", "nFeature_RNA","nCount_RNA"), pt.size = 0.2)


Early_Merge  <- subset(Early_Merge , subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 35)
VlnPlot(Early_Merge, features = c("percent.mt", "nFeature_RNA","nCount_RNA"), pt.size = 0.2)


Early_Merge  <-subset(Early_Merge , subset = nCount_RNA > 500 & nCount_RNA < 20000)
VlnPlot(Early_Merge, features = c("percent.mt", "nFeature_RNA","nCount_RNA"), pt.size = 0.2)

# Section 4. Trimming data annotations

head(Early_Merge@meta.data)
columns.to.remove <- c("replicate")

for(i in columns.to.remove) {
  Early_Merge[[i]] <- NULL
}


head(Early_Merge@meta.data)



# From here is downsampling
