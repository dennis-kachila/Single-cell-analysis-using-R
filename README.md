# Single-Cell RNA Sequencing Analysis with Seurat

This project processes single-cell RNA sequencing data using the Seurat package in R. The analysis includes data preprocessing, quality control, and merging of datasets from two samples.

## Requirements

- R (version 4.0 or higher)
- Seurat package

Install the Seurat package if not already installed:

```r
install.packages("Seurat")



Collecting workspace informationHere is a suggested README.md content based on the script:

```md
# Single-Cell RNA Sequencing Analysis with Seurat

This project processes single-cell RNA sequencing data using the Seurat package in R. The analysis includes data preprocessing, quality control, and merging of datasets from two samples.

## Requirements

- R (version 4.0 or higher)
- Seurat package

Install the Seurat package if not already installed:

```r
install.packages("Seurat")
```

## Data

The analysis uses filtered feature-barcode matrices for two samples:
- **Female1**: Located in `F1/filtered_feature_bc_matrix/`
- **Male9**: Located in `M9/filtered_feature_bc_matrix/`

## Workflow

### 1. Data Loading
The script reads the filtered feature-barcode matrices for the two samples using the `Read10X` function from Seurat.

### 2. Data Inspection
- Checks the class of the loaded data objects (`Female1` and `Male9`).
- Extracts specific gene expression data for genes `Gapdh`, `Jag1`, `Notch2`, and `Wnt7b` from the first 50 cells of the `Female1` dataset.

### 3. Create Seurat Objects
- Creates Seurat objects for both datasets with a minimum of 3 cells and 300 features.
- Assigns project names "F1" and "M9" to the respective datasets.

### 4. Merging Datasets
- Merges the `Female1` and `Male9` datasets into a single Seurat object named `Early_Merge`.
- Adds metadata columns for `replicate` and `sex` to distinguish between the datasets.

### 5. Quality Control
- Calculates the percentage of mitochondrial gene expression (`percent.mt`).
- Visualizes quality metrics (`percent.mt`, `nFeature_RNA`, `nCount_RNA`) using violin plots.
- Filters out low-quality cells based on thresholds:
  - `nFeature_RNA` between 500 and 5000
  - `nCount_RNA` between 500 and 20000
  - `percent.mt` less than 35%

### 6. Data Annotation
- Removes unnecessary metadata columns (e.g., `replicate`).

## Outputs
- Merged and filtered Seurat object (`Early_Merge`) ready for downstream analysis.

## Visualization
The script includes violin plots to visualize the distribution of quality metrics before and after filtering.

## File Structure
- `F1/filtered_feature_bc_matrix/`: Contains the Female1 sample data.
- `M9/filtered_feature_bc_matrix/`: Contains the Male9 sample data.
- `assignment script.R`: The main R script for the analysis.

## Notes
- Ensure the file paths in the script match the location of your data files.
- Modify thresholds for quality control as needed based on your dataset.

## License
This project is licensed under the MIT License.