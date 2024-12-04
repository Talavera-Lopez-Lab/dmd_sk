import os
import numpy as np
import anndata
import datetime
import scanpy as sc
import pandas as pd
import sc_toolbox
import scipy.io
import matplotlib.pyplot as plt
import scipy.sparse as sparse
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import r
import logging

# Activate pandas interface for R objects
pandas2ri.activate()
MAST = importr('MAST')

# Configure logging
timestamp = datetime.datetime.now().strftime("%d_%m_%y,%H:%M")
# Set up log directory and filename
log_dir = "/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results"
log_filename = f"{log_dir}/de_analysis_{timestamp}.log"

# Ensure the log directory exists; if not, create it
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
    print(f"Created log directory: {log_dir}")

# Ensure the log file exists; if not, create it
if not os.path.exists(log_filename):
    open(log_filename, 'w').close()  # This will create an empty log file
    print(f"Created log file: {log_filename}")

logging.basicConfig(
    filename=log_filename,
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

logging.info("Starting DE analysis script.")

# Set Scanpy settings
sc.settings.verbosity = 3
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300, color_map='RdPu', dpi_save=300, vector_friendly=True, format='svg')

# Load data files into AnnData objects and concatenate
files = [
    '/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/data/dmd_annotated_human_wt_cmc_1k_hvg_25_10_24,11:40.h5ad',
    '/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/data/dmd_annotated_wt_others_cmc_5k_hvg_25_10_24,11:49.h5ad'
]
logging.info("Loading AnnData files.")
adatas = [anndata.read_h5ad(fp) for fp in files]
adata = anndata.concat(adatas, join='inner')

# Filter for desired cell states
desired_cell_state = ['vCM1', 'vCM2', 'vCM3', 'vCM4']
adata = adata[adata.obs['cell_state'].isin(desired_cell_state)]
logging.info("Filtered AnnData for desired cell states.")

# Function to create SingleCellAssay object in R
def create_single_cell_assay():
    logging.info("Creating SingleCellAssay object...")
    robjects.r('''
    library(MAST)
    library(Matrix)
    
    exprsMatrix <- readMM("/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results/de_matrix_all_geno_08_11_24,10:58.mtx")
    genes <- read.table("/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results/de_genes_all_geno_08_11_24,10:58.tsv", header=FALSE, stringsAsFactors=FALSE)
    barcodes <- read.table("/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results/de_barcodes_all_geno_08_11_24,10:58.tsv", header=FALSE, stringsAsFactors=FALSE)
    metadata <- read.table("/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results/de_metadata_all_geno_08_11_24,10:58.tsv", header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)

    colnames(exprsMatrix) <- barcodes$V1
    rownames(exprsMatrix) <- genes$V1
    exprsMatrix <- as.matrix(exprsMatrix)

    common_cells <- intersect(colnames(exprsMatrix), rownames(metadata))
    exprsMatrix <- exprsMatrix[, common_cells]
    metadata <- metadata[common_cells, ]

    sca <- FromMatrix(
        exprsArray=exprsMatrix,
        cData=metadata,
        fData=data.frame(gene_id=rownames(exprsMatrix)),
        class="SingleCellAssay",
        check_sanity=FALSE
    )

    cdr2 <- colSums(assay(sca) > 0)
    colData(sca)$n_genes_per_cell <- scale(cdr2)
    colData(sca)$genotype <- factor(colData(sca)$genotype)
    #colData(sca)$cell_state <- factor(colData(sca)$cell_state)
    
    sca
    ''')
    logging.info("SingleCellAssay object created.")
    return robjects.r['sca']

# Create SingleCellAssay object
sca = create_single_cell_assay()

# Function to perform differential expression analysis based on genotype only
# Function to perform differential expression analysis based only on genotype
def find_de_MAST(sca):
    logging.info("Finding differentially expressed genes based on genotype only...")
    robjects.r('''
        find_de <- function(sca) {
            # Fit model with genotype as the main effect, excluding cell_state from the model
            zlmCond <- zlm(~ n_genes_per_cell + genotype, sca)
            summaryCond <- summary(zlmCond, doLRT = TRUE)
            summaryDt <- summaryCond$datatable
            
            results <- list()
            # Extract unique genotype comparisons
            contrasts <- unique(summaryDt[component == 'H', .(contrast)])
            
            for (contrast in contrasts$contrast) {
                if (contrast != 'n_genes_per_cell' && grepl("genotype", contrast)) {
                    # Extract log fold change and p-value for genotype contrasts
                    contrast_lfc = summaryDt[contrast == contrast & component == 'logFC', .(primerid, coef)]
                    contrast_p = summaryDt[contrast == contrast & component == 'H', .(primerid, `Pr(>Chisq)`)]
                    
                    # Merge results and calculate log fold change in base 2
                    tmp <- merge(contrast_lfc, contrast_p, by='primerid', allow.cartesian=TRUE)
                    tmp$log_fold_change <- tmp$coef / log(2)
                    tmp$FDR <- p.adjust(tmp$`Pr(>Chisq)`, 'fdr')
                    
                    # Label genotype
                    genotype <- sub(".*genotype", "genotype", contrast)
                    tmp$genotype <- genotype
                    
                    # Rename columns
                    colnames(tmp) <- c('gene_id', 'log_fold_change', 'p_value', 'FDR', 'genotype')
                    
                    # Save result
                    results[[genotype]] <- tmp
                }
            }
            results <- lapply(results, na.omit)
            return(results)
        }
    ''')
    
    logging.info("Differentially expressed genes based on genotype identified.")
    result = robjects.r['find_de'](sca)
    return result

# Run the DE analysis
de_results = find_de_MAST(sca)

# Convert results to DataFrames and prepare for merging
DE_results_df = {genotype: pandas2ri.rpy2py_dataframe(df) for genotype, df in de_results.items()}

# Concatenate all results and add `cell_state` as an additional column (without considering it in DE analysis)
all_results = pd.concat(
    [df.assign(genotype=genotype, cell_state=adata.obs.loc[adata.obs['genotype'] == genotype, 'cell_state'].values[0]) for genotype, df in DE_results_df.items()],
    ignore_index=True
)

# Save combined DE results to a single file
combined_file_name = f"/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results/DE_combined_all_genotypes_{timestamp}.txt"
all_results.to_csv(combined_file_name, sep='\t', index=False)
logging.info(f"Saved combined DE results to {combined_file_name}")

# Add DE results to AnnData
try:
    sc_toolbox.tools.de_res_to_anndata(
        adata,
        all_results,
        groupby="genotype",
        gene_id_col='gene_id',
        score_col='log_fold_change',
        pval_col='p_value',
        pval_adj_col="FDR",
        lfc_col='log_fold_change',
        key_added='MAST_results'
    )
except ValueError as e:
    logging.error(f"Error updating AnnData: {e}")

# Ensure all obs columns in AnnData are strings before saving
for col in adata.obs.columns:
    adata.obs[col] = adata.obs[col].astype(str)

# Save the updated AnnData object
adata_path = f"/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/data/dmd_annotated_with_MAST_{timestamp}.h5ad"
adata.write_h5ad(adata_path)
logging.info(f"Annotated AnnData saved to {adata_path}")

combined_file_name = f"/mnt/LaCIE/skolla/Github/DMD_mouse_Lopez-Hoffman_Kolla/differential_gene_expression/DE_results/DE_combined_all_genotypes_{timestamp}.txt"
all_results.to_csv(combined_file_name, sep='\t', index=False)
logging.info(f"Saved combined DE results to {combined_file_name}")