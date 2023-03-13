# Query the oneK1K dataset

This tool allows you to query the OneK1K single cell RNA-seq dataset. Given a list of SNPs of interest and eventually specific cell types, it will generate violin plots of gene expression across each SNP genotypes for each gene. Expression data from raw counts and normalized counts is saved as well for further analysis.

The tools will also generate a table of the top genes with the highest differential expression across genotypes for each SNP.

## Usage

You first need to activate the `single-cell-basic` conda environment:

```bash
conda activate single-cell-basic
```

Then the essential command is:

```bash
./onek1k_query.py --snp_id chr8_172477_T_C --genes ENSG00000198899,ENSG00000198786 --output test_out2
```

And you can specify cell types of interest using `--cell_types`:

```bash
./onek1k_query.py --snp_id chr8_172477_T_C --genes ENSG00000198899,ENSG00000198786 --output test_out2 --cell_types "CD4 TCM, CD8 TCM, CD8 TEM"
```

**NB.** Quotes are required when specifying cell types with spaces in their names.

### Inputs

All the three inputs `snp__id`, `genes` and `cell_types` can be specificed either as a comma-separated list or as a file containing one entry per line.

Cell types are tags from azymuth as annotated in the OneK1K dataset, and can be listed using the `--list_cell_types` option.

## Full list of options

```bash
usage: onek1k_query.py [-h] --snp_id SNP_ID --genes GENES [--cell_types CELL_TYPES] [--output_dir OUTPUT_DIR]
                       [--scdata SCDATA] [--genos_path GENOS_PATH] [--separate_genes] [--separate_cell_types]
                       [--n_top_genes N_TOP_GENES] [--rank_method {wilcoxon,t-test,t-test_overestim_var,logreg}]
                       [--list_cell_types]

Query the OneK1K datasets for a list of SNPs

options:
  -h, --help            show this help message and exit
  --snp_id SNP_ID       SNP ID in the format chr1_123456_A_T. Accept a comma-separated list or a file with 1 SNP per line
  --genes GENES         Ensembl IDs for genes of interest. Accept a comma-separated list or a file with 1 gene per line
  --cell_types CELL_TYPES
                        Cell types of interest. Accept a comma-separated list or a file with 1 cell type per line. Use
                        quotes if you have spaces
  --output_dir OUTPUT_DIR
                        Output directory
  --scdata SCDATA       Path to OneK1K dataset in h5ad formta
  --genos_path GENOS_PATH
                        Path to OneK1K genotypes folder
  --separate_genes      If set, plot a separate file for each gene
  --separate_cell_types
                        If set, and cell types are defined, analyse each cell type separately
  --n_top_genes N_TOP_GENES
                        Number of top genes to compute comparing genotypes
  --rank_method {wilcoxon,t-test,t-test_overestim_var,logreg}
                        Method to use for ranking genes
  --list_cell_types     List the available cell types and exit
```
