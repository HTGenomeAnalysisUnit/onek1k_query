#! /usr/bin/env python3

import scanpy as sc
import pandas as pd
import subprocess
import re
import os
import io
import sys
from datetime import datetime
import argparse
from shutil import which

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
ONEK1K_ADATA_FILE=os.path.join(SCRIPT_DIR, "onek1k_normalized.h5ad")
ONEK1K_GENOS_PATH="/project/alfredo/reference_data/OneK1K/genotypes/filter_vcf_r08/"
BCFTOOLS_BASE_CMD = ["bcftools", "query", "-f", "'[%CHROM %POS %REF %ALT %SAMPLE %GT\n]'"]
CELL_TYPES_FILE=os.path.join(SCRIPT_DIR, "cell_types.txt")

SNP_PARSE_RE = re.compile("chr([0-9XY]+)_([0-9]+)_([ATGC]+)_([ATGC]+)")

# parse the command line arguments
parser = argparse.ArgumentParser(description="Query the OneK1K datasets for a list of SNPs")
parser.add_argument("--snp_id", help="SNP ID in the format chr1_123456_A_T. Accept a comma-separated list or a file with 1 SNP per line", required='--list_cell_types' not in sys.argv)
parser.add_argument("--genes", help="Ensembl IDs for genes of interest. Accept a comma-separated list or a file with 1 gene per line", required='--list_cell_types' not in sys.argv)
parser.add_argument("--cell_types", help="Cell types of interest. Accept a comma-separated list or a file with 1 cell type per line. Use quotes if you have spaces", required=False, default="all")
parser.add_argument("--output_dir", help="Output directory", required=False, default="output")
parser.add_argument("--scdata", help="Path to OneK1K dataset in h5ad formta", required=False, default=ONEK1K_ADATA_FILE)
parser.add_argument("--genos_path", help="Path to OneK1K genotypes folder", required=False, default=ONEK1K_GENOS_PATH)
parser.add_argument("--separate_genes", help="If set, plot a separate file for each gene", required=False, action="store_true")
parser.add_argument("--separate_cell_types", help="If set, and cell types are defined, analyse each cell type separately", required=False, action="store_true")
parser.add_argument("--n_top_genes", help="Number of top genes to compute comparing genotypes", required=False, default=20, type=int)
parser.add_argument("--rank_method", help="Method to use for ranking genes", required=False, default="t-test", choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"])
parser.add_argument("--list_cell_types", help="List the available cell types and exit", required=False, action="store_true")
args = parser.parse_args()

# verify if a command is available
def program_is_available(fpath):
	return which(fpath) is not None

# given a start time return the time elapsed nicely formatted
def get_elapsed_time(start_time):
	elapsed_time = datetime.now() - start_time
	return f"{elapsed_time.seconds // 60}m {elapsed_time.seconds % 60}s"

# given a string, if it is a file path, return lines as a list, otherwise split the string on commas and return as a list
def get_list(s):
	if os.path.isfile(s):
		with open(s, 'r') as f:
			res = [x.strip() for x in f.readlines()]
	else:
		res = [x.strip() for x in s.split(",")]
	
	if len(res) == 0:
		raise ValueError(f"Invalid input: {s}")
	
	return res

# check all SNP IDs string are valid
def validate_snp_ids(snp_ids):
	for snp_id in snp_ids:
		if not re.match(SNP_PARSE_RE, snp_id):
			raise ValueError(f"Invalid snp_id: {snp_id}")
	log(f"Imported {len(snp_ids)} snp_ids")

# parse the snp_id to get region, ref, alt, and genotype file path
def get_vcf_search_params(snp_id, geno_dir):
	m = re.match(SNP_PARSE_RE, snp_id)
	chrom, pos, ref, alt = m.groups()
	pos = int(pos)
	geno_path = os.path.join(geno_dir, f"chr{chrom}.dose.filtered.R2_0.8.vcf.gz")
	region = f"{chrom}:{pos}-{pos}"
	return region, ref, alt, geno_path

# capture stdout from external command
def run_cmd(cmd):
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	if proc.returncode != 0:
		raise RuntimeError(f"Command failed: {cmd}")
	return stdout.decode('utf-8')

def annotate_genotypes(df):
	# add a column to tag the genotypes as 0, 1, or 2
	df.loc[df['genotype'] == '0|0', 'geno_tag'] = "0"
	df.loc[df['genotype'] == '0|1', 'geno_tag'] = "1"
	df.loc[df['genotype'] == '0|1', 'geno_tag'] = "1"
	df.loc[df['genotype'] == '1|1', 'geno_tag'] = "2"
	# convert genotypes to actual alleles
	df.loc[df['geno_tag'] == '0', 'geno_alleles'] = df['ref'].astype(str) + '/' + df['ref'].astype(str)
	df.loc[df['geno_tag'] == '1', 'geno_alleles'] = df['ref'].astype(str) + '/' + df['alt'].astype(str)
	df.loc[df['geno_tag'] == '2', 'geno_alleles'] = df['alt'].astype(str) + '/' + df['alt'].astype(str)
	# count of each genotype
	geno_counts = df['geno_tag'].value_counts()
	df.loc[df['geno_tag'] == '0', 'geno_label'] = df['geno_alleles'].astype(str) + f"({geno_counts[0]})"
	df.loc[df['geno_tag'] == '1', 'geno_label'] = df['geno_alleles'].astype(str) + f"({geno_counts[1]})"
	df.loc[df['geno_tag'] == '2', 'geno_label'] = df['geno_alleles'].astype(str) + f"({geno_counts[2]})"
	return df

def df_from_bcftools(stdout, ref, alt):
	# read the stdout into a dataframe
	df = pd.read_csv(io.StringIO(stdout), sep=" ", header=0, names=["chrom", "pos", "ref", "alt", "sample_id", "genotype"])
	df.set_index("sample_id", inplace=True)
	df = filter_alleles(df, ref, alt)
	if df.shape[0] == 0: 
		print("No data found for this variant")
		raise ValueError
	# translate genotypes to actual alleles
	df = annotate_genotypes(df)
	# convert all columns with dtype object to type category	
	df[['chrom', 'ref', 'alt', 'genotype', 'geno_tag', 'geno_alleles', 'geno_label']] = df[['chrom', 'ref', 'alt', 'genotype', 'geno_tag', 'geno_alleles', 'geno_label']].astype('category')
	return df[df.index.notnull()]

def filter_alleles(df, ref, alt):
	# filter the dataframe to only the alleles we want
	df = df[(df.ref == ref) & (df.alt == alt)]
	return df

def log(msg):
	print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ## {msg}")

def save_rank_results(uns_data, output_file):
	# save the gene rank results to a csv file, assuming 2 comparisons were made
	# first get the column names from uns keys
	col_names = [k for k in uns_data.keys() if k not in ['params','pts']]
	col_names={k : (f"{k}_1", f"{k}_2") for k in col_names}
	# then we make a dictionary of the results from uns using the column names
	res_dict = {}
	for key, names in col_names.items():
		for i in range(len(names)):
			res_dict[names[i]] = [x[i] for x in uns_data[key]]
	rnk_results = pd.DataFrame.from_dict(res_dict)
	rnk_results.to_csv(output_file, index=False)

# Print list of available cell types if requested
if args.list_cell_types:
	with open(CELL_TYPES_FILE) as f:
		cell_types = f.readlines()
	print("## Available cell types ##")
	for cell_type in cell_types:
		print(f"{cell_type.strip()}")
	quit()

# Check if we have bcftools installed
if not program_is_available("bcftools"):
	raise RuntimeError("bcftools is not installed")

# Variables
snp_ids = get_list(args.snp_id) #chr8_172477_T_C
validate_snp_ids(snp_ids)
output_main_folder = args.output_dir
os.makedirs(output_main_folder, exist_ok=True)
genes = get_list(args.genes) #['ENSG00000198899', 'ENSG00000198786']
log(f"Imported {len(genes)} genes")

cell_types = []
if args.cell_types != "all":
	cell_types = get_list(args.cell_types) #['CD4 T cells', 'CD8 T cells']
log(f"Imported {len(cell_types)} cell types")

separate_genes = args.separate_genes
separate_cell_types = args.separate_cell_types
n_top_genes = args.n_top_genes
rank_method = args.rank_method
log(f"Find the top {n_top_genes} marker genes using {rank_method} method")

genotype_path = args.genos_path

# load single-cell data
log("loading single-cell data")
scdata = sc.read_h5ad(args.scdata)

# open a file to store ids of not found SNPs
not_found_file = open(os.path.join(output_main_folder, "SNPs_not_found.txt"), "w")

start_time = datetime.now()
for snp_id in snp_ids:
	start_time_snp = datetime.now()
	log(f"processing snp {snp_id}")
	out_prefix = f"_{snp_id}"
	output_data_folder = os.path.join(output_main_folder, snp_id, "data")
	output_fig_folder = os.path.join(output_main_folder, snp_id, "figs")
	sc.settings.figdir = output_fig_folder

	region, ref, alt, geno_path = get_vcf_search_params(snp_id, genotype_path)

	# search for the snp in the vcf
	log(f"searching for snp {snp_id} in genotype data")
	cmd = BCFTOOLS_BASE_CMD + ["-r", region, geno_path]
	stdout = run_cmd(cmd)
	df = df_from_bcftools(stdout, ref, alt)
	if df.shape[0] == 0: 
		log("WARN - No data found for this variant")
		not_found_file.write(f"{snp_id}\n")
		continue

	# create output folder if it doesn't exist
	os.makedirs(output_data_folder, exist_ok=True)
	os.makedirs(output_fig_folder, exist_ok=True)

	# merge the genotypes with the scdata
	scdata.obs = scdata.obs.merge(df, left_on='donor_id', right_index=True, how='left')

	# now add the snp_id to the obs
	scdata.obs['snp_id'] = snp_id

	# make violin plot stratified by genotype for a given gene
	log(f"making violin plots for {len(genes)} genes")
	if len(cell_types) > 0:
		if separate_cell_types:
			for cell_type in cell_types:
				if separate_genes:
					for gene in genes:
						sc.pl.violin(scdata[scdata.obs['predicted.celltype.l2'] == cell_type], keys=gene, groupby="geno_label", save=f"{out_prefix}-{cell_type.replace(' ','_')}-{gene}.png", inner='box', show=False, jitter=False, log=True, layer='PFlog1pPF_normalization')
				else:
					sc.pl.violin(scdata[scdata.obs['predicted.celltype.l2'] == cell_type], keys=genes, groupby="geno_label", save=f"{out_prefix}-{cell_type.replace(' ','_')}.png", inner='box', show=False, jitter=False, log=True, layer='PFlog1pPF_normalization')
		else:
			if separate_genes:
				for gene in genes:
					sc.pl.violin(scdata[scdata.obs['predicted.celltype.l2'].isin(cell_types)], keys=gene, groupby="geno_label", save=f"{out_prefix}-{gene}.png", inner='box', show=False, jitter=False, log=True, layer='PFlog1pPF_normalization')
			else:
				sc.pl.violin(scdata[scdata.obs['predicted.celltype.l2'].isin(cell_types)], keys=genes, groupby="geno_label", save=f"{out_prefix}.png", inner='box', show=False, jitter=False, log=True, layer='PFlog1pPF_normalization')
	else:
		if separate_genes:
			for gene in genes:
				sc.pl.violin(scdata, keys=gene, groupby="geno_label", save=f"{out_prefix}_{gene}.png", inner='box', show=False, jitter=False, log=True, layer='PFlog1pPF_normalization')
		else:
			sc.pl.violin(scdata, keys=genes, groupby="geno_label", save=f"{out_prefix}.png", inner='box', show=False, jitter=False, log=True, layer='PFlog1pPF_normalization')

	# save obs for the relevant cell types and expression of the genes of interest to a file
	log("saving expression data to file")
	if len(cell_types) > 0:
		output_scdata = scdata[scdata.obs['predicted.celltype.l2'].isin(cell_types)].copy()
	else:
		output_scdata = scdata

	output_df = scdata.obs 

	output_df[genes] = scdata[:, genes].X.toarray()
	output_df.to_csv(os.path.join(output_data_folder, f"{snp_id}_expression_raw.csv"))

	output_df[genes] = scdata[:, genes].layers.get('PFlog1pPF_normalization').toarray()
	output_df.to_csv(os.path.join(output_data_folder, f"{snp_id}_expression_norm.csv"))

	# identify marker genes across the genotypes
	if len(cell_types) > 0:
		if separate_cell_types:
			for cell_type in cell_types:
				cell_scdata = scdata[scdata.obs['predicted.celltype.l2'] == cell_type].copy()
				log(f"finding the {n_top_genes} top marker genes for cell type {cell_type}")
				log("computing for raw counts")
				sc.tl.rank_genes_groups(cell_scdata, groupby="geno_tag", groups=['0', '1', '2'], reference='0', method=rank_method, n_genes=n_top_genes, use_raw=True, key_added=f"{rank_method}_raw", pts=True)
				save_rank_results(cell_scdata.uns[f"{rank_method}_raw"], os.path.join(output_data_folder, f"{snp_id}_{cell_type}_top_{n_top_genes}_genes.raw_counts.csv"))
				log("computing for normalized counts")				
				sc.tl.rank_genes_groups(cell_scdata, groupby="geno_tag", groups=['0', '1', '2'], reference='0', method=rank_method, n_genes=n_top_genes, use_raw=False, layer='PFlog1pPF_normalization', key_added=f"{rank_method}_norm", pts=True)
				save_rank_results(cell_scdata.uns[f"{rank_method}_norm"], os.path.join(output_data_folder, f"{snp_id}_{cell_type}_top_{n_top_genes}_genes.norm_counts.csv"))
		else:
			cell_scdata = scdata[scdata.obs['predicted.celltype.l2'].isin(cell_types)].copy()
			log(f"finding the {n_top_genes} top marker genes for cell types {cell_types}")
			log("computing for raw counts")
			sc.tl.rank_genes_groups(cell_scdata, groupby="geno_tag", groups=['0', '1', '2'], reference='0', method=rank_method, n_genes=n_top_genes, use_raw=True, key_added=f"{rank_method}_raw", pts=True)
			save_rank_results(cell_scdata.uns[f"{rank_method}_raw"], os.path.join(output_data_folder, f"{snp_id}_top_{n_top_genes}_genes.raw_counts.csv"))
			log("computing for normalized counts")
			sc.tl.rank_genes_groups(cell_scdata, groupby="geno_tag", groups=['0', '1', '2'], reference='0', method=rank_method, n_genes=n_top_genes, use_raw=False, layer='PFlog1pPF_normalization', key_added=f"{rank_method}_norm", pts=True)
			save_rank_results(cell_scdata.uns[f"{rank_method}_norm"], os.path.join(output_data_folder, f"{snp_id}_top_{n_top_genes}_genes.norm_counts.csv"))
		params_raw = cell_scdata.uns[f"{rank_method}_raw"]["params"]
		params_norm = cell_scdata.uns[f"{rank_method}_norm"]["params"]
	else:
		log(f"find the {n_top_genes} top marker genes using {rank_method} method")
		log("computing for raw counts")
		sc.tl.rank_genes_groups(scdata, groupby="geno_tag", groups=['0', '1', '2'], reference='0', method=rank_method, n_genes=n_top_genes, use_raw=True, key_added=f"{rank_method}_raw", pts=True)
		save_rank_results(scdata.uns[f"{rank_method}_raw"], os.path.join(output_data_folder, f"{snp_id}_raw_top_{n_top_genes}_genes.csv"))
		# or we can use layer='PFlog1pPF_normalization' to use the normalized expression
		log("computing for normalized counts")
		sc.tl.rank_genes_groups(scdata, groupby="geno_tag", groups=['0', '1', '2'], reference='0', method=rank_method, n_genes=n_top_genes, use_raw=False, layer='PFlog1pPF_normalization', key_added=f"{rank_method}_norm", pts=True)
		save_rank_results(scdata.uns[f"{rank_method}_norm"], os.path.join(output_data_folder, f"{snp_id}_norm_top_{n_top_genes}_genes.csv"))
		params_raw = scdata.uns[f"{rank_method}_raw"]["params"]
		params_norm = scdata.uns[f"{rank_method}_norm"]["params"]

	# save params of rank_genes_groups to file
	with open(os.path.join(output_data_folder, f"rank_genes_groups_params.txt"), 'w') as f:
		f.write("=== RAW DATA ===\n")
		f.write(f"{params_raw}\n")
		f.write("=== NORM DATA ===\n")
		f.write(f"{params_norm}\n")
	
	log(f"{snp_id} completed in {get_elapsed_time(start_time_snp)}")

not_found_file.close()

log("All done!")
log(f"Total time: {get_elapsed_time(start_time)}")