#!/usr/bin/env python

import argparse
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

#Parse arguments
parser = argparse.ArgumentParser(description='Make high quality plots from regenie results')
parser.add_argument('-s','--sumstats', type=str, required=True, help='Path to summary stats table')
parser.add_argument('-f','--format', type=str, default='regenie', help='Summary stats format among those accepted by gwaslab')
parser.add_argument('-o','--outdir', type=str, default='plots', help='Folder where to save the plots')
parser.add_argument('-a','--annotation_type', type=str, default='genes', choices=['genes','snpid'], help='Type of annotation to use for the Manhattan plot')
parser.add_argument('-m','--additional_log10p', type=float, help='Additional log10p threshold to use for the Manhattan plot')
parser.add_argument('-t','--annotated_tophits', type=str, help='Filename of the annotated top hits')
parser.add_argument('-l', '--toploci', type=str, help='Tab-separated file with header listing loci to plot. Columns: CHR, START, END, SNP')
parser.add_argument('-s', '--single_locus', type=str, help='Generate a regional plot for this single locus. Format CHR,START,END,SNP')
parser.add_argument('-w', '--window_kb', type=int, default=300, help='Window size in kb to use for the regional plots')
parser.add_argument('-g','--genome_build', type=str, default='hg19', choices=['hg19', 'hg38'], help='Genome build to use for the plots')
args = parser.parse_args()

regenie_merged = args.sumstats
manhattan_annotation_type = args.annotation_type
annotation_min_log10p = args.additional_log10p
annotated_tophits_filename = args.annotated_tophits
annotated_toploci_filename = args.toploci
regional_plot_window = args.window_kb * 1000
genome_build = args.genome_build.replace("hg", "")
single_locus = args.single_locus.split(",") if args.single_locus else []

#Import main libraries
import gwaslab as gl
import pandas as pd

#Load regenie results
logging.info(f"Loading regenie summary stat {regenie_merged}")
regenieTable = gl.Sumstats(regenie_merged, fmt="regenie", build=genome_build, sep="\t", verbose=False)

if 'MLOG10P' in regenieTable.data.columns:
  logging.info("Compute P values")
  regenieTable.data['P'] = 10 ** (-regenieTable.data['MLOG10P'])

#Load top hits
annotated_tophits=pd.DataFrame()
if annotated_tophits_filename:
	logging.info(f"Loading annotated top hits {annotated_tophits_filename}")
	annotated_tophits = pd.read_csv(annotated_tophits_filename, sep="\t")
	annotated_tophits.rename(columns={'ID': 'SNPID'}, inplace=True)

logging.info(f"Loaded tophits: {annotated_tophits.shape[0]}")

if manhattan_annotation_type == "genes" and not annotated_tophits.empty:
  logging.info("Using genes from tophits table for Manhattan plot")
  #Group by gene name and select the top hit
  top_snp_per_gene = annotated_tophits.groupby('CLOSEST_GENE_NAME').apply(lambda x: x.loc[x['LOG10P'].idxmax()]).loc[:, ['CHROM','GENPOS','SNPID','LOG10P']].reset_index()

  #To avoid overlapping SNPs, group by 500kb windows
  top_snp_per_gene['group'] = (top_snp_per_gene['GENPOS'] // 500000).astype(int)
  top_snp_per_gene = top_snp_per_gene.groupby(['CHROM','group']).apply(lambda x: x.loc[x['LOG10P'].idxmax()])

  #Reset index and select only the SNPID column
  top_snp_per_gene = top_snp_per_gene.loc[:, ['SNPID','CLOSEST_GENE_NAME']].reset_index(drop=True)

  #Merge top hits with regenie results
  regenieTable.data = pd.merge(regenieTable.data, top_snp_per_gene, on='SNPID', how='left')

  anno_value = "CLOSEST_GENE_NAME"
  anno_set = top_snp_per_gene['SNPID'].tolist()

else:
  #Get lead SNPs in each locus
  logging.info("Computing lead SNPs from regenie results for Manhattan plot")
  logging.info(f"sig_level: {annotation_min_log10p}")
  lead_snps = regenieTable.get_lead(
     windowsizekb=500,
     sig_level=annotation_min_log10p,
     anno=False,
     build=genome_build)

  anno_value = "SNPID"
  anno_set = lead_snps['SNPID'].tolist()

#Plot Manhattan and QQ
plot_filename = f"{args.outdir}/manhattan_qq.png"
additional_line = None
if annotation_min_log10p:
  additional_line = [10 ** -annotation_min_log10p]
logging.info(f"Plotting Manhattan and QQ to {plot_filename}")
figure, log = regenieTable.plot_mqq(
  mode="qqm",
  stratified=True, skip=2, 
  build=genome_build,
  anno=anno_value, anno_set=anno_set,
  highlight=anno_set,
  highlight_windowkb = 500,
  additional_line=additional_line, 
  additional_line_color=['orange'],
  suggestive_sig_line=True, sig_line=True,
  save=plot_filename,
  saveargs={"dpi":300, "height":6, "width":12},
  verbose=False)

#Plot regional plots
if len(single_locus) > 0:
  chrom = int(single_locus[0].replace("chr", ""))
  start_value = int(single_locus[1]) - regional_plot_window
  end_value = int(single_locus[2]) + regional_plot_window
  snp = single_locus[3]
  logging.info("Plotting regional plot for single locus")
  logging.info(f"Locus - {chrom}:{single_locus[1]}-{single_locus[2]}")
  plot_filename = f"{args.outdir}/${snp}_regionalplot.png"
  figure, log = regenieTable.plot_mqq(
    mode="r",
    build=genome_build,
    region_grid=True,
    region_ref = snp,
    region = (chrom, start_value, end_value),
    vcf_path=gl.get_path(f"1kg_eur_hg{genome_build}"),
    save=plot_filename,
    saveargs={"dpi":300, "height":6, "width":12},
    verbose=False)    
if annotated_toploci_filename:
	logging.info("Plotting regional plots")
	toploci_df = pd.read_csv(annotated_toploci_filename, sep="\t")
	for index, row in toploci_df.iterrows():
		locus_id = f"locus_{index+1}"
		logging.info(f"Locus {locus_id} - {row['CHR']}:{row['START']}-{row['END']}")
		plot_filename = f"{args.outdir}/{locus_id}_regionalplot.png"
		start_value = int(row['START']) - regional_plot_window
		end_value = int(row['END']) + regional_plot_window
		if start_value < 0: start_value = 0

		figure, log = regenieTable.plot_mqq(
			mode="r",
			build=genome_build,
			region_grid=True,
			region_ref = row['SNP'],
			region=(int(row['CHR']),start_value,end_value),
			vcf_path=gl.get_path(f"1kg_eur_hg{genome_build}"),
			verbose=False,
			save=plot_filename,
  		saveargs={"dpi":300, "height":5, "width":5},
			)
                
logging.info("All done")