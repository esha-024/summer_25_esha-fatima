import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
from collections import defaultdict

# -------------PARSING GFF FILE------------

def parse_gff(file):
    print("\n Parsing GFF File for Exon Coordinates...")
    chrom_map = {'NC_004354.4': '2L', 'NC_004353.4': '2R', 'NC_004351.4': '3L',
                 'NC_004350.4': '3R', 'NC_004352.4': 'X', 'NC_004355.4': '4',
                 'NC_004356.4': 'Y', 'NC_024511.2': 'X', 'NC_024512.1': 'Y'}
    exons = []
    with open(file) as f:
        for line in f:
            if line.startswith('#') or '\texon\t' not in line:
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            seqid = chrom_map.get(parts[0], parts[0])
            start, end, strand = int(parts[3]), int(parts[4]), parts[6]
            attr = {k: v for item in parts[8].split(';') if '=' in item for k, v in [item.split('=')]}
            gene_id = attr.get('Parent') or attr.get('gene_id')
            if gene_id:
                exons.append({'seqid': seqid, 'start': start, 'end': end, 'strand': strand, 'gene_id': gene_id})
    exons_df = pd.DataFrame(exons)
    print(f"Total exons parsed: {len(exons_df)}")
    print("Preview of Exons:")
    print(exons_df.head(), "\n")
    return exons_df

#------------PARSING VCF FILE----------------

def parse_vcf(file):
    print(" Parsing VCF File for SNPs...")
    snps = []
    with open(file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            try:
                snps.append({'chrom': parts[0], 'pos': int(parts[1])})
            except:
                continue
    snps_df = pd.DataFrame(snps)
    print(f"Total SNPs parsed: {len(snps_df)}")
    print("Preview of SNPs:")
    print(snps_df.head(), "\n")
    return snps_df

#-----------------MAPPING SNPS TO EXONS-----------------------

def map_snps_to_exons(snps, exons):
    print(" Mapping SNPs to Exons...")
    snp_counts = defaultdict(int)
    for chrom in set(exons['seqid']) & set(snps['chrom']):
        snp_chr = snps[snps['chrom'] == chrom]
        for _, exon in exons[exons['seqid'] == chrom].iterrows():
            match = snp_chr[(snp_chr['pos'] >= exon['start']) & (snp_chr['pos'] <= exon['end'])]
            if not match.empty:
                snp_counts[exon['gene_id']] += len(match)
    mapped_df = pd.DataFrame(list(snp_counts.items()), columns=['gene_id', 'snp_count'])
    print(f"SNPs mapped to {len(mapped_df)} genes")
    print("Preview of Mapped SNP Counts:")
    print(mapped_df.head(), "\n")
    return mapped_df

#------------NORMALIZING GENE IDs-------------------------

def normalize_id(gid):
    if pd.isna(gid): return None
    gid = str(gid).strip()
    gid = re.sub(r'^(rna-|gene-|Dmel_)|(\.\d+)$', '', gid)
    return gid.split('-')[0] if gid.startswith("FBtr") else gid

#-------------BUIL MAPPING DICTIONARY-------------------

def build_mapping(mapping_file):
    print(" Building Mapping Dictionary...")
    mapping = {}
    df = pd.read_csv(mapping_file, sep='\t')
    for _, row in df.iterrows():
        fbgn = str(row.get('Gene stable ID', '')).strip()
        if not fbgn: continue
        mapping[fbgn] = fbgn
        for col in ['RefSeq mRNA ID', 'Transcript stable ID', 'FlyBase gene ID', 'Gene name']:
            val = str(row.get(col, '')).strip()
            if val:
                norm = normalize_id(val)
                mapping[norm] = fbgn
                for match in re.findall(r'(NM_\d+|NR_\d+)', val):
                    mapping[match] = fbgn
    print(f"Mapping dictionary with {len(mapping)} entries created\n")
    return mapping

#-------------MAIN ANALYSIS-----------------
def main(gff, vcf, expr_file, map_file=None):
    if not all(os.path.exists(f) for f in [gff, vcf, expr_file]):
        print("One or more input files not found.")
        return

    exons = parse_gff(gff)
    snps = parse_vcf(vcf)
    if exons.empty or snps.empty:
        print("Parsed data is empty.")
        return

    snp_df = map_snps_to_exons(snps, exons)

#----------------CALCULATING SNP DENSITY-----------------------

    print(" Calculating Exon Lengths and SNP Density...")
    exon_len = exons.assign(length=exons.end - exons.start + 1).groupby('gene_id')['length'].sum().to_dict()
    snp_df['exon_length'] = snp_df['gene_id'].map(exon_len)
    snp_df['snp_density'] = snp_df['snp_count'] / snp_df['exon_length']
    print("Preview of SNP Density Data:")
    print(snp_df.head(), "\n")

#---------------READING EXPRESSION FILE-----------------

    print(" Reading Expression File...")
    exp = pd.read_csv(expr_file)
    exp = exp.rename(columns={'primary_FBid': 'gene_id'})
    exp['expression'] = exp[[c for c in exp.columns if c not in ['gene_id', 'current_symbol']]].mean(axis=1)
    exp = exp[['gene_id', 'expression']]
    print(f"Expression values parsed: {len(exp)}")
    print("Preview of Expression Data:")
    print(exp.head(), "\n")

    snp_df['gene_id_norm'] = snp_df['gene_id'].apply(normalize_id)
    exp['gene_id_norm'] = exp['gene_id'].apply(normalize_id)

    if map_file and os.path.exists(map_file):
        mapping = build_mapping(map_file)
        snp_df['gene_id_mapped'] = snp_df['gene_id_norm'].map(mapping).fillna(snp_df['gene_id_norm'])
        exp['gene_id_mapped'] = exp['gene_id_norm'].map(mapping).fillna(exp['gene_id_norm'])
    else:
        snp_df['gene_id_mapped'] = snp_df['gene_id_norm']
        exp['gene_id_mapped'] = exp['gene_id_norm']

#---------------MERGING SNPS AND EXPRESSION DATA IN CSV-----------------

    print("Merging SNP and Expression Data...")
    merged = pd.merge(snp_df, exp, on='gene_id_mapped')
    merged = merged.dropna(subset=['expression', 'snp_density'])
    print(f"Merged dataset contains {len(merged)} genes")
    print("Preview of Merged Data:")
    print(merged.head(), "\n")

   
    output_csv = "snp_expression_merged_results.csv"
    merged.to_csv(output_csv, index=False)
    print(f"Merged results saved to: {output_csv}\n")

    if merged.empty:
        print("No common genes after merging.")
        return

#----------------------CALCULATING CORRELATION-------------------
    print("Calculating Pearson Correlation...")
    r = merged['expression'].corr(merged['snp_density'])
    print(f"Pearson Correlation (original data): r = {r:.4f}")

#-------------VISUALIZATION-----------------    

    print(" Plotting Scatter Plot ...")
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=merged, x='expression', y='snp_density', alpha=0.6, s=50)
    plt.title(f"SNP Density vs Gene Expression\nPearson r = {r:.3f}")
    plt.xlabel("Gene Expression")
    plt.ylabel("SNP Density")
    plt.tight_layout()
    plt.savefig("snp_vs_expression.png", dpi=300)
    plt.show()

    print(" Generating Heatmap of Top 20 Variable Genes...")
    exp_raw = pd.read_csv(expr_file)
    exp_raw.set_index("primary_FBid", inplace=True)
    exp_raw.iloc[:, 1:] = exp_raw.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
    filter_mask = (exp_raw.iloc[:, 1:] > 10).sum(axis=1) >= 3
    df_filtered = exp_raw[filter_mask]
    df_log2 = df_filtered.copy()
    df_log2.iloc[:, 1:] = np.log2(df_log2.iloc[:, 1:] + 1)
    expression_only = df_log2.iloc[:, 1:]
    gene_variance = expression_only.var(axis=1)
    top_genes = gene_variance.sort_values(ascending=False).head(20).index
    plt.figure(figsize=(14, 10))
    sns.heatmap(expression_only.loc[top_genes], cmap="viridis", xticklabels=True, yticklabels=True)
    plt.title("Top 20 Most Variable Genes (Log2 Normalized)", fontsize=14)
    plt.xlabel("Samples")
    plt.ylabel("Genes")
    plt.tight_layout()
    plt.savefig("top_variable_genes_heatmap.png", dpi=300)
    plt.close()
    print("Heatmap saved as: top_variable_genes_heatmap.png\n")


# --- Run the script ---
if __name__ == "__main__":
    main(
        gff="Drosophila_Melanogaster.gff",
        vcf="snps.vcf",
        expr_file="expression.csv",
        map_file="map_file.txt"
    )