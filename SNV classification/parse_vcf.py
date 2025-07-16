import pandas as pd

def parse_info_field(info_str):
    """
    Convert INFO column string into a dictionary.
    """
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def vcf_to_csv(vcf_file, output_csv="clinvar_parsed.csv"):
    data = []

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.strip().lstrip('#').split('\t')
                continue

            cols = line.strip().split('\t')
            row = dict(zip(header, cols))
            info = parse_info_field(row["INFO"])

            # Interpret CLNSIG (pathogenic or benign)
            clnsig = info.get('CLNSIG', '').lower()
            if 'pathogenic' in clnsig:
                label = 1
            elif 'benign' in clnsig:
                label = 0
            else:
                continue  # skip uncertain/conflicting

            gene = info.get("GENEINFO", "unknown").split(":")[0]
            consequence = info.get("MC", "unknown|unknown").split('|')[1] if '|' in info.get("MC", "") else "unknown"

            try:
                af_exac = float(info.get("AF_EXAC", 0))
            except:
                af_exac = 0

            data.append({
                "chrom": row["CHROM"],
                "pos": int(row["POS"]),
                "ref": row["REF"],
                "alt": row["ALT"],
                "gene": gene,
                "consequence": consequence,
                "af_exac": af_exac,
                "variant_class": info.get("CLNVC", "unknown"),
                "review_status": info.get("CLNREVSTAT", "unknown"),
                "label": label
            })

    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f" Saved {len(df)} variants to {output_csv}")

if __name__ == "__main__":
    vcf_to_csv("clinvar.vcf")  
