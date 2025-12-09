import pandas as pd
import sys

# Usage: python clean_report.py raw_snpeff_output.txt > clean_report.csv

# 1. Map Locus Tags to Gene Names
gene_map = {
    "B9J08_03698": "ERG11",
    "B9J08_02922": "FKS1",
    "B9J08_05449": "FUR1"
}

def clean_snp_report(file_path):
    # Read the pipe-delimited file
    df = pd.read_csv(file_path, sep='|')

    # Remove the empty last column if it exists (pandas often reads trailing pipes as empty cols)
    df = df.dropna(axis=1, how='all')

    # Rename Genes
    df['GENENAME'] = df['GENEID'].map(gene_map).fillna(df['GENEID'])

    # Clean Column Headers (Remove _sorted_marked.bam)
    df.columns = [c.replace('_sorted_marked.bam', '') for c in df.columns]

    # Define metadata columns (everything before the first sample)
    # Based on your header, samples start at column index 22 (0-based)
    meta_cols = df.columns[:22]
    sample_cols = df.columns[22:]

    # Prepare output list
    results = []

    for _, row in df.iterrows():
        gene = row['GENENAME']
        # Use HGVS.p (Protein) if available, otherwise HGVS.c (DNA)
        mutation = row['HGVS.p'] if pd.notna(row['HGVS.p']) else row['HGVS.c']
        dna_change = row['HGVS.c']
        impact = row['PUTATIVE_IMPACT']
        
        # Find which samples have '1' (Present)
        present_in = []
        for sample in sample_cols:
            if row[sample] == 1:
                present_in.append(sample)
        
        # Only add if at least one sample has it
        if present_in:
            results.append({
                "Gene": gene,
                "Protein_Change": mutation,
                "DNA_Change": dna_change,
                "Impact": impact,
                "Samples": ", ".join(present_in)
            })

    # Create cleaned DataFrame
    clean_df = pd.DataFrame(results)
    
    # Print or Save
    print(clean_df.to_markdown(index=False))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        clean_snp_report(sys.argv[1])
    else:
        print("Please provide the input filename.")
