import pandas as pd
import sys
import os

# Usage: python clean_report_final.py targeted_final.csv > summary_report.csv

# 1. Exact Gene Mapping Provided by User
gene_map = {
    "B9J08_01595": "ERG3",
    "B9J08_03698": "ERG11",
    "B9J08_02359": "FLO8",
    "B9J08_00960": "MEC3",
    "B9J08_04780": "TAC1b",
    "B9J08_02922": "FKS1",
    "B9J08_03606": "PEA2",
    "B9J08_01093": "CIS2",
    "B9J08_01933": "FUR1",
    "B9J08_02123": "CDR1",
    "B9J08_01918": "MRR1",
    "B9J08_05449": "rpsU"
}

def clean_snp_report(file_path):
    # Read the pipe-delimited file
    try:
        df = pd.read_csv(file_path, sep='|')
    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    # Remove empty trailing column (common with pipe delimiters)
    df = df.dropna(axis=1, how='all')

    # Clean Column Headers (Remove _sorted_marked.bam)
    df.columns = [c.replace('_sorted_marked.bam', '').strip() for c in df.columns]

    # Map Locus Tags to Gene Names
    # If a tag is NOT in the map, it remains as the original ID
    df['GENENAME'] = df['GENEID'].map(gene_map).fillna(df['GENEID'])

    # Define metadata vs sample columns
    # Based on standard output, samples start after 'ERRORS' column
    try:
        start_idx = df.columns.get_loc("ERRORS") + 1
    except KeyError:
        # Fallback if ERRORS column missing
        start_idx = 22 
    
    sample_cols = df.columns[start_idx:]

    # Prepare output list
    results = []

    for _, row in df.iterrows():
        gene = row['GENENAME']
        locus = row['GENEID'] # Original locus tag
        mutation = row['HGVS.p'] if pd.notna(row['HGVS.p']) and str(row['HGVS.p']).strip() != '' else row['HGVS.c']
        dna_change = row['HGVS.c']
        impact = row['PUTATIVE_IMPACT']
        
        # Find which samples have '1' (Present)
        present_in = []
        for sample in sample_cols:
            val = row[sample]
            # Check for 1 (integer or string representation)
            if val == 1 or val == '1':
                present_in.append(sample)
        
        # Only add if at least one sample has mutation
        if present_in:
            results.append({
                "Locus_tag": locus,  # First column as requested
                "Gene": gene,        # Second column
                "Protein Change": mutation,
                "DNA Change": dna_change,
                "Impact": impact,
                "Samples": ", ".join(present_in)
            })

    # Create cleaned DataFrame
    clean_df = pd.DataFrame(results)
    
    # Sort by Gene for readability
    if not clean_df.empty:
        clean_df = clean_df.sort_values(by="Gene")
        # Use to_csv instead of to_markdown to avoid dependencies
        print(clean_df.to_csv(index=False))
    else:
        print("No mutations found in target genes.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        clean_snp_report(sys.argv[1])
    else:
        print("Usage: python clean_report_final.py <input_file>")
