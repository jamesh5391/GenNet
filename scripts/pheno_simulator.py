import pandas as pd
import numpy as np
import argparse

def select_causal_snps(bim_file_path: str, polygenicity: int, output_file: str):
    """
    Reads a PLINK .bim file and randomly selects a number of SNPs to be causal.

    Args:
        bim_file_path (str): Path to the .bim file.
        polygenicity (int): The number of causal SNPs to select.
        output_file (str): Path to save the causal SNP list.
    """
    print(f"Reading SNPs from {bim_file_path}...")
    bim_df = pd.read_csv(bim_file_path, sep='\s+', header=None, names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
    
    # Randomly select SNPs
    causal_snps_df = bim_df.sample(n=polygenicity)
    
    # Assign random effect sizes (betas) between -0.5 and 0.5
    causal_snps_df['beta'] = np.random.uniform(-0.5, 0.5, size=polygenicity)
    
    # Save to a format PLINK can read
    causal_snps_df[['snp', 'beta']].to_csv(output_file, sep='\t', index=False, header=False)
    print(f"Saved {polygenicity} causal SNPs and their effects to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select causal SNPs for phenotype simulation.")
    parser.add_argument("--bim", required=True, help="Path to the input .bim file.")
    parser.add_argument("--polygenicity", type=int, default=1000, help="Number of causal SNPs.")
    parser.add_argument("--out", default="causal_snps.txt", help="Output file name.")
    args = parser.parse_args()
    
    select_causal_snps(args.bim, args.polygenicity, args.out)
