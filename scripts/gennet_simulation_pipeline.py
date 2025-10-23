#!/usr/bin/env python
"""
gennet_simulation_pipeline.py - Fixed GenNet output directory issue
"""

import subprocess
import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path
import shutil

class GenNetSimulation:
    def __init__(self, vcf_file, output_dir="gennet_simulation"):
        self.vcf_file = vcf_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
    def prepare_plink_data(self, n_samples=1000):
        """Step 1: Convert VCF to PLINK format"""
        print("=== Step 1: Preparing PLINK data ===")
        
        # First, convert VCF to PLINK without filtering
        temp_prefix = self.output_dir / "temp_all"
        cmd = f"""
        plink2 --vcf {self.vcf_file} \
               --make-bed \
               --out {temp_prefix}
        """
        subprocess.run(cmd, shell=True, check=True)
        
        # Read FAM file to get sample IDs in correct format
        fam = pd.read_csv(f"{temp_prefix}.fam", sep=r'\s+', header=None,
                         names=['fid', 'iid', 'father', 'mother', 'sex', 'pheno'])
        
        # Select random samples
        n_available = len(fam)
        n_to_select = min(n_samples, n_available)
        selected_indices = np.random.choice(n_available, n_to_select, replace=False)
        selected_samples = fam.iloc[selected_indices]
        
        # Create keep file with FID and IID columns
        samples_file = self.output_dir / "samples_keep.txt"
        selected_samples[['fid', 'iid']].to_csv(samples_file, sep='\t', header=False, index=False)
        
        print(f"Selected {n_to_select} samples from {n_available} available")
        
        # Now filter with the keep file - output to a clean directory
        plink_clean_dir = self.output_dir / "plink_clean"
        plink_clean_dir.mkdir(exist_ok=True)
        plink_prefix = plink_clean_dir / "chr22_data"
        
        cmd = f"""
        plink2 --bfile {temp_prefix} \
               --keep {samples_file} \
               --make-bed \
               --maf 0.01 \
               --geno 0.05 \
               --out {plink_prefix}
        """
        subprocess.run(cmd, shell=True, check=True)
        
        # Clean up temp files
        for ext in ['.bed', '.bim', '.fam', '.log']:
            temp_file = Path(f"{temp_prefix}{ext}")
            if temp_file.exists():
                temp_file.unlink()
        
        return str(plink_prefix)
    
    def convert_to_gennet(self, plink_prefix):
        """Step 2: Use GenNet's converter with correct output directory"""
        print("=== Step 2: Converting to GenNet format ===")
        
        # Create a temporary clean directory with only PLINK files
        convert_dir = self.output_dir / "gennet_convert_temp"
        if convert_dir.exists():
            shutil.rmtree(convert_dir)
        convert_dir.mkdir()
        
        # Copy only the PLINK files to the clean directory
        base_name = Path(plink_prefix).name
        for ext in ['.bed', '.bim', '.fam']:
            src = Path(f"{plink_prefix}{ext}")
            dst = convert_dir / f"{base_name}{ext}"
            if src.exists():
                shutil.copy2(src, dst)
        
        # GenNet wants an output DIRECTORY, not a file
        # Clean up any existing h5 output directory
        h5_output_dir = self.output_dir / "h5_output"
        if h5_output_dir.exists():
            shutil.rmtree(h5_output_dir)
        h5_output_dir.mkdir()
        
        # GenNet convert command with output directory
        gennet_script = Path("../GenNet.py")
        if not gennet_script.exists():
            gennet_script = Path("GenNet.py")
        
        if gennet_script.exists():
            # Use output directory, not file
            cmd = f"""python {gennet_script} convert \
                -g {convert_dir} \
                -study_name {base_name} \
                -o {h5_output_dir}"""
            
            print(f"Running GenNet convert command:\n{cmd}")
            
            try:
                subprocess.run(cmd, shell=True, check=True)
                print(f"Successfully converted. Output in {h5_output_dir}")
                
                # Find the actual h5 file created
                h5_files = list(h5_output_dir.glob("**/*.h5"))
                if h5_files:
                    print(f"Found HDF5 file: {h5_files[0]}")
                    return str(h5_files[0])
                else:
                    print("Warning: No .h5 file found after conversion")
                    return str(h5_output_dir)
                    
            except subprocess.CalledProcessError as e:
                print(f"Warning: GenNet conversion failed: {e}")
                print("Continuing with simulation setup...")
                return str(h5_output_dir)
            finally:
                # Clean up temporary convert directory
                if convert_dir.exists():
                    shutil.rmtree(convert_dir)
        else:
            print("Warning: GenNet.py not found, skipping HDF5 conversion")
            return str(h5_output_dir)
    
    def create_topology(self, plink_prefix, n_genes=500):
        """Step 3: Create topology file with correct GenNet format"""
        print("=== Step 3: Creating topology ===")

        # Read BIM file
        bim = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None,
                         names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])

        print(f"Creating topology for {len(bim)} SNPs")

        # Create gene to pathway mapping
        n_pathways = max(1, n_genes // 20)  # ~20 genes per pathway
        gene_to_pathway = {}
        pathway_idx = 0
        for gene_idx in range(n_genes):
            if gene_idx > 0 and gene_idx % 20 == 0:
                pathway_idx += 1
            gene_to_pathway[gene_idx] = pathway_idx

        topology_rows = []
        for snp_idx, row in bim.iterrows():
            # Assign SNP to gene (distribute evenly)
            gene_idx = snp_idx % n_genes
            pathway_idx = gene_to_pathway[gene_idx]

            topology_rows.append({
                'chr': row['chr'],
                'layer0_node': snp_idx,  # SNP index
                'layer0_name': row['snp'],  # SNP name
                'layer1_node': gene_idx,  # Gene index
                'layer1_name': f"GENE_{gene_idx}",  # Gene name
                'layer2_node': pathway_idx,  # Pathway index
                'layer2_name': f"PATHWAY_{pathway_idx}"  # Pathway name
            })

        topology_file = self.output_dir / "topology.csv"
        pd.DataFrame(topology_rows).to_csv(topology_file, index=False)

        print(f"Topology created: {len(topology_rows)} SNPs -> {n_genes} genes -> {len(set(gene_to_pathway.values()))} pathways")

        return str(topology_file), topology_rows
    
    def simulate_phenotype(self, plink_prefix, topology_data, n_causal_genes=50, h2=0.5, seed=None):
        """Step 4: Simulate phenotype using plink2 - Fixed subject file format"""
        if seed:
            np.random.seed(seed)
        
        print(f"=== Simulating phenotype (n_genes={n_causal_genes}, h2={h2}) ===")
        
        # Get FAM file for sample IDs
        fam = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None,
                        names=['fid', 'iid', 'father', 'mother', 'sex', 'pheno'])
        
        # Get unique genes from topology
        topology_df = pd.DataFrame(topology_data)
        genes = topology_df['layer1_name'].unique()

        # Select causal genes
        causal_genes = np.random.choice(genes, min(n_causal_genes, len(genes)), replace=False)

        # Get causal SNPs
        causal_snps = topology_df[topology_df['layer1_name'].isin(causal_genes)]['layer0_name'].values
        
        print(f"Selected {len(causal_snps)} causal SNPs from {n_causal_genes} genes")
        
        # Create effect file for PLINK scoring
        bim = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None,
                        names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
        
        effect_size = 1.0 / np.sqrt(len(causal_snps))
        effects = []
        for snp in causal_snps:
            snp_info = bim[bim['snp'] == snp]
            if len(snp_info) > 0:
                effects.append({
                    'SNP': snp,
                    'A1': snp_info.iloc[0]['a1'],
                    'BETA': effect_size
                })
        
        effects_file = self.output_dir / f"effects_h{h2}_g{n_causal_genes}.txt"
        pd.DataFrame(effects).to_csv(effects_file, sep='\t', header=False, index=False)
        
        # Calculate GRS using plink2
        grs_prefix = self.output_dir / f"grs_h{h2}_g{n_causal_genes}"
        
        cmd = f"""plink2 --bfile {plink_prefix} \
                --score {effects_file} 1 2 3 \
                --out {grs_prefix}"""
        
        subprocess.run(cmd, shell=True, check=True)
        
        # Read the score output
        score_file = f"{grs_prefix}.sscore"
        if not Path(score_file).exists():
            score_file = f"{grs_prefix}.profile"
        
        grs_df = pd.read_csv(score_file, sep=r'\s+')
        
        # Get the score column
        score_col = None
        for col in ['SCORE1_SUM', 'SCORE1_AVG', 'SCORESUM', 'SCORE']:
            if col in grs_df.columns:
                score_col = col
                break
        
        if score_col is None:
            print(f"Warning: Could not find score column. Available columns: {grs_df.columns.tolist()}")
            grs = np.random.normal(0, 1, len(fam))
        else:
            grs = grs_df[score_col].values
        
        # Standardize GRS
        if grs.std() > 0:
            grs_std = (grs - grs.mean()) / grs.std()
        else:
            grs_std = grs - grs.mean()
        
        # Add environmental component
        genetic_component = np.sqrt(h2) * grs_std
        environmental_component = np.sqrt(1 - h2) * np.random.normal(0, 1, len(grs))
        liability = genetic_component + environmental_component
        
        # Binary phenotype
        threshold = np.percentile(liability, 50)  # 50% prevalence
        phenotype = (liability > threshold).astype(int)
        
        print(f"Created phenotype: {phenotype.sum()} cases, {len(phenotype) - phenotype.sum()} controls")
        
        # Create subject file with CORRECT GenNet format
        n_samples = len(fam)
        
        # Split into train/val/test (60/20/20 as per GenNet paper)
        np.random.seed(seed)
        indices = np.arange(n_samples)
        np.random.shuffle(indices)
        
        train_idx = indices[:int(0.6 * n_samples)]
        val_idx = indices[int(0.6 * n_samples):int(0.8 * n_samples)]
        test_idx = indices[int(0.8 * n_samples):]
        
        # Create set assignments (GenNet expects 1=train, 2=val, 3=test)
        set_assignment = np.zeros(n_samples, dtype=int)
        set_assignment[train_idx] = 1  # train
        set_assignment[val_idx] = 2    # validation
        set_assignment[test_idx] = 3   # test
        
        # Create GenNet-format subject file
        subject_df = pd.DataFrame({
            'patient_id': fam['iid'].values,  # Use IID as patient ID
            'labels': phenotype,
            'genotype_row': np.arange(n_samples),  # Row index in genotype matrix
            'set': set_assignment
        })
        
        subject_file = self.output_dir / f"subjects_h{h2}_g{n_causal_genes}_gennet.csv"
        subject_df.to_csv(subject_file, index=False)
        
        print(f"Train: {sum(set_assignment==1)}, Val: {sum(set_assignment==2)}, Test: {sum(set_assignment==3)}")
        
        # Save metadata
        metadata = {
            'n_causal_genes': int(n_causal_genes),
            'n_causal_snps': int(len(causal_snps)),
            'heritability': float(h2),
            'n_train': int(sum(set_assignment==1)),
            'n_val': int(sum(set_assignment==2)),
            'n_test': int(sum(set_assignment==3)),
            'causal_genes': [str(gene) for gene in causal_genes.tolist()[:10]]
        }
        
        metadata_file = self.output_dir / f"metadata_h{h2}_g{n_causal_genes}.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
    
        return str(subject_file)
    
    def run_experiments(self):
        """Run complete experiment suite"""
        
        # Prepare data once
        plink_prefix = self.prepare_plink_data()
        
        # Check the data
        bim = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None)
        fam = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None)
        print(f"Data ready: {len(bim)} SNPs, {len(fam)} samples")
        
        # Create topology
        topology_file, topology_data = self.create_topology(plink_prefix)
        
        # Convert to GenNet format
        h5_file = self.convert_to_gennet(plink_prefix)
        
        # Parameter grid
        experiments = [
            {'h2': 0.2, 'n_genes': 50},
            {'h2': 0.5, 'n_genes': 50},
            {'h2': 0.8, 'n_genes': 50},
            {'h2': 0.5, 'n_genes': 10},
            {'h2': 0.5, 'n_genes': 100},
        ]
        
        results = []
        for i, params in enumerate(experiments):
            print(f"\n=== Experiment {i+1}/{len(experiments)} ===")
            print(f"Parameters: {params}")
            
            # Simulate phenotype
            subject_file = self.simulate_phenotype(
                plink_prefix, 
                topology_data,
                n_causal_genes=params['n_genes'],
                h2=params['h2'],
                seed=42 + i
            )
            
            results.append({
                'experiment': i,
                'h2': params['h2'],
                'n_genes': params['n_genes'],
                'subject_file': subject_file,
                'topology_file': topology_file,
                'genotype_file': h5_file,
                'plink_prefix': plink_prefix
            })
        
        # Save experiment summary
        summary_df = pd.DataFrame(results)
        summary_df.to_csv(self.output_dir / "experiment_summary.csv", index=False)
        
        print(f"\n=== Simulation setup complete! ===")
        print(f"Results in: {self.output_dir}")
        print(f"\nNext steps:")
        print(f"1. Check HDF5 output: {h5_file}")
        print(f"2. Review topology: {topology_file}")
        print(f"3. Subject files: subjects_h*_g*.csv")

def main():
    parser = argparse.ArgumentParser(description='Run GenNet simulation experiments')
    parser.add_argument('--vcf', required=True, help='1000 Genomes VCF file')
    parser.add_argument('--output-dir', default='gennet_simulation', help='Output directory')
    parser.add_argument('--n-samples', type=int, default=1000, help='Number of samples')
    
    args = parser.parse_args()
    
    # Run simulation
    sim = GenNetSimulation(args.vcf, args.output_dir)
    sim.run_experiments()

if __name__ == '__main__':
    main()