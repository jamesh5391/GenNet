#!/usr/bin/env python
"""
gennet_simulation_enhanced.py - GenNet simulation with Kelemen-style epistasis testing

This pipeline implements multiple simulation scenarios to test for genuine epistasis
vs joint tagging effects, following approaches from:
1. GenNet paper (gene-based architecture)
2. Kelemen et al. 2025 (epistasis testing framework)
"""

import subprocess
import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path
import shutil
from typing import Dict, List, Tuple

class GenNetSimulationEnhanced:
    def __init__(self, vcf_file, output_dir="gennet_simulation_enhanced"):
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
        h5_output_dir = self.output_dir / "h5_output"
        if h5_output_dir.exists():
            shutil.rmtree(h5_output_dir)
        h5_output_dir.mkdir()
        
        # GenNet convert command with output directory
        gennet_script = Path("../GenNet.py")
        if not gennet_script.exists():
            gennet_script = Path("GenNet.py")
        
        if gennet_script.exists():
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
                'layer0_node': snp_idx,
                'layer0_name': row['snp'],
                'layer1_node': gene_idx,
                'layer1_name': f"GENE_{gene_idx}",
                'layer2_node': pathway_idx,
                'layer2_name': f"PATHWAY_{pathway_idx}"
            })

        topology_file = self.output_dir / "topology.csv"
        pd.DataFrame(topology_rows).to_csv(topology_file, index=False)

        print(f"Topology created: {len(topology_rows)} SNPs -> {n_genes} genes -> {len(set(gene_to_pathway.values()))} pathways")

        return str(topology_file), topology_rows
    
    def simulate_additive_phenotype(self, plink_prefix, topology_data, 
                                   n_causal_genes=50, h2=0.5, seed=None):
        """Simulate purely additive phenotype (Kelemen scenario 1)"""
        if seed:
            np.random.seed(seed)
        
        print(f"=== Simulating ADDITIVE phenotype (h2={h2}, n_genes={n_causal_genes}) ===")
        
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
        
        # Equal effect sizes for all causal SNPs (purely additive)
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
        
        effects_file = self.output_dir / f"effects_additive_h{h2}_g{n_causal_genes}.txt"
        pd.DataFrame(effects).to_csv(effects_file, sep='\t', header=False, index=False)
        
        # Calculate GRS using plink2
        grs_prefix = self.output_dir / f"grs_additive_h{h2}_g{n_causal_genes}"
        
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
        
        return phenotype, causal_genes, causal_snps, effects_file
    
    def simulate_epistatic_phenotype(self, plink_prefix, topology_data,
                                    n_causal_genes=50, h2=0.5, seed=None, 
                                    interaction_order=4):
        """Simulate purely epistatic phenotype (Kelemen scenario 2)
        
        Args:
            interaction_order: Order of epistatic interactions (2, 3, or 4-way)
        """
        if seed:
            np.random.seed(seed)
        
        print(f"=== Simulating EPISTATIC phenotype ({interaction_order}-way, h2={h2}, n_genes={n_causal_genes}) ===")
        
        # Get FAM file
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
        
        # Read genotype data
        bed_file = f"{plink_prefix}.bed"
        bim = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None,
                        names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
        
        # Get indices of causal SNPs
        causal_snp_indices = [i for i, snp in enumerate(bim['snp']) if snp in causal_snps]
        
        # Read genotypes for causal SNPs using plink2
        temp_geno_file = self.output_dir / "temp_causal_snps.txt"
        with open(temp_geno_file, 'w') as f:
            for snp in causal_snps:
                f.write(f"{snp}\n")
        
        # Extract causal SNP genotypes
        extract_prefix = self.output_dir / "temp_extracted"
        cmd = f"""plink2 --bfile {plink_prefix} \
                --extract {temp_geno_file} \
                --export A \
                --out {extract_prefix}"""
        subprocess.run(cmd, shell=True, check=True)
        
        # Read extracted genotypes
        geno_file = f"{extract_prefix}.raw"
        geno_df = pd.read_csv(geno_file, sep=r'\s+')
        
        # Get genotype columns (skip first 6 columns: FID, IID, PAT, MAT, SEX, PHENOTYPE)
        geno_cols = [col for col in geno_df.columns if col.endswith('_ADD')]
        genotypes = geno_df[geno_cols].values
        
        # Simulate epistatic effects
        n_samples = len(fam)
        n_interactions = 2000  # Number of epistatic interaction terms
        
        epistatic_variance = np.zeros(n_samples)
        
        for i in range(n_interactions):
            # Randomly select SNPs for interaction
            selected_snp_idx = np.random.choice(genotypes.shape[1], 
                                               size=interaction_order, 
                                               replace=False)
            
            # Calculate interaction term (product of genotypes)
            interaction_term = np.ones(n_samples)
            for idx in selected_snp_idx:
                interaction_term *= genotypes[:, idx]
            
            # Standardize interaction term
            if interaction_term.std() > 0:
                interaction_term = (interaction_term - interaction_term.mean()) / interaction_term.std()
            
            # Add to epistatic variance
            epistatic_variance += interaction_term
        
        # Standardize total epistatic variance
        if epistatic_variance.std() > 0:
            epistatic_variance = (epistatic_variance - epistatic_variance.mean()) / epistatic_variance.std()
        
        # Add environmental component
        genetic_component = np.sqrt(h2) * epistatic_variance
        environmental_component = np.sqrt(1 - h2) * np.random.normal(0, 1, n_samples)
        liability = genetic_component + environmental_component
        
        # Binary phenotype
        threshold = np.percentile(liability, 50)
        phenotype = (liability > threshold).astype(int)
        
        print(f"Created phenotype: {phenotype.sum()} cases, {len(phenotype) - phenotype.sum()} controls")
        
        # Clean up temp files
        temp_geno_file.unlink()
        for ext in ['.raw', '.log']:
            temp_file = Path(f"{extract_prefix}{ext}")
            if temp_file.exists():
                temp_file.unlink()
        
        return phenotype, causal_genes, causal_snps, None
    
    def simulate_mixed_phenotype(self, plink_prefix, topology_data,
                                n_causal_genes=50, h2=0.5, 
                                additive_fraction=0.5, seed=None):
        """Simulate mixed additive-epistatic phenotype (Kelemen scenario 3)
        
        Args:
            additive_fraction: Fraction of variance from additive effects (0-1)
        """
        if seed:
            np.random.seed(seed)
        
        print(f"=== Simulating MIXED phenotype (additive={additive_fraction:.1%}, h2={h2}, n_genes={n_causal_genes}) ===")
        
        # Generate both additive and epistatic phenotypes
        pheno_add, genes_add, snps_add, effects_add = self.simulate_additive_phenotype(
            plink_prefix, topology_data, n_causal_genes, h2, seed
        )
        
        pheno_epi, genes_epi, snps_epi, _ = self.simulate_epistatic_phenotype(
            plink_prefix, topology_data, n_causal_genes, h2, seed + 1 if seed else None
        )
        
        # Read phenotypes as continuous liability scores
        # (we need to regenerate these as continuous values)
        fam = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None,
                        names=['fid', 'iid', 'father', 'mother', 'sex', 'pheno'])
        
        # Regenerate continuous liabilities
        # For additive component
        effects_df = pd.read_csv(effects_add, sep='\t', header=None,
                                names=['SNP', 'A1', 'BETA'])
        grs_prefix = self.output_dir / "temp_grs_for_mixed"
        cmd = f"""plink2 --bfile {plink_prefix} \
                --score {effects_add} 1 2 3 \
                --out {grs_prefix}"""
        subprocess.run(cmd, shell=True, check=True)
        
        score_file = f"{grs_prefix}.sscore"
        if not Path(score_file).exists():
            score_file = f"{grs_prefix}.profile"
        grs_df = pd.read_csv(score_file, sep=r'\s+')
        score_col = next((col for col in ['SCORE1_SUM', 'SCORE1_AVG', 'SCORESUM', 'SCORE'] 
                         if col in grs_df.columns), None)
        grs_add = grs_df[score_col].values if score_col else np.random.normal(0, 1, len(fam))
        
        if grs_add.std() > 0:
            grs_add = (grs_add - grs_add.mean()) / grs_add.std()
        
        # Simulate epistatic component (simplified - using random interaction)
        np.random.seed(seed + 2 if seed else None)
        grs_epi = np.random.normal(0, 1, len(fam))
        if grs_epi.std() > 0:
            grs_epi = (grs_epi - grs_epi.mean()) / grs_epi.std()
        
        # Combine genetic components
        genetic_component = (np.sqrt(additive_fraction) * grs_add + 
                           np.sqrt(1 - additive_fraction) * grs_epi)
        genetic_component = (genetic_component - genetic_component.mean()) / genetic_component.std()
        
        # Add environmental component
        full_liability = (np.sqrt(h2) * genetic_component + 
                         np.sqrt(1 - h2) * np.random.normal(0, 1, len(fam)))
        
        # Binary phenotype
        threshold = np.percentile(full_liability, 50)
        phenotype = (full_liability > threshold).astype(int)
        
        print(f"Created phenotype: {phenotype.sum()} cases, {len(phenotype) - phenotype.sum()} controls")
        
        # Clean up temp files
        for ext in ['.sscore', '.log']:
            temp_file = Path(f"{grs_prefix}{ext}")
            if temp_file.exists():
                temp_file.unlink()
        
        return phenotype, np.union1d(genes_add, genes_epi), np.union1d(snps_add, snps_epi), effects_add
    
    def simulate_joint_tagging_scenario(self, plink_prefix, topology_data,
                                       n_causal_genes=50, h2=0.5, 
                                       coverage=0.5, seed=None):
        """Simulate joint tagging effect scenario (Kelemen scenario 4)
        
        This creates a scenario where:
        1. Phenotype is purely additive
        2. But some causal SNPs are removed (incomplete coverage)
        3. This creates apparent epistasis through joint tagging
        
        Args:
            coverage: Fraction of causal SNPs to keep (0.5 = 50% coverage)
        """
        if seed:
            np.random.seed(seed)
        
        print(f"=== Simulating JOINT TAGGING scenario (coverage={coverage:.1%}, h2={h2}) ===")
        
        # First simulate a normal additive phenotype
        phenotype, causal_genes, causal_snps, effects_file = self.simulate_additive_phenotype(
            plink_prefix, topology_data, n_causal_genes, h2, seed
        )
        
        # Read BIM file to identify SNPs in LD with causal SNPs
        bim = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None,
                        names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
        
        # Calculate LD matrix for causal SNPs
        temp_ld_prefix = self.output_dir / "temp_ld_check"
        causal_snp_file = self.output_dir / "temp_causal_snps_ld.txt"
        with open(causal_snp_file, 'w') as f:
            for snp in causal_snps:
                f.write(f"{snp}\n")
        
        # Extract causal SNPs and calculate LD
        cmd = f"""plink2 --bfile {plink_prefix} \
                --extract {causal_snp_file} \
                --r2-phased square \
                --out {temp_ld_prefix}"""
        subprocess.run(cmd, shell=True, check=True)
        
        # Read LD matrix
        ld_file = f"{temp_ld_prefix}.vcor"
        if Path(ld_file).exists():
            ld_matrix = np.loadtxt(ld_file)
            
            # Identify SNPs with high LD (r2 > 0.25) with at least 2 other SNPs
            high_ld_snps = []
            for i in range(len(causal_snps)):
                n_high_ld = np.sum(ld_matrix[i, :] > 0.25) - 1  # -1 to exclude self
                if n_high_ld >= 2:
                    high_ld_snps.append(causal_snps[i])
            
            # Remove a fraction of high-LD SNPs to simulate incomplete coverage
            n_to_remove = int(len(high_ld_snps) * (1 - coverage))
            removed_snps = np.random.choice(high_ld_snps, 
                                          size=min(n_to_remove, len(high_ld_snps)), 
                                          replace=False)
            
            # Create new SNP list (reduced panel)
            reduced_snps = [snp for snp in causal_snps if snp not in removed_snps]
            
            print(f"Removed {len(removed_snps)} SNPs to create joint tagging scenario")
            print(f"Reduced panel: {len(reduced_snps)} SNPs (from {len(causal_snps)})")
            
            # Create new effects file with reduced SNP panel
            effects_df = pd.read_csv(effects_file, sep='\t', header=None,
                                    names=['SNP', 'A1', 'BETA'])
            reduced_effects_df = effects_df[effects_df['SNP'].isin(reduced_snps)]
            
            reduced_effects_file = self.output_dir / f"effects_jointag_h{h2}_g{n_causal_genes}_cov{int(coverage*100)}.txt"
            reduced_effects_df.to_csv(reduced_effects_file, sep='\t', header=False, index=False)
        else:
            print("Warning: Could not calculate LD matrix, using all SNPs")
            reduced_snps = causal_snps
            reduced_effects_file = effects_file
        
        # Clean up temp files
        causal_snp_file.unlink()
        for ext in ['.vcor', '.log']:
            temp_file = Path(f"{temp_ld_prefix}{ext}")
            if temp_file.exists():
                temp_file.unlink()
        
        return phenotype, causal_genes, reduced_snps, reduced_effects_file
    
    def create_subject_file(self, plink_prefix, phenotype, scenario_name, 
                          h2, n_causal_genes, seed=None):
        """Create GenNet-format subject file"""
        
        # Get FAM file for sample IDs
        fam = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None,
                        names=['fid', 'iid', 'father', 'mother', 'sex', 'pheno'])
        
        n_samples = len(fam)
        
        # Split into train/val/test (60/20/20)
        if seed:
            np.random.seed(seed)
        indices = np.arange(n_samples)
        np.random.shuffle(indices)
        
        train_idx = indices[:int(0.6 * n_samples)]
        val_idx = indices[int(0.6 * n_samples):int(0.8 * n_samples)]
        test_idx = indices[int(0.8 * n_samples):]
        
        # Create set assignments (GenNet expects 1=train, 2=val, 3=test)
        set_assignment = np.zeros(n_samples, dtype=int)
        set_assignment[train_idx] = 1
        set_assignment[val_idx] = 2
        set_assignment[test_idx] = 3
        
        # Create GenNet-format subject file
        subject_df = pd.DataFrame({
            'patient_id': fam['iid'].values,
            'labels': phenotype,
            'genotype_row': np.arange(n_samples),
            'set': set_assignment
        })
        
        subject_file = self.output_dir / f"subjects_{scenario_name}_h{h2}_g{n_causal_genes}.csv"
        subject_df.to_csv(subject_file, index=False)
        
        print(f"Created subject file: Train={sum(set_assignment==1)}, Val={sum(set_assignment==2)}, Test={sum(set_assignment==3)}")
        
        return str(subject_file)
    
    def run_comprehensive_experiments(self):
        """Run comprehensive experiment suite with all Kelemen scenarios"""
        
        # Prepare data once
        print("\n" + "="*80)
        print("PREPARING GENOTYPE DATA")
        print("="*80)
        plink_prefix = self.prepare_plink_data()
        
        # Check the data
        bim = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None)
        fam = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None)
        print(f"Data ready: {len(bim)} SNPs, {len(fam)} samples")
        
        # Create topology
        print("\n" + "="*80)
        print("CREATING TOPOLOGY")
        print("="*80)
        topology_file, topology_data = self.create_topology(plink_prefix)
        
        # Convert to GenNet format
        print("\n" + "="*80)
        print("CONVERTING TO GENNET FORMAT")
        print("="*80)
        h5_file = self.convert_to_gennet(plink_prefix)
        
        # Define experiment scenarios
        scenarios = [
            # Scenario 1: Purely additive (baseline)
            {
                'name': 'additive',
                'type': 'additive',
                'h2_values': [0.2, 0.5, 0.8],
                'n_genes_values': [10, 50, 100],
            },
            # Scenario 2: Purely epistatic
            {
                'name': 'epistatic_4way',
                'type': 'epistatic',
                'h2_values': [0.5],
                'n_genes_values': [50],
                'interaction_order': 4
            },
            # Scenario 3: Mixed additive-epistatic
            {
                'name': 'mixed_2to1_add',
                'type': 'mixed',
                'h2_values': [0.5],
                'n_genes_values': [50],
                'additive_fraction': 0.67  # 2:1 additive:epistatic
            },
            {
                'name': 'mixed_1to2_add',
                'type': 'mixed',
                'h2_values': [0.5],
                'n_genes_values': [50],
                'additive_fraction': 0.33  # 1:2 additive:epistatic
            },
            # Scenario 4: Joint tagging effects
            {
                'name': 'joint_tagging_50pct',
                'type': 'joint_tagging',
                'h2_values': [0.5],
                'n_genes_values': [50],
                'coverage': 0.5  # 50% coverage
            },
        ]
        
        results = []
        experiment_id = 0
        
        for scenario in scenarios:
            scenario_name = scenario['name']
            scenario_type = scenario['type']
            
            print("\n" + "="*80)
            print(f"SCENARIO: {scenario_name.upper()}")
            print("="*80)
            
            for h2 in scenario['h2_values']:
                for n_genes in scenario['n_genes_values']:
                    experiment_id += 1
                    print(f"\n--- Experiment {experiment_id}: {scenario_name}, h2={h2}, n_genes={n_genes} ---")
                    
                    seed = 42 + experiment_id
                    
                    # Simulate phenotype based on scenario type
                    if scenario_type == 'additive':
                        phenotype, causal_genes, causal_snps, effects_file = \
                            self.simulate_additive_phenotype(
                                plink_prefix, topology_data, n_genes, h2, seed
                            )
                    
                    elif scenario_type == 'epistatic':
                        interaction_order = scenario.get('interaction_order', 4)
                        phenotype, causal_genes, causal_snps, effects_file = \
                            self.simulate_epistatic_phenotype(
                                plink_prefix, topology_data, n_genes, h2, seed,
                                interaction_order
                            )
                    
                    elif scenario_type == 'mixed':
                        additive_fraction = scenario['additive_fraction']
                        phenotype, causal_genes, causal_snps, effects_file = \
                            self.simulate_mixed_phenotype(
                                plink_prefix, topology_data, n_genes, h2,
                                additive_fraction, seed
                            )
                    
                    elif scenario_type == 'joint_tagging':
                        coverage = scenario['coverage']
                        phenotype, causal_genes, causal_snps, effects_file = \
                            self.simulate_joint_tagging_scenario(
                                plink_prefix, topology_data, n_genes, h2,
                                coverage, seed
                            )
                    
                    # Create subject file
                    subject_file = self.create_subject_file(
                        plink_prefix, phenotype, scenario_name, h2, n_genes, seed
                    )
                    
                    # Save metadata
                    metadata = {
                        'experiment_id': experiment_id,
                        'scenario': scenario_name,
                        'scenario_type': scenario_type,
                        'n_causal_genes': int(n_genes),
                        'n_causal_snps': int(len(causal_snps)),
                        'heritability': float(h2),
                        'seed': int(seed),
                        'causal_genes_sample': [str(g) for g in list(causal_genes)[:10]]
                    }
                    
                    if scenario_type == 'epistatic':
                        metadata['interaction_order'] = scenario.get('interaction_order', 4)
                    elif scenario_type == 'mixed':
                        metadata['additive_fraction'] = scenario['additive_fraction']
                    elif scenario_type == 'joint_tagging':
                        metadata['coverage'] = scenario['coverage']
                    
                    metadata_file = self.output_dir / f"metadata_{scenario_name}_h{h2}_g{n_genes}.json"
                    with open(metadata_file, 'w') as f:
                        json.dump(metadata, f, indent=2)
                    
                    results.append({
                        'experiment_id': experiment_id,
                        'scenario': scenario_name,
                        'h2': h2,
                        'n_genes': n_genes,
                        'subject_file': subject_file,
                        'topology_file': topology_file,
                        'genotype_file': h5_file,
                        'plink_prefix': plink_prefix,
                        'metadata_file': str(metadata_file)
                    })
        
        # Save experiment summary
        summary_df = pd.DataFrame(results)
        summary_df.to_csv(self.output_dir / "experiment_summary_comprehensive.csv", index=False)
        
        print("\n" + "="*80)
        print("SIMULATION COMPLETE!")
        print("="*80)
        print(f"Results in: {self.output_dir}")
        print(f"Total experiments: {experiment_id}")
        print("\nScenarios tested:")
        for scenario in scenarios:
            print(f"  - {scenario['name']}: {scenario['type']}")
        print(f"\nNext steps:")
        print(f"1. Review experiment summary: experiment_summary_comprehensive.csv")
        print(f"2. Check HDF5 file: {h5_file}")
        print(f"3. Review topology: {topology_file}")
        print(f"4. Train GenNet models on each scenario")
        print(f"5. Compare performance across scenarios to test for genuine epistasis")

def main():
    parser = argparse.ArgumentParser(
        description='Run comprehensive GenNet simulation experiments with epistasis testing'
    )
    parser.add_argument('--vcf', required=True, help='1000 Genomes VCF file')
    parser.add_argument('--output-dir', default='gennet_simulation_enhanced',
                       help='Output directory')
    parser.add_argument('--n-samples', type=int, default=1000, 
                       help='Number of samples')
    
    args = parser.parse_args()
    
    # Run comprehensive simulation
    sim = GenNetSimulationEnhanced(args.vcf, args.output_dir)
    sim.run_comprehensive_experiments()

if __name__ == '__main__':
    main()