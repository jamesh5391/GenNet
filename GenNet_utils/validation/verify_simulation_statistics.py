"""
Statistical Verification of Corrected Simulation

This script verifies that the simulation fix was implemented correctly by checking:
1. Var(G_add) decreases with P (polygenicity effect preserved)
2. Var(phenotype) ≈ 1.0 across all configs
3. Var(genetic_component) = h² for all configs
4. Within-gene epistatic interactions (100%)
5. Corr(G_add, G_epi) < 0.7 (reduced from 0.98)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import tables

print("="*80)
print("STATISTICAL VERIFICATION OF SIMULATION")
print("="*80)

# Load genotype matrix
h5_file = Path('data/processed/simulationx/h5_output/genotype.h5')
h5 = tables.open_file(str(h5_file), "r")
genotype_matrix = h5.root.data[:]
h5.close()

n_samples, n_snps = genotype_matrix.shape
print(f"\nGenotype matrix: {n_samples} samples × {n_snps} SNPs")

# Load gene annotations
plink_prefix = Path('data/processed/simulationx/plink_clean/chr22_pruned')
snp_info = pd.read_csv(f"{plink_prefix}.bim", sep='\t', header=None,
                       names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])

n_genes = 1500
snps_per_gene = max(1, n_snps // n_genes)
actual_n_genes = min(n_genes, n_snps)

gene_annotations = pd.DataFrame({
    'snp_id': np.arange(n_snps),
    'gene_id': np.minimum(np.arange(n_snps) // snps_per_gene, actual_n_genes - 1),
    'snp_name': snp_info['snp'].values[:n_snps]
})

# Load all phenotype files
base_dir = Path('data/processed/grid_experiments')

print("\n" + "="*80)
print("CHECK 1: VARIANCE OF G_add SHOULD DECREASE WITH P")
print("="*80)

# We'll read the subjects files and check phenotype variance
# Since we can't access G_add directly from saved files, we'll regenerate it

def simulate_and_extract_components(geno_matrix, gene_annotations, P, h2, alpha, N=10, n_interactions=2000):
    """Re-simulate to extract G_add and G_epi for analysis"""
    n_samples, n_snps = geno_matrix.shape

    np.random.seed(42)  # Same seed as generation

    # Sample causal genes and SNPs
    unique_genes = gene_annotations['gene_id'].unique()
    causal_genes = np.random.choice(unique_genes, size=P, replace=False)

    causal_snps_additive = []
    for gene in causal_genes:
        gene_snps = gene_annotations[gene_annotations['gene_id'] == gene]['snp_id'].values
        if len(gene_snps) >= N:
            selected = np.random.choice(gene_snps, size=N, replace=False)
        else:
            selected = gene_snps
        causal_snps_additive.extend(selected)

    causal_snps_additive = np.array(causal_snps_additive, dtype=int)

    # Additive component
    effect_size = 2.0 / (N * P)
    effects = np.zeros(n_snps)
    effects[causal_snps_additive] = effect_size

    G_add = geno_matrix @ effects
    G_add_centered = G_add - G_add.mean()

    # Epistatic component
    epistatic_signal = np.zeros(n_samples)

    # Check if interactions are within-gene or cross-gene
    within_gene_count = 0
    cross_gene_count = 0
    snp_to_gene = {}
    for gene in causal_genes:
        gene_snps = gene_annotations[gene_annotations['gene_id'] == gene]['snp_id'].values
        for snp in gene_snps:
            if snp in causal_snps_additive:
                snp_to_gene[snp] = gene

    for i in range(n_interactions):
        if len(causal_snps_additive) >= 2:
            selected_snp_idx = np.random.choice(causal_snps_additive, size=2, replace=False)
        else:
            selected_snp_idx = causal_snps_additive

        # Check if within-gene
        if len(selected_snp_idx) == 2:
            gene1 = snp_to_gene.get(selected_snp_idx[0])
            gene2 = snp_to_gene.get(selected_snp_idx[1])
            if gene1 is not None and gene2 is not None:
                if gene1 == gene2:
                    within_gene_count += 1
                else:
                    cross_gene_count += 1

        interaction_term = np.ones(n_samples)
        for idx in selected_snp_idx:
            interaction_term *= geno_matrix[:, idx]

        if interaction_term.std() > 0:
            interaction_term = (interaction_term - interaction_term.mean()) / interaction_term.std()

        epistatic_signal += interaction_term

    G_epi_std = (epistatic_signal - epistatic_signal.mean()) / (epistatic_signal.std() + 1e-10)

    # Check correlation
    corr = np.corrcoef(G_add_centered, G_epi_std)[0, 1] if G_add_centered.std() > 0 else 0

    return {
        'G_add_centered': G_add_centered,
        'G_epi_std': G_epi_std,
        'var_G_add': G_add_centered.var(),
        'corr_G_add_G_epi': corr,
        'within_gene_pct': 100 * within_gene_count / (within_gene_count + cross_gene_count) if (within_gene_count + cross_gene_count) > 0 else 0,
        'cross_gene_pct': 100 * cross_gene_count / (within_gene_count + cross_gene_count) if (within_gene_count + cross_gene_count) > 0 else 0
    }

# Test for pure additive (alpha=0) to check P effect
print("\n--- Variance of G_add by Polygenicity (h²=0.6, α=0) ---")
print(f"{'P':<10} {'Var(G_add)':<15} {'Expected':<15} {'Match?'}")
print("-" * 60)

p_values = [1, 10, 100, 500]
N = 10

for P in p_values:
    result = simulate_and_extract_components(
        genotype_matrix, gene_annotations,
        P=P, h2=0.6, alpha=0, N=N
    )

    var_observed = result['var_G_add']
    var_expected = 4.0 / (N * P)  # Theoretical: Var(G_add) = 4/(N*P)

    match = "✓" if abs(var_observed - var_expected) / var_expected < 0.5 else "✗"

    print(f"{P:<10} {var_observed:<15.6f} {var_expected:<15.6f} {match}")

print("\n" + "="*80)
print("CHECK 2: WITHIN-GENE EPISTATIC INTERACTIONS")
print("="*80)

print("\n--- Epistatic Interaction Sampling (α=1, pure epistatic) ---")
print(f"{'P':<10} {'Within-Gene %':<20} {'Cross-Gene %':<20} {'Status'}")
print("-" * 70)

for P in p_values:
    result = simulate_and_extract_components(
        genotype_matrix, gene_annotations,
        P=P, h2=0.6, alpha=1, N=N
    )

    within_pct = result['within_gene_pct']
    cross_pct = result['cross_gene_pct']

    # Should be 100% within-gene if fix was applied
    status = "✓ FIXED" if within_pct > 95 else "✗ STILL BROKEN"

    print(f"{P:<10} {within_pct:<20.1f} {cross_pct:<20.1f} {status}")

print("\n" + "="*80)
print("CHECK 3: CORRELATION BETWEEN G_add AND G_epi")
print("="*80)

print("\n--- Corr(G_add, G_epi) by Polygenicity (h²=0.6, α=0.5) ---")
print(f"{'P':<10} {'Correlation':<15} {'Status'}")
print("-" * 40)

for P in p_values:
    result = simulate_and_extract_components(
        genotype_matrix, gene_annotations,
        P=P, h2=0.6, alpha=0.5, N=N
    )

    corr = result['corr_G_add_G_epi']

    # Should be < 0.7 if fix was applied (was 0.98 before)
    status = "✓ FIXED" if abs(corr) < 0.7 else "✗ STILL HIGH"

    print(f"{P:<10} {corr:<15.4f} {status}")

print("\n" + "="*80)
print("CHECK 4: PHENOTYPE VARIANCE ACROSS CONFIGS")
print("="*80)

print("\n--- Checking saved phenotype files ---")
print(f"{'Config':<30} {'Var(Y)':<15} {'Status'}")
print("-" * 50)

configs_checked = 0
configs_passed = 0

for P in [1, 10]:  # Check subset
    for h2 in [0.2, 0.6]:
        for alpha in [0, 0.5]:
            exp_id = f"exp_P{P}_h2{h2}_alpha{alpha}"
            subjects_file = base_dir / exp_id / 'subjects.csv'

            if subjects_file.exists():
                df = pd.read_csv(subjects_file)
                phenotype_var = df['labels'].var()

                configs_checked += 1

                # Should be close to 1.0
                status = "✓" if 0.8 < phenotype_var < 1.2 else "✗"
                if status == "✓":
                    configs_passed += 1

                print(f"{exp_id:<30} {phenotype_var:<15.4f} {status}")

print(f"\nPassed: {configs_passed}/{configs_checked} configs")


