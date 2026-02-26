"""
Deep dive into why additive and epistatic components are highly correlated.

This demonstrates:
1. Mathematical relationship between additive and epistatic terms
2. Why sampling from the same SNP set creates correlation
3. How LD structure amplifies the correlation
4. Concrete numerical examples
"""

import numpy as np
import pandas as pd
from pathlib import Path
import tables
import matplotlib.pyplot as plt

# Load data
base_dir = Path(__file__).parent.parent
h5_file = base_dir / 'data/processed/simulationx/h5_output/genotype.h5'
h5 = tables.open_file(str(h5_file), "r")
genotype_matrix = h5.root.data[:]
h5.close()

n_samples, n_snps = genotype_matrix.shape

# Load gene annotations
plink_prefix = base_dir / 'data/processed/simulationx/plink_clean/chr22_data'
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


# Now show what happens in your actual simulation
print("\n" + "="*80)
print("PART 2: YOUR ACTUAL SIMULATION")
print("="*80)

def simulate_with_analysis(P=10, N=10, n_interactions=2000):
    """Simulate and analyze correlation sources."""
    np.random.seed(42)

    # Sample causal genes and SNPs
    unique_genes = gene_annotations['gene_id'].unique()
    causal_genes = np.random.choice(unique_genes, size=P, replace=False)

    causal_snps = []
    snp_to_gene = {}
    for gene in causal_genes:
        gene_snps = gene_annotations[gene_annotations['gene_id'] == gene]['snp_id'].values
        if len(gene_snps) >= N:
            selected = np.random.choice(gene_snps, size=N, replace=False)
        else:
            selected = gene_snps

        for snp in selected:
            snp_to_gene[snp] = gene

        causal_snps.extend(selected)

    causal_snps = np.array(causal_snps, dtype=int)

    # Compute additive component
    effect_size = 2.0 / (N * P)
    effects = np.zeros(n_snps)
    effects[causal_snps] = effect_size

    G_add = genotype_matrix @ effects
    G_add_std = (G_add - G_add.mean()) / (G_add.std() + 1e-10)

    # Compute epistatic component (with tracking)
    epistatic_signal = np.zeros(n_samples)

    within_gene_interactions = []
    cross_gene_interactions = []

    for i in range(n_interactions):
        selected_snp_idx = np.random.choice(causal_snps, size=2, replace=False)

        # Calculate interaction
        interaction = genotype_matrix[:, selected_snp_idx[0]] * genotype_matrix[:, selected_snp_idx[1]]
        interaction_std = (interaction - interaction.mean()) / (interaction.std() + 1e-10)

        # Track within vs cross gene
        genes = set([snp_to_gene[idx] for idx in selected_snp_idx])
        if len(genes) == 1:
            within_gene_interactions.append(interaction_std)
        else:
            cross_gene_interactions.append(interaction_std)

        epistatic_signal += interaction_std

    G_epi_std = (epistatic_signal - epistatic_signal.mean()) / (epistatic_signal.std() + 1e-10)

    # Analyze correlation
    total_corr = np.corrcoef(G_add_std, G_epi_std)[0, 1]

    # Correlation from within-gene interactions
    if len(within_gene_interactions) > 0:
        within_gene_sum = np.sum(within_gene_interactions, axis=0)
        within_gene_sum_std = (within_gene_sum - within_gene_sum.mean()) / (within_gene_sum.std() + 1e-10)
        within_corr = np.corrcoef(G_add_std, within_gene_sum_std)[0, 1]
    else:
        within_corr = np.nan

    # Correlation from cross-gene interactions
    if len(cross_gene_interactions) > 0:
        cross_gene_sum = np.sum(cross_gene_interactions, axis=0)
        cross_gene_sum_std = (cross_gene_sum - cross_gene_sum.mean()) / (cross_gene_sum.std() + 1e-10)
        cross_corr = np.corrcoef(G_add_std, cross_gene_sum_std)[0, 1]
    else:
        cross_corr = np.nan

    # Compute average LD among causal SNPs
    causal_geno = genotype_matrix[:, causal_snps]
    # Standardize
    for i in range(causal_geno.shape[1]):
        col = causal_geno[:, i]
        causal_geno[:, i] = (col - col.mean()) / (col.std() + 1e-10)

    ld_matrix = np.corrcoef(causal_geno.T)
    # Get upper triangle (exclude diagonal)
    upper_tri_indices = np.triu_indices_from(ld_matrix, k=1)
    ld_values = ld_matrix[upper_tri_indices]
    mean_ld = np.mean(np.abs(ld_values))

    return {
        'P': P,
        'n_causal_snps': len(causal_snps),
        'n_interactions': n_interactions,
        'n_within_gene': len(within_gene_interactions),
        'n_cross_gene': len(cross_gene_interactions),
        'total_corr': total_corr,
        'within_corr': within_corr,
        'cross_corr': cross_corr,
        'mean_ld': mean_ld
    }


# Run analysis for different P values
print("\nAnalyzing correlation breakdown by polygenicity (P):\n")

configs = [
    {'P': 1, 'N': 10},
    {'P': 10, 'N': 10},
    {'P': 100, 'N': 10},
    {'P': 500, 'N': 10}
]

results = []
for config in configs:
    result = simulate_with_analysis(P=config['P'], N=config['N'])
    results.append(result)

    print(f"P = {result['P']:3d} ({result['n_causal_snps']:4d} causal SNPs):")
    print(f"  Total correlation:        {result['total_corr']:.4f}")
    print(f"  Within-gene contribution: {result['within_corr']:.4f} ({result['n_within_gene']} interactions)")
    print(f"  Cross-gene contribution:  {result['cross_corr']:.4f} ({result['n_cross_gene']} interactions)")
    print(f"  Mean |LD| among SNPs:     {result['mean_ld']:.4f}")
    print()
