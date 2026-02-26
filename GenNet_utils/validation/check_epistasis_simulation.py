"""
Diagnostic script to verify epistatic simulation is working correctly.

This checks:
1. Whether epistatic interactions are truly cross-gene or within-gene
2. Whether the epistatic component has meaningful signal
3. Whether cross-gene interactions can be learned by GenNet architecture
"""

import numpy as np
import pandas as pd
from pathlib import Path
import tables

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


def simulate_phenotype_diagnostic(geno_matrix, gene_annotations, P, h2, alpha,
                                   N=10, interaction_order=2, n_interactions=2000,
                                   seed=42):
    """
    Modified simulation function with diagnostic output.
    """
    np.random.seed(seed)
    n_samples, n_snps = geno_matrix.shape

    # ========== ADDITIVE COMPONENT ==========
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

    effect_size = 2.0 / (N * P)
    effects = np.zeros(n_snps)
    effects[causal_snps_additive] = effect_size

    G_add = geno_matrix @ effects
    G_add_std = (G_add - G_add.mean()) / (G_add.std() + 1e-10)

    # ========== EPISTATIC COMPONENT (WITH DIAGNOSTICS) ==========
    epistatic_signal = np.zeros(n_samples)

    # Track which interactions are within-gene vs cross-gene
    within_gene_count = 0
    cross_gene_count = 0

    # Map SNPs to genes for diagnostic
    snp_to_gene = {row['snp_id']: row['gene_id'] for _, row in gene_annotations.iterrows()}

    for i in range(n_interactions):
        if len(causal_snps_additive) >= interaction_order:
            selected_snp_idx = np.random.choice(
                causal_snps_additive,
                size=interaction_order,
                replace=False
            )
        else:
            selected_snp_idx = causal_snps_additive

        # Check if interaction is within-gene or cross-gene
        genes_in_interaction = set([snp_to_gene[idx] for idx in selected_snp_idx])
        if len(genes_in_interaction) == 1:
            within_gene_count += 1
        else:
            cross_gene_count += 1

        # Calculate interaction term
        interaction_term = np.ones(n_samples)
        for idx in selected_snp_idx:
            interaction_term *= geno_matrix[:, idx]

        if interaction_term.std() > 0:
            interaction_term = (interaction_term - interaction_term.mean()) / interaction_term.std()

        epistatic_signal += interaction_term

    G_epi_std = (epistatic_signal - epistatic_signal.mean()) / (epistatic_signal.std() + 1e-10)

    # ========== COMBINE ==========
    G_combined = np.sqrt(alpha) * G_epi_std + np.sqrt(1 - alpha) * G_add_std
    G_combined_std = (G_combined - G_combined.mean()) / (G_combined.std() + 1e-10)

    genetic_component = np.sqrt(h2) * G_combined_std
    environmental_component = np.sqrt(1 - h2) * np.random.randn(n_samples)
    phenotype = genetic_component + environmental_component

    # ========== DIAGNOSTICS ==========
    diagnostics = {
        'P': P,
        'h2': h2,
        'alpha': alpha,
        'n_causal_genes': len(causal_genes),
        'n_causal_snps': len(causal_snps_additive),
        'n_interactions_total': n_interactions,
        'n_interactions_within_gene': within_gene_count,
        'n_interactions_cross_gene': cross_gene_count,
        'pct_cross_gene': 100 * cross_gene_count / n_interactions,
        'G_add_var': G_add_std.var(),
        'G_epi_var': G_epi_std.var(),
        'G_combined_var': G_combined_std.var(),
        'correlation_add_epi': np.corrcoef(G_add_std, G_epi_std)[0, 1],
        'phenotype_mean': phenotype.mean(),
        'phenotype_var': phenotype.var(),
        'genetic_component_var': genetic_component.var(),
        'environmental_component_var': environmental_component.var()
    }

    return phenotype, diagnostics


# Run diagnostics for a few key configurations
print("="*80)
print("EPISTATIC SIMULATION DIAGNOSTICS")
print("="*80)

test_configs = [
    {'P': 10, 'h2': 0.6, 'alpha': 0.0},
    {'P': 10, 'h2': 0.6, 'alpha': 0.6},
    {'P': 10, 'h2': 0.6, 'alpha': 1.0},
    {'P': 100, 'h2': 0.6, 'alpha': 0.0},
    {'P': 100, 'h2': 0.6, 'alpha': 1.0},
]

for config in test_configs:
    phenotype, diag = simulate_phenotype_diagnostic(
        genotype_matrix, gene_annotations,
        P=config['P'], h2=config['h2'], alpha=config['alpha'],
        seed=42
    )

    print(f"\n--- Config: P={diag['P']}, h²={diag['h2']}, α={diag['alpha']} ---")
    print(f"Causal genes: {diag['n_causal_genes']}")
    print(f"Causal SNPs: {diag['n_causal_snps']}")
    print(f"\nInteraction breakdown:")
    print(f"  Within-gene interactions: {diag['n_interactions_within_gene']} ({100-diag['pct_cross_gene']:.1f}%)")
    print(f"  Cross-gene interactions:  {diag['n_interactions_cross_gene']} ({diag['pct_cross_gene']:.1f}%)")
    print(f"\nGenetic component variances:")
    print(f"  Var(G_add):      {diag['G_add_var']:.4f}")
    print(f"  Var(G_epi):      {diag['G_epi_var']:.4f}")
    print(f"  Var(G_combined): {diag['G_combined_var']:.4f}")
    print(f"  Corr(G_add, G_epi): {diag['correlation_add_epi']:.4f}")
    print(f"\nPhenotype statistics:")
    print(f"  Mean: {diag['phenotype_mean']:.4f}")
    print(f"  Var:  {diag['phenotype_var']:.4f}")
    print(f"  Var(genetic):      {diag['genetic_component_var']:.4f}")
    print(f"  Var(environment):  {diag['environmental_component_var']:.4f}")

print("\n" + "="*80)
print("KEY FINDINGS:")
print("="*80)
print("\n1. CROSS-GENE INTERACTION PERCENTAGE")
print("   - If ~90-99% are cross-gene → GenNet CANNOT learn them!")
print("   - This explains why ReLU ≈ Linear")

print("\n2. CORRELATION BETWEEN ADDITIVE AND EPISTATIC")
print("   - If correlation is high → epistatic signal is redundant with additive")
print("   - This could explain why α has no effect")

print("\n3. EXPECTED PROBABILITY OF WITHIN-GENE INTERACTION")
# For P=10 genes out of 1500 total genes:
# Probability that 2 randomly sampled SNPs are from same gene:
# = (P * N * (N-1)) / (P*N * (P*N-1))
# ≈ (N-1) / (P*N - 1)
P = 10
N = 10
prob_within = (N - 1) / (P * N - 1)
print(f"   - For P={P}, N={N}: Expected within-gene rate ≈ {100*prob_within:.1f}%")
print(f"   - Remaining {100*(1-prob_within):.1f}% are cross-gene (unlearnable!)")

print("\n" + "="*80)
