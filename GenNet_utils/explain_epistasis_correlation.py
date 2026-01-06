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


print("="*80)
print("WHY IS THE CORRELATION BETWEEN ADDITIVE AND EPISTATIC SO HIGH?")
print("="*80)

# Simple example with toy data
print("\n" + "="*80)
print("PART 1: TOY EXAMPLE - Why Products Correlate with Sums")
print("="*80)

print("\nConsider a simple case with 4 SNPs and 100 individuals:")
print("We'll use REAL genotype data to show realistic LD patterns\n")

np.random.seed(42)

# Select 4 SNPs from real data
example_snp_indices = np.random.choice(n_snps, size=4, replace=False)
example_geno = genotype_matrix[:100, example_snp_indices]  # 100 individuals, 4 SNPs

# Standardize
for i in range(4):
    col = example_geno[:, i]
    example_geno[:, i] = (col - col.mean()) / (col.std() + 1e-10)

# Create additive and epistatic signals
additive = example_geno.sum(axis=1)  # Sum of all 4 SNPs
additive_std = (additive - additive.mean()) / (additive.std() + 1e-10)

# Create epistatic: all pairwise products
interactions = []
interaction_pairs = []
for i in range(4):
    for j in range(i+1, 4):
        interaction = example_geno[:, i] * example_geno[:, j]
        interaction_std = (interaction - interaction.mean()) / (interaction.std() + 1e-10)
        interactions.append(interaction_std)
        interaction_pairs.append((i, j))

epistatic = np.sum(interactions, axis=0)
epistatic_std = (epistatic - epistatic.mean()) / (epistatic.std() + 1e-10)

corr_toy = np.corrcoef(additive_std, epistatic_std)[0, 1]

print(f"Additive component:   G_add = SNP₀ + SNP₁ + SNP₂ + SNP₃")
print(f"Epistatic component:  G_epi = (SNP₀×SNP₁) + (SNP₀×SNP₂) + ... + (SNP₂×SNP₃)")
print(f"\nCorrelation: {corr_toy:.4f}")

# Show why this happens mathematically
print("\n--- Mathematical Explanation ---")
print("When SNPs are in LD (correlated with each other):")
print("  SNP_i × SNP_j ≈ SNP_i + SNP_j  (approximately!)")
print("\nWhy? Because standardized genotypes are centered around 0:")
print("  - When SNP_i = 1, SNP_j = 1 → product = +1")
print("  - When SNP_i = -1, SNP_j = -1 → product = +1")
print("  - When SNP_i = 1, SNP_j = -1 → product = -1")
print("  - When SNP_i = -1, SNP_j = 1 → product = -1")
print("\nIf SNP_i and SNP_j are positively correlated (LD):")
print("  → First two cases are more common → products tend to be positive")
print("  → This mimics additive effects!")

# Demonstrate with actual LD
print("\n--- Real LD in Example SNPs ---")
print("Pairwise correlations between SNPs:")
snp_corr_matrix = np.corrcoef(example_geno.T)
for i in range(4):
    for j in range(i+1, 4):
        print(f"  Corr(SNP_{i}, SNP_{j}) = {snp_corr_matrix[i, j]:6.3f}")

print("\nCorrelation between each interaction and additive component:")
for idx, (i, j) in enumerate(interaction_pairs):
    corr_with_add = np.corrcoef(interactions[idx], additive_std)[0, 1]
    snp_ld = snp_corr_matrix[i, j]
    print(f"  Corr(SNP_{i}×SNP_{j}, G_add) = {corr_with_add:6.3f} | LD(SNP_{i},SNP_{j}) = {snp_ld:6.3f}")


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


print("="*80)
print("PART 3: THE FUNDAMENTAL PROBLEM")
print("="*80)

print("\nYour simulation creates epistatic signal from products of causal SNPs:")
print("  G_epi = Σ (SNP_i × SNP_j)  where SNP_i, SNP_j are causal")
print("\nYour additive signal is the sum of the SAME causal SNPs:")
print("  G_add = Σ SNP_k  where SNP_k are the same causal SNPs")

print("\n--- Why This Creates High Correlation ---")
print("\n1. **Same SNP Set:**")
print("   - Both components use the EXACT same causal SNPs")
print("   - Products of SNPs are algebraically related to sums of SNPs")

print("\n2. **LD Structure:**")
print(f"   - Mean |LD| among causal SNPs: {results[1]['mean_ld']:.3f}")
print("   - When SNP_i correlates with SNP_j, their product SNP_i×SNP_j")
print("     correlates with both SNP_i and SNP_j individually")
print("   - Sum of products ≈ weighted sum of individual SNPs")

print("\n3. **Mathematical Identity:**")
print("   - For centered/standardized variables:")
print("     E[SNP_i × SNP_j] = Cov(SNP_i, SNP_j)")
print("   - If LD > 0, then E[SNP_i × SNP_j] > 0")
print("   - This creates positive contribution to both G_add and G_epi")

print("\n4. **Central Limit Theorem:**")
print("   - You're summing 2000 interaction terms")
print("   - Each term is a function of the same underlying SNPs")
print("   - By CLT, the sum converges to a similar distribution as G_add")

print("\n--- Concrete Example ---")
print("\nSuppose SNP_1 and SNP_2 have LD = 0.8:")
print("  Individual with SNP_1 = 1 likely has SNP_2 = 1 (LD > 0)")
print("  → Additive contribution: 1 + 1 = 2")
print("  → Epistatic contribution: 1 × 1 = 1")
print("  Both are positive!")
print("\n  Individual with SNP_1 = -1 likely has SNP_2 = -1")
print("  → Additive contribution: -1 + -1 = -2")
print("  → Epistatic contribution: -1 × -1 = 1  (still positive!)")
print("  → But the SIGN of genetic risk is the same (both negative)")

print("\n--- Why Changing α Doesn't Matter ---")
print("\nBecause Corr(G_add, G_epi) ≈ 0.98, changing α is like:")
print("  G_combined = √α × G_epi + √(1-α) × G_add")
print("             ≈ √α × 0.98×G_add + √(1-α) × G_add")
print("             = [√α × 0.98 + √(1-α)] × G_add")
print("             ≈ constant × G_add")
print("\nNo matter what α is, G_combined is just a rescaled version of G_add!")
print("This is why performance doesn't change with α.")


print("\n" + "="*80)
print("PART 4: WHAT WOULD ACTUALLY CREATE INDEPENDENT EPISTASIS?")
print("="*80)

print("\n**Option 1: Use DIFFERENT SNPs for epistasis**")
print("  - Additive: SNPs 1-100")
print("  - Epistatic: SNPs 101-200")
print("  → But this isn't biologically realistic")

print("\n**Option 2: Use WITHIN-GENE interactions only**")
print("  - Ensures GenNet can learn the interactions")
print("  - Creates less correlation because you're sampling from")
print("    a restricted subset (within each gene)")
print("  - More biologically plausible (gene-gene interactions")
print("    happen via protein products, not random SNP pairs)")

print("\n**Option 3: Introduce NEGATIVE interactions**")
print("  - Current: all interactions have positive effects")
print("  - Alternative: randomly assign positive/negative effects")
print("  - This would decorrelate G_epi from G_add")

print("\n**Option 4: Use HIGHER-ORDER interactions**")
print("  - Current: 2-way (SNP_i × SNP_j)")
print("  - Alternative: 3-way (SNP_i × SNP_j × SNP_k)")
print("  - Higher-order products are less correlated with sums")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)
print("\n✅ Implement within-gene epistasis (Option 2)")
print("   This is:")
print("   - Biologically motivated")
print("   - Learnable by GenNet architecture")
print("   - Will reduce correlation")
print("   - Will show clear ReLU > Linear when α > 0")

print("\n❌ Current cross-gene approach is fundamentally broken because:")
print("   1. GenNet cannot learn cross-gene interactions (architecture mismatch)")
print("   2. High correlation makes α parameter meaningless")
print("   3. Doesn't test what you want to test (epistasis detection)")

print("\n" + "="*80)
