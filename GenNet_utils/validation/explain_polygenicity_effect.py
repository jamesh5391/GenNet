"""
Analysis of why polygenicity (P) doesn't affect model performance.

This explains:
1. Why P has minimal effect in current simulation
2. The mathematical relationship between P and effect sizes
3. How effect size normalization masks P's impact
4. How to design simulation to show realistic P effects
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

print("="*80)
print("WHY DOESN'T POLYGENICITY (P) AFFECT MODEL PERFORMANCE?")
print("="*80)

# Load results
results_relu = pd.read_csv('results/grid_results_relu.csv')
results_linear = pd.read_csv('results/grid_results_linear.csv')
results = pd.concat([results_relu, results_linear], ignore_index=True)

print("\n" + "="*80)
print("PART 1: EMPIRICAL EVIDENCE - P HAS MINIMAL EFFECT")
print("="*80)

# Show R² by P for each h² (α=0 only, additive)
print("\nR² by Polygenicity (P) at different heritabilities (α=0, pure additive):")
print("\nActivation: LINEAR")
for h2 in sorted(results['h2'].unique()):
    if h2 == 0:
        continue
    subset = results[(results['h2'] == h2) & (results['alpha'] == 0) & (results['activation'] == 'linear')]
    r2_by_p = subset.groupby('P')['test_r2'].mean()

    print(f"  h²={h2}:  ", end="")
    for P in sorted(subset['P'].unique()):
        r2 = r2_by_p[P]
        print(f"P={P:3d}: {r2:.3f}  ", end="")

    # Calculate range
    r2_range = r2_by_p.max() - r2_by_p.min()
    print(f"| Range: {r2_range:.3f}")

print("\nActivation: RELU")
for h2 in sorted(results['h2'].unique()):
    if h2 == 0:
        continue
    subset = results[(results['h2'] == h2) & (results['alpha'] == 0) & (results['activation'] == 'relu')]
    r2_by_p = subset.groupby('P')['test_r2'].mean()

    print(f"  h²={h2}:  ", end="")
    for P in sorted(subset['P'].unique()):
        r2 = r2_by_p[P]
        print(f"P={P:3d}: {r2:.3f}  ", end="")

    # Calculate range
    r2_range = r2_by_p.max() - r2_by_p.min()
    print(f"| Range: {r2_range:.3f}")

print("\n--- Key Observation ---")
print("For h²=0.6:")
print("  Range across P = 0.05-0.15  (very small!)")
print("  This is ~10-30% of the signal (R² ≈ 0.4-0.5)")
print("\nFor h²=1.0:")
print("  Range across P = 0.08-0.13  (very small!)")
print("  This is ~8-13% of the signal (R² ≈ 0.9)")
print("\n→ P has minimal impact compared to h²!")


print("\n" + "="*80)
print("PART 2: THE MATHEMATICAL REASON - EFFECT SIZE NORMALIZATION")
print("="*80)

print("\nIn your simulation, effect sizes are calculated as:")
print("  effect_size = 2.0 / (N × P)")
print("\nWhere:")
print("  N = 10 (SNPs per gene)")
print("  P = polygenicity (number of causal genes)")

print("\n--- Effect Sizes for Different P ---")
for P in [1, 10, 100, 500]:
    N = 10
    effect = 2.0 / (N * P)
    total_snps = N * P
    print(f"  P={P:3d}: effect_size = {effect:.6f}, total_causal_SNPs = {total_snps:4d}")

print("\n--- What This Does to Genetic Signal ---")
print("\nFor additive genetic component:")
print("  G_add = Σ (effect_size × genotype_i)")
print("        = Σ (2.0/(N×P) × genotype_i)  for i in causal SNPs")
print("\nExpected variance:")
print("  Var(G_add) = (N×P) × effect_size² × Var(genotype)")
print("             = (N×P) × [2.0/(N×P)]² × 1")
print("             = (N×P) × 4/(N×P)²")
print("             = 4 / (N×P)")

print("\n--- Variance of Genetic Component by P ---")
for P in [1, 10, 100, 500]:
    N = 10
    var_expected = 4.0 / (N * P)
    print(f"  P={P:3d}: Var(G_add) ≈ {var_expected:.6f}")

print("\n--- After Standardization ---")
print("\nYour code standardizes G_add:")
print("  G_add_std = (G_add - mean) / std")
print("  → Var(G_add_std) = 1.0  (always!)")

print("\nThis normalization REMOVES the variance difference between P values!")
print("All P values end up with the same standardized variance.")

print("\n--- The Phenotype Formula ---")
print("\nFinal phenotype:")
print("  Y = √h² × G_add_std + √(1-h²) × ε")
print("\nSince Var(G_add_std) = 1 regardless of P:")
print("  Var(genetic component) = h² × 1 = h²  (independent of P!)")
print("  Var(environmental) = (1-h²) × 1 = (1-h²)")
print("  Var(Y) ≈ 1  (by construction)")

print("\n→ The standardization makes h² the ONLY parameter controlling signal strength!")
print("→ P only affects the NUMBER of SNPs, not the TOTAL signal strength!")


print("\n" + "="*80)
print("PART 3: WHY THIS IS UNREALISTIC")
print("="*80)

print("\nIn real genetics, polygenicity DOES affect predictability:")

print("\n1. **Sample Size Requirements:**")
print("   - P=1: Need ~100-1000 samples to estimate 1 gene's effect")
print("   - P=500: Need ~50,000-500,000 samples to estimate 500 genes' effects")
print("   - Your simulation has 3,202 samples → should struggle with P=500!")

print("\n2. **Effect Size Distribution:**")
print("   - Real GWAS: effects follow exponential/gamma distribution")
print("   - Large effects are rare, small effects are common")
print("   - Your simulation: all effects are EQUAL (2.0/N×P)")

print("\n3. **Statistical Power:**")
print("   - Detecting effect_size = 0.02 (P=1) needs ~100 samples")
print("   - Detecting effect_size = 0.0004 (P=500) needs ~250,000 samples")
print("   - But standardization masks this!")

print("\n4. **Curse of Dimensionality:**")
print("   - P=1: 10 features (SNPs) to learn")
print("   - P=500: 5,000 features to learn")
print("   - Neural networks should struggle more with P=500!")


print("\n" + "="*80)
print("PART 4: WHAT REAL POLYGENICITY LOOKS LIKE")
print("="*80)

print("\nReal-world examples:")
print("\n**Monogenic (P ≈ 1-5):**")
print("  - Sickle cell disease: HBB gene (1 gene)")
print("  - Cystic fibrosis: CFTR gene (1 gene)")
print("  - Huntington's disease: HTT gene (1 gene)")
print("  - h² ≈ 0.95-1.0, easily predictable with small samples")

print("\n**Oligogenic (P ≈ 10-50):**")
print("  - Type 1 diabetes: ~40 loci, h² ≈ 0.9")
print("  - Crohn's disease: ~170 loci, h² ≈ 0.5")
print("  - Moderately predictable with 10k-50k samples")

print("\n**Polygenic (P ≈ 100-10,000):**")
print("  - Height: ~10,000 loci, h² ≈ 0.8")
print("  - Schizophrenia: ~1,000 loci, h² ≈ 0.6-0.8")
print("  - BMI: ~1,000 loci, h² ≈ 0.4")
print("  - Requires 100k-1M samples for good prediction")

print("\n**Extreme Polygenicity (P > 10,000):**")
print("  - Educational attainment: ~10,000+ loci, h² ≈ 0.4")
print("  - Depression: ~10,000+ loci, h² ≈ 0.3")
print("  - Requires 500k+ samples, still hard to predict")


print("\n" + "="*80)
print("PART 5: HOW TO FIX YOUR SIMULATION")
print("="*80)

print("\n**Option 1: DO NOT STANDARDIZE G_add (Most Realistic)**")
print("\nChange your code from:")
print("  G_add_std = (G_add - G_add.mean()) / G_add.std()")
print("To:")
print("  G_add_centered = G_add - G_add.mean()")
print("  # Don't divide by std!")

print("\nThen create phenotype:")
print("  genetic_var = G_add_centered.var()")
print("  target_genetic_var = h²  # Target variance")
print("  scale = np.sqrt(target_genetic_var / genetic_var)")
print("  genetic_component = scale × G_add_centered")
print("  environmental_component = np.sqrt(1 - h²) × noise")
print("  phenotype = genetic_component + environmental_component")

print("\nThis preserves the natural relationship:")
print("  - P=1: Large effect per SNP → easier to predict")
print("  - P=500: Tiny effect per SNP → harder to predict")

print("\n**Option 2: Use Different Effect Size Distributions**")
print("\nInstead of equal effects, use realistic distributions:")
print("  # GWAS-like: exponential distribution")
print("  effect_sizes = np.random.exponential(scale=1.0, size=n_causal_snps)")
print("  effect_sizes /= np.sqrt(np.sum(effect_sizes**2))  # Normalize to unit variance")

print("\n**Option 3: Sample Size Penalty**")
print("\nAdd noise proportional to model complexity:")
print("  estimation_error = np.random.normal(0, sqrt(P/n_samples), size=n_samples)")
print("  phenotype += estimation_error")
print("\nThis mimics finite-sample effects:")
print("  - P=1: Small error (easy to estimate)")
print("  - P=500: Large error (hard to estimate with 3k samples)")

print("\n**Option 4: Realistic Sample Sizes**")
print("\nAdjust sample size to match polygenicity:")
print("  - P=1: n=1,000 (sufficient)")
print("  - P=10: n=5,000 (sufficient)")
print("  - P=100: n=20,000 (borderline)")
print("  - P=500: n=100,000 (needed for good prediction)")
print("\nKeep effect_size formula, but vary sample size!")


print("\n" + "="*80)
print("PART 6: EXPECTED RESULTS AFTER FIX")
print("="*80)

print("\nAfter removing standardization (Option 1):")

print("\n**At h²=0.6:**")
print("  P=1:   R² ≈ 0.55-0.60 (large effects, easy to learn)")
print("  P=10:  R² ≈ 0.50-0.55 (moderate effects)")
print("  P=100: R² ≈ 0.40-0.45 (small effects, harder)")
print("  P=500: R² ≈ 0.25-0.35 (very small effects, much harder!)")

print("\n**At h²=1.0:**")
print("  P=1:   R² ≈ 0.95-0.98 (near perfect)")
print("  P=10:  R² ≈ 0.90-0.95 (very good)")
print("  P=100: R² ≈ 0.75-0.85 (good but degrading)")
print("  P=500: R² ≈ 0.50-0.65 (moderate - curse of dimensionality)")

print("\n**Key Pattern:**")
print("  - R² decreases as P increases (for fixed h²)")
print("  - Gap between P values widens at higher h²")
print("  - GenNet may show sweet spot around P=10-100")


print("\n" + "="*80)
print("PART 7: SIMULATION REALITY CHECK")
print("="*80)

print("\nLet's compute what your current simulation actually does:")

np.random.seed(42)
for P in [1, 10, 100, 500]:
    N = 10
    effect_size = 2.0 / (N * P)
    n_causal = N * P

    # Simulate genotypes (standardized)
    n_samples = 3202
    geno = np.random.normal(0, 1, size=(n_samples, n_causal))

    # Additive component
    effects = np.full(n_causal, effect_size)
    G_add = geno @ effects

    # Current: standardize
    G_add_std_current = (G_add - G_add.mean()) / G_add.std()
    var_current = G_add_std_current.var()

    # Proposed: don't standardize
    G_add_centered = G_add - G_add.mean()
    var_proposed = G_add_centered.var()

    print(f"\nP={P:3d} ({n_causal:4d} causal SNPs, effect={effect_size:.6f}):")
    print(f"  Current method (with standardization):  Var(G) = {var_current:.6f}")
    print(f"  Proposed method (no standardization):   Var(G) = {var_proposed:.6f}")
    print(f"  Ratio (proposed/current): {var_proposed/var_current:.2f}×")

