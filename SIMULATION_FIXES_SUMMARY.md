# GenNet Simulation Framework: Critical Fixes Required

**Date:** 2025-12-31
**Status:** Issues identified, fixes specified, ready for implementation

---

## Executive Summary

Your GenNet grid experiments revealed that **neither epistasis (α) nor polygenicity (P) affect model performance**. Deep analysis identified **two fundamental flaws** in the simulation framework that must be fixed before re-running experiments.

---

## Problem 1: Epistatic Interactions Cannot Be Learned

### **The Issue**

**Cross-gene interaction sampling:**
- 91-99% of epistatic interactions are between SNPs in **different genes**
- GenNet topology: `SNP → Gene → Output`
- **No pathway exists** to learn interactions between SNPs in different genes
- Both Linear and ReLU fail equally → explains why activations don't matter

**Perfect correlation between additive and epistatic components:**
- Corr(G_add, G_epi) = **0.98** (nearly perfect!)
- Changing α from 0→1 just rescales the **same signal**
- No meaningful difference in phenotypes across α values

### **Root Causes**

1. **Architecture mismatch:**
   ```python
   # Current: samples pairs from ALL causal SNPs
   selected_snp_idx = np.random.choice(
       causal_snps_additive,  # All SNPs across all genes!
       size=interaction_order,
       replace=False
   )
   ```

   - For P=10: Only 8% of interactions are within-gene
   - For P=100: Only 0.8% are within-gene
   - GenNet cannot learn the other 92-99%!

2. **Mathematical redundancy:**
   - Both G_add and G_epi use the **same causal SNPs**
   - Products of SNPs are algebraically related to sums of SNPs
   - With 2000 interaction terms, Central Limit Theorem makes them converge to similar distributions
   - Expected value: E[SNP_i × SNP_j] = Cov(SNP_i, SNP_j) = LD
   - Sum of products ≈ weighted sum of individual SNPs

### **The Fix**

**Sample interactions ONLY within genes:**

```python
# ========== EPISTATIC COMPONENT (FIXED) ==========
epistatic_signal = np.zeros(n_samples)

# Create SNP-to-gene mapping
snp_to_gene = {}
gene_to_snps = {}
for gene in causal_genes:
    gene_snps = gene_annotations[gene_annotations['gene_id'] == gene]['snp_id'].values
    for snp in gene_snps:
        if snp in causal_snps_additive:
            snp_to_gene[snp] = gene
            if gene not in gene_to_snps:
                gene_to_snps[gene] = []
            gene_to_snps[gene].append(snp)

# Generate interactions WITHIN genes only
for i in range(n_interactions):
    # Pick a random causal gene
    gene = np.random.choice(causal_genes)
    gene_snps = gene_to_snps[gene]

    # Sample interaction from SNPs within this gene
    if len(gene_snps) >= interaction_order:
        selected_snp_idx = np.random.choice(
            gene_snps,
            size=interaction_order,
            replace=False
        )
    else:
        # Use all available SNPs if gene has fewer than interaction_order
        selected_snp_idx = gene_snps

    # Calculate interaction term
    interaction_term = np.ones(n_samples)
    for idx in selected_snp_idx:
        interaction_term *= geno_matrix[:, idx]

    # Standardize and accumulate
    if interaction_term.std() > 0:
        interaction_term = (interaction_term - interaction_term.mean()) / interaction_term.std()

    epistatic_signal += interaction_term

# Standardize total epistatic signal
G_epi_std = (epistatic_signal - epistatic_signal.mean()) / (epistatic_signal.std() + 1e-10)
```

### **Expected Improvements**

- ✅ 100% within-gene interactions (was 8-0.8%)
- ✅ GenNet CAN learn all interactions
- ✅ Corr(G_add, G_epi) reduced from 0.98 to ~0.4-0.7
- ✅ ReLU will outperform Linear when α > 0
- ✅ R² will decrease as α increases

---

## Problem 2: Polygenicity Has No Effect

### **The Issue**

**Standardization removes P's natural effect:**
- Current results show R² range across P values: only 0.05-0.13
- This is tiny compared to the h² effect (range: 0.0 to 0.9)
- In reality, P **should** strongly affect predictability

### **Root Cause**

```python
# Current code standardizes G_add
G_add_std = (G_add - G_add.mean()) / G_add.std()
# → Var(G_add_std) = 1.0 regardless of P!
```

**Why this is wrong:**

1. **Effect size calculation:**
   ```
   effect_size = 2.0 / (N × P)

   P=1:   effect = 0.20,   Var(G_add) ≈ 0.40
   P=10:  effect = 0.02,   Var(G_add) ≈ 0.04
   P=100: effect = 0.002,  Var(G_add) ≈ 0.004
   P=500: effect = 0.0004, Var(G_add) ≈ 0.0008
   ```

2. **Standardization erases this:**
   ```python
   # After standardization
   Var(G_add_std) = 1.0  # Always! For all P values!
   ```

3. **Result:** h² is the ONLY parameter controlling signal strength
   ```python
   phenotype = √h² × G_add_std + √(1-h²) × noise
   # Var(genetic) = h² × 1 = h²  (independent of P!)
   ```

### **Why This Is Unrealistic**

In real genetics:

**Monogenic (P ≈ 1-5):**
- Examples: Sickle cell (HBB), Cystic fibrosis (CFTR)
- h² ≈ 0.95-1.0
- Easily predictable with 100-1,000 samples

**Oligogenic (P ≈ 10-50):**
- Examples: Type 1 diabetes (~40 loci), Crohn's (~170 loci)
- h² ≈ 0.5-0.9
- Needs 10k-50k samples

**Polygenic (P ≈ 100-10,000):**
- Examples: Height (~10k loci), Schizophrenia (~1k loci)
- h² ≈ 0.4-0.8
- Requires 100k-1M samples for good prediction

**Your simulation:**
- 3,202 samples (constant)
- Should easily handle P=1-10
- Should struggle significantly with P=500
- But standardization masks this completely!

### **The Fix**

**Remove standardization of G_add:**

```python
# ========== ADDITIVE COMPONENT (FIXED) ==========
# Calculate additive effects (same as before)
effect_size = 2.0 / (N * P)
effects = np.zeros(n_snps)
effects[causal_snps_additive] = effect_size

# Additive GRS
G_add = geno_matrix @ effects

# OLD (WRONG): Standardize
# G_add_std = (G_add - G_add.mean()) / G_add.std()

# NEW (CORRECT): Only center, don't divide by std
G_add_centered = G_add - G_add.mean()

# ========== COMBINE COMPONENTS (FIXED) ==========
# For epistasis, still standardize G_epi (it's aggregated, needs standardization)
G_epi_std = (epistatic_signal - epistatic_signal.mean()) / (epistatic_signal.std() + 1e-10)

# But for G_add, use centered (not standardized) version
# Weight by α
G_combined = np.sqrt(alpha) * G_epi_std + np.sqrt(1 - alpha) * G_add_centered

# Scale combined signal to match target heritability
genetic_var = G_combined.var()
if genetic_var > 0:
    scale = np.sqrt(h2) / np.sqrt(genetic_var)
else:
    scale = 0

genetic_component = scale * G_combined

# Environmental component
environmental_component = np.sqrt(1 - h2) * np.random.randn(n_samples)

# Final phenotype (DO NOT re-standardize!)
phenotype = genetic_component + environmental_component
```

### **Expected Results After Fix**

**At h²=0.6:**
- P=1:   R² ≈ 0.55-0.60 (large effects, easy)
- P=10:  R² ≈ 0.50-0.55 (moderate)
- P=100: R² ≈ 0.40-0.45 (small effects, harder)
- P=500: R² ≈ 0.25-0.35 (tiny effects, much harder!)

**At h²=1.0:**
- P=1:   R² ≈ 0.95-0.98 (near perfect)
- P=10:  R² ≈ 0.90-0.95 (very good)
- P=100: R² ≈ 0.75-0.85 (good but degrading)
- P=500: R² ≈ 0.50-0.65 (curse of dimensionality)

**Key pattern:**
- R² decreases as P increases (for fixed h²)
- Gap widens at higher h²
- Shows realistic sample size limitations

---

## Implementation Checklist

### **1. Update `simulate_phenotype` function**
Location: `GenNet_utils/simulation_utils.py` (or notebook version)

- [ ] **Fix epistatic sampling**: Change to within-gene only
- [ ] **Remove G_add standardization**: Use centered only
- [ ] **Update variance scaling**: Scale to match h² properly
- [ ] **DO NOT standardize final phenotype**

### **2. Verify fixes work**

```bash
# Run diagnostic to check correlation reduction
python GenNet_utils/check_epistasis_simulation.py
```

**What to verify:**
- [ ] Within-gene interactions: 100% (was 8-0.8%)
- [ ] Corr(G_add, G_epi): <0.7 (was 0.98)
- [ ] Var(G_add) decreases with P
- [ ] Var(phenotype) ≈ 1.0 across all configs
- [ ] Var(genetic_component) = h² across all configs

### **3. Run pilot experiment**

Small test before full grid:
- [ ] P=10, h²=0.6, α∈{0, 0.6, 1.0}
- [ ] Both Linear and ReLU activations
- [ ] Check: ReLU > Linear for α > 0
- [ ] Check: R² decreases as α increases

### **4. Re-run full grid**

- [ ] 64 configs (4 P × 4 h² × 4 α)
- [ ] 2 activations (Linear, ReLU)
- [ ] Total: 128 experiments
- [ ] Estimated time: ~20 hours

### **5. Optional: Add independent replicates**

For robust confidence intervals:
- [ ] 5 replicates per config
- [ ] Different causal genes each replicate
- [ ] Total: 640 experiments
- [ ] Estimated time: ~100 hours

---

## Expected Results After Both Fixes

### **Effect of Heritability (h²)**

All models should show strong h² dependence:
- h²=0.0: R² ≈ -0.1 to 0.0 (no signal)
- h²=0.2: R² ≈ 0.0 to 0.15 (weak)
- h²=0.6: R² ≈ 0.3 to 0.55 (moderate)
- h²=1.0: R² ≈ 0.5 to 0.98 (strong)

**This already works in your current results!**

### **Effect of Polygenicity (P)** — NEW!

At h²=0.6, α=0:
- P=1:   R² ≈ 0.55 (easy)
- P=10:  R² ≈ 0.50 (moderate)
- P=100: R² ≈ 0.42 (harder)
- P=500: R² ≈ 0.30 (much harder)

**Clear downward trend!**

### **Effect of Epistasis (α)** — NEW!

At h²=0.6, P=10:

**Linear activation:**
- α=0.0: R² ≈ 0.50 (pure additive, good)
- α=0.2: R² ≈ 0.45 (mostly additive, slight drop)
- α=0.6: R² ≈ 0.30 (mostly epistatic, large drop)
- α=1.0: R² ≈ 0.05 (pure epistatic, very poor)

**ReLU activation:**
- α=0.0: R² ≈ 0.50 (pure additive, same as Linear)
- α=0.2: R² ≈ 0.48 (mostly additive, slight drop)
- α=0.6: R² ≈ 0.45 (mostly epistatic, moderate performance)
- α=1.0: R² ≈ 0.42 (pure epistatic, still good!)

**Clear ReLU advantage for epistasis!**

### **Interaction Effects**

**P × α interaction:**
- Higher P + higher α = hardest (curse of dimensionality + epistasis)
- Lower P + lower α = easiest (few genes + additive)

**Activation × α interaction:**
- ReLU advantage grows with α
- Gap is minimal at α=0 (pure additive)
- Gap is maximal at α=1 (pure epistatic)

---

## Files to Modify

### Primary file
- **`GenNet_utils/simulation_utils.py`** (if using centralized utils)
- OR **`jupyter_notebooks/simulate_data.ipynb`** (if using notebook)

Function to update: `simulate_phenotype()`

### Diagnostic scripts (already created)
- `GenNet_utils/check_epistasis_simulation.py` — verify epistasis fix
- `GenNet_utils/explain_epistasis_correlation.py` — understand correlation
- `GenNet_utils/explain_polygenicity_effect.py` — understand P effect

---

## What NOT to Change

Keep these the same for comparability:
- ✅ P values: [1, 10, 100, 500]
- ✅ h² values: [0, 0.2, 0.6, 1]
- ✅ α values: [0, 0.2, 0.6, 1]
- ✅ Sample size: 3,202
- ✅ interaction_order: 2 (pairwise)
- ✅ n_interactions: 2,000
- ✅ N_snps_per_gene: 10
- ✅ Total genes: 1,500

---

## Validation Tests

After implementing fixes, verify:

1. **Epistasis fix:**
   ```python
   # Run diagnostic
   python GenNet_utils/check_epistasis_simulation.py

   # Expected output:
   # - Within-gene: 100% (was 8%)
   # - Corr(G_add, G_epi): 0.4-0.7 (was 0.98)
   ```

2. **Polygenicity fix:**
   ```python
   # Check variance scales correctly
   for P in [1, 10, 100, 500]:
       G_add_var = simulate_and_get_variance(P, h2=0.6, alpha=0)
       print(f"P={P}: Var(G_add) = {G_add_var:.4f}")

   # Expected: decreasing trend
   # P=1: ~0.4, P=10: ~0.04, P=100: ~0.004, P=500: ~0.0008
   ```

3. **Overall variance check:**
   ```python
   # All phenotypes should have var ≈ 1.0
   for config in all_configs:
       phenotype = simulate_phenotype(**config)
       print(f"{config}: Var(Y) = {phenotype.var():.4f}")

   # Expected: all close to 1.0
   ```

---

## Timeline

1. **Implement fixes**: 1-2 hours
2. **Run diagnostics**: 15 minutes
3. **Pilot test (6 experiments)**: 1-1.5 hours
4. **Review pilot results**: 30 minutes
5. **Full grid (128 experiments)**: ~20 hours
6. **Analysis and visualization**: 2-3 hours

**Total: ~1-2 days**

---

## Questions or Clarifications?

Contact the analysis team or review:
- `/Users/jhu/Documents/broad/GenNet/GenNet_utils/check_epistasis_simulation.py`
- `/Users/jhu/Documents/broad/GenNet/GenNet_utils/explain_epistasis_correlation.py`
- `/Users/jhu/Documents/broad/GenNet/GenNet_utils/explain_polygenicity_effect.py`

These diagnostic scripts provide detailed explanations and can be re-run after fixes to verify success.

---

**Last Updated:** 2025-12-31
**Status:** Ready for implementation
