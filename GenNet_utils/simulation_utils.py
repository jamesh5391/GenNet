"""
Simulation utility functions for GenNet phenotype simulation.

This module provides reusable functions for:
- Calculating genetic risk scores (GRS)
- Applying liability threshold models
- Creating GenNet subject files
- Creating effect files for PLINK scoring

These utilities are shared between the Jupyter notebook simulations
and the comprehensive simulation pipeline script.
"""

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Dict


def create_effect_file_from_snps(
    causal_snps: np.ndarray,
    bim_file: str,
    effect_size: Optional[float] = None,
    output_file: str = None,
    equal_effects: bool = True
) -> str:
    """
    Create a PLINK effect file from a list of causal SNPs.

    Args:
        causal_snps: Array of causal SNP IDs
        bim_file: Path to PLINK .bim file
        effect_size: Effect size per SNP (if None, calculated as 1/sqrt(n_snps))
        output_file: Path to output effect file
        equal_effects: If True, all SNPs get equal effect sizes

    Returns:
        Path to the created effect file
    """
    # Read BIM file
    bim = pd.read_csv(bim_file, sep='\t', header=None,
                     names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])

    # Calculate effect size if not provided
    if effect_size is None:
        effect_size = 1.0 / np.sqrt(len(causal_snps))

    # Filter to causal SNPs (vectorized for performance)
    causal_snps_set = set(causal_snps)
    causal_bim = bim[bim['snp'].isin(causal_snps_set)].copy()

    # Assign effect sizes
    if equal_effects:
        causal_bim['BETA'] = effect_size
    else:
        # Could implement variable effect sizes here if needed
        causal_bim['BETA'] = effect_size

    # Create effect file in PLINK format (SNP, A1, BETA)
    effects_df = causal_bim[['snp', 'a1', 'BETA']].copy()
    effects_df.columns = ['SNP', 'A1', 'BETA']

    # Save effect file
    if output_file is None:
        output_file = "effects.txt"

    effects_df.to_csv(output_file, sep='\t', header=False, index=False)

    print(f"Created effect file with {len(effects_df)} causal SNPs")
    print(f"Effect size per SNP: {effect_size:.6f}")

    return str(output_file)


def calculate_grs_from_effects(
    plink_prefix: str,
    effects_file: str,
    output_prefix: str
) -> np.ndarray:
    """
    Calculate genetic risk scores (GRS) using PLINK2 scoring.

    Args:
        plink_prefix: Path to PLINK files (without extension)
        effects_file: Path to effect file (SNP, A1, BETA)
        output_prefix: Prefix for PLINK output files

    Returns:
        Array of GRS values for each sample
    """
    # Run PLINK2 scoring command
    cmd = f"""plink2 --bfile {plink_prefix} \
            --score {effects_file} 1 2 3 \
            --out {output_prefix}"""

    subprocess.run(cmd, shell=True, check=True)

    # Read the score output file
    score_file = f"{output_prefix}.sscore"
    if not Path(score_file).exists():
        # Try alternative naming convention
        score_file = f"{output_prefix}.profile"

    if not Path(score_file).exists():
        raise FileNotFoundError(f"Could not find score file: {score_file}")

    grs_df = pd.read_csv(score_file, sep=r'\s+')

    # Find the score column (different PLINK versions use different names)
    score_col = None
    for col in ['SCORE1_SUM', 'SCORE1_AVG', 'SCORESUM', 'SCORE']:
        if col in grs_df.columns:
            score_col = col
            break

    if score_col is None:
        available_cols = grs_df.columns.tolist()
        raise ValueError(
            f"Could not find score column. Available columns: {available_cols}"
        )

    grs = grs_df[score_col].values

    print(f"Calculated GRS for {len(grs)} samples")
    print(f"  GRS mean: {grs.mean():.6f}")
    print(f"  GRS std: {grs.std():.6f}")
    print(f"  GRS range: [{grs.min():.6f}, {grs.max():.6f}]")

    return grs


def apply_liability_threshold_model(
    grs: np.ndarray,
    h2: float = 0.5,
    prevalence: float = 0.5,
    seed: Optional[int] = None,
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Apply liability threshold model to convert GRS to binary phenotype.

    This implements the classic liability threshold model where:
    - Liability = sqrt(h2) * GRS + sqrt(1-h2) * Environmental_noise
    - Phenotype = 1 if Liability > threshold, else 0

    Args:
        grs: Array of genetic risk scores
        h2: Heritability (proportion of variance explained by genetics)
        prevalence: Disease prevalence (determines threshold percentile)
        seed: Random seed for environmental component
        verbose: Whether to print statistics

    Returns:
        Tuple of (phenotype, liability, genetic_component, environmental_component, threshold)
    """
    if seed is not None:
        np.random.seed(seed)

    # Standardize GRS
    if grs.std() > 0:
        grs_std = (grs - grs.mean()) / grs.std()
    else:
        grs_std = grs - grs.mean()

    # Calculate genetic and environmental components
    genetic_component = np.sqrt(h2) * grs_std
    environmental_component = np.sqrt(1 - h2) * np.random.normal(0, 1, len(grs))

    # Total liability
    liability = genetic_component + environmental_component

    # Apply threshold based on prevalence
    threshold_percentile = 100 * (1 - prevalence)
    threshold = np.percentile(liability, threshold_percentile)

    # Binary phenotype
    phenotype = (liability > threshold).astype(int)

    if verbose:
        n_cases = phenotype.sum()
        n_controls = len(phenotype) - n_cases

        print(f"\nLiability Threshold Model Applied:")
        print(f"  Heritability (h²): {h2}")
        print(f"  Target prevalence: {prevalence:.1%}")
        print(f"  Genetic component: mean={genetic_component.mean():.4f}, std={genetic_component.std():.4f}")
        print(f"  Environmental component: mean={environmental_component.mean():.4f}, std={environmental_component.std():.4f}")
        print(f"  Liability threshold: {threshold:.4f}")
        print(f"  Phenotype distribution:")
        print(f"    Cases: {n_cases} ({100*phenotype.mean():.1f}%)")
        print(f"    Controls: {n_controls} ({100*(1-phenotype.mean()):.1f}%)")

    return phenotype, liability, genetic_component, environmental_component, threshold


def create_gennet_subject_file(
    plink_prefix: str,
    phenotype: np.ndarray,
    output_file: str,
    train_frac: float = 0.6,
    val_frac: float = 0.2,
    test_frac: float = 0.2,
    seed: Optional[int] = None,
    verbose: bool = True
) -> Tuple[str, Dict[str, int]]:
    """
    Create a GenNet-format subject file with train/val/test splits.

    GenNet expects a CSV file with columns:
    - patient_id: Sample identifier
    - labels: Phenotype (0/1 for binary)
    - genotype_row: Row index in genotype matrix
    - set: 1=train, 2=validation, 3=test

    Args:
        plink_prefix: Path to PLINK files (to get sample IDs from .fam)
        phenotype: Array of phenotype values
        output_file: Path to output subject file
        train_frac: Fraction of samples for training
        val_frac: Fraction of samples for validation
        test_frac: Fraction of samples for testing
        seed: Random seed for splitting
        verbose: Whether to print split statistics

    Returns:
        Tuple of (subject_file_path, split_info_dict)
    """
    # Validate split fractions
    if not np.isclose(train_frac + val_frac + test_frac, 1.0):
        raise ValueError(
            f"Split fractions must sum to 1.0, got {train_frac + val_frac + test_frac}"
        )

    # Read FAM file to get sample IDs
    fam = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None,
                     names=['fid', 'iid', 'father', 'mother', 'sex', 'pheno'])

    n_samples = len(fam)

    if len(phenotype) != n_samples:
        raise ValueError(
            f"Phenotype length ({len(phenotype)}) doesn't match FAM file ({n_samples})"
        )

    # Create random train/val/test split
    if seed is not None:
        np.random.seed(seed)

    indices = np.arange(n_samples)
    np.random.shuffle(indices)

    # Calculate split points
    train_end = int(train_frac * n_samples)
    val_end = train_end + int(val_frac * n_samples)

    train_idx = indices[:train_end]
    val_idx = indices[train_end:val_end]
    test_idx = indices[val_end:]

    # Create set assignments (GenNet convention: 1=train, 2=val, 3=test)
    set_assignment = np.zeros(n_samples, dtype=int)
    set_assignment[train_idx] = 1
    set_assignment[val_idx] = 2
    set_assignment[test_idx] = 3

    # Create GenNet-format DataFrame
    subject_df = pd.DataFrame({
        'patient_id': fam['iid'].values,
        'labels': phenotype,
        'genotype_row': np.arange(n_samples),
        'set': set_assignment
    })

    # Save to CSV
    subject_df.to_csv(output_file, index=False)

    # Calculate split statistics
    split_info = {
        'train': int(sum(set_assignment == 1)),
        'val': int(sum(set_assignment == 2)),
        'test': int(sum(set_assignment == 3)),
        'train_cases': int(phenotype[train_idx].sum()),
        'val_cases': int(phenotype[val_idx].sum()),
        'test_cases': int(phenotype[test_idx].sum())
    }

    if verbose:
        print(f"\nCreated GenNet subject file: {output_file}")
        print(f"  Train: {split_info['train']} samples ({split_info['train_cases']} cases)")
        print(f"  Validation: {split_info['val']} samples ({split_info['val_cases']} cases)")
        print(f"  Test: {split_info['test']} samples ({split_info['test_cases']} cases)")

    return str(output_file), split_info


def simulate_pathway_phenotype(
    plink_prefix: str,
    topology_data: list,
    causal_pathways: np.ndarray,
    h2: float = 0.5,
    prevalence: float = 0.5,
    seed: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """
    Simulate phenotype based on pathway-level effects.

    This is a convenience function that combines several steps:
    1. Get genes in causal pathways
    2. Get SNPs in those genes
    3. Create effect file
    4. Calculate GRS
    5. Apply liability threshold model

    Args:
        plink_prefix: Path to PLINK files
        topology_data: List of topology rows (dicts with layer info)
        causal_pathways: Array of causal pathway names
        h2: Heritability
        prevalence: Disease prevalence
        seed: Random seed

    Returns:
        Tuple of (phenotype, causal_genes, causal_snps, effects_file)
    """
    if seed is not None:
        np.random.seed(seed)

    print(f"\n=== Simulating Pathway-based Phenotype ===")
    print(f"Causal pathways: {causal_pathways}")

    # Convert topology to DataFrame for easier querying
    topology_df = pd.DataFrame(topology_data)

    # Get genes in causal pathways
    causal_genes = topology_df[
        topology_df['layer2_name'].isin(causal_pathways)
    ]['layer1_name'].unique()

    print(f"Number of genes in causal pathways: {len(causal_genes)}")

    # Get SNPs in those genes
    causal_snps = topology_df[
        topology_df['layer1_name'].isin(causal_genes)
    ]['layer0_name'].unique()

    print(f"Number of causal SNPs: {len(causal_snps)}")

    # Print pathway composition
    print(f"\n--- Causal Pathway Composition ---")
    for pathway in causal_pathways:
        genes_in_pathway = topology_df[
            topology_df['layer2_name'] == pathway
        ]['layer1_name'].nunique()
        snps_in_pathway = topology_df[
            topology_df['layer2_name'] == pathway
        ]['layer0_name'].nunique()
        print(f"{pathway}: {genes_in_pathway} genes, {snps_in_pathway} SNPs")

    # Create effect file
    bim_file = f"{plink_prefix}.bim"
    effects_file = f"{Path(plink_prefix).parent}/effects_pathway.txt"

    create_effect_file_from_snps(
        causal_snps=causal_snps,
        bim_file=bim_file,
        output_file=effects_file
    )

    # Calculate GRS
    grs_prefix = f"{Path(plink_prefix).parent}/grs_pathway"
    grs = calculate_grs_from_effects(
        plink_prefix=plink_prefix,
        effects_file=effects_file,
        output_prefix=grs_prefix
    )

    # Apply liability threshold model
    phenotype, liability, genetic_comp, env_comp, threshold = apply_liability_threshold_model(
        grs=grs,
        h2=h2,
        prevalence=prevalence,
        seed=seed,
        verbose=True
    )

    return phenotype, causal_genes, causal_snps, effects_file
