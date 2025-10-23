#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error.
set -u
# The return value of a pipeline is the status of the last command to exit with a non-zero status.
set -o pipefail

# --- Default Parameters ---
VCF_FILE=""
GENNET_SCRIPT_PATH="GenNet.py"
ANNOVAR_PATH=""
ANNOVAR_DB=""
GENOME_BUILD="hg38"
SUBSET_N=""
POLYGENICITY=1000
HERITABILITY=0.5
CONTINUOUS=false
SEED=42
OUT_DIR="./gennet_pipeline_output"

# --- Helper Functions ---
print_usage() {
    echo "Usage: $0 --vcf <path> --annovar-path <path> --annovar-db <path> [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  --vcf <path>              Path to the input VCF file (.vcf.gz)."
    echo "  --annovar-path <path>     Path to ANNOVAR's table_annovar.pl script."
    echo "  --annovar-db <path>       Path to ANNOVAR's human database directory (humandb)."
    echo ""
    echo "Optional arguments:"
    echo "  --gennet-script-path <path> Path to the main GenNet.py script. [Default: GenNet.py]"
    echo "  --genome-build <string>   Genome build version for ANNOVAR. [Default: hg38]"
    echo "  --subset-n <int>          Number of samples to subset for testing. [Default: all samples]"
    echo "  --polygenicity <int>      Number of causal SNPs for the simulation. [Default: 1000]"
    echo "  --heritability <float>    Heritability (h^2) of the simulated trait. [Default: 0.5]"
    echo "  --continuous              Flag to simulate a continuous trait instead of binary. [Default: false]"
    echo "  --seed <int>              Random seed for reproducibility. [Default: 42]"
    echo "  --out-dir <path>          Main output directory for all generated files. [Default: ./gennet_pipeline_output]"
    echo "  -h, --help                Show this help message."
}

# --- Parse Command-Line Arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --vcf) VCF_FILE="$2"; shift ;;
        --gennet-script-path) GENNET_SCRIPT_PATH="$2"; shift ;;
        --annovar-path) ANNOVAR_PATH="$2"; shift ;;
        --annovar-db) ANNOVAR_DB="$2"; shift ;;
        --genome-build) GENOME_BUILD="$2"; shift ;;
        --subset-n) SUBSET_N="$2"; shift ;;
        --polygenicity) POLYGENICITY="$2"; shift ;;
        --heritability) HERITABILITY="$2"; shift ;;
        --continuous) CONTINUOUS=true ;;
        --seed) SEED="$2"; shift ;;
        --out-dir) OUT_DIR="$2"; shift ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; print_usage; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [ -z "$VCF_FILE" ] || [ -z "$ANNOVAR_PATH" ] || [ -z "$ANNOVAR_DB" ]; then
    echo "Error: Missing required arguments."
    print_usage
    exit 1
fi

# --- 0. Setup Directories ---
echo "--- Setting up directory structure in $OUT_DIR ---"
PLINK_FULL_DIR="$OUT_DIR/plink_full"
PLINK_SUBSET_DIR="$OUT_DIR/plink_subset"
ANNOVAR_DIR="$OUT_DIR/annovar"
GENNET_INPUT_DIR="$OUT_DIR/gennet_input"
CLEAN_PLINK_DIR="$OUT_DIR/clean_plink_for_gennet"

mkdir -p "$PLINK_FULL_DIR" "$PLINK_SUBSET_DIR" "$ANNOVAR_DIR" "$GENNET_INPUT_DIR" "$CLEAN_PLINK_DIR"

# --- 1. Convert VCF to PLINK ---
echo -e "\n--- STEP 1: Convert VCF to PLINK format ---"
FULL_PLINK_PREFIX="$PLINK_FULL_DIR/all_samples"
plink2 --vcf "$VCF_FILE" --make-bed --out "$FULL_PLINK_PREFIX"
CURRENT_PLINK_PREFIX=$FULL_PLINK_PREFIX

# --- 2. Subset to N samples (if requested) ---
if [ ! -z "$SUBSET_N" ]; then
    echo -e "\n--- STEP 2: Subsetting to $SUBSET_N samples ---"
    KEEP_FILE="$PLINK_SUBSET_DIR/subset_ids.txt"
    head -n "$SUBSET_N" "$FULL_PLINK_PREFIX.fam" | awk '{print $1, $2}' > "$KEEP_FILE"
    
    SUBSET_PLINK_PREFIX="$PLINK_SUBSET_DIR/subset_data"
    plink2 --bfile "$FULL_PLINK_PREFIX" --keep "$KEEP_FILE" --make-bed --out "$SUBSET_PLINK_PREFIX"
    CURRENT_PLINK_PREFIX=$SUBSET_PLINK_PREFIX
fi

# --- 3. Simulate Phenotype ---
echo -e "\n--- STEP 3: Simulate Phenotype ---"
BIM_FILE="$CURRENT_PLINK_PREFIX.bim"
FAM_FILE="$CURRENT_PLINK_PREFIX.fam"
CAUSAL_FILE="$(dirname "$CURRENT_PLINK_PREFIX")/causal_snps.txt"

echo "Using Python for reproducible random sampling of causal SNPs."
python3 -c "
import pandas as pd
import numpy as np
bim_df = pd.read_csv('$BIM_FILE', sep='\s+', header=None, names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
causal_snps_df = bim_df.sample(n=$POLYGENICITY, random_state=$SEED)
rng = np.random.RandomState($SEED)
causal_snps_df['beta'] = rng.uniform(-0.5, 0.5, size=$POLYGENICITY)
causal_snps_df[['snp', 'beta']].to_csv('$CAUSAL_FILE', sep='\t', index=False, header=False)
"
echo "Saved $POLYGENICITY causal SNPs to $CAUSAL_FILE"

PHENO_PREFIX="$(dirname "$CURRENT_PLINK_PREFIX")/simulated_pheno"

# This command now uses PLINK 2.0 flags
if [ "$CONTINUOUS" = true ]; then
    # Simulate a continuous trait
    SIMU_CMD="plink2 --bfile $CURRENT_PLINK_PREFIX --simu-causal-snps $CAUSAL_FILE --simu-h2 $HERITABILITY --out $PHENO_PREFIX"
else
    # Simulate a binary (case/control) trait
    NUM_SAMPLES=$(wc -l < "$FAM_FILE")
    NUM_CASES=$((NUM_SAMPLES / 2))
    NUM_CONTROLS=$((NUM_SAMPLES - NUM_CASES))
    SIMU_CMD="plink2 --bfile $CURRENT_PLINK_PREFIX --simu-causal-snps $CAUSAL_FILE --simu-h2 $HERITABILITY --simu-cc $NUM_CASES $NUM_CONTROLS --out $PHENO_PREFIX"
fi
eval "$SIMU_CMD"
# Rename PLINK2 output to match expected .pheno extension
mv "$PHENO_PREFIX.sscore.pheno" "$PHENO_PREFIX.pheno"


# --- 4. Format Data for GenNet ---
echo -e "\n--- STEP 4: Format Data for GenNet ---"
PHENO_FILE="$PHENO_PREFIX.pheno"
SUBJECT_CSV="$GENNET_INPUT_DIR/subject.csv"

python3 -c "
import pandas as pd; import numpy as np;
fam=pd.read_csv('$FAM_FILE', sep='\s+', header=None, usecols=[0,1], names=['patient_id','IID']);
pheno=pd.read_csv('$PHENO_FILE', sep='\s+', header=None, usecols=[2], names=['labels']);
# PLINK2 uses 1 for control, 2 for case. GenNet expects 0 and 1.
if not $CONTINUOUS:
    pheno['labels'] = pheno['labels'] - 1
df=pd.concat([fam, pheno], axis=1);
df['genotype_row']=np.arange(len(df));
np.random.seed($SEED);
df['set']=np.random.choice([1,2,3], size=len(df), p=[0.6,0.2,0.2]);
df[['patient_id','labels','genotype_row','set']].to_csv('$SUBJECT_CSV', index=False);
"
echo "Created GenNet subject file at: $SUBJECT_CSV"

echo "Converting PLINK to HDF5..."
rm -rf "$CLEAN_PLINK_DIR"/*
cp "$CURRENT_PLINK_PREFIX".{bed,bim,fam} "$CLEAN_PLINK_DIR"/
python3 "$GENNET_SCRIPT_PATH" convert -g "$CLEAN_PLINK_DIR" -study_name "$(basename "$CURRENT_PLINK_PREFIX")" -o "$GENNET_INPUT_DIR"
mv "$GENNET_INPUT_DIR/$(basename "$CURRENT_PLINK_PREFIX").h5" "$GENNET_INPUT_DIR/genotype.h5"
echo "Created GenNet genotype file at: $GENNET_INPUT_DIR/genotype.h5"

# --- 5. Generate Topology File ---
echo -e "\n--- STEP 5: Generate Topology File ---"
ANNOVAR_INPUT="$ANNOVAR_DIR/annovar_input.txt"
awk '{print $1, $4, $4, $5, $6, $2}' "$BIM_FILE" > "$ANNOVAR_INPUT"
ANNOVAR_PREFIX="$ANNOVAR_DIR/snp_gene_map"
perl "$ANNOVAR_PATH" "$ANNOVAR_INPUT" "$ANNOVAR_DB" -buildver "$GENOME_BUILD" -out "$ANNOVAR_PREFIX" -protocol refGene -operation g -nastring .
SNP_GENE_FILE="$ANNOVAR_PREFIX.refGene.variant_function"

GENE_PATHWAY_FILE="$ANNOVAR_DIR/gene_pathway_map.txt"
python3 pathway_mapper.py --annovar "$SNP_GENE_FILE" --out "$GENE_PATHWAY_FILE"

TOPOLOGY_FILE="$GENNET_INPUT_DIR/topology.tsv"
python3 topology_builder.py --snp_gene "$SNP_GENE_FILE" --gene_pathway "$GENE_PATHWAY_FILE" --out "$TOPOLOGY_FILE"

echo -e "\n--- Pipeline Complete! ---"
echo "Your GenNet input files are ready in: $GENNET_INPUT_DIR"
echo "  - $GENNET_INPUT_DIR/genotype.h5"
echo "  - $GENNET_INPUT_DIR/subject.csv"
echo "  - $GENNET_INPUT_DIR/topology.tsv"

# Example usage:
# ./prepare_gennet_data.sh \
#   --vcf ./data/raw/geno_phased_fixed.vcf.gz \
#   --annovar-path /path/to/annovar/table_annovar.pl \
#   --annovar-db /path/to/annovar/humandb/ \
#   --subset-n 1000 \
#   --polygenicity 500 \
#   --heritability 0.6 \
#   --out-dir ./my_simulation_run

