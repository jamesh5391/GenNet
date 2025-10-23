python pheno_simulator.py --bim your_data.bim --polygenicity 1000 --out ../data/processed/causal_snps.txt

#### **Step 1.3: Simulate the Phenotype using PLINK**

Now, use the `causal_snps.txt` file to generate the phenotype. This command calculates a genetic risk score and adds random noise to achieve the desired heritability.

```bash
# This command simulates a binary phenotype with 50% prevalence and a heritability of 0.5
plink \
    --bfile your_data \
    --make-pheno ../data/processed/causal_snps.txt '*' \
    --simu-qt \
    --simu-h2 0.5 \
    --out ../data/processed/simulated_pheno
* `--simu-h2 0.5`: Adjust this value to change the **heritability**.
* `--simu-prev 0.5`: Adjust for case/control prevalence. Remove for a continuous trait.

This will create a `simulated_pheno.pheno` file.

#### **Step 1.4: Format Data for GenNet**

The following script converts the PLINK output into the `subject.csv` file and provides the command to create `genotype.h5`.


http://googleusercontent.com/immersive_entry_chip/1

**How to run it:**
```bash
# 1. Format the phenotype into subject.csv
python formatter.py --fam your_data.fam --pheno ../data/processed/simulated_pheno.pheno --out ../data/processed/subject.csv

# 2. Convert the PLINK genotype data to HDF5 format using the GenNet tool
python ../GenNet.py convert --plink your_data --output ../data/processed/genotype.h5

***
### ## Part 2: Annotation and Topology Generation

This part uses standard bioinformatics tools to create the network structure.

#### **Step 2.1: Prerequisites**

* You need to have **ANNOVAR** installed.
* You need Python with `gseapy` installed (`pip install gseapy`).

#### **Step 2.2: Map SNPs to Genes with ANNOVAR**

First, create an ANNOVAR input file from your `.bim` file.
```bash
awk '{print $1, $4, $4, $5, $6, $2}' your_data.bim > annovar_input.txt

Now, run ANNOVAR. You will need to download the RefSeq gene database for your genome build first.
```bash
# Example for hg19 build
table_annovar.pl annovar_input.txt /path/to/humandb/ -buildver hg19 -out snp_gene_map -protocol refGene -operation g -nastring .
This produces `snp_gene_map.txt.variant_function`, which contains the SNP-to-gene mappings.

#### **Step 2.3: Map Genes to Pathways**

This script takes the list of genes from the ANNOVAR output and maps them to KEGG pathways.


http://googleusercontent.com/immersive_entry_chip/2

**How to run it:**
```bash
python pathway_mapper.py --annovar snp_gene_map.txt.variant_function --out gene_pathway_map.txt

#### **Step 2.4: Build the Final Topology File**

Finally, this script combines the two mapping files into the `topology.tsv` file required by GenNet.


http://googleusercontent.com/immersive_entry_chip/3

**How to run it:**
```bash
python topology_builder.py \
    --snp_gene snp_gene_map.txt.variant_function \
    --gene_pathway gene_pathway_map.txt \
    --bim your_data.bim \
    --out topology.tsv

After completing all these steps, you will have the three required files (`genotype.h5`, `subject.csv`, `topology.tsv`) to run your `SNP -> Gene -> Pathway` simulation in the GenNet command-line environment.