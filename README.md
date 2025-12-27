# gnomAD-population-af-analysis
Scrapes population-specific allele/genotype frequencies from the gnomAD website for a list of SNPs, produces allele and genotype frequency and count tables, and runs statistical comparisons between chosen populations and a reference population. This pipeline can also compute a cross-SNP “risk divergence” summary (unweighted or beta-weighted).

# Overview
1. Input file instructions
2. "One-and-done" script instructions
3. Scrape allele count data from the genome frequency site
4. Calculate aggregated allele frequencies & genotype counts
5. Reformat into analysis-ready wide tables (AF / allele counts / genotype counts)
6. Run statistical tests vs gnomAD_EUR per SNP + write per-SNP result tables
7. Caclculate a cumulative risk summary across SNPs (optional beta weighting)

## 1. Input file instructions
### Populations file (csv)
This file defines how subpopulations should be grouped into larger comparison populations. 
Reference for subpopulation availability: https://gnomad.broadinstitute.org/variant/4-99318162-T-C?dataset=gnomad_r4
- **dataset**: Dataset of the subpopulation as listed on website.
- **genetic_ancestry_group**: Genetic ancestry group corresponding to that subpopulation as listed on website.
- **group_name**: The name of the overarching group that each population All subpopulations that share the same group label will be summed together (allele counts and genotype counts).
- **subpop_label**: Name of subpopulation label as listed on website.
  
**Important requirements**
1. Spelling and spacing must match the scraped labels exactly.
2. Any trailing spaces in population labels will break matching and yield missing values.

### SNP file (csv)
- **snp_id**: identifier used throughout the pipeline (often an rsID).
- **chromosome, position**: GRCh38 coordinates.
  - Note: different databases sometimes report slightly different positions for the same rsID. If scraping fails for a SNP, confirm GRCh38 location and alleles using NCBI dbSNP:https://www.ncbi.nlm.nih.gov/snp/rs{rsid}
- **reference_allele, alternate_allele**: alleles for this variant.
- **EA**: whether the alternate allele is the effect allele (“y” = alt is effect allele; “n” = ref is effect allele).
- **beta**: effect size magnitude for the effect allele (your current code treats this as a float; interpretation is yours—GWAS beta, meta-analysis beta, etc.).
- **directionality**: whether the effect allele increases (+) or decreases (-) the trait being studied.

## 2. "One-and-done" script instructions
WIP--coming soon(ish)

## 3. Scrape allele count data from the genome frequency site (1_scrape_gnomad_af.py)
Pull per-population allele count information for every SNP from the online resource, without downloading huge public datasets.

1. First you'll need to install the required packages:

`pip install pandas selenium webdriver-manager`

(You’ll also need a compatible Chrome/Chromedriver; webdriver-manager usually handles this.)

2. Run the script
`python scrape_gnomad_af.py populations.csv snps.csv intermediate.csv`

This will: 
- Open an automated Chrome instance
- Visit SNP pages and scrape population-level counts
- Produce an intermediate csv with the raw scraping output

**Troubleshooting**: 
If the scrape fails for a SNP, you’ll usually see missing/empty values in the later columns for that SNP. This is likely a mismatch in:
- chrom/position
- ref/alt alleles
- build version

## 4. Calculate aggregated allele frequencies & genotype counts (2_calc_allele_freqs.py)
This will convert the file from the previous step into: aggregated allele and genotype counts/frequencies

`python calc_allele_freqs.py populations.csv snps.csv intermediate.csv final_output.csv`


## 5. Reformat into wide tables (3_reframe_allele_freqs.py)
This will create wide tables containing allele/genotype frequencies/counts.

`python reframe_allele_freqs.py final_output.csv`

Produces **3** output files:
1. **prefix_AF_reformatted_output.csv**:Alternate allele frequency per population/group (wide format).
2. **prefix_allele_count_reformatted_output.csv**: Allele counts per SNP per population/group (wide format). This file is also the key input for allele-count statistical testing.
3. **prefix_genotype_count_reformatted_output.csv**: Genotype counts per SNP per population/group (wide format). This file is also the key input for genotype-count statistical testing.

## 6. Run statistical tests vs reference population per SNP + write per-SNP result tables + caclculate a cumulative risk summary across SNPs (with optional beta weighting) (4_count_stats.py)
This script performs per-SNP statistical comparisons vs a reference population (default: gnomAD_EUR) and writes multiple outputs into a folder.

`python 4_count_stats.py snps.csv prefix output_prefix output_folder [beta_weighting yes/no]`

Input notes: 
1. "prefix" must match the prefix used for the _allele_count_reformatted_output.csv and _final_output.csv inputs.
2. output_folder is created automatically.
3. beta_weighting toggles whether you compute the beta-weighted cumulative score.

Sub-functions: 
1. **LD-ish pruning by proximity (snp_loci())**: Before testing, the script prunes SNPs to reduce redundancy from physically nearby variants--SNPs are sorted by chromosome and position, clustered by 50kb (the width of this window can be modified within the script), and within each cluster, the SNP with the largest absolute beta is retained.
2. **Allele-count analysis output (per SNP) allele_count_analysis()**:
- Produces one csv file (in the new output_folder) per SNP ({snp_id}_allele_count_analysis_output.csv), providing 2x2 allele-count comparisons, where “effect allele” is defined by EA (y/n) (not directionality).
- For each population vs reference:
  - Odds ratio (OR) (note: (Haldane) correction = 0.5 is added to each cell when computing OR and its CI to avoid division by zero--this can be modified within the script)
  - SE, logOR, 95% CI (lowerCL, upperCL)
  - z statistic for log(OR)
  - p_fisher: Fisher’s exact test (two-sided)
  - chi2, p_chi2: chi-square test (2×2)
  - min_expected: minimum expected count in chi-square expected table
  - primary_p / primary_test: chooses Fisher if min expected < 5 else chi-square
- How to interpet the allele-count file:
  - **OR_vs_gnomAD_EUR**: OR > 1 → effect allele is enriched in that population relative to EUR
  - **primary_p**: the “main” p-value to report for pop vs EUR for this SNP (depending on expected counts)
  - **primary_test**: indicates whether Fisher or chi-square was used
- **Important note**: this analysis is about the effect allele, not necessarily the risk allele.
Risk allele depends on directionality and is used in the risk summary functions, not in allele_count_analysis().
3. **Genotype-count analysis output (per SNP) genotype_count_analysis_long()**:
- Produces one csv file (in the new output_folder) per SNP ({snp_id}_genotype_count_analysis_output.csv), providing 2x3 genotype-count comparisons
- For each population vs reference:
  - chi2_2x3, p_chi2_2x3: chi-square test on the 2×3 table
  - min_expected_2x3
  - p_fisher_carrier_2x2: Fisher exact on collapsed carrier table:
    - Carrier = het + alt_hom
    - Noncarrier = ref_hom
  - primary_p / primary_test: Fisher if min expected < 5 else chi-square
- How to interpret genotype outputs
  - The 2×3 test asks: “Does the genotype distribution differ from EUR?”
  - The carrier Fisher fallback asks: “Is the carrier rate different from EUR?”
4. **Cumulative risk summaries (across SNPs) (cumulative_risk_allele_summary() and cumulative_risk_allele_summary_beta_weighting())**:
- Will incorporate beta weighting when the [yes] flag is used.
- Produces one {prefix}_risk_analysis_output.csv or {prefix}_risk_analysis_output_beta_weighted.csv file
- Includes:
  - **Ns_sig_vs_ref**: number of SNPs with primary_p < alpha
  - **Pct_sig_vs_EUR**: Ns / Ntotal
  - **Risk_allele_index**: among significant SNPs, count where OR > 1 (Here OR is the OR for the risk allele (via risk_stats_for_snp()))
  - **Pct_risk_among_sig**: Risk_allele_index / Ns
  - **Risk_among_sig_fraction**: formatted R/Ns 
- With beta weighting, there will be two additional output columns giving the weighted absolute divergence from reference in risk allele frequency (weighted_abs_score_sig) + the normalized version (weighted_abs_score_sig_norm)
  - How to interpret beta-weighted results:
    - Weighted_abs_score_sig increases iff the population differs more from the reference population in risk allele frequency and/or those differences occur at SNPs with larger effect sizes
    - Weighted_abs_score_sig_norm is a “per significant SNP” scaled score; approximately interpretable on a 0–1-ish scale
    - **Note**: because we're using absolute values for both beta and AF difference, this is a divergence metric, not a “higher risk burden” metric.
