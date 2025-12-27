
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("populations_file")
parser.add_argument("snp_file")
parser.add_argument("prefix")
args = parser.parse_args()

populations_file = args.populations_file
snp_file = args.snp_file
pD= Path(__file__).parent

PREFIX = args.prefix
allele_count_file= f"{pD}/{PREFIX}_raw_scrape_data.csv"
output_filename = f"{pD}/{PREFIX}_grouped_AF.csv"

"""
example usage: 
python /Users/ainemacdermott/PersonalProjects/dyslipidemia/calc_allele_freqs.py /Users/ainemacdermott/PersonalProjects/dyslipidemia/populations_list.csv /Users/ainemacdermott/PersonalProjects/dyslipidemia/example_snp.csv /Users/ainemacdermott/PersonalProjects/dyslipidemia/example_intermediate_output.csv /Users/ainemacdermott/PersonalProjects/dyslipidemia/example_final_output.csv

"""

import pandas as pd

def main(
    snp_file: str = snp_file,
    populations_file: str = populations_file,
    allele_count_file: str = allele_count_file,
    output_filename: str = output_filename,
):

    # Read population grouping info and allele counts
    pops_df = pd.read_csv(populations_file)
    counts_df = pd.read_csv(allele_count_file)

    # Make sure these columns exist with the expected names
    required_cols_pops = ["dataset", "genetic_ancestry_group", "group_name", "subpop_label"]
    required_cols_counts = [
        "snp_id", "chrom", "pos", "ref", "alt",
        "dataset", "genetic_ancestry_group", "group_name", "subpop_label",
        "alternate_allele_count", "total_count", "homozygote_count"
    ]

    for col in required_cols_pops:
        if col not in pops_df.columns:
            raise ValueError(f"Populations file is missing required column: {col}")

    for col in required_cols_counts:
        if col not in counts_df.columns:
            raise ValueError(f"Allele count file is missing required column: {col}")

    # Ensure counts are numeric
    for col in ["alternate_allele_count", "total_count", "homozygote_count"]:
        counts_df[col] = pd.to_numeric(counts_df[col], errors="raise")

    # Keep only rows in counts_df that appear in populations_file
    merged = counts_df.merge(
        pops_df[required_cols_pops],  # just these columns are needed for the join
        on=["dataset", "genetic_ancestry_group", "group_name", "subpop_label"],
        how="inner"
    )

    # Group by SNP + dataset + ancestry + group_name (collapsing subpop_label)
    group_cols = [
        "snp_id", "chrom", "pos", "ref", "alt",
        "dataset", "genetic_ancestry_group", "group_name"
    ]

    agg = (
        merged
        .groupby(group_cols, as_index=False)
        .agg(
            aa_count=("alternate_allele_count", "sum"),
            total_count=("total_count", "sum"),
            alt_homozygote_count=("homozygote_count", "sum"),
        )
    )

    # Compute frequencies from the summed counts
    agg["allele_frequency"] = agg["aa_count"] / agg["total_count"]
    agg["alt_homozygote_frequency"] = agg["alt_homozygote_count"] / (agg["total_count"] / 2)
    agg["heterozygote_count"] = agg["aa_count"] - (2 * agg["alt_homozygote_count"])
    agg["ref_homozygote_count"] = (agg["total_count"] / 2) - agg["alt_homozygote_count"] - agg["heterozygote_count"]


    # Write final output
    agg.to_csv(output_filename, index=False)
    print ("Output written to", output_filename)

if __name__ == "__main__":
    main()
