import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("prefix")
args = parser.parse_args()
pD = Path(__file__).parent

prefix = args.prefix
input_file = f"{pD}/{prefix}_grouped_AF.csv"
af_output_file = f"{pD}/{prefix}_AF_reformatted.csv"
count_output_file = f"{pD}/{prefix}_allele_count_reformatted.csv"
genotype_output_file = f"{pD}/{prefix}_genotype_frequency_reformatted.csv"

input = pd.read_csv(input_file)
print(input.columns.tolist())
input.columns = input.columns.str.strip()

### basic allele freqs
new_df=(input.groupby(["snp_id", "dataset","group_name"])
       .agg(A=("aa_count", "sum"),
            B=("total_count", "sum")  
       ).reset_index()      
        .assign(total_allele_frequency=lambda d: d["A"] / d["B"])
        .pivot(index=["snp_id"], columns="group_name", values="total_allele_frequency")
).reset_index()
#print(new_df)
new_df.to_csv(af_output_file, index=False)

### basic allele counts
agg_counts = (input.groupby(["snp_id", "group_name"])
         .agg(ALT=("aa_count", "sum"),
              TOTAL=("total_count", "sum"),
              REF_ALLELE=("ref", "first"),
              ALT_ALLELE=("alt", "first"))
         .reset_index()
         .assign(REF=lambda d: d["TOTAL"] - d["ALT"])
)

count_long = pd.concat(
    [
        agg_counts.assign(Allele=lambda d: d["REF_ALLELE"], Count=lambda d: d["REF"]),
        agg_counts.assign(Allele=lambda d: d["ALT_ALLELE"], Count=lambda d: d["ALT"]),
    ],
    ignore_index=True,
)[["snp_id", "Allele", "group_name", "Count"]]

count_wide = (count_long.pivot_table(
        index=["snp_id", "Allele"],
        columns="group_name",
        values="Count",
        aggfunc="sum",
        fill_value=0
    ).reset_index()
)

pop_cols = [c for c in count_wide.columns if c not in ["snp_id", "Allele"]]
total_rows = (
    count_wide.groupby("snp_id", as_index=False)[pop_cols].sum()
             .assign(Allele="Total:")
)

count_df = (
    pd.concat([count_wide, total_rows], ignore_index=True)
      .sort_values(["snp_id", "Allele"], kind="stable")  # you can customize ordering later
      .reset_index(drop=True)
)

print(count_df)
count_df.to_csv(count_output_file, index=False)


### genotype counts
geno_agg = (
    input.groupby(["snp_id", "group_name"])
         .agg(
             REF=("ref", "first"),
             ALT=("alt", "first"),
             ref_hom=("ref_homozygote_count", "sum"),
             het=("heterozygote_count", "sum"),
             alt_hom=("alt_homozygote_count", "sum"),
             total_alleles=("total_count", "sum")
         )
         .reset_index()
         .assign(n_individuals=lambda d: d["total_alleles"] / 2)
)

geno_long = pd.concat(
    [geno_agg.assign(
            genotype=lambda d: d["REF"] + d["REF"],
            Frequency=lambda d: d["ref_hom"] / d["n_individuals"]
        ),
        geno_agg.assign(
            genotype=lambda d: d["REF"] + d["ALT"],
            Frequency=lambda d: d["het"] / d["n_individuals"]
        ),
        geno_agg.assign(
            genotype=lambda d: d["ALT"] + d["ALT"],
            Frequency=lambda d: d["alt_hom"] / d["n_individuals"]
        ),
    ],
    ignore_index=True,
)[["snp_id", "genotype", "group_name", "Frequency"]]

geno_wide = (geno_long.pivot_table(
            index=["snp_id", "genotype"],
            columns="group_name",
            values="Frequency",
            aggfunc="mean"
        ) .reset_index()
)

ref_alt = (
    input.groupby("snp_id")[["ref", "alt"]]
         .first()
)

def genotype_order(row):
    r, a = ref_alt.loc[row["snp_id"], ["ref", "alt"]]
    if row["genotype"] == r + r:
        return 0
    if row["genotype"] == r + a:
        return 1
    return 2

geno_wide["__order"] = geno_wide.apply(genotype_order, axis=1)

geno_wide = (
    geno_wide
        .sort_values(["snp_id", "__order"])
        .drop(columns="__order")
        .reset_index(drop=True)
)

geno_wide.to_csv(genotype_output_file, index=False)
print(geno_wide)