import argparse
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import csv
import math
import os
from collections import defaultdict


#settings to arrange:
reference_pop = "gnomAD_EUR"
correction = 0.5 # haldane correction
alpha = 0.05   # significance cutoff
nearby_window_bp = 50000

parser = argparse.ArgumentParser()
parser.add_argument("snp_file")
parser.add_argument("prefix", help=("Same prefix from 3_reframe_allele_freqs.py\n" 
                                                    "ex: 'final_output' will create 'final_output_AF_reformatted.csv' \n"
                                                    "'final_output_allele_count_reformatted.csv', \n"
                                                    " and 'final_output_genotype_frequency_reformatted.csv'"))
parser.add_argument("new_output_foldername", help=("Folder name for these output files; will be automatically created in the same folder as this script exists."))
parser.add_argument("beta_weighting", help=("Whether to apply beta-weighting to risk allele counts: 'yes' or 'no'; default 'no'"), nargs='?', default='no')
args = parser.parse_args()

pD = os.path.dirname(os.path.abspath(__file__))
snp_file = args.snp_file
actual_output_file_folder = (pD + "/" + args.new_output_foldername) if args.new_output_foldername else (pD + "stats")
input_intermeds_filename = args.prefix
input_file = (pD + "/" + input_intermeds_filename + "_AF_reformatted.csv")
input_file2 = (pD + "/" + input_intermeds_filename + "_grouped_AF.csv")
beta_weighting = False if args.beta_weighting.lower() != 'yes' and args.beta_weighting.lower() != 'y' and args.beta_weighting.lower() == None else True

input = pd.read_csv(input_file)
snp = pd.read_csv(snp_file)
af_output_file = pd.read_csv(pD + "/" + input_intermeds_filename + "_AF_reformatted.csv")
count_output_filepath = (pD + "/" + input_intermeds_filename + "_allele_count_reformatted.csv")
count_output_file = pd.read_csv(count_output_filepath)
genotype_output_file = pd.read_csv(pD + "/" + input_intermeds_filename + "_genotype_frequency_reformatted.csv")
print(input.columns.tolist())


def snp_loci(snp_file):#finish this later
    snps = []
    with open(snp_file, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            snps.append({
                "snp_id": row["snp_id"].strip(),
                "chrom": row["chromosome"].strip(),
                "pos": int(row["position"]),
                "beta": float(row["beta"]),
            })

    # sort by chrom and then pos
    snps.sort(key=lambda x: (x["chrom"], x["pos"]))

    retained = []
    i = 0
    while i < len(snps):
        cluster = [snps[i]]
        j = i + 1
        while j < len(snps) and snps[j]["chrom"] == snps[i]["chrom"]:
            if snps[j]["pos"] - snps[j - 1]["pos"] <= nearby_window_bp:
                cluster.append(snps[j])
                j += 1
            else:
                break

        # pick snp with highest abs value beta in the cluster --> greatest effect regardless of direction
        best = max(cluster, key=lambda x: abs(x["beta"]))
        retained.append(best["snp_id"])

        # next unclustered snp
        i = j

    return retained

# mathy fxs
def norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def z_pvalue_two_sided(z: float) -> float:
    if z is None or math.isnan(z):
        return float("nan")
    return 2.0 * (1.0 - norm_cdf(abs(z)))

def or_and_stats(a, b, c, d):
    """
    for a 2x2 table:
        comparison popn: effect=a, other=b
        reference popn: effect=c, other=d

    applies correction (for haldane, add 0.5 to each)
    + return OR, SE, logOR, lowerCL, upperCL, z, p for z test
    """
    a = a + correction
    b = b + correction
    c = c + correction
    d = d + correction

    or_val = (a * d) / (b * c)
    log_or = math.log(or_val)
    se = math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    lower = math.exp(log_or - 1.96 * se)
    upper = math.exp(log_or + 1.96 * se)
    z = log_or / se if se else float("nan")
    p = z_pvalue_two_sided(z)
    return or_val, se, lower, upper, log_or, z, p

def fisher_exact_two_sided(a, b, c, d):
    """
    two-sided Fisher exact p-value for a 2x2 table.
    """
    try:
        _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        return p
    except Exception:
        return float("nan")
    
def chi_square_2x2(a, b, c, d):
    chi2, p, dof, expected = chi2_contingency([[a, b], [c, d]], correction=False)
    # no correction because of haldane already applied
    return chi2, p, dof, expected

def chi_square_2x3(row_pop, row_ref):
    chi2, p, dof, expected = chi2_contingency([row_pop, row_ref], correction=False)
    return chi2, p, dof, expected

# input parsing fx
def load_snp_metadata(snp_file):
    meta = {}
    with open(snp_file, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            snp_id = row["snp_id"].strip()
            ref = row["reference_allele"].strip()
            alt = row["alternate_allele"].strip()
            ea = row["EA"].strip()
            meta[snp_id] = {"ref": ref, "alt": alt, "effect": ea}
    return meta

def load_snp_metadata_with_risk(snp_file):
    """
      directionality = effect of EA on trait of interest
        '+' => EA increases trait => EA is the RISK allele
        '-' => EA decreases trait => the OTHER allele is the RISK allele
    """
    meta = {}
    with open(snp_file, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            snp_id = row["snp_id"].strip()
            ref = row["reference_allele"].strip()
            alt = row["alternate_allele"].strip()

            ea_flag = row["EA"].strip().lower()
            if ea_flag == "y":
                ea = alt
                other = ref
            elif ea_flag == "n":
                ea = ref
                other = alt
            else:
                raise ValueError(f"{snp_id}: EA must be y/n, got '{row['EA']}'")

            direction = row["directionality"].strip()
            if direction not in {"+", "-"}:
                raise ValueError(f"{snp_id}: directionality must be +/-, got '{direction}'")

            risk = ea if direction == "+" else other

            meta[snp_id] = {
                "ref": ref,
                "alt": alt,
                "ea": ea,
                "risk": risk,
                "directionality": direction,
                "beta": float(row["beta"]) if row.get("beta", "").strip() != "" else float("nan"),
            }
    return meta

def parse_count_output(count_output_file):
    """
    where
    popn_cols: ordered list of population columns
    snp_rows: dict snp_id -> dict allele_label -> {popn: count}
        where allele_label is "C", "T", or "Total:"
    """
    with open(count_output_file, newline="") as f:
        r = csv.DictReader(f)
        fieldnames = r.fieldnames
        popn_cols = fieldnames[2:]

        snp_rows = defaultdict(dict)
        for row in r:
            snp_id = row["snp_id"].strip()
            allele_label = row["Allele"].strip() 
            popn_counts = {}
            for popn in popn_cols:
                raw = row[popn].strip()
                popn_counts[popn] = int(raw) if raw != "" else 0
            snp_rows[snp_id][allele_label] = popn_counts
    return popn_cols, snp_rows

def parse_genotype_long_counts(input_file2):
    """
      pops: sorted list of population labels found (group_name)
      geno_counts: dict snp_id -> dict pop -> dict with keys:
          ref_hom, het, alt_hom
    """
    geno_counts = defaultdict(dict)
    pops = set()

    with open(input_file2, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            snp_id = row["snp_id"].strip()

            pop = row["group_name"].strip()  # important: strip spaces like " gnomAD_SAS"
            pops.add(pop)

            # counts (sometimes appear as floats like "3.0")
            ref_hom = int(round(float(row["ref_homozygote_count"])))
            het = int(round(float(row["heterozygote_count"])))
            alt_hom = int(round(float(row["alt_homozygote_count"])))

            geno_counts[snp_id][pop] = {
                "ref_hom": ref_hom,
                "het": het,
                "alt_hom": alt_hom,
            }

    return sorted(pops), geno_counts

# put together
def allele_count_analysis(snp_id, popn_cols, snp_rows_for_id, snp_meta, out_dir):
    ref_allele = snp_meta[snp_id]["ref"]
    effect_allele = snp_meta[snp_id]["effect"]

    if "Total:" not in snp_rows_for_id:
        raise ValueError(f"{snp_id}: missing 'Total:' row in count file.")

    total_counts = snp_rows_for_id["Total:"]

    if effect_allele =="y":
        effect_allele = snp_meta[snp_id]["alt"]
        other_allele = snp_meta[snp_id]["ref"]
    elif effect_allele =="n":
        effect_allele = snp_meta[snp_id]["ref"]
        other_allele = snp_meta[snp_id]["alt"]

    eff_counts = snp_rows_for_id.get(effect_allele)
    other_counts = snp_rows_for_id.get(other_allele)

    if reference_pop not in popn_cols:
        raise ValueError(f"Reference population column '{reference_pop}' not found in count file header.")

    eff_ref = eff_counts[reference_pop]
    other_ref = total_counts[reference_pop] - eff_ref
    # prep output rows
    rows = {
        f"other_allele_count ({other_allele})": [],
        f"effect_allele_count ({effect_allele})": [],
        "total_allele_count": [],
        "OR_vs_gnomAD_ref": [],
        "SE": [],
        "logOR": [],
        "lowerCL": [],
        "upperCL": [],
        "z": [],
        "p_fisher": [],
        "chi2": [],
        "p_chi2": [],
        "min_expected": [],
        "primary_p": [],
        "primary_test": [],
    }

    for pop in popn_cols:
        other_ct = other_counts[pop]
        eff_ct = eff_counts[pop]
        tot_ct = total_counts[pop]

        rows[f"other_allele_count ({other_allele})"].append(other_ct)
        rows[f"effect_allele_count ({effect_allele})"].append(eff_ct)
        rows["total_allele_count"].append(tot_ct)

        if pop == reference_pop:
            rows["OR_vs_gnomAD_ref"].append(1.0)
            rows["SE"].append(float("nan"))
            rows["logOR"].append(0.0)
            rows["lowerCL"].append(float("nan"))
            rows["upperCL"].append(float("nan"))
            rows["z"].append(float("nan"))
            rows["p_fisher"].append(float("nan"))
            rows["chi2"].append(float("nan"))
            rows["p_chi2"].append(float("nan"))
            rows["min_expected"].append(float("nan"))
            rows["primary_p"].append(float("nan"))
            rows["primary_test"].append("NA")
            continue

        other_pop = tot_ct - eff_ct

        chi2, p_chi, _, expected = chi_square_2x2(eff_ct, other_pop, eff_ref, other_ref)
        min_exp = float(expected.min())
        pf = fisher_exact_two_sided(eff_ct, other_pop, eff_ref, other_ref)

        use_fisher = min_exp < 5
        primary_p = pf if use_fisher else p_chi
        primary_test = "fisher" if use_fisher else "chi2"


        OR, SE, logOR, lower, upper, z, pz = or_and_stats(
            a=eff_ct, b=other_pop, c=eff_ref, d=other_ref
        )
        pf = fisher_exact_two_sided(
            a=eff_ct, b=other_pop, c=eff_ref, d=other_ref
        )

        rows["OR_vs_gnomAD_ref"].append(OR)
        rows["SE"].append(SE)
        rows["logOR"].append(logOR)
        rows["lowerCL"].append(lower)
        rows["upperCL"].append(upper)
        rows["z"].append(z)
        rows["p_fisher"].append(pf)
        rows["chi2"].append(chi2)
        rows["p_chi2"].append(p_chi)
        rows["min_expected"].append(min_exp)
        rows["primary_p"].append(primary_p)
        rows["primary_test"].append(primary_test)

    # write output
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{snp_id}_allele_count_analysis_output.csv")

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric"] + popn_cols)
        for metric, values in rows.items():
            w.writerow([metric] + values)

    return out_path

def genotype_count_analysis_long(
    snp_id,
    popn_cols,            
    geno_counts_for_id,   # geno_counts[snp_id]
    out_dir
):
    """
    2x3 chi2 most of the time; backup = fisher on carrier vs non-carrier 2x2
    """
    if reference_pop not in popn_cols:
        raise ValueError(f"Reference population '{reference_pop}' must be in popn_cols.")

    if reference_pop not in geno_counts_for_id:
        raise ValueError(f"{snp_id}: missing genotype counts for reference pop '{reference_pop}' in input_file2.")

    ref_refhom = geno_counts_for_id[reference_pop]["ref_hom"]
    ref_het = geno_counts_for_id[reference_pop]["het"]
    ref_althom = geno_counts_for_id[reference_pop]["alt_hom"]

    rows = {
        "ref_homozygote_count": [],
        "heterozygote_count": [],
        "alt_homozygote_count": [],
        "chi2_2x3": [],
        "p_chi2_2x3": [],
        "min_expected_2x3": [],
        "p_fisher_carrier_2x2": [],
        "primary_p": [],
        "primary_test": [],
    }

    for pop in popn_cols:
        if pop not in geno_counts_for_id:
            rows["ref_homozygote_count"].append(float("nan"))
            rows["heterozygote_count"].append(float("nan"))
            rows["alt_homozygote_count"].append(float("nan"))
            rows["chi2_2x3"].append(float("nan"))
            rows["p_chi2_2x3"].append(float("nan"))
            rows["min_expected_2x3"].append(float("nan"))
            rows["p_fisher_carrier_2x2"].append(float("nan"))
            rows["primary_p"].append(float("nan"))
            rows["primary_test"].append("MISSING")
            continue

        pop_refhom = geno_counts_for_id[pop]["ref_hom"]
        pop_het = geno_counts_for_id[pop]["het"]
        pop_althom = geno_counts_for_id[pop]["alt_hom"]

        rows["ref_homozygote_count"].append(pop_refhom)
        rows["heterozygote_count"].append(pop_het)
        rows["alt_homozygote_count"].append(pop_althom)

        if pop == reference_pop:
            rows["chi2_2x3"].append(float("nan"))
            rows["p_chi2_2x3"].append(float("nan"))
            rows["min_expected_2x3"].append(float("nan"))
            rows["p_fisher_carrier_2x2"].append(float("nan"))
            rows["primary_p"].append(float("nan"))
            rows["primary_test"].append("NA")
            continue

        chi2, p_chi, _, expected = chi_square_2x3(
            [pop_refhom, pop_het, pop_althom],
            [ref_refhom, ref_het, ref_althom],
        )
        min_exp = float(expected.min())

        # fisher fallback
        carrier_pop = pop_het + pop_althom
        noncarrier_pop = pop_refhom
        carrier_ref = ref_het + ref_althom
        noncarrier_ref = ref_refhom

        p_fisher_carrier = fisher_exact_two_sided(
            carrier_pop, noncarrier_pop,
            carrier_ref, noncarrier_ref
        )

        use_fisher = min_exp < 5
        primary_p = p_fisher_carrier if use_fisher else p_chi
        primary_test = "fisher_2x2_carrier" if use_fisher else "chi2_2x3"

        rows["chi2_2x3"].append(chi2)
        rows["p_chi2_2x3"].append(p_chi)
        rows["min_expected_2x3"].append(min_exp)
        rows["p_fisher_carrier_2x2"].append(p_fisher_carrier)
        rows["primary_p"].append(primary_p)
        rows["primary_test"].append(primary_test)

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{snp_id}_genotype_count_analysis_output.csv")

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric"] + popn_cols)
        for metric, values in rows.items():
            w.writerow([metric] + values)

    return out_path

# risk stats
def risk_stats_for_snp(snp_id, popn_cols, snp_rows_for_id, snp_meta):
    """
    Returns dict pop -> dict with:
      OR, primary_p, primary_test, risk_af, risk_af_ref
    """
    if reference_pop not in popn_cols:
        raise ValueError(f"Reference pop '{reference_pop}' not in count file columns.")

    if "Total:" not in snp_rows_for_id:
        raise ValueError(f"{snp_id}: missing Total: row.")

    total = snp_rows_for_id["Total:"]  # pop -> total alleles

    risk_allele = snp_meta[snp_id]["risk"]
    if risk_allele not in snp_rows_for_id:
        raise ValueError(f"{snp_id}: missing risk allele row '{risk_allele}' in count file.")

    risk_counts = snp_rows_for_id[risk_allele]  # pop -> risk allele count

    # reference cells
    a_ref = risk_counts[reference_pop]
    b_ref = total[reference_pop] - a_ref

    out = {}

    for pop in popn_cols:
        a = risk_counts[pop]
        b = total[pop] - a

        # risk allele frequency
        risk_af = a / total[pop] if total[pop] else float("nan")
        risk_af_ref = a_ref / total[reference_pop] if total[reference_pop] else float("nan")

        if pop == reference_pop:
            out[pop] = {
                "OR": 1.0,
                "primary_p": float("nan"),
                "primary_test": "NA",
                "risk_af": risk_af,
                "risk_af_ref": risk_af_ref,
            }
            continue

        chi2, p_chi, _, expected = chi_square_2x2(a, b, a_ref, b_ref)
        min_exp = float(expected.min())
        p_fisher = fisher_exact_two_sided(a, b, a_ref, b_ref)

        use_fisher = (min_exp < 5)
        primary_p = p_fisher if use_fisher else p_chi
        primary_test = "fisher_2x2" if use_fisher else "chi2_2x2"

        OR, _, _, _, _, _, _ = or_and_stats(a, b, a_ref, b_ref)

        out[pop] = {
            "OR": OR,
            "primary_p": primary_p,
            "primary_test": primary_test,
            "risk_af": risk_af,
            "risk_af_ref": risk_af_ref,
        }

    return out

def cumulative_risk_allele_summary(
    count_output_filepath,
    snp_file,
    populations_of_interest=None,
    alpha=0.05,
    out_csv_path=None,
):
    snp_meta = load_snp_metadata_with_risk(snp_file)
    popn_cols, snp_rows = parse_count_output(count_output_filepath)

    if populations_of_interest is None:
        populations_of_interest = popn_cols
    populations_of_interest = [p for p in populations_of_interest if p in popn_cols]

    snps = [s for s in snp_meta.keys() if s in snp_rows]

    sig_counts = {p: 0 for p in populations_of_interest}
    risk_enriched_counts = {p: 0 for p in populations_of_interest}

    for snp_id in snps:
        stats = risk_stats_for_snp(snp_id, popn_cols, snp_rows[snp_id], snp_meta)
        for pop in populations_of_interest:
            if pop == reference_pop:
                continue
            pval = stats[pop]["primary_p"]
            OR = stats[pop]["OR"]
            if pval is not None and not math.isnan(pval) and pval < alpha:
                sig_counts[pop] += 1
                if OR > 1:
                    risk_enriched_counts[pop] += 1

    n_total = len(snps)

    rows = []
    for pop in populations_of_interest:
        if pop == reference_pop:
            rows.append({
                "population": pop,
                "Ns_sig_vs_ref": f"Reference ({n_total} SNPs)",
                "pct_sig_vs_ref": "",
                "risk_allele_index": "",
                "pct_risk_among_sig": "",
                "risk_among_sig_fraction": "",
            })
            continue

        Ns = sig_counts[pop]
        R = risk_enriched_counts[pop]

        pct_sig = (Ns / n_total * 100) if n_total else float("nan")
        pct_risk = (R / Ns * 100) if Ns else float("nan")

        rows.append({
            "population": pop,
            "Ns_sig_vs_ref": f"{Ns}/{n_total}",
            "pct_sig_vs_ref": f"{pct_sig:.1f}%",
            "risk_allele_index": R,
            "pct_risk_among_sig": f"{pct_risk:.1f}%" if Ns else "0.0%",
            "risk_among_sig_fraction": f"{R}/{Ns}" if Ns else "0/0",
        })

    if out_csv_path:
        fieldnames = [
            "population",
            "Ns_sig_vs_ref",
            "pct_sig_vs_ref",
            "risk_allele_index",
            "pct_risk_among_sig",
            "risk_among_sig_fraction",
        ]
        with open(out_csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

    return rows

def cumulative_risk_allele_summary_beta_weighting(
    count_output_filepath,
    snp_file,
    populations_of_interest=None,
    alpha=0.05,
    out_csv_path=None,
):
    """
    Significant-only, abs(beta)-weighted divergence score.

    For each population pop != reference_pop:
      - Ns = number of SNPs with primary_p < alpha vs reference_pop
      - risk_allele_index = among significant SNPs, count with OR > 1
      - weighted_abs_score_sig = sum_{sig SNPs} abs(beta) * abs(AF_pop - AF_refref)
      - weighted_abs_score_sig_norm = weighted_abs_score_sig / (beta_max * Ns)
        (per-pop normalization; avoids undefined sig_snps/raw_score)
    """
    snp_meta = load_snp_metadata_with_risk(snp_file)
    popn_cols, snp_rows = parse_count_output(count_output_filepath)

    if populations_of_interest is None:
        populations_of_interest = popn_cols
    populations_of_interest = [p for p in populations_of_interest if p in popn_cols]

    snps = [s for s in snp_meta.keys() if s in snp_rows]

    sig_counts = {p: 0 for p in populations_of_interest}
    risk_enriched_counts = {p: 0 for p in populations_of_interest}
    weighted_abs_score_sig = {p: 0.0 for p in populations_of_interest}

    # Compute beta_max over SNPs that actually have valid betas (and are in snps list)
    betas = []
    for s in snps:
        b = snp_meta[s].get("beta", float("nan"))
        if b is not None and not math.isnan(b):
            betas.append(abs(b))
    beta_max = max(betas) if betas else float("nan")

    for snp_id in snps:
        beta = snp_meta[snp_id].get("beta", float("nan"))
        abs_beta = None if (beta is None or math.isnan(beta)) else abs(beta)

        stats = risk_stats_for_snp(snp_id, popn_cols, snp_rows[snp_id], snp_meta)

        for pop in populations_of_interest:
            if pop == reference_pop:
                continue

            pval = stats[pop]["primary_p"]
            OR = stats[pop]["OR"]

            is_sig = (pval is not None and not math.isnan(pval) and pval < alpha)
            if not is_sig:
                continue

            sig_counts[pop] += 1
            if OR > 1:
                risk_enriched_counts[pop] += 1

            if abs_beta is None:
                continue

            af_pop = stats[pop].get("risk_af")
            af_refref = stats[pop].get("risk_af_ref")
            if (af_pop is None or af_refref is None or
                math.isnan(af_pop) or math.isnan(af_refref)):
                continue

            delta = abs(af_pop - af_refref)
            weighted_abs_score_sig[pop] += abs_beta * delta

    n_total = len(snps)

    rows = []
    for pop in populations_of_interest:
        if pop == reference_pop:
            rows.append({
                "population": pop,
                "Ns_sig_vs_ref": f"Reference ({n_total} SNPs)",
                "pct_sig_vs_ref": "",
                "risk_allele_index": "",
                "pct_risk_among_sig": "",
                "risk_among_sig_fraction": "",
                "weighted_abs_score_sig": "",
                "weighted_abs_score_sig_norm": "",
            })
            continue

        Ns = sig_counts[pop]
        R = risk_enriched_counts[pop]

        pct_sig = (Ns / n_total * 100) if n_total else float("nan")
        pct_risk = (R / Ns * 100) if Ns else 0.0

        raw = weighted_abs_score_sig[pop]

        # per-pop normalization; avoid divide-by-zero
        norm = (raw / (beta_max * Ns)) if (Ns > 0 and beta_max and not math.isnan(beta_max)) else float("nan")

        rows.append({
            "population": pop,
            "Ns_sig_vs_ref": f"{Ns}/{n_total}",
            "pct_sig_vs_ref": f"{pct_sig:.1f}%",
            "risk_allele_index": R,
            "pct_risk_among_sig": f"{pct_risk:.1f}%" if Ns else "0.0%",
            "risk_among_sig_fraction": f"{R}/{Ns}" if Ns else "0/0",
            "weighted_abs_score_sig": f"{raw:.6f}",
            "weighted_abs_score_sig_norm": f"{norm:.6f}" if not math.isnan(norm) else "nan",
        })

    if out_csv_path:
        fieldnames = [
            "population",
            "Ns_sig_vs_ref",
            "pct_sig_vs_ref",
            "risk_allele_index",
            "pct_risk_among_sig",
            "risk_among_sig_fraction",
            "weighted_abs_score_sig",
            "weighted_abs_score_sig_norm",
        ]
        with open(out_csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

    return rows

def main(count_output_file, snp_file, actual_output_file_folder):
    snp_list = snp_loci(snp_file)                 
    snp_meta = load_snp_metadata(snp_file)       
    popn_cols, snp_rows = parse_count_output(count_output_file)
    pops_from_geno, geno_counts = parse_genotype_long_counts(input_file2)

    for snp_id in snp_list:
        if snp_id not in snp_rows:
            print(f"Skipping {snp_id}: not found in count_output_file")
            continue

        # alleles
        out_path = allele_count_analysis(
            snp_id=snp_id,
            popn_cols=popn_cols,
            snp_rows_for_id=snp_rows[snp_id],
            snp_meta=snp_meta,
            out_dir=actual_output_file_folder
        )
        print(f"Wrote: {out_path}")

        # genos
        geno_out = genotype_count_analysis_long(
        snp_id=snp_id,
        popn_cols=popn_cols,              # from allele count file
        geno_counts_for_id=geno_counts[snp_id],
        out_dir=actual_output_file_folder
        )
        print(f"Wrote: {geno_out}")

    print(f"Risk analysis files located at {actual_output_file_folder}.")

    risk_analysis_output_file = (actual_output_file_folder + "/" + input_intermeds_filename + "_risk_analysis_output.csv")
    beta_risk_analysis_output_file = (actual_output_file_folder + "/" + input_intermeds_filename + "_risk_analysis_output_beta_weighted.csv")

    if beta_weighting == False:
        cumulative_risk_allele_summary(
        count_output_filepath,
        snp_file,
        populations_of_interest=None,   # e.g. ["1KG_EUR", "1KG_AFR", ...] or ["ASW","CHS"...] if present in columns
        alpha=0.05, 
        out_csv_path=risk_analysis_output_file)
    if beta_weighting == True:
        cumulative_risk_allele_summary_beta_weighting(
        count_output_filepath,
        snp_file,
        populations_of_interest=None,   # e.g. ["1KG_EUR", "1KG_AFR", ...] or ["ASW","CHS"...] if present in columns
        alpha=0.05, 
        out_csv_path=beta_risk_analysis_output_file)
    print(f"Wrote cumulative risk analysis to {risk_analysis_output_file}.")

main(count_output_filepath, snp_file, actual_output_file_folder)
