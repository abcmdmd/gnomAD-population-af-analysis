"""
This script requires TWO input files. I recommend writing out all of the information into an Excel file and then downloading
as CSVs to ensure proper formatting. Please refer to "example_snp.csv" and "example_populations_list.csv" for formatting examples.

1. snp info: 
    - snp_id
    - chromosome
    - position
    - reference allele
    - alternate allele

2. populations list:
    - datasets: (choose between gnomAD, HGDP, and 1KG)
    - format: "dataset", "genetic ancestry group, "name of group", "subpopn as listed on website"
    - HGDP,East Asian,OSEA,Cambodian
    - HGDP,East Asian,OSEA,Lahu

    Note: make sure there are no trailing spaces/tabs. This will prevent proper matching + data scraping.

** gnomAD URL pattern for variant lookup:
    template: https://gnomad.broadinstitute.org/variant/4-99318162-T-C?dataset=gnomad_r4
    example:  https://gnomad.broadinstitute.org/variant/{chr}-{position}-{ref}-{alt}?dataset=gnomad_r4

** resource to look up SNP info:
    If you're having any trouble matching SNP info or finding the right alleles/positions, use NCBI's SNP database:
    template: https://www.ncbi.nlm.nih.gov/snp/rs1229984
    example: https://www.ncbi.nlm.nih.gov/snp

"""

import argparse
import time
import csv
from pathlib import Path

import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

parser = argparse.ArgumentParser()
parser.add_argument("populations_file")
parser.add_argument("snp_file")
#parser.add_argument("output_filename")
parser.add_argument("prefix")
pD = Path(__file__).parent
args = parser.parse_args()

POPS_CSV = args.populations_file
SNPS_CSV = args.snp_file
#OUTPUT_CSV = args.output_filename
PREFIX = args.prefix
OUTPUT_CSV= f"{pD}/{PREFIX}_raw_scrape_data.csv"

# config
GNOMAD_BASE_URL = "https://gnomad.broadinstitute.org/variant/{chrom}-{pos}-{ref}-{alt}?dataset=gnomad_r4"

# map dataset label to the text that appears on the tabs in the UI
DATASET_TAB_LABELS = {
    "gnomAD": "gnomAD",           
    "HGDP": "HGDP",              
    "1KG": "1KG",               
}

# OR, you can hardcode paths here:
"""
POPS_CSV = "/Users/ainemacdermott/PersonalProjects/dyslipidemia/test_populations_list.csv"
SNPS_CSV = "/Users/ainemacdermott/PersonalProjects/dyslipidemia/test_snp.csv"
OUTPUT_CSV = "/Users/ainemacdermott/PersonalProjects/dyslipidemia/test_output.csv
"""

# don't overhwlem w requests
PER_VARIANT_DELAY = 2

# selenium setup
def make_driver(headless: bool = True):
    options = Options()
    if headless:
        options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    driver = webdriver.Chrome(options=options)
    driver.set_page_load_timeout(60)
    return driver

# functions
def open_genetic_ancestry_table(driver, dataset_label: str):
    """
     already-loaded gnomAD variant page:
      - find the 'Genetic Ancestry Group Frequencies' section
      - click into dataset tab (gnomAD / HGDP / 1KG)
      - return <table> inside dataset tab panel
    """

    wait = WebDriverWait(driver, 20)

    #  right H2
    section = wait.until(
        EC.presence_of_element_located(
            (
                By.XPATH,
                "//section[h2[contains(normalize-space(), 'Genetic Ancestry Group Frequencies')]]"
            )
        )
    )

    # find dataset_label
    tab = section.find_element(
        By.XPATH,
        (
            ".//li[@role='tab'][.//div[normalize-space()="
            f"'{dataset_label}']]"
        ),
    )
    tab.click()
    panel_id = tab.get_attribute("aria-controls")
    panel = wait.until(
        EC.presence_of_element_located(
            (By.XPATH, f".//div[@id='{panel_id}']")
        )
    )

    # find ancestry table
    table = panel.find_element(By.XPATH, ".//table")

    # Debug dump for just this table, if needed:
    # debug_dump_table_html(table, filename=f"debug_genetic_ancestry_{dataset_label}.html")

    return table

def safe_filename(name: str) -> str:
    # replace anything scary (slashes, spaces) with underscores
    return name.replace("/", "_").replace(" ", "_")


def debug_dump_table_html(table_element, filename="debug_table.html"):
    html = table_element.get_attribute("outerHTML")
    with open(filename, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[DEBUG] Wrote table HTML to {filename}")

def parse_table_row_counts(table_element, subpop_label: str, debug: bool = True):
    """
    Given a <table> element and a desired label somewhere in a row,
    return (alternate_allele_count, total_count, homozygote_count).

    We:
      - look for any row where *any* cell (th or td) == subpop_label
      - then take the next three numeric cells in that row as:
          alternate_allele_count, total_count, homozygote_count
    """
    rows = table_element.find_elements(By.XPATH, ".//tbody/tr")

    #if debug:
        #print(f"[DEBUG] Looking for subpop_label {subpop_label!r}")
        #print("[DEBUG] All row texts in this table:")

    for row in rows:
        cells = row.find_elements(By.XPATH, ".//th | .//td")
        if not cells:
            continue

        texts = [c.text.strip() for c in cells]

        if debug and any(texts):
            print("  ROW:", texts)

        # is our label anywhere in this row?
        if any(t == subpop_label for t in texts):
            label_idx = texts.index(subpop_label)

            def parse_int_safe(x):
                try:
                    return int(x.replace(",", ""))
                except ValueError:
                    return None

            alt = total = hom = None
            count_texts = texts[label_idx + 1: label_idx + 4]

            if len(count_texts) > 0:
                alt = parse_int_safe(count_texts[0])
            if len(count_texts) > 1:
                total = parse_int_safe(count_texts[1])
            if len(count_texts) > 2:
                hom = parse_int_safe(count_texts[2])

            return alt, total, hom

    if debug:
        print(f"[DEBUG] subpop_label {subpop_label!r} not found in this table.")
    return None

def get_group_overall_counts(table_element, ancestry_group: str, debug: bool = True):
    """
    For gnomAD dataset: return counts for the ancestry_group row itself
    (e.g., 'East Asian', 'African', 'Central/South Asian'), ignoring any 'Overall' labels.
    """
    target = ancestry_group.strip().lower()
    rows = table_element.find_elements(By.XPATH, ".//tbody/tr")

    if debug:
        print(f"[DEBUG] Looking for ancestry_group row {ancestry_group!r} in gnomAD table")

    for row in rows:
        try:
            label_el = row.find_element(By.XPATH, ".//th[1]")
            label_text = label_el.text.strip()
        except Exception:
            continue

        if debug and label_text:
            print("  ROW LABEL:", repr(label_text))

        label_norm = label_text.lower()

        # primary: exact match after strip/lower
        if label_norm == target or target in label_norm:
            tds = row.find_elements(By.TAG_NAME, "td")

            def parse_int_safe(text):
                try:
                    return int(text.replace(",", ""))
                except ValueError:
                    return None

            alt = parse_int_safe(tds[0].text) if len(tds) > 0 else None
            total   = parse_int_safe(tds[1].text) if len(tds) > 1 else None
            hom     = parse_int_safe(tds[2].text) if len(tds) > 2 else None

            return alt, total, hom

    if debug:
        print(f"[DEBUG] ancestry_group row {ancestry_group!r} not found in gnomAD table.")
    return None

            

    if debug:
        print(f"[DEBUG] ancestry_group row {ancestry_group!r} not found in gnomAD table.")
    return None


def get_counts_for_population(driver, dataset_label: str, ancestry_group: str, subpop_label: str, debug: bool = True):
    """
    For a given dataset (gnomAD / HGDP / 1KG), ancestry_group ('East Asian', etc),
    and subpop_label ('Cambodian', 'Lahu', 'Overall', etc), return
      (alternate_allele_count, total_count, homozygote_count).

    Logic:
      1) Open the dataset-level Genetic Ancestry Group Frequencies table.
      2) Try to parse subpop_label directly from that table.
         - If it’s already expanded (from a previous call), this should succeed.
      3) If not found, click the ancestry_group row to expand it.
      4) Re-open the table and try parsing again.
    """

    if debug:
        print(f"[DEBUG] Getting counts for dataset={dataset_label!r}, "
              f"ancestry_group={ancestry_group!r}, subpop_label={subpop_label!r}")

    # Step 1: open the Genetic Ancestry Group Frequencies table for this dataset
    table = open_genetic_ancestry_table(driver, dataset_label=dataset_label)

    # Step 2: first try to parse directly (works if group is already expanded)
    counts = parse_table_row_counts(table, subpop_label=subpop_label, debug=debug)
    if counts is not None:
        if debug:
            print(f"[DEBUG] Found {subpop_label!r} without expanding ancestry_group {ancestry_group!r}.")
        return counts

    if debug:
        print(f"[DEBUG] {subpop_label!r} not found yet. Expanding ancestry_group {ancestry_group!r}...")

    # Step 3: click ancestry_group row (e.g. 'East Asian') to expand
    try:
        group_row = table.find_element(
            By.XPATH,
            f".//tr[.//th[normalize-space()='{ancestry_group}']]"
        )
    except Exception as e:
        if debug:
            print(f"[DEBUG] Could not find ancestry_group row {ancestry_group!r}: {e}")
        return None

    try:
        button = group_row.find_element(By.XPATH, ".//button")
        button.click()
    except Exception:
        # if no button, try clicking the whole row
        group_row.click()

    # Let the DOM update
    time.sleep(1)

    # Step 4: re-open the table (the same panel, now with the group expanded)
    table = open_genetic_ancestry_table(driver, dataset_label=dataset_label)

    # Optional: dump this expanded table for debugging
    fname = f"debug_nested_{dataset_label}_{ancestry_group}_{subpop_label}.html"
    newname = safe_filename(fname)
    debug_dump_table_html(table, filename=newname)


    # Step 5: try parsing again; now the subpop rows (Lahu, Cambodian, etc.) should be visible
    return parse_table_row_counts(table, subpop_label=subpop_label, debug=debug)

#### MAIN MAIN MAIN

def main(snps_csv=SNPS_CSV, pops_csv=POPS_CSV, output_csv=OUTPUT_CSV):
    snps_df = pd.read_csv(snps_csv)
    pops_df = pd.read_csv(pops_csv)

    # sanity check: required columns
    for col in ["chromosome", "position", "reference_allele", "alternate_allele"]:
        if col not in snps_df.columns:
            raise ValueError(f"SNV file is missing required column: {col}")

    for col in ["dataset", "genetic_ancestry_group", "group_name", "subpop_label"]:
        if col not in pops_df.columns:
            raise ValueError(f"Population file is missing required column: {col}")

    driver = make_driver(headless=True)
    output_rows = []

    try:
        for _, snp in snps_df.iterrows():
            chrom = str(snp["chromosome"])
            pos = str(snp["position"])
            ref = str(snp["reference_allele"])
            alt = str(snp["alternate_allele"])
            snp_id = snp.get("snp_id", "")

            variant_url = GNOMAD_BASE_URL.format(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
            )

            print(f"\n--- {snp_id}: {chrom}-{pos}-{ref}-{alt} ---")
            print(f"URL: {variant_url}")

            # Load variant page
            driver.get(variant_url)
            time.sleep(3)  # give time for JS to hydrate

            # For each requested pop/group for this SNP
            for _, row in pops_df.iterrows():
                dataset = str(row["dataset"])                    # 'gnomAD', 'HGDP', or '1KG'
                ancestry_group = str(row["genetic_ancestry_group"])
                group_name = str(row["group_name"])
                subpop_label = str(row["subpop_label"])

                try:
                    if dataset in ["gnomAD", "1KG"] and subpop_label == "Overall":
                        # Use group-level row from gnomAD, no nested subpops
                        table = open_genetic_ancestry_table(driver, dataset_label=dataset)
                        counts = get_group_overall_counts(
                            table_element=table,
                            ancestry_group=ancestry_group,
                            debug=True,
                        )
                    else:
                        # Use HGDP/1KG-style logic: select dataset, expand ancestry_group, then find subpop_label
                        counts = get_counts_for_population(
                            driver=driver,
                            dataset_label=dataset,
                            ancestry_group=ancestry_group,
                            subpop_label=subpop_label,
                            debug=True,
                        )

                except Exception as e:
                    print(f"  [WARN] Failed for {dataset}, {ancestry_group}, {subpop_label}: {e}")
                    counts = None

                if counts is None:
                    print(f"  No row found for subpop '{subpop_label}' in {dataset}/{ancestry_group}")
                    alternate_allele_count = total_count = homozygote_count = None
                else:
                    alternate_allele_count, total_count, homozygote_count = counts
                    print(
                        f"  {dataset}, {ancestry_group}, {subpop_label} → "
                        f"alt={alternate_allele_count}, total={total_count}, hom={homozygote_count}"
                    )

                output_rows.append(
                    {
                        "snp_id": snp_id,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "dataset": dataset,
                        "genetic_ancestry_group": ancestry_group,
                        "group_name": group_name,
                        "subpop_label": subpop_label,
                        "alternate_allele_count": alternate_allele_count,
                        "total_count": total_count,
                        "homozygote_count": homozygote_count,
                    }
                )

            # be gentle to gnomAD
            time.sleep(PER_VARIANT_DELAY)

    finally:
        driver.quit()

    out_df = pd.DataFrame(output_rows)
    out_df.to_csv(output_csv, index=False)
    print(f"\nWrote {len(out_df)} rows to {output_csv}")


if __name__ == "__main__":
    main()


