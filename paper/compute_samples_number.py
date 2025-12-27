""" """

from typing import Set, Dict
from time import time
from tqdm import tqdm
from glob import glob

import pandas as pd

import pysam
import os

RESULTSDIR = "results"
PAMS = {"CAS9": "NGG", "CPF1": "TTTV"}
SAMPLESNUM = {"1000G": 2504, "HGDP": 929, "GNOMAD": 76215}

# download from gnomad data webserver (tbi required), and replace values with new paths
GNOMADVCFS_ORIGINAL = {
    "chr2": "/storage_2/mtognon/datasets/variants/GNOMAD/v4.1/genomes/bcfs/gnomad.genomes.v4.1.sites.chr2.vcf.bgz",
    "chr3": "/storage_2/mtognon/datasets/variants/GNOMAD/v4.1/genomes/bcfs/gnomad.genomes.v4.1.sites.chr3.vcf.bgz",
    "chr7": "/storage_2/mtognon/datasets/variants/GNOMAD/v4.1/genomes/bcfs/gnomad.genomes.v4.1.sites.chr7.vcf.bgz",
    "chr11": "/storage_2/mtognon/datasets/variants/GNOMAD/v4.1/genomes/bcfs/gnomad.genomes.v4.1.sites.chr11.vcf.bgz",
}

def _compute_id(chrom: str, start: int, stop: int, strand: str) -> str:
    return f"{chrom}_{start}_{stop}_{strand}"

def compute_samples_num(df: pd.DataFrame, samplesnum: int) -> pd.DataFrame:
    assert "samples" in df.columns.tolist()
    # count samples carrying alternative grnas
    df["n_samples"] = df["samples"].apply(
        lambda x: len(str(x).split(",")) if pd.notna(x) and x!= "REF" else None
    )  
    # compute guide id for each grna
    df["guide_id"] = df.apply(
        lambda x: _compute_id(x["chr"], x["start"], x["stop"], x["strand"]), axis=1
    )
    # sum samples carrying alternative guides
    sum_alternative = df[df["samples"] != "REF"].groupby("guide_id")["n_samples"].sum()
    # retrieve number of samples for reference grnas
    refmask = df["samples"] == "REF"
    df.loc[refmask, "n_samples"] = samplesnum - df.loc[refmask, "guide_id"].map(sum_alternative).fillna(0)
    return df

# gnomad-only post-processing (NB: may require hours to run)
def _retrieve_position(variant_id: str) -> int:
    return int(variant_id.split("-")[1])

def extract_variant_positions(df: pd.DataFrame) -> Set[int]:
    # retrieve position of each variant
    variant_ids = df["variant_id"].dropna().tolist()
    return {_retrieve_position(v) for vs in variant_ids for v in vs.split(",")}

def get_relevant_variants(vcf_fname: str, chrom: str, positions: Set[int]) -> pd.DataFrame:
    print(f"Scanning VCF for {len(positions)} positions on {chrom}")
    vcf = pysam.VariantFile(vcf_fname)  # load vcf
    variants, matches = [], set()
    # requires bgzipped vcf with tbi index
    for r in tqdm(vcf.fetch(chrom), desc="Scanning VCF", unit="variants"):
        if r.pos not in positions:
            continue
        ac = r.info.get("AC")  # read allele count
        if ac is None:
            continue
        nhomalt = r.info.get("nhomalt", 0)  # count number of homozygous for alternative
        # iterate on each alternative allele
        for i, alt in enumerate(r.alts):
            # scalar vs per-allele fields for AC and nhomalt
            if isinstance(ac, (list, tuple)):
                ac_ = ac[i] if i < len(ac) else None
            else:
                ac_ = ac
            if ac_ is None:
                continue
            if isinstance(nhomalt, (list, tuple)):
                nhomalt_ = nhomalt[i] if i < len(nhomalt) else 0
            else:
                nhomalt_ = nhomalt
            # compute number of heterozygous samples
            n_het = ac_ - 2 * nhomalt_
            n_total = nhomalt_ + n_het
            # add variants and matched position
            variant_id = f"{r.chrom}-{r.pos}-{r.ref}/{alt}"
            variants.append(
                {
                    "variant_id": variant_id,
                    "chrom": r.chrom,
                    "pos": r.pos,
                    "ref": r.ref,
                    "alt": alt,
                    "AC": ac_,
                    "n_hom": nhomalt_,
                    "n_het": n_het,
                    "n_samples": n_total
                }
            )
            matches.add(r.pos)
    print(f"Matched {len(matches)} of {len(positions)} positions")
    return pd.DataFrame(variants)

def build_samples_dict(df_variants: pd.DataFrame) -> Dict:
    # map variant IDs to total sample carrier counts
    df_sums = df_variants.groupby("variant_id", sort=False)["n_samples"].sum().reset_index()
    return df_sums.set_index("variant_id")["n_samples"].to_dict()

def get_samples_total(variant_ids, samples_dict):
    if pd.isna(variant_ids):
        return "NA"
    total = sum(
        int(samples_dict[v]) 
        for v in variant_ids.split(",") 
        if v in samples_dict
    )
    return str(total) if total else "NA"

def annotate_gnomad_report(vcf_path: str, chrom: str, report_df: pd.DataFrame) -> pd.DataFrame:
    print("Extracting variant positions...")
    positions = extract_variant_positions(report_df)
    if not positions:
        print("No valid variant positions found. Exiting.")
        return report_df
    start = time()
    df_variants = get_relevant_variants(vcf_path, chrom, positions)
    # Write carriers file as COMMA-separated 
    variant_out = f"variants_carriers_{chrom}_pysam.tsv"
    df_variants.to_csv(variant_out, index=False)
    print(f"Saved variant carrier data to {variant_out} in {time() - start:.2f} seconds")
    # Build samples dict from grouped sums
    samples_dict = build_samples_dict(df_variants)
    # Annotate the report
    tqdm.pandas(desc="Annotating report")
    report_df["n_samples"] = report_df["variant_id"].progress_apply(lambda x: get_samples_total(x, samples_dict))
    print(f"\nDone")
    return report_df


def main():
    outdir = f"{RESULTSDIR}_nsamples"
    for cas in PAMS:
        results_cas_dir = os.path.join(RESULTSDIR, cas)
        outdir_cas = os.path.join(outdir, cas)
        for dataset in SAMPLESNUM:
            results_dataset_dir = os.path.join(results_cas_dir, dataset)
            outdir_dataset = os.path.join(outdir_cas, dataset)
            os.makedirs(outdir_dataset, exist_ok=True)
            reports = glob(os.path.join(results_dataset_dir, "crisprhawk_guides__*"))
            for report in reports:
                df = pd.read_csv(report, sep="\t")
                if dataset == "GNOMAD":  # this may take a long time
                    chrom = os.path.basename(report).replace("crisprhawk_guides__", "").split("_")[0]
                    df = annotate_gnomad_report(GNOMADVCFS_ORIGINAL[chrom], chrom, df)
                else:
                    df = compute_samples_num(df, SAMPLESNUM[dataset])
                outfname = os.path.join(outdir_dataset, os.path.basename(report).replace(".tsv", "_nsamples.tsv"))
                df.to_csv(outfname, sep="\t", index=False)

if __name__ == "__main__":
    main()