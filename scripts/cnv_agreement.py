"""
Caclulate CNV agreement between two SEG files.

Agreement is calculated as the fraction of the genome where the CNVs are the same.

OBS! Only for female samples (no chrY).
"""
import argparse
import logging
import sys

import pandas as pd
import pyranges as pr

from seg2vcf import parse_fai

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def _drop_invalid_intervals(df, label):
    """Drop malformed intervals where End <= Start.

    Negative or zero interval lengths can otherwise create negative weighted
    counts (FP/FN/TP) when length is used as a multiplier.
    """
    invalid = df["End"] <= df["Start"]
    n_invalid = int(invalid.sum())
    if n_invalid:
        logger.warning(
            "Dropping %d malformed interval(s) from %s where End <= Start.",
            n_invalid,
            label,
        )
        logger.debug("Malformed intervals:\n%s", df.loc[invalid])
    return df.loc[~invalid].copy()



def main(args):
    logger.info("Running cnv_agreement.py")
    for k,v in vars(args).items():
        logger.info(f" {k}: {v}")
    logger.info("-"*30)

    run(
        fai=args.fai,
        truth_seg=args.truth,
        query_seg=args.query,
        ploidy=args.ploidy,
        prefix=args.output_prefix
    )

def run(
        fai,
        truth_seg,
        query_seg,
        ploidy=2,
        prefix="cnv_agreement",
):
    logger.info("Reading reference genome lengths...")
    chr_lengths = {chrom: length for chrom, length in parse_fai(fai)}
    chr1toX = set(f"chr{i}" for i in list(range(1, 23)) + ["X"])

    # Background for filling out neutral regions not covered by CNVs
    background = pd.DataFrame([
        {"Chromosome": chrom, "Start": 0, "End": length, "CopyNumber": ploidy}
        for chrom, length in chr_lengths.items()
    ])
    background_pr = pr.PyRanges(background)

    # Load truth set 
    truth_seg_df = pd.read_csv(
        truth_seg, 
        sep="\t", 
        header=None, 
        skiprows=2, 
        names=["Name", "Chromosome", "Start", "End", "CopyNumber"],
        dtype={"Name": str, "Chromosome": str, "Start": int, "End": int, "CopyNumber": int}
    )

    # Check that truth file has only one sample
    if len(truth_seg_df["Name"].unique()) != 1:
        logger.error("Truth file should only have one sample.")
        sys.exit(1)

    # Create truth PyRanges
    truth_cnvs_df = truth_seg_df[["Chromosome", "Start", "End", "CopyNumber"]]
    truth_cnvs_df = truth_cnvs_df[truth_cnvs_df["Chromosome"].isin(chr1toX)]
    truth_cnvs_df = _drop_invalid_intervals(truth_cnvs_df, "truth")
    truth_cnvs_pr = pr.PyRanges(truth_cnvs_df)

    # Load query set    
    query_seg_df = pd.read_csv(
        query_seg, 
        sep="\t", 
        skiprows=2, 
        header=None, 
        names=["Name", "Chromosome", "Start", "End", "CopyNumber"],
        dtype={"Name": str, "Chromosome": str, "Start": int, "End": int, "CopyNumber": int}
    )

    seg = []
    output_stats = prefix + ".stats.tsv"
    with open(output_stats, "w") as f:
        print(
            "Name",
            "TotalBases",
            "TPs",
            "FPs",
            "FNs",
            "Agreement",
            "Precision",
            "Recall",
            "F1Score",
            sep="\t",
            file=f
        )       

        # Loop over queries and calculate agreement    
        for name, df in query_seg_df.groupby("Name"):
            # Create query PyRanges
            query_cnvs_df = df[["Chromosome", "Start", "End", "CopyNumber"]]
            query_cnvs_df = _drop_invalid_intervals(query_cnvs_df, f"query sample {name}")
            query_cnvs_pr = pr.PyRanges(query_cnvs_df)
            
            # Add missing neutral regions if any
            missing = background_pr.subtract(query_cnvs_pr)
            query_cnvs_df = pd.concat([query_cnvs_pr.df,missing.df])
            query_cnvs_pr = pr.PyRanges(query_cnvs_df)

            # Intersect truth and query to get query CNVs
            query_cnvs_pr = query_cnvs_pr.intersect(truth_cnvs_pr)

            # Intersect truth and query to get truth CNVs
            truth_cnvs_pr_inter = truth_cnvs_pr.intersect(query_cnvs_pr)

            # Merge truth and query
            query_cnvs_pr = query_cnvs_pr.df.merge(truth_cnvs_pr_inter.df, how="left", on=["Chromosome", "Start", "End"], suffixes=("", "Truth"))

            # Evaluate agreement / True positives
            query_cnvs_pr["Agree"] = (query_cnvs_pr["CopyNumber"] == query_cnvs_pr["CopyNumberTruth"]).astype(int)
            disagree = query_cnvs_pr["Agree"] == 0

            # For the disagreements, calculate false positives and false negatives
            # FP = CNV called in query
            # FN = CNV called in truth
            query_cnvs_pr["FP"] = ((query_cnvs_pr["CopyNumberTruth"] == ploidy) & disagree).astype(int)
            query_cnvs_pr["FN"] = ((query_cnvs_pr["CopyNumberTruth"] != ploidy) & disagree).astype(int)
            
            query_cnvs_pr = pr.PyRanges(query_cnvs_pr)

            # Calculate weight and length
            query_cnvs_pr.Length = query_cnvs_pr.End - query_cnvs_pr.Start
            assert (query_cnvs_pr.Length >= 0).all(), "Negative interval lengths found. Check input SEG files."
            query_cnvs_pr.TPs = query_cnvs_pr.Length * query_cnvs_pr.Agree
            assert (query_cnvs_pr.TPs >= 0).all(), "Negative TP counts found. Check input SEG files."
            query_cnvs_pr.FPs = query_cnvs_pr.Length * query_cnvs_pr.FP
            assert (query_cnvs_pr.FPs >= 0).all(), "Negative FP counts found. Check input SEG files."
            query_cnvs_pr.FNs = query_cnvs_pr.Length * query_cnvs_pr.FN
            assert (query_cnvs_pr.FNs >= 0).all(), "Negative FN counts found. Check input SEG files."

            # Calculate agreement
            total_bases = query_cnvs_pr.Length.sum()
            total_TP = query_cnvs_pr.TPs.sum()
            total_FP = query_cnvs_pr.FPs.sum()
            total_FN = query_cnvs_pr.FNs.sum()
            agreement = total_TP / total_bases if total_bases > 0 else 0.0  # this is the jaccard similarity/index
            precision_denom = total_TP + total_FP
            recall_denom = total_TP + total_FN
            precision = total_TP / precision_denom if precision_denom > 0 else 0.0
            recall = total_TP / recall_denom if recall_denom > 0 else 0.0
            f1_denom = precision + recall
            f1_score = 2 * (precision * recall) / f1_denom if f1_denom > 0 else 0.0

            # Write stats
            print(
                name, 
                total_bases, 
                total_TP, 
                total_FP, 
                total_FN, 
                agreement,
                precision,
                recall,
                f1_score,
                sep="\t", 
                file=f
            )

            # Save data for visualization of agreement
            s = query_cnvs_pr[["Chromosome", "Start", "End", "Agree"]].df
            s["Agree"] = s["Agree"].replace({0: 4, 1: 0}) # --> Agree=0 (blue), NotAgree=4 (red)
            s["Id"] = name
            seg.append(s)

    # Write SEG for visualization of agreement
    seg = pd.concat(seg)
    seg = seg[["Id", "Chromosome", "Start", "End", "Agree"]]
    output_seg = prefix + ".agreements.seg"
    with open(output_seg, "w") as f:
        print("#type=COPY_NUMBER", file=f)
        print("id\tchromosome\tstart\tend\tcn", file=f)

    seg.to_csv(output_seg, mode="a", sep="\t", index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query", help="Query SEG file")
    parser.add_argument("truth", help="Truth SEG file. Should only have one sample.")
    parser.add_argument("-f", "--fai", help="Reference fasta index file", required=True)
    parser.add_argument("-p", "--ploidy", help="Ploidy. Default: %(default)s", type=int, default=2)
    parser.add_argument("-o", "--output-prefix", help="Output prefix. Default: %(default)s", default="cnv_agreement")
    args = parser.parse_args()
    main(args)