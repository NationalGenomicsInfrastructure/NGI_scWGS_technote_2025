"""
Caclulate CNV agreement between two SEG files.

Agreement is calculated as the fraction of the genome where the CNVs are the same.

OBS! Only for female samples (no chrY).
"""
import argparse
from collections import defaultdict
import dataclasses
import logging
import sys

import numpy as np

from seg2vcf import parse_fai, parse_seg

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclasses.dataclass
class CNV:
    chrom: str
    start: int
    end: int
    cn: int
    
    def gain_loss(self, ploidy=2):
        """
        Return 0 for copy neutral, 1 for gain, -1 for loss.
        """
        if self.cn == ploidy:
            return 0
        elif self.cn > ploidy:
            return 1
        return -1


def get_cnvs(seg, ploidy=2):
    """
    Get CNVs from a SEG file.
    """
    chr1toX = set(f"chr{i}" for i in list(range(1, 23)) + ["X"])

    data = defaultdict(list)
    for name, chrom, start, end, cn in parse_seg(seg):
        if cn == ploidy:
            continue

        if chrom not in chr1toX:
            continue

        data[name].append(CNV(chrom, start, end, cn))

    return data


def normalize_cnvs(cnvs, chr_lengths, ploidy=2):
    """
    Normalize CNVs by filling in gaps between CNVs with copy neutral segments.
    """
    normalized = []
    current_pos = 0
    current_chrom = None
    for cnv in cnvs:
        # New chromosome
        if cnv.chrom != current_chrom:
            # Add the copy neutral CNV if gap at end of the previous chromosome
            if current_chrom is not None and current_pos < chr_lengths[current_chrom]:
                normalized.append(
                    CNV(current_chrom, current_pos, chr_lengths[current_chrom], ploidy)
                )
            
            # Update current position
            current_chrom = cnv.chrom
            current_pos = 0
        
        # Add copy neutral segment if there is a gap between CNVs (or at start)
        if cnv.start > current_pos:
            normalized.append(CNV(cnv.chrom, current_pos, cnv.start, ploidy))

        if cnv.end > chr_lengths[cnv.chrom]:
            logger.warning(f"End position {cnv.end} is greater than chromosome length {chr_lengths[cnv.chrom]}.")
            logger.warning(f"Truncating CNV to chromosome length.")
            cnv.end = chr_lengths[cnv.chrom]

        current_pos = cnv.end 
        normalized.append(cnv)
 
    # Add the copy neutral CNV if gap at end of the last chromosome
    if current_pos < chr_lengths[current_chrom]:
        normalized.append(
            CNV(current_chrom, current_pos, chr_lengths[current_chrom], ploidy)
        )

    return normalized


def group_by_chrom(cnvs):
    """
    Group CNVs by chromosome.
    """
    data = defaultdict(list)
    for cnv in cnvs:
        data[cnv.chrom].append(cnv)
    
    return data


def assign_cn(intervals, cnvs, ploidy=2, gain_loss=False):
    """
    Assign CNV to intervals.
    """
    if gain_loss:
        cn_array = np.zeros(len(intervals))
    else:
        cn_array = np.ones(len(intervals)) * ploidy
    
    for cnv in cnvs:
        # Find intervals that are fully contained in the CNV
        idx = ((intervals[:,0] >= cnv.start) & (intervals[:,1] <= cnv.end)).nonzero()[0]
        
        # Update CN for these intervals
        cn = cnv.cn if not gain_loss else cnv.gain_loss(ploidy)
        cn_array[idx] = cn

    return cn_array


def count_types(cnvs, ploidy=2, gain_loss=False):
    """
    Count the number of CNVs of each type.
    """
    types = defaultdict(int)
    for cnv in cnvs:
        cntype = cnv.cn if not gain_loss else cnv.gain_loss(ploidy)
        types[cntype] += 1
    return dict(types)


def process_truth(truth_seg, chr_lengths, ploidy=2, gain_loss=False):
    """
    Process truth SEG file to get CNVs grouped by chromosome.
    """
    truth_cnvs = get_cnvs(truth_seg, ploidy)
    
    if len(truth_cnvs) != 1:
        logger.error("Truth file should only have one sample.")
        sys.exit(1)
    
    truth_cnvs = list(truth_cnvs.values())[0]

    truth_types = count_types(truth_cnvs, ploidy, gain_loss)
    logger.info(f"Truth CNVs: {len(truth_cnvs)} ({truth_types})")

    truth_cnvs = normalize_cnvs(truth_cnvs, chr_lengths, ploidy)
    truth_cnvs = group_by_chrom(truth_cnvs)
    return truth_cnvs


def main(args):
    logger.info("Running cnv_agreement.py")
    for k,v in vars(args).items():
        logger.info(f" {k}: {v}")
    logger.info("-"*30)

    logger.info("Reading reference genome lengths...")
    chr_lengths = {chrom: length for chrom, length in parse_fai(args.fai)}
    
    truth_cnvs = process_truth(args.truth, chr_lengths, args.ploidy, args.gl)

    all_querys = get_cnvs(args.query, args.ploidy)

    for name, query_cnvs in all_querys.items():
        logger.info(f"Processing {name}...")
        query_types = count_types(query_cnvs, args.ploidy, args.gl)
        logger.info(f"Query CNVs: {len(query_cnvs)} ({query_types})")
        
        query_cnvs = normalize_cnvs(query_cnvs, chr_lengths, args.ploidy)
        query_cnvs = group_by_chrom(query_cnvs)

        chroms_union = set(truth_cnvs.keys()) | set(query_cnvs.keys())
        logger.info(f"Chromosomes: {','.join(sorted(chroms_union))}")
        agree_length = 0
        total_length = 0
        for chrom in chroms_union:
            truth = truth_cnvs.get(chrom, [])
            query = query_cnvs.get(chrom, [])

            # Get all positions
            truth_positions = set(cnv.start for cnv in truth) | set(cnv.end for cnv in truth)
            query_positions = set(cnv.start for cnv in query) | set(cnv.end for cnv in query)
            positions = sorted(truth_positions | query_positions)

            # Get intervals
            intervals = np.array([(positions[i], positions[i+1]) for i in range(len(positions)-1)])

            # Calculate CN for each interval 
            truth_cn = assign_cn(intervals, truth, args.ploidy, args.gl)
            query_cn = assign_cn(intervals, query, args.ploidy, args.gl)

            # Calculate agreement
            agree = intervals[np.where(truth_cn == query_cn)]
            agree_length += np.sum(agree[:,1] - agree[:,0])

            length = np.sum(intervals[:,1] - intervals[:,0])
            assert length == chr_lengths[chrom], f"Length mismatch for chromosome {chrom}: {length} != {chr_lengths[chrom]}"

            total_length += length

        agree_fraction = agree_length / total_length

        logger.info(f"Agreement: {agree_length:,} / {total_length:,} ({agree_fraction:.2%})")

        logger.info("-"*30)

        print(f"{name}\t{agree_fraction}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query", help="Query SEG file")
    parser.add_argument("truth", help="Truth SEG file. Should only have one sample.")
    parser.add_argument("-f", "--fai", help="Reference fasta index file", required=True)
    parser.add_argument("-p", "--ploidy", help="Ploidy. Default: %(default)s", type=int, default=2)
    parser.add_argument("--gl", "--gain-loss", help="Only consider gains and losses", action="store_true")
    args = parser.parse_args()
    main(args)