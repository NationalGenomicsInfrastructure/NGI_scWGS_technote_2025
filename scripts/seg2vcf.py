"""
Convert SEG file to VCF files
"""

import argparse
from pathlib import Path
from collections import defaultdict


def parse_fai(fai):
    with open(fai) as f:
        for line in f:
            chrom, length, *_ = line.strip().split("\t")
            yield chrom, int(length)


def parse_seg(seg):
    with open(seg) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            # The third column should be an integer
            if not line.strip().split("\t")[2].isdigit():
                continue

            els = line.strip().split("\t")
            name = els[0]
            chrom = els[1]
            start = int(els[2])
            end = int(els[3])
            cn = int(els[-1])

            yield name, chrom, start, end, cn


def write_header(chr_lengths, writer):
    print("##fileformat=VCFv4.2", file=writer)
    print("##source=\"seg2vcf.py\"", file=writer)

    for chrom, length in chr_lengths.items():
        print(f"##contig=<ID={chrom},length={length}>", file=writer)

    print("""\
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=TARGET_SIZE,Description="Call is smaller than the target size of 100kbp">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\
""", file=writer)
    
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE", file=writer)

def group_by_name(seg):
    data = defaultdict(list)
    for name, chrom, start, end, cn in parse_seg(seg):
        data[name].append((chrom, start, end, cn))
    
    return data

def main(args):
    print("Running seg2vcf.py")
    for k, v in vars(args).items():
        print(f" {k}: {v}")
    print("-"*30)
    
    print("Reading reference genome lengths...")
    chr_lengths = {chrom: length for chrom, length in parse_fai(args.fai)}
    if args.add_chr_ref:
        print("Adding 'chr' to chromosome names in reference genome...")
        chr_lengths = {f"chr{chrom}": length for chrom, length in chr_lengths.items()}
    
    print("Creating output directory...")
    outpath = Path(args.outdir)
    outpath.mkdir(exist_ok=True)

    for name, segs in group_by_name(args.seg).items():
        print(f"Processing {name}...")
        n_cnvs = 0
        with open(outpath / f"{name}.vcf", "w") as f:
            write_header(chr_lengths, f)
            for chrom, start, end, cn in segs:
                # Skip copy neutral segments
                if cn == args.ploidy:
                    continue

                if args.add_chr:
                    if chrom == "23":
                        chrom = "X"
                    elif chrom == "24":
                        chrom = "Y"
                    
                    chrom = f"chr{chrom}"
                
                assert chrom in chr_lengths, f"Chromosome '{chrom}' not found in reference genome."
                assert end <= chr_lengths[chrom], f"End position {end} is greater than chromosome length {chr_lengths[chrom]}."

                cnv_type = "DEL" if cn < args.ploidy else "DUP"
                svlen = end - start

                print(
                    chrom,
                    start,
                    ".",
                    "N",
                    f"<{cnv_type}>",
                    ".",
                    "PASS",
                    f"IMPRECISE;END={end};SVLEN={svlen};SVTYPE={cnv_type}",
                    "GT:CN",
                    f"0/1:{cn}",
                    sep="\t",
                    file=f
                )
                n_cnvs += 1
        
        print(f"Found {n_cnvs} CNVs in {name}.")
    
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seg", help="Segment file")
    parser.add_argument("-o", "--outdir", help="Output dir. Default: %(default)s (CWD).", default=Path(".").cwd())
    parser.add_argument("-f", "--fai", help="Reference genome FAI", required=True)
    parser.add_argument("-p", "--ploidy", help="Ploidy of the sample", type=int, default=2)
    parser.add_argument("--add-chr", help="Add 'chr' to chromosome names in SEG", action="store_true")
    parser.add_argument("--add-chr-ref", help="Add 'chr' to chromosome names in reference genome", action="store_true")
    args = parser.parse_args()
    main(args)