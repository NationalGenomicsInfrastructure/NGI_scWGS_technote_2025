


rule extract_sample_values:
    # https://penncnv.openbioinformatics.org/en/latest/user-guide/input/#penncnv-input-signal-intensity-files
    input:
        vcf = "samples.rename.vcf.gz",
    output:
        txt = "penncnv/{sample}/1_signal_intensity.txt"
    singularity: config["containers"]["bcftools"]
    shell:
        # Write header
        "echo 'Name\tChr\tPosition\t{wildcards.sample}.B Allele Freq\t{wildcards.sample}.Log R Ratio'"
        " > {output.txt}"
        " && "
        # Append values
        "bcftools query"
        " -s {wildcards.sample}"
        # Removing low GS variants did not improve CNV calling
        # " -i 'GS > 0.2'" # https://www.illumina.com/Documents/products/technotes/technote_gencall_data_analysis_software.pdf 
        " -f '%ID\t%CHROM\t%POS\t[%BAF\t%LRR]\n'"
        " {input.vcf}"
        " | "
        "awk '$4!=\".\"'" # Filter missing values
        " >> {output.txt}"


rule create_pfb_file:
    # https://penncnv.openbioinformatics.org/en/latest/user-guide/input/#pfb-population-frequency-of-b-allele-file
    # file from: https://emea.support.illumina.com/downloads/infinium-global-screening-array-v3-0-support-files.html
    input:
        txt = "infinium-global-screening-array-24-v3-0-a1-population-reports-maf-copy-numbers/GSA-24v3-0_A1_PopulationReport_MAF.txt"
    output:
        txt = "penncnv/name_to_pfb.txt"
    shell:
        "echo 'Name\tChr\tPosition\tPFB' > {output.txt}"
        " && "
        "cut -f 1-4 {input.txt} | tail -n+2 >> {output.txt}"


rule map_id_to_pfb_file:
    # map file from: https://emea.support.illumina.com/downloads/infinium-global-screening-array-v3-0-support-files.html
    # "Infinium Global Screening Array v3.0 Loci Name to rsID Conversion File"
    input:
        txt = "penncnv/name_to_pfb.txt",
        map = "GSA-24v3-0_A1_b151_rsids.txt"
    output:
        txt = "penncnv/rsid_to_pfb.txt"
    run:
        name_rsid = {}
        with open(input.map) as f:
            for line in f:
                name, rsid  = line.strip().split("\t")
                if rsid == ".":
                    continue
                
                name_rsid[name] = rsid 
        
        n = 0
        n_miss = 0
        with open(input.txt) as fin, open(output.txt, "w") as fout:
            next(fin) # skip header

            # Write header
            print("Name", "Chr", "Position", "PFB", sep="\t", file=fout)
            for line in fin:
                if line.startswith("Name"):
                    continue
                
                n += 1
                name, chrom, pos, pfb = line.strip().split("\t")

                if name not in name_rsid:
                    n_miss += 1
                    continue
                
                rsid = name_rsid[name]
                
                print(rsid, chrom, pos, pfb, sep="\t", file=fout)
        
        print(f"Map size: {len(name_rsid):,}")
        print(f"PFB entries: {n:,}")
        print(f"Names without RsID: {n_miss:,}")


rule get_snp_positions:
    input:
        vcf = "samples.rename.vcf.gz",
    output:
        bed = "penncnv/rsid.bed"
    singularity: config["containers"]["bcftools"]
    shell:
        "bcftools query"
        " -f '%CHROM\t%POS\t%POS\t%ID\n'"
        " {input.vcf}"
        " > {output.bed}"

rule genome_no_chr:
    input:
        genome = "GRCh38.genome"
    output:
        genome = "GRCh38.nochr.genome"
    shell:
        "sed 's|^chr||g' {input.genome} > {output.genome}"
    
rule expand_snp_positions:
    # https://penncnv.openbioinformatics.org/en/latest/user-guide/input/#hmm-file
    input:
        bed = "penncnv/rsid.bed",
        genome = "GRCh38.nochr.genome"
    output:
        bed = "penncnv/rsid_slop_1Mb.bed"
    singularity: config["containers"]["bedtools"]
    shell:
        "bedtools slop"
        " -g {input.genome}"
        " -i {input.bed}"
        " -b 500000" # 500 kb on each side
        " > {output.bed}"


import string
pieces = list(string.ascii_lowercase[:20])


rule split_bed:
    input:
        bed = "penncnv/rsid_slop_1Mb.bed"
    output:
        beds = temp(expand("penncnv/rsid_slop_1Mb.bed.tmp{suf}", suf=pieces))
    params:
        n = len(pieces),
        prefix = "penncnv/rsid_slop_1Mb.bed.tmp"
    shell:
        "split -l $((1 + `wc -l < {input.bed}`/{params.n})) -a 1 {input.bed} {params.prefix}"


rule calc_gc:
    input:
        bed = "penncnv/rsid_slop_1Mb.bed.tmp{suf}"
    output:
        txt = temp("penncnv/rsid_slop_1Mb.bed.tmp{suf}.gc")
    singularity: config["containers"]["bedtools"]
    params:
        fasta = config["fasta"]
    shell:
        "bedtools nuc"
        " -fi {params.fasta}"
        " -bed {input.bed}"
        " | "
        # Remove header
        "tail -n+2"
        " | "
        "awk '{{OFS=\"\\t\"; print $4,$1,$2,$6*100}}'"
        " >> {output.txt}"
    

rule merge_gcmodel_file:
    input:
        beds = expand("penncnv/rsid_slop_1Mb.bed.tmp{suf}.gc", suf=pieces)
    output:
        txt = "penncnv/gcmodel.txt"
    shell:
        "echo 'Name\tChr\tPosition\tGC' > {output.txt}"
        " && "
        "cat {input.beds}"
        " >> {output.txt}"


rule call_cnv:
    input:
        pfb = "penncnv/name_to_pfb.txt",
        gcmodel = "penncnv/gcmodel.txt",
        signal = "penncnv/{sample}/1_signal_intensity.txt"
    output:
        cnvs = "penncnv/{sample}/2_raw_cnvs.txt"
    log: "penncnv/{sample}/2_raw_cnvs.txt.log"
    singularity: config["containers"]["penncnv"]
    shell:
        "/home/user/PennCNV/detect_cnv.pl"
        " --test"
        " --confidence"
        # Same HMM file as used in https://doi.org/10.1016/j.ygeno.2024.110962
        " --hmmfile /home/user/PennCNV/lib/hhall.hmm"
        #" --hmmfile penncnv/hhall_uf_10pct.hmm" # 10% of markers expected as outliers (default is 1%)
        #" --hmmfile penncnv/hhall_uf_5pct.hmm" # 5% of markers expected as outliers (default is 1%)
        #" --hmmfile penncnv/hhall_uf_.1pct.hmm" # 0.1% of markers expected as outliers (default is 1%)
        " --pfbfile {input.pfb}"
        " --output {output.cnvs}"
        " --minsnp 10"
        " --minlength 100000"
        #" --sexfile <(echo '{input.signal}\tf'})"
        " --gcmodelfile {input.gcmodel}"
        " --coordinate_from_input" # Coordinates from pfbfile are from hg19
        " {input.signal}"
        " &> {log}"

rule call_cnv_joint:
    input:
        pfb = "penncnv/name_to_pfb.txt",
        gcmodel = "penncnv/gcmodel.txt",
        signals = expand("penncnv/{sample}/1_signal_intensity.txt", sample=[s["name"] for s in samples]),
    output:
        cnvs = "penncnv/cnvs.raw.txt"
    log: "penncnv/cnvs.raw.txt.log"
    singularity: config["containers"]["penncnv"]
    shell:
        "/home/user/PennCNV/detect_cnv.pl"
        " --test"
        " --confidence"
        # Same HMM file as used in https://doi.org/10.1016/j.ygeno.2024.110962
        #" --hmmfile /home/user/PennCNV/lib/hhall.hmm"
        #" --hmmfile penncnv/hhall_uf_10pct.hmm" # 10% of markers expected as outliers (default is 1%)
        " --hmmfile penncnv/hhall_uf_5pct.hmm" # 5% of markers expected as outliers (default is 1%)
        " --pfbfile {input.pfb}"
        " --output {output.cnvs}"
        " --minsnp 10"
        " --minlength 100000"
        #" --sexfile <(echo '{input.signal}\tf'})"
        " --gcmodelfile {input.gcmodel}"
        " --coordinate_from_input" # Coordinates from pfbfile are from hg19
        " {input.signals}"
        " &> {log}"


rule get_rsid_positions:
    input:
        vcf = "samples.rename.vcf.gz",
    output:
        txt = "penncnv/rsid.txt"
    singularity: config["containers"]["bcftools"]
    shell:
        "echo 'Name\tChr\tPosition' > {output.txt}"
        " && "
        "bcftools query"
        " -f '%ID\t%CHROM\t%POS\n'"
        " {input.vcf}"
        " >> {output.txt}"


rule clean_cnvs:
    input:
        signal = "penncnv/{sample}/1_signal_intensity.txt",
        cnvs = "penncnv/{sample}/2_raw_cnvs.txt"
    output:
        cnvs = "penncnv/{sample}/3_clean_cnvs.txt"
    log: "penncnv/{sample}/3_clean_cnvs.txt.log",
    singularity: config["containers"]["penncnv"]
    shell:
        "/home/user/PennCNV/clean_cnv.pl"
        " combineseg"
        " {input.cnvs}"
        " --signalfile {input.signal}"
        " --output {output.cnvs}"
        " --fraction 0.2"
        " 2> {log}"


rule convert_cnv_format:
    input:
        cnvs = "penncnv/{sample}/3_clean_cnvs.txt"
    output:
        cnvs = "penncnv/{sample}/4_clean_cnvs.tsv"
    singularity: config["containers"]["penncnv"]
    shell:
        "/home/user/PennCNV/convert_cnv.pl"
        " --intype penncnv"
        " --outtype tab"
        " {input.cnvs}"
        " > {output.cnvs}"


rule merge_tsvs:
    input:
        cnvs = expand("penncnv/{sample}/4_clean_cnvs.tsv", sample=[s["name"] for s in samples]),
    output:
        cnvs = "penncnv/cnvs.all.tsv"
    shell:
        "cat {input} > {output}"


rule convert_to_seg:
    input:
        cnvs = "{file}.tsv"
    output:
        seg = "{file}.seg"
    run:
        with open(input.cnvs) as fin, open(output.seg, "w") as fout:
            print("#type=COPY_NUMBER", file=fout)
            print("id", "chromosome", "start", "end", "cn", sep="\t", file=fout)

            for line in fin:
                chrom, start, end, cn, file, *_ = line.strip().split("\t")
                fileid = file.split("/")[1]
                print(fileid, "chr" + chrom, start, end, cn, sep="\t", file=fout)
