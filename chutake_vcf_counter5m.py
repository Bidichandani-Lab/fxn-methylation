#!/usr/bin/env python3
import gzip
import subprocess
from statistics import mean
from argparse import ArgumentParser, Namespace

"""
This script was translated to python from the previous manual pipeline. 
It generates potentially useful statistics.
"""


def set_from_lines(filepath: str):
    with open(filepath, "r") as file:
        return set(line.strip() for line in file)


def main(args: Namespace):
    # prep a set to confirm bisulfite conversion
    confirm_hash = set_from_lines(args.confirm_file)
    # prep a set to check certain bases
    check_hash = set_from_lines(args.check_file)

    filtered_lines = []
    with gzip.open(args.vcf_file, "rt") as vcf:
        for line in vcf:
            stripped_line = str(line.strip())
            if stripped_line and not stripped_line.startswith("#"):
                filtered_lines.append(stripped_line)

    conversion_rates = []
    targets = []
    for line in filtered_lines:
        vcf_line = line.split("\t")
        info1 = vcf_line[-1].split(":")
        info2 = info1[-1].split(",")

        if len(info2) > 1:
            output = f"{vcf_line[1]}\t{info2[0]}\t{info2[1]}"
            rate = int(info2[0]) / (int(info2[0]) + int(info2[1]))
        else:
            converted = int(info1[1]) - int(info2[0])
            output = f"{vcf_line[1]}\t{info1[2]}\t{converted}"
            rate = 100 - (converted / int(info2[0]))

        if vcf_line[1] in confirm_hash:
            conversion_rates.append(rate)
        if vcf_line[1] in check_hash:
            targets.append(output)

    bs_conversion_rate = mean(conversion_rates) if conversion_rates else 0

    bam_read_count = subprocess.getoutput(f"samtools view -c {args.bam_file}").strip()

    outfile = f"{args.vcf_file.replace('.vcf.gz', '')}.counts-py.tsv"
    with open(outfile, "w") as out:
        out.write("Total_Reads\tConversion_Rate\n")
        out.write(f"{bam_read_count}\t{float(bs_conversion_rate):.2}\n")
        out.write("Base\tNum_UnMethylated_(G->A)\tNum_Methylated_(G)\n")
        out.writelines(targets)


if __name__ == "__main__":
    parer = ArgumentParser(
        description="Count the number of methylated and unmethylated reads in a VCF file.")
    parer.add_argument("vcf_file",
                       help="The VCF file from which to count methylated and unmethylated reads.")
    parer.add_argument("confirm_file",
                       help="The file containing the target bases to confirm BS conversion.")
    parer.add_argument("check_file",
                       help="The file containing the target bases to check methylation status.")
    parer.add_argument("bam_file",
                       help="The BAM file from which the VCF was generated.")
    args_ = parer.parse_args()
    main(args_)
