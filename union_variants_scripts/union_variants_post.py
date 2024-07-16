import argparse
import pandas as pd
import json
import os
import pickle
from collections import defaultdict
from union_variants_pre import Variant

def read_bamcounts(bamcounts_file, variants_dict, sample, report_file):
    bamcounts = dict()
    with open(bamcounts_file, 'r') as infile, open(report_file, 'a') as outfile:
        for line in infile:
            fields = line.split('\t')
            chrom, pos = fields[0], fields[1]

            try: 
                variants = variants_dict[f"{chrom}_{pos}"] # insertion or snv
            except KeyError:
                variants = variants_dict[f"{chrom}_{int(pos)-1}"] # deletion
            
            all_counts = fields[4:]
            base_counts = {counts.split(":")[0]: counts.split(":")[1] for counts in all_counts}
            bam_depth = fields[3]

            for variant in variants:
                ref = variant.ref
                alt = variant.alt

                bam_alt_count = ""

                if len(ref) > len(alt): #deletion
                    assert(len(alt) == 1)
                    alt = "-"+ref[1:]
                
                elif len(ref) < len(alt): #insertion
                    assert(len(alt) > len(ref))
                    assert(len(ref) == 1)
                    alt = "+"+alt[1:]
                
                try:
                    bam_alt_count = base_counts[alt]
                except KeyError:
                    bam_alt_count = "0" # how should i handle this case?
                

                bamcounts[variant.id] =f"{bam_depth}:{bam_alt_count}"

                # TODO: option to report discrepancy between vcf and bamcounts 
                if sample in variant.samples:
                    vcf_alt_count, vcf_depth = variant.get_sample_counts(sample)
                    #print(vcf_alt_count, vcf_depth)
                    if not (vcf_alt_count == bam_alt_count and vcf_depth == bam_depth): 
                        outfile.write(f"{variant.id}\t{sample}\t{vcf_depth}:{vcf_alt_count}\t{bam_depth.strip()}:{bam_alt_count.strip()}\n")
    
    return bamcounts


def write_vcf_file(file, variant_to_counts):
    def sort_key(key):
        prefix, suffix = key.split('_', 1)
        if prefix.isdigit():
            return (0, int(prefix), int(suffix.split('_')[0]))  # Numeric chromosomes first
        elif prefix == 'X':
            return (1, float('inf'), int(suffix.split('_')[0]))  # 'X' after numeric
        elif prefix == 'Y':
            return (2, float('inf'), int(suffix.split('_')[0]))  # 'Y' after 'X'
        return (3, float('inf'), float('inf'))  # Any other non-numeric values
    
    sorted_variants = sorted(list(variant_to_counts.keys()), key=sort_key)

    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']
    
    # reporting format as DR:RD:AP <==> total depth : reference depth : alternate depth
    
    with open(file, 'wt') as f:
        f.write("\t".join(columns) + "\n")
        for variant in sorted_variants:
            sample_counts = variant_to_counts[variant]
            chrom, pos, ref, alt = variant.split("_")
            f.write(f"{chrom}\t{pos}\t{variant}\t{ref}\t{alt}\t.\tPASS\tGENE_PLACEHOLDER\tDP:AP\t0:0\t{sample_counts}\n")
    f.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-patient_id", required=True)
    parser.add_argument("-file_info", required=True)
    parser.add_argument("--ignore_samples", nargs="+")
    parser.add_argument("--single_file", action="store_true")

    args = parser.parse_args()

    file_info = pd.read_csv(args.file_info, sep='\t')
    
    patient_info = file_info.loc[file_info['Patient'] == args.patient_id]

    out_dir = patient_info['Directory'].iloc[0]
    Samples = {}

    ignore_samples = args.ignore_samples if args.ignore_samples is not None else []
    Samples = {}

    for index, row in patient_info.iterrows():
        sample = row['Sample']
        if sample not in ignore_samples:
            Samples[sample] = {'bam': row['BAM'],'vcf': row['VCF'].split(', '), 'bamcounts': row['BAM_COUNTS']}
    
    variants_pkl = os.path.join(out_dir, patient_info['VARIANTS'][0])

    all_variants = defaultdict(list)
    with open(variants_pkl, 'rb') as f:
        loaded_variants = pickle.load(f)

    for variant in loaded_variants:
        all_variants[f"{variant.chrom}_{variant.pos}"].append(variant)
    
    # TODO: add case for writing to single vcf file with all sample columns
    # TODO: handle 'NORMAL' column data (currently just setting DP:AP = 0:0)

    report_file = os.path.join(out_dir, "vcf_bamcounts_report.txt")

    with open(report_file, 'w') as f:
        f.write("variant_id\tsample_origin\tvcf_counts\tbam_counts\n")
    
    for sample in Samples.keys():
        print(f"{sample}: Reading bam-readcount output and writing counts to vcf file...", end="")
        bamcounts_file = os.path.join(out_dir, sample, Samples[sample]['bamcounts'])
        variant_to_counts = read_bamcounts(bamcounts_file, all_variants, sample, report_file)
        vcf_outfile = os.path.join(out_dir, sample, f"{sample}_union.vcf")
        write_vcf_file(vcf_outfile, variant_to_counts)
        print("Done.")



if __name__ == "__main__":
    main()