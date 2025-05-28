import argparse
import pandas as pd
import json
import os
import pickle
from collections import defaultdict
from union_variants_pre import Variant


def impute_readcounts(bamcounts_dict, variants_dict):

    gene_regions = {"1": [(148602086,148688393)],                                                #ACVR2A
                    "7": [(151832010, 152133090)],                                               #KMT2C
                    "6": [(80713604, 80752244)],                                                 #TTK
                    "3": [(50712672, 51421629), (53880607, 53899615), (114070533, 114863756)],   #DOCK3, IL17RB, and ZBTB20
                    "4": [(83769979, 83821722)],                                                 #SEC31A
                    "11": [(31811484, 31833725)],                                                #PAX6
                    "X": [(47420516, 47431307), (105937024, 106040223)],                         #ARAF and RNF128
                    "20": [(30780306, 30795594)],                                                #PLAGL2
                    } 

    for variant_id, variant in variants_dict.items():
        if variant_id in bamcounts_dict.keys():
            continue
        if variant.chrom not in gene_regions.keys():
            continue
        variant_in_gene = False
        for region in gene_regions[variant.chrom]:
            if int(variant.pos) >= region[0] and int(variant.pos) <= region[1]:
                variant_in_gene = True
                break
        if variant_in_gene:
            print(f"Imputing bam readcounts for variant {variant.id}")
            bamcounts_dict[variant_id] = "0:0"



# TODO: handle strelka and mutect vcf counts
def read_bamcounts(bamcounts_file, variants_dict, sample, report_file=None):
    bamcounts = dict()
    with open(bamcounts_file, 'r') as infile:
        outfile = open(report_file, 'w') if report_file is not None else None

        for line in infile.readlines():
            fields = line.split('\t')
            
            if len(fields) < 4: 
                print(f"Warning: unable to process bamcounts line: \n {line} \n continuing...")
                continue
            
            chrom, pos = fields[0], fields[1]

            try:
                int(pos)
            except (ValueError, TypeError):
                print(f"Warning: unable to parse bamcounts line:\n {line}.\n Continuing.")
                continue
    
            variants = variants_dict[f"{chrom}_{pos}"]  # insertion or SNV
            
            if not variants:
                variants = variants_dict[f"{chrom}_{int(pos)-1}"]  # deletion
            
            if not variants: continue

            all_counts = fields[4:]
            base_counts = dict()

            try:
                base_counts = {counts.split(":")[0]: counts.split(":")[1] for counts in all_counts}
            except IndexError:
                print(f"Warning: unable to parse bamcounts line:\n {line}.\n Continuing.")
                continue

            bam_depth = fields[3].strip()

            for variant in variants:
                ref = variant.ref
                alt = variant.alt

                bam_alt_count, bam_ref_count = "", ""

                if len(ref) > len(alt): #deletion
                    assert(len(alt) == 1)
                    alt = "-"+ref[1:]
                
                elif len(ref) < len(alt): #insertion
                    assert(len(alt) > len(ref))
                    assert(len(ref) == 1)
                    alt = "+"+alt[1:]
                
                # DP:AP <-> total:ref

                try:
                    bam_alt_count = base_counts[alt].strip()
                    # sometimes, total != ref + alt, but want alt count to be accurate --> ref = total - alt
                    bamcounts[variant.id] =f"{bam_depth}:{int(bam_depth) - int(bam_alt_count)}" 

                except KeyError:
                    try:
                        # sometimes (esp for indels), alt allele not in bamcounts --> use ref
                        bam_ref_count = base_counts[ref].strip()
                        bamcounts[variant.id] =f"{bam_depth}:{int(bam_ref_count)}" 
                    except KeyError:
                        # shouldn't really get to this case
                        print(f"Warning: ref {ref} and alt {alt} alleles not present in bam readcounts. Using depth count")
                        bamcounts[variant.id] =f"{bam_depth}:{int(bam_depth)}" 


                if outfile is not None:
                    # TODO: command line option to report discrepancy between vcf and bamcounts 
                    strelka_vcf_alt_count, strelka_vcf_depth = "", ""
                    mutect_vcf_alt_count, mutect_vcf_depth = "", ""

                    if sample in variant.samples['strelka']:
                        strelka_vcf_alt_count, strelka_vcf_depth = variant.get_sample_counts('strelka', sample)
                
                    if sample in variant.samples['mutect']:
                        mutect_vcf_alt_count, mutect_vcf_depth = variant.get_sample_counts('mutect', sample)
                        
                    if (not (strelka_vcf_alt_count == bam_alt_count and strelka_vcf_depth == bam_depth)) or (not (mutect_vcf_alt_count == bam_alt_count and mutect_vcf_depth == bam_depth)): 
                            if (mutect_vcf_alt_count, mutect_vcf_depth) == ("", ""):
                                if ((strelka_vcf_alt_count == bam_alt_count and strelka_vcf_depth == bam_depth)):
                                    continue
                            if (strelka_vcf_alt_count, strelka_vcf_depth) == ("", ""):
                                if ((mutect_vcf_alt_count == bam_alt_count and mutect_vcf_depth == bam_depth)):
                                    continue
                            outfile.write(f"{variant.id}\t{sample}\t{strelka_vcf_depth.strip()}:{strelka_vcf_alt_count.strip()}\t{mutect_vcf_depth.strip()}:{mutect_vcf_alt_count.strip()}\t{bam_depth.strip()}:{bam_alt_count.strip()}\n")

    return bamcounts


def write_vcf_file(sample_to_variants, final_variants, out_dir, single_file=False):
    print("Writing output vcf file(s)...", end="")
    
    def sort_key(key):
        prefix, suffix = key.split('_', 1)
        if prefix.isdigit():
            return (0, int(prefix), int(suffix.split('_')[0]))  # Numeric chromosomes first
        elif prefix == 'X':
            return (1, float('inf'), int(suffix.split('_')[0]))  # 'X' after numeric
        elif prefix == 'Y':
            return (2, float('inf'), int(suffix.split('_')[0]))  # 'Y' after 'X'
        return (3, float('inf'), float('inf'))  # Any other non-numeric values
    
    sorted_variants = sorted(list(final_variants), key=sort_key)

    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL']
    
    # reporting format as DP:AP <==> total depth : ref depth
    
    if not single_file:
        for sample, variant_counts in sample_to_variants.items():
            print(f"{sample}...", end="")
            vcf_file = os.path.join(out_dir, f"{sample}.vcf")
            with open(vcf_file, 'wt') as f:
                f.write("\t".join(columns) + "\tTUMOR\n")
                for variant in sorted_variants:
                    try:    
                        sample_counts = variant_counts[variant]
                    except KeyError:
                        sample_counts = "0:0"
                    chrom, pos, ref, alt = variant.split("_")
                    f.write(f"{chrom}\t{pos}\t{variant}\t{ref}\t{alt}\t.\tPASS\tGENE_PLACEHOLDER\tDP:AP\t0:0\t{sample_counts}\n")
            f.close()
    
    else:
        print("all samples in single file...", end="")
        sample_names =[s for s in sample_to_variants.keys()]
        sample_ids = [s.replace("Sample", "") for s in sample_names]
        columns += sample_ids
        vcf_file = os.path.join(out_dir, "_".join(sample_ids)+".vcf")
        with open(vcf_file, 'wt') as f:
            f.write("\t".join(columns))
            for variant in sorted_variants:
                chrom, pos, ref, alt = variant.split("_")
                counts = []
                for s in sample_names:
                    try:    
                        sample_counts = variant_counts[variant]
                    except KeyError:
                        sample_counts = "0:0"
                    counts.append(sample_counts)
                all_counts = "\t".join(counts)
                f.write(f"{chrom}\t{pos}\t{variant}\t{ref}\t{alt}\t.\tPASS\tGENE_PLACEHOLDER\tDP:AP\t0:0\t{all_counts}\n")
            f.close()
    
    print("Done.")




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-patient_id", required=True)
    parser.add_argument("-samples", nargs='+', required=True)
    parser.add_argument("-data_dir", required=True)
    parser.add_argument("--zero_coverage_ok", default=False, type=lambda x: (str(x).lower() == '1'))
    parser.add_argument("--single_combined_vcf", default=False, type=lambda x: (str(x).lower() == '1'))
    #parser.add_argument("--use_vcf_counts", default=False, type=lambda x: (str(x).lower() == '1'))
    
    args = parser.parse_args()

    variants_pkl = os.path.join(args.data_dir, f"{args.patient_id}_variants.pkl")

    all_variants = defaultdict(list)
    with open(variants_pkl, 'rb') as f:
        loaded_variants = pickle.load(f)

    for variant in loaded_variants:
        all_variants[f"{variant.chrom}_{variant.pos}"].append(variant)
    

    report_file = os.path.join(args.data_dir, "vcf_bamcounts_report.txt")
    with open(report_file, 'w') as f:
        f.write("variant\tsample\tstrelka_vcf_depth:strelka_vcf_alt_count\tmutect_vcf_depth:mutect_vcf_alt_count\tbam_depth:bam_alt_count\n")
    
    # TODO: consider incorporating 'NORMAL' column data from original vcf file (currently just setting 'NORMAL' with DP:AP = 0:0 in output vcf)

    sample_to_variants = dict()

    for sample in args.samples:
        if sample == 'Normal': continue
        print(f"{sample}: Reading bam-readcount output...")
        bamcounts_file = os.path.join(args.data_dir, sample, f'{sample}_bamcounts.txt')
        variant_to_counts = read_bamcounts(bamcounts_file, all_variants, sample, report_file)
        sample_to_variants[sample] = variant_to_counts
        print("Done.\n")


    variant_sets = (set(sample_to_variants[s].keys()) for s in args.samples if s != 'Normal')
    final_variants = set.union(*variant_sets) if args.zero_coverage_ok else set.intersection(*variant_sets)
    
    out_dir = os.path.join(args.data_dir, "..", "..","VCF", args.patient_id)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    print(f"Output Directory: {out_dir}\n")

    # TODO: implement option to use vcf counts instead of bamcounts for original sample variants 
    # same with vcf and bam count average?

    write_vcf_file(sample_to_variants, final_variants, out_dir, args.single_combined_vcf)


if __name__ == "__main__":
    main()