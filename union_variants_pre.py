import argparse
import os
import gzip
import pickle
from collections import defaultdict
import pandas as pd
import sys


# TODO: add caller attribute to variant
class Variant:
    def __init__(self, chrom, pos, ref, alt, sample, alt_count, total_count, caller):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.mut_type = 'snv' if len(ref) == len(alt) else 'indel'
        self.counts = defaultdict(dict)
        self.counts[caller]['alt_counts'] = [alt_count]
        self.counts[caller]['total_counts'] = [total_count]
        self.samples = defaultdict(list)
        self.samples[caller] = [sample]
        self.id = f"{chrom}_{pos}_{ref}_{alt}"
    
    def add_caller_sample_origin(self, caller, sample_name, alt_count, total_count):
        self.samples[caller].append(sample_name)
        if caller in self.counts.keys():
            self.counts[caller]['alt_counts'].append(alt_count)
            self.counts[caller]['total_counts'].append(total_count)
        else:
            self.counts[caller]['alt_counts'] = [alt_count]
            self.counts[caller]['total_counts'] = [total_count]
        
    def get_sample_counts(self, caller, sample_name):
        assert(sample_name in self.samples[caller])
        for i, s in enumerate(self.samples[caller]):
            if sample_name == s:
                return (self.counts[caller]['alt_counts'][i], self.counts[caller]['total_counts'][i])


def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')


def parse_strelka_vcf_line(line):
    tumor_total, tumor_alt, normal_total, normal_alt = "","","",""
    line_fields = line.split('\t')
    format = line_fields[8]
    tumor_counts = line_fields[10].split(":") 
    normal_counts = line_fields[9].split(":")
    ref, alt = line_fields[3], line_fields[4]

    if len(alt) != len(ref):
        if format == "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50":
            tumor_alt = tumor_counts[3].split(",")[1] # use TIR tier2
            tumor_total = tumor_counts[0]
            normal_alt = normal_counts[3].split(",")[1]
            normal_total = normal_counts[0]
    else:
        if format == "DP:FDP:SDP:SUBDP:AU:CU:GU:TU":
            BASE_TO_FORMATID = {"A": 4, "C": 5, "G": 6, "T": 7}
            tumor_total, tumor_alt = tumor_counts[0], tumor_counts[BASE_TO_FORMATID[alt]].split(',')[1]
            normal_total, normal_alt = normal_counts[0], normal_counts[BASE_TO_FORMATID[alt]].split(',')[1]
    
    return tumor_total, tumor_alt, normal_total, normal_alt 


def parse_mutect_vcf_line(line):
    line_fields = line.split('\t')
    format = line_fields[8]
    if format != "GT:AD:AF:DP:F1R2:F2R1:FAD:SB": 
        #print("Warning: mutect snv format not recognized.")
        return "", "", "", ""
    tumor_counts = line_fields[10].split(":") if line_fields[10].split(":")[0] == '0/1' else line_fields[9].split(":")
    normal_counts = line_fields[9].split(":") if line_fields[9].split(":")[0] == '0/0' else line_fields[10].split(":")
    tumor_alt = tumor_counts[1][1]
    tumor_total = tumor_counts[3]
    normal_alt = normal_counts[1][1]
    normal_total = normal_counts[3]
    
    return tumor_total, tumor_alt, normal_total, normal_alt 


def pass_variant_filter(tumor_total, tumor_alt, normal_total, normal_alt):
    """Filtering criteria are 
        1) total coverage for tumor ≥10, 
        2) variant allele frequency (VAF) for tumor ≥4%, 
        3) number of reads with alternative allele ≥9 for tumor, 
        4) total coverage for normal ≥7, and 
        5) VAF for normal ≤1% at a given mutation."""
    
    if int(tumor_total) < 10: return False
    if int(tumor_alt)/int(tumor_total) < 0.04: return False
    if int(tumor_alt) < 9: return False
    if int(normal_total) < 7: return False
    if int(normal_alt)/int(normal_total) > 0.01: return False

    return True



def read_vcf(file, sample_name, variants_dict, pass_filter):

    caller = None  #assume vcf file name has caller info
    if 'strelka' in file:
        caller = 'strelka'
    elif 'mutect' in file:
        caller = 'mutect' 

    if caller is None:
        print(f"Error: {file} unkown variant caller")
        return None 
    
    regions = []
    f = open_file(file)
    for line in f.readlines():
        if line.startswith('#'):
            continue
        if 'PASS' not in line:
            continue

        # VCF file columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NORMAL, TUMOR
        
        # Ignoring 'NORMAL' counts --> #TODO: add option to not ignore
        # TODO: add option to not overwrite existing vcf 'TUMOR' counts / report discrepancies between vcf and bam counts

        chrom, pos, _, ref, alt, _, _, _, _, _, _ = line.split('\t')
        
        tumor_total, tumor_alt, normal_total, normal_alt = "","","","" # Note: total_counts = depth, so not necessarily ref_counts+alt_counts = total_counts
        
        if caller == 'strelka':
            tumor_total, tumor_alt, normal_total, normal_alt = parse_strelka_vcf_line(line) 
        elif caller == 'mutect':
            tumor_total, tumor_alt, normal_total, normal_alt = parse_mutect_vcf_line(line)
        
        if tumor_total == "" or tumor_alt == "" or normal_total == "" or normal_alt == "":
            print(f"Warning: unable to parse vcf line. Continuing.")
            continue
        
        if pass_filter and not pass_variant_filter(tumor_total, tumor_alt, normal_total, normal_alt):
            print(f"Variant {chrom}_{pos}_{ref}_{alt} did not pass filtering")
            # TODO: check if variant in shared antigen gene, and keep anyway
            continue

        variant_id = f"{chrom}_{pos}_{ref}_{alt}"

        if variant_id in variants_dict.keys(): 
            variant = variants_dict[variant_id]
            variant.add_caller_sample_origin(caller, sample_name, tumor_alt.strip(), tumor_total.strip())

        else:
            variants_dict[variant_id] = Variant(chrom, pos, ref, alt, sample_name, tumor_alt.strip(), tumor_total.strip(), caller)
        
        if len(ref) > len(alt): # need to handle deletion region
            pos = str(int(pos)+len(alt))

        regions.append(f"{chrom}\t{pos}\t{pos}")
        

    return set(regions)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-patient_id", required=True)
    parser.add_argument("-samples", nargs='+', required=True)
    parser.add_argument("-data_dir", required=True)
    parser.add_argument("--additional_filter", default=False, type=lambda x: (str(x).lower() == '1'))
    parser.add_argument("--strelka_mutect_snv_intersect", default=False, type=lambda x: (str(x).lower() == '1'))
    parser.add_argument("--strelka_mutect_indel_intersect", default=False, type=lambda x: (str(x).lower() == '1'))


    args = parser.parse_args()

    sample_to_vcfs = dict()

    for sample in args.samples:
        if sample == 'Normal': continue
        samplesheet = pd.read_csv(os.path.join(args.data_dir, f"samplesheet_{sample}.csv"))
        file_info = samplesheet[samplesheet['sample'] == sample]['vcf'].values[0]
        sample_to_vcfs[sample] = file_info.split('|')
    
    out_dir = args.data_dir

    
    all_regions = set()
    all_variants = dict()

    parse_error = False
    problem_vcf = ""
    for sample, vcf_files in sample_to_vcfs.items():
        print(f"{sample}: Reading and parsing vcf file(s)...")
        for vcf in vcf_files:
            vcf_file = os.path.join(out_dir, sample, vcf)
            regions = read_vcf(vcf_file, sample, all_variants, args.additional_filter)
            if regions == None:
                parse_error = True
                problem_vcf = vcf_file
                break
            else:
                all_regions.update(regions)
                print(f"{os.path.basename(vcf_file)}: Done.")
        print("\n")
    


    if not parse_error:

        regions_file = os.path.join(out_dir, f"{args.patient_id}_regions.txt")
        print(f"Writing regions to {os.path.basename(regions_file)} ...", end="")
        with open(regions_file, 'w') as f:
            f.write("\n".join(all_regions))
        print("Done.")
        
        if args.strelka_mutect_snv_intersect: # assume both caller vcf files processed --> #TODO: add extra check
            for variant_id, variant in all_variants.items():
                if variant.mut_type == 'indel': continue
                if 'mutect' not in variant.counts.keys() or 'strelka' not in variant.counts.keys():
                    del all_variants[variant_id]
                    continue
                if set(variant.samples['mutect']).isdisjoint(set(variant.samples['strelka'])):
                    del all_variants[variant_id]
                    continue
        
        if args.strelka_mutect_indel_intersect: # assume both caller vcf files processed --> #TODO: add extra check
            for variant_id, variant in all_variants.items():
                if variant.mut_type == 'snv': continue
                if 'mutect' not in variant.counts.keys() or 'strelka' not in variant.counts.keys():
                    del all_variants[variant_id]
                    continue
                if set(variant.samples['mutect']).isdisjoint(set(variant.samples['strelka'])):
                    del all_variants[variant_id]
                    continue

        variants_file = os.path.join(out_dir, f"{args.patient_id}_variants.pkl")
        print(f"Writing variant objects to {os.path.basename(variants_file)} ...", end="")
        variants_list = []
        for _, variant in all_variants.items():
            variants_list.append(variant)
        with open(variants_file, 'wb') as f:
            pickle.dump(variants_list, f)
        print("Done.\n")

        sys.exit(0)
    
    
    else:
        print(f"Error: Unable to parse {problem_vcf}. Aborting.")
        sys.exit(1)


if __name__ == "__main__":
    main()
        
