import argparse
import json
import os
import gzip
import pickle
from collections import defaultdict
import pandas as pd

class Variant:
    def __init__(self, chrom, pos, ref, alt, sample, alt_count=None, total_count=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.alt_counts = [alt_count]
        self.total_counts = [total_count]
        self.samples = [sample]
        self.id = f"{chrom}_{pos}_{ref}_{alt}"
    
    def add_sample_origin(self, sample_name, alt_count, total_count):
        self.samples.append(sample_name)
        self.alt_counts.append(alt_count)
        self.total_counts.append(total_count)
    
    def get_sample_counts(self, sample):
        assert(sample in self.samples)
        for i, s in enumerate(self.samples):
            if sample == s:
                return (self.alt_counts[i], self.total_counts[i])


def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')


def parse_strelka_vcf_line(line):
    alt_count, total_count = "", ""
    line_fields = line.split('\t')
    format = line_fields[8]
    counts = line_fields[10].split(":") #only care about 'TUMOR' column1
    ref, alt = line_fields[3], line_fields[4]

    if len(alt) != len(ref):
        if format != "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50":
            #print(f"Warning: strelka indel format not recognized.")
            return None, None
        alt_count = counts[3].split(",")[1] # use TIR tier2
        total_count = counts[0]
    else:
        if format != "DP:FDP:SDP:SUBDP:AU:CU:GU:TU":
            #print(f"Warning: strelka snv format not recognized.")
            return None, None
        BASE_TO_FORMATID = {"A": 4, "C": 5, "G": 6, "T": 7}
        total_count, alt_count = counts[0], counts[BASE_TO_FORMATID[alt]].split(',')[1]
    
    return total_count, alt_count


def parse_mutect_vcf_line(line):
    line_fields = line.split('\t')
    format = line_fields[8]
    if format != "GT:AD:AF:DP:F1R2:F2R1:FAD:SB":
        #print("Warning: mutect snv format not recognized.")
        return None, None
    counts = line_fields[10].split(":")
    alt_count = counts[1][1]
    total_count = counts[3]
    return total_count, alt_count




def read_vcf(file, sample_name, variants_dict, pass_filter=True):

    caller = None  #assume vcf file name has caller info
    if 'strelka' in file:
        caller = 'strelka'
    elif 'mutect' in file:
        caller = 'mutect'  
    
    regions = []
    f = open_file(file)
    for line in f.readlines():
        if line.startswith('#'):
            continue
        if pass_filter and ('PASS' not in line):
            continue

        # VCF file columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NORMAL, TUMOR
        
        # Ignoring 'NORMAL' counts --> #TODO: add option to not ignore
        # TODO: add option to not overwrite existing vcf 'TUMOR' counts / report discrepancies between vcf and bam counts

        chrom, pos, id, ref, alt, _, _, _, _, _, tumor_counts = line.split('\t')
        
        total_counts, alt_counts = "",""
        
        if caller == 'strelka':
            total_counts, alt_counts = parse_strelka_vcf_line(line) 
        elif caller == 'strelka':
            total_counts, alt_counts = parse_mutect_vcf_line(line)
        else:
            print(f"Warning: {file} unkown variant caller")
            return None
        
        if total_counts == None and alt_counts == None:
            print(f"Warning: unable to parse vcf line")
            continue

        variant_id = f"{chrom}_{pos}_{ref}_{alt}"

        if variant_id in variants_dict.keys():
            variant = variants_dict[variant_id]
            variant.add_sample_origin(sample_name, alt_counts, total_counts)

        else:
            variants_dict[variant_id] = Variant(chrom, pos, ref, alt, sample_name, alt_counts, total_counts)
        
        if len(ref) > len(alt): # need to handle deletion region
            pos = str(int(pos)+len(alt))

        regions.append(f"{chrom}\t{pos}\t{pos}")
        

    return set(regions)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-patient_id", required=True)
    parser.add_argument("-file_info", required=True, help="config file with sample file infos")
    parser.add_argument("-samples", required=True, nargs="+", help="list of sample ids to ignore")

    args = parser.parse_args()

    file_info = pd.read_csv(args.file_info, sep='\t')
    
    patient_info = file_info.loc[file_info['Patient'] == args.patient_id]

    out_dir = patient_info['Directory'].iloc[0]
    Samples = {}

    for index, row in patient_info.iterrows():
        sample = row['Sample']
        if sample in args.samples:
            Samples[sample] = {'bam': row['BAM'],'vcf': row['VCF'].split(', ')}
    
    all_regions = set()
    all_variants = dict()

    parse_error = False
    for sample in Samples.keys():
        print(f"{sample}: Reading and parsing vcf file(s)...", end="")
        for vcf in Samples[sample]["vcf"]:
            vcf_file = os.path.join(out_dir, sample, vcf)
            regions = read_vcf(vcf_file, sample, all_variants)
            if regions == None:
                parse_error = True
                print("Error.")
            else:
                all_regions.update(regions)
                print("Done.")
    

    if not parse_error:
        regions_file = os.path.join(out_dir, f"{args.patient_id}_regions.txt")
        print(f"Writing regions to {regions_file} ...", end="")
        with open(regions_file, 'w') as f:
            f.write("\n".join(all_regions))
        print("Done.")
        
        variants_file = os.path.join(out_dir, f"{args.patient_id}_variants.pkl")
        print(f"Writing variant objects to {variants_file} ...", end="")
        variants_list = []
        for _, variant in all_variants.items():
            variants_list.append(variant)
        with open(variants_file, 'wb') as f:
            pickle.dump(variants_list, f)
        print("Done.")


if __name__ == "__main__":
    main()
        
