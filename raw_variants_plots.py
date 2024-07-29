import argparse
import matplotlib.pyplot as plt
import gzip
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import itertools
from matplotlib_venn import venn2


def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def intersect_dicts(dicts):
    if not dicts:
        return {}
    # Get intersection of keys
    keys = set(dicts[0]).intersection(*dicts[1:])
    return {k: [d[k] for d in dicts if k in d] for k in keys}
    

class Sample():

    def __init__(self, name):
        self.name = name
        self.strelka_snv = dict()
        self.strelka_indel = dict()
        self.mutect_snv = dict()
        self.mutect_indel = dict()
    
    def add_variant(self, variant_id, variant_type, total_count, alt_count):
        if int(total_count) == 0: 
            return 
        if int(alt_count) / int(total_count) > 1:
            print(variant_id)
        if variant_type == 'strelka_snv':
            self.strelka_snv[variant_id] = int(alt_count) / int(total_count)
        elif variant_type == 'strelka_indel':
            self.strelka_indel[variant_id] = int(alt_count) / int(total_count)
        elif variant_type == 'mutect_snv':
            self.mutect_snv[variant_id] = int(alt_count) / int(total_count)
        elif variant_type == 'mutect_indel':
            self.mutect_indel[variant_id] = int(alt_count) / int(total_count)
    
    def get_variants(self, variant_type, thr=0):
        variants = dict()
        if variant_type == 'strelka_snv':
            variants = self.strelka_snv
        elif variant_type == 'strelka_indel':
            variants = self.strelka_indel
        elif variant_type == 'mutect_snv':
            variants = self.mutect_snv
        elif variant_type == 'mutect_indel':
            variants = self.mutect_indel

        if thr == 0:
            return variants
        result = dict()
        for variant, vaf in variants.items():
            if vaf >= thr:
                result[variant] = vaf
        return result



def parse_strelka_vcf_line(line):
    tumor_total, tumor_alt, normal_total, normal_alt = "","","",""
    line_fields = line.split('\t')
    format = line_fields[8]
    tumor_counts = line_fields[10].split(":") 
    normal_counts = line_fields[9].split(":")
    ref, alt = line_fields[3], line_fields[4]

    if len(alt) != len(ref):
        if format == "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50":
            tumor_alt = tumor_counts[3].split(",")[0] # use TIR tier1
            tumor_total = tumor_counts[0]
            normal_alt = normal_counts[3].split(",")[0]
            normal_total = normal_counts[0]
        else:
            print("strelka indel format not recognized")
    else:
        if format == "DP:FDP:SDP:SUBDP:AU:CU:GU:TU":
            BASE_TO_FORMATID = {"A": 4, "C": 5, "G": 6, "T": 7}
            tumor_total, tumor_alt = tumor_counts[0], tumor_counts[BASE_TO_FORMATID[alt]].split(',')[0]
            normal_total, normal_alt = normal_counts[0], normal_counts[BASE_TO_FORMATID[alt]].split(',')[0]
        else:
            print("strelka snv format not recognized")
    return tumor_total, tumor_alt, normal_total, normal_alt 


def parse_mutect_vcf_line(line):
    line_fields = line.split('\t')
    format = line_fields[8]
    if format != "GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB" and format != "GT:AD:AF:DP:F1R2:F2R1:FAD:SB":
        print("Warning: mutect format not recognized.")
        print(format)
        return "", "", "", ""
    tumor_counts = line_fields[10].split(":") if line_fields[10].split(":")[0] in ['0/1', '0|1'] else line_fields[9].split(":")
    normal_counts = line_fields[9].split(":") if line_fields[9].split(":")[0] in ['0/0', '0|0'] else line_fields[10].split(":")
    tumor_alt = tumor_counts[1].split(',')[1]
    tumor_total = tumor_counts[3]
    normal_alt = normal_counts[1].split(',')[1]
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


def read_vcf(file, sample_obj, pass_filter):

    caller = None  #assume vcf file name has caller info
    if 'strelka' in file:
        caller = 'strelka'
    elif 'mutect' in file:
        caller = 'mutect' 

    if caller is None:
        print(f"Error: {file} unkown variant caller")
        return None 
    
    f = open_file(file)
    for line in f.readlines():
        if line.startswith('#'):
            continue
        if 'PASS' not in line:
            continue

        chrom, pos, _, ref, alt, _, _, _, _, _, _ = line.split('\t')
        
        tumor_total, tumor_alt, normal_total, normal_alt = "","","","" # Note: total_counts = depth, so not necessarily ref_counts+alt_counts = total_counts
        
        if caller == 'strelka':
            tumor_total, tumor_alt, normal_total, normal_alt = parse_strelka_vcf_line(line) 
        elif caller == 'mutect':
            tumor_total, tumor_alt, normal_total, normal_alt = parse_mutect_vcf_line(line)
        
        if tumor_total == "" or tumor_alt == "" or normal_total == "" or normal_alt == "":
            #print(f"Warning: unable to parse vcf line. Continuing.")
            continue
        
        if pass_filter and not pass_variant_filter(tumor_total, tumor_alt, normal_total, normal_alt):
            #print(f"Variant {chrom}_{pos}_{ref}_{alt} did not pass filtering")
            continue

        variant_id = f"{chrom}_{pos}_{ref}_{alt}"

        mut_type = f'{caller}_snv' if len(ref) == len(alt) else f'{caller}_indel'
        
        sample_obj.add_variant(variant_id, mut_type, tumor_total, tumor_alt)


def plot_single_caller_vaf(Samples, caller, outdir, thr=0, separate=False):

    pairs = list(itertools.combinations(list(Samples.keys()), 2))
    
    for pair in pairs:
        
        snv_dicts = []
        indel_dicts =[]

        for sample_name in pair:
            sample_obj = Samples[sample_name]
            snv, indel = sample_obj.get_variants(f"{caller}_snv", thr), sample_obj.get_variants(f"{caller}_indel", thr)
            snv_dicts.append(snv)
            indel_dicts.append(indel)
        
        snv_intersection = intersect_dicts(snv_dicts)
        indel_intersection = intersect_dicts(indel_dicts)

        n_snv = len(snv_intersection)

        #plot snv vaf
        S1, S2 = [], []

        for variant, vaf in snv_intersection.items():
            S1.append(vaf[0])
            S2.append(vaf[1])

        for variant, vaf in indel_intersection.items():
            S1.append(vaf[0])
            S2.append(vaf[1])
            
    
        plt.scatter(S1[:n_snv], S2[:n_snv], c='red', label='SNV')
        plt.scatter(S1[n_snv:], S2[n_snv:], c='green', label='Indel')
        plt.xlabel(f'Sample {pair[0]} vaf')
        plt.ylabel(f'Sample {pair[1]} vaf')
        x = np.linspace(0, 1, 100)
        plt.plot(x, x, linestyle=':', color='gray')
        plt.title(f"{caller} VAF of all shared mutations (vaf threshold {thr})")
        plt.legend()
        plt.savefig(os.path.join(outdir, f"{pair[0]}_vs_{pair[1]}_vaf_{caller}.png"))
        plt.close()



def plot_strelka_mutect_vaf(Samples, mut, outdir, thr=0):
        
    strelka_variants = {s: Samples[s].get_variants(f"strelka_{mut}",thr) for s in Samples.keys()}
    mutect_variants = {s: Samples[s].get_variants(f"mutect_{mut}",thr) for s in Samples.keys()}
    
    assert(strelka_variants.keys() == mutect_variants.keys())
    
    colors = ['orange', 'cyan', 'pink', 'olive', 'purple', 'magenta']  # enough colors for all samples

    sample_to_color = {sample: colors[i] for i, sample in enumerate(list(strelka_variants.keys()))}

    strelka_vaf, mutect_vaf, C, labels = [], [], [], []

    for sample in Samples.keys():
        strelka = strelka_variants[sample]
        mutect = mutect_variants[sample]

        for variant in strelka.keys():
            if variant not in mutect.keys():
                continue

            strelka_vaf.append(strelka[variant])
            mutect_vaf.append(mutect[variant])

            C.append(sample_to_color[sample])
            labels.append(sample)
    

    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(strelka_vaf, mutect_vaf, c=C, alpha=0.8)
    unique_labels = list(set(labels))
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=sample_to_color[label], markersize=10) for label in unique_labels]
    x = np.linspace(0, 1, 100)
    plt.plot(x, x, linestyle=':', color='gray')
    plt.legend(handles, unique_labels, title="Samples")
    plt.xlabel('Strelka VAF')
    plt.ylabel('Mutect VAF')
    plt.title(f"Strelka vs Mutect2 VAF for {mut} mutations (vaf threshold {thr})")
    plt.savefig(os.path.join(outdir, f"strelka_vs_mutect_{mut}_vaf.png"))
    plt.close()


def plot_single_caller_venn(Samples, caller, mut, outdir, thr=0):

    variant_sets = [set(Samples[s].get_variants(f"{caller}_{mut}", thr).keys()) for s in Samples.keys()]
    set_names = [s for s in Samples.keys()]

    # Compute all intersections
    intersections = {}
    for r in range(1, len(variant_sets) + 1):
        for combination in itertools.combinations(range(len(variant_sets)), r):
            intersected_set = set.intersection(*[variant_sets[i] for i in combination])
            intersections[combination] = intersected_set

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))

    # Calculate the number of sets and their combinations
    bars = []
    labels = []

    for combination, intersected_set in intersections.items():
        bars.append(len(intersected_set))
        labels.append(', '.join([set_names[i] for i in combination]))

    # Create bar plot
    y_pos = range(len(bars))
    ax.barh(y_pos, bars, align='center', color='olive')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Number of elements')
    ax.set_title(f'{caller} {mut} sets (vaf threshold {thr})')

    # Add annotations
    for i, v in enumerate(bars):
        ax.text(v + 0.2, i, str(v), color='black', va='center')  # Adjust position for better visibility

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{caller}_{mut}_venn.png"))
    plt.close()


def plot_strelka_mutect_venn(Samples, sample, mut, outdir, thr=0):
    strelka_variants = set(Samples[sample].get_variants(f"strelka_{mut}", thr).keys())
    mutect_variants = set(Samples[sample].get_variants(f"mutect_{mut}", thr).keys())
    venn2([strelka_variants, mutect_variants], (f'strelka', f'mutect2'))
    plt.title(f'Venn diagram for sample {sample} {mut} (vaf threshold {thr})')
    plt.savefig(os.path.join(outdir, f"{sample}_{mut}_venn.png"))
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-patient_id', required=True, help="Patient id.")
    parser.add_argument('-samplesheet', required=True, help="CSV file with vcf information.")
    parser.add_argument('-hdir', required=True, help="Pipeline data home directory.")
    parser.add_argument('-additional_filter', action='store_true', help="Flag to pass variant through additional filtering.")
    
    args = parser.parse_args()
    
    samplesheet = pd.read_csv(args.samplesheet)
    samplesheet = samplesheet.loc[samplesheet['patient'] == args.patient_id]

    sample_to_vcfs= dict()

    for index, row in samplesheet.iterrows():
        sample = row['sample']
        if sample == 'Normal': continue
        sample_to_vcfs[sample] = row['vcf'].split('|')
    
    out_dir = os.path.join(args.hdir, 'Plots', args.patient_id)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    Samples = dict()

    for sample, vcf_files in sample_to_vcfs.items():
        print(f"{sample}: Parsing vcf files...", end="")
        sample_obj = Sample(sample)
        for vcf in vcf_files:
            print(f"{os.path.basename(vcf)}...", end="")
            read_vcf(vcf, sample_obj, args.additional_filter)
        Samples[sample] = sample_obj
        print("Done.")
    
    plot_single_caller_vaf(Samples, 'strelka', out_dir)
    plot_single_caller_vaf(Samples, 'mutect', out_dir)
    plot_strelka_mutect_vaf(Samples, 'snv', out_dir)
    plot_strelka_mutect_vaf(Samples, 'indel', out_dir)
    plot_single_caller_venn(Samples, 'strelka', 'snv', out_dir)
    plot_single_caller_venn(Samples, 'mutect', 'snv', out_dir)
    plot_single_caller_venn(Samples, 'strelka', 'indel', out_dir)
    plot_single_caller_venn(Samples, 'mutect', 'indel', out_dir)
    
    for sample in Samples.keys():
        plot_strelka_mutect_venn(Samples, sample, 'snv', out_dir)
        plot_strelka_mutect_venn(Samples, sample, 'indel', out_dir)

if __name__ == "__main__":
    main()

    
    



