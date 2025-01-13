import os
import argparse
from collections import defaultdict
import itertools
import pandas as pd
import gzip
import numpy as np

"""
5) For each  variant in <500nM set, report NMD escape status (from SnpEff) and CCF (from CFIT or pairtree)
"""

# assumes output file {lesion}_sub500_variants.txt from lesion_loads.py present in HDIR/LesionVariantsComparisons


def get_variants(hdir, lesion, patient):
    variants = dict()
    max_vaf = 0
    with open(os.path.join(hdir, 'VCF', patient, lesion+'_ann.vcf'), 'r') as f:
        for line in f.readlines():
            if line.startswith('#'): continue
            var = line.split('\t')[2]
            nmd = 1 if 'nonsense_mediated_decay' in line.split('\t')[7] else 0
            tumor_total, tumor_ref = line.split('\t')[10].split(':') # DP:AP format <=> total:ref
            vaf = (int(tumor_total) - int(tumor_ref)) / int(tumor_total)
            if vaf > max_vaf:
                max_vaf = vaf
            variants[var] = {'nmd': nmd, 'vaf': vaf}
    return variants, max_vaf





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="path to home directory with /VCF folder")
    parser.add_argument('-lesions', required=True, help="comma-separated list of lesion ids")
    parser.add_argument('-patients', required=True, help="comma-separated list of patient ids, corresponding to lesion ids list")
   
    args = parser.parse_args()

    lesions = args.lesions.split(',')
    patients = args.patients.split(',')

    if len(lesions) != len(patients):
        print("Error: length of lesions and patients list do not match.")
        return 
    
    outdir = os.path.join(args.hdir, 'LesionVariantsComparisons')

    for lesion, patient in zip(lesions, patients):
        sub500_variants_file = os.path.join(outdir,f'{lesion}_sub500_variants.txt')
        if not os.path.exists(sub500_variants_file):
            print(f"Error: file for variants that encode at least one <500nM peptide does not exist for lesion {lesion}.")
            print(f"Please run lesions_loads.py for lesion {lesion}")
        variants, max_vaf = get_variants(args.hdir, lesion, patient)
        

        sub500_variants = dict()
        with open(sub500_variants_file, 'r') as f:
            f.readline()
            for line in f.readlines():
                var, effect, neos, scores = line.split('\t')
                vcf_data = variants.get(var, None)
                if vcf_data is None:
                    print(f'Warning: variant {var} present in {lesion}_sub500_variants.txt but not found in {lesion}_ann.vcf')
                sub500_variants[var] = f'{var}\t{effect}\t{neos}\t{scores}\t{vcf_data['nmd']}\t{vcf_data['vaf']/max_vaf}\n'

        with open(sub500_variants_file, 'w') as f:
            f.write('variant\teffect\tneoantigens\tnetmhc_scores\tnmd_status\tccf_vcf\tccf_cfit\n')
            for _, line in sub500_variants.items():
                f.write(line)


