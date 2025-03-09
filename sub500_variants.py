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
    variants_vcf_data = dict()
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
            variants_vcf_data[var] = {'nmd': nmd, 'vaf': vaf}
    return variants_vcf_data, max_vaf





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

    patient_to_sub500variants = dict()

    for patient in patients:
        sub500_variants = set()
        with open(sub500_variants_file, 'r') as f:
            f.readline()
            for line in f.readlines():
                var, effect, neos, scores = line.split('\t')
                vcf_data = variants_vcf_data.get(var, None)
                if vcf_data is None:
                    print(f'Warning: variant {var} present in {patient}_sub500_variants.txt but not found in {lesion}_ann.vcf')
                sub500_variants.add(var)
        patient_to_sub500variants[patient] = sub500_variants

    for lesion, patient in zip(lesions, patients):
        sub500_variants_file = os.path.join(outdir,f'{patient}_sub500_variants.txt')
        if not os.path.exists(sub500_variants_file):
            print(f"Error: file for variants that encode at least one <500nM peptide does not exist for patient {patient}.")
            print(f"Please run lesions_loads.py for lesion {lesion} and patient {patient}.")
        
        variants_vcf_data, max_vaf = get_variants(args.hdir, lesion, patient)
        
        sub500_variants_data = dict()
        for var in patient_to_sub500variants[patient]:
            sub500_variants_data[var] = {'nmd': vcf_data['nmd'], 'ccf_vcf': vcf_data['vaf'] / max_vaf}
        
        out_df = pd.DataFrame()
        out_df.columns = ['variant', 'nmd', 'ccf_vcf', 'ccf_cfit']

 


