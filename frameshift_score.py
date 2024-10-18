import os
import argparse
from collections import defaultdict
import itertools
import pandas as pd
import gzip
from tqdm import tqdm 
import logging
import math
from openpyxl import load_workbook

# For each frameshift peptide, report total number of <50nm 9mers and total number of potential 9mers encoded by the peptide

def get_raw_variants(lesion, patient, hdir):
    variants = set()
    raw_vcf_files = [os.path.join(hdir, 'Raw', patient, lesion, f'{lesion}_vs_Normal.mutect2.filtered.vcf.gz'),
                     os.path.join(hdir, 'Raw', patient, lesion, f'{lesion}_vs_Normal.strelka.somatic_indels.vcf.gz'),
                     os.path.join(hdir, 'Raw', patient, lesion, f'{lesion}_vs_Normal.strelka.somatic_snvs.vcf.gz')]

    #TODO: cache set of raw variants to pkl file

    for file in raw_vcf_files:
        f = gzip.open(file, 'rt') if file.endswith('.gz') else open(file, 'r')
        for line in f.readlines():
            if line.startswith('#'): continue
            if 'PASS' in line:
                chrom, pos, ref, alt = line.split('\t')[0], line.split('\t')[1], line.split('\t')[3], line.split('\t')[4]
                variant = f"{chrom}_{pos}_{ref}_{alt}"
                variants.add(variant)
    return variants

def get_fs_variants(patients, lesions, hdir, check_raw, fs_annotation):
    result = defaultdict(list)
    for lesion, patient in zip(lesions, patients):

        raw_variants = get_raw_variants(lesion, patient, hdir) if check_raw else None
        
        with open(os.path.join(hdir, 'VCF', patient, lesion+'_ann.vcf'), 'r') as snpeff_file, open(os.path.join(hdir, 'VCF', patient, lesion+'_varcode.vcf'), 'r') as varcode_file:
            snpeff_lines = snpeff_file.readlines()
            varcode_lines = varcode_file.readlines()

            for snpeff_line, varcode_line in zip(snpeff_lines, varcode_lines):
            
                assert(snpeff_line.split('\t')[i] == varcode_line.split('\t')[i] for i in range(len(snpeff_line.split('\t'))) if i != 7)
                
                if snpeff_line.startswith('#'): continue
                
                variant = snpeff_line.split('\t')[2]
                count = snpeff_line.split('\t')[10].strip()
                
                if not check_raw and (count == '0:0' or count.split(':')[-1].strip() == '0'):
                    print(f"Warning: variant {variant} has count {count} in {lesion}")
                    continue
                
                elif check_raw and not variant in raw_variants:
                    print(f"Warning: variant {variant} not present with 'PASS' filter in raw (strelka/mutect) vcf file for {lesion}")
                    continue

                if fs_annotation:
                    snpeff_ann, varcode_ann = snpeff_line.split('\t')[7], varcode_line.split('\t')[7]
                    if 'frameshift' in snpeff_ann or 'FrameShift' in varcode_ann:
                        result[variant].append(lesion)

                else:
                    ref, alt = variant.split('_')[-2], variant.split('_')[-1]
                    if (len(ref) == 1 and len(alt) != 1) or (len(alt) == 1 and len(ref) != 1) :
                        result[variant].append(lesion)

    return result



def get_neoantigens(patients, hdir, kd=500):
    variant_to_neos = dict()
    passed_variants = defaultdict(set) # variants included in list passed to netMHC --> some NeoPipe filters (see SnpEff.py)
    for patient in patients:
        var_to_neo = defaultdict(list)
        pass_var = defaultdict(int)
        noe_files =[os.path.join(hdir, 'Neoantigens', 'pan41', f'neoantigens_{patient}.txt'), os.path.join(hdir, 'Neoantigens', 'pan41', f'neoantigens_other_{patient}.txt')]
        for file in noe_files:
            f = open(file, 'r')
            f.readline()
            for line in f.readlines():
                variant = line.split('\t')[1] if 'other' not in file else ('_').join(line.split('\t')[1].split('_')[:-1])
                pass_var[variant] += 1
                neo = line.split('\t')[3]
                score = line.split('\t')[6]
                if float(score) <= float(kd): # TODO: check this with Matt
                    var_to_neo[variant].append(neo)
        variant_to_neos[patient] = var_to_neo
        passed_variants[patient] = pass_var
    return passed_variants, variant_to_neos

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="path to home directory with /VCF folder")
    parser.add_argument('-lesions', required=True, help="comma-separated list of lesion ids")
    parser.add_argument('-patients', required=True, help="comma-separated list of patient ids, corresponding to lesion ids list")
    parser.add_argument('-fs_annotation', default=False, action='store_true', help='count frameshift only through annotation (i.e. not just looking if indel is multiple of three)')
    parser.add_argument('-check_raw', default=False, action='store_true', help="check for variant in lesion raw (mutect/strelka) vcf file before counting it in a lesion")
    args = parser.parse_args()

    lesions = args.lesions.split(',')
    patients = args.patients.split(',')

    if len(lesions) != len(patients):
        print("Error: length of lesions and patients list do not match.")
        return 

    if len(lesions) < 2:
        print("Error: please specify at least two lesion ids.")
        return

    fsvariant_to_lesions = get_fs_variants(patients, lesions, args.hdir, args.check_raw, args.fs_annotation)
    passed_variants, variant_to_neos = get_neoantigens(patients, args.hdir, kd=50)

    data = dict()

    for patient in patients:
        for variant in passed_variants[patient]:
            total_nmers = passed_variants[patient][variant]
            sub50_nmers = len(variant_to_neos[patient][variant])
            if variant not in fsvariant_to_lesions.keys():
                print(f"variant {variant} not found in variant->lesion dict")
            lesion_source = (', ').join(fsvariant_to_lesions[variant])
            if variant in data.keys() and (data[variant]['total_nmers'] != total_nmers or data[variant]['sub50_nmers'] != sub50_nmers):
                print(f"data mismatch: curr entry: {data[variant]}, new entry: total={total_nmers}, sub50={sub50_nmers}")
            elif variant not in data.keys():
                data[variant] = dict()
                data[variant]['total_nmers'] = total_nmers
                data[variant]['sub50_nmers'] = sub50_nmers
                data[variant]['lesions'] = set(lesion_source)
            else:
                data[variant]['lesions'].add(lesion_source)

    outdir = os.path.join(args.hdir, 'LesionVariantsComparisons')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    
    out_file = os.path.join(outdir, 'frameshift_scores.xlsx')
    out_df = pd.DataFrame(columns=['frameshift variant', 'total nmers', 'sub 50nm nmers', 'lesions'])
    if os.path.exists(out_file):
        out_df = pd.read_excel(out_file, index_col=0) 
    out_df.index = out_df.index.map(str)

    for variant in data.keys():
        out_df.loc[variant] = pd.Series({'frameshift variant': variant,
                                         'total nmers': data[variant]['total_nmers'],
                                         'sub 50nm nmers': data[variant]['sub50_nmers'],
                                         'lesions': data[variant]['lesions']})
    
    with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
       out_df.to_excel(writer)
 
    

if __name__ == "__main__":
    main()