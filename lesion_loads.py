import os
import argparse
from collections import defaultdict
import itertools
import pandas as pd
import gzip
import numpy as np

"""
3) For each lesion, output total load of: 
    a) frameshifts 
    b) nonsynonymous substitutions (aka missense_variant in SnpEff) 
    c) in-frame insertion/deletion 
    d) frameshift truncation (aka stop_gained in SnpEff)
4) For each lesion, output total load (subset of which encode at least one peptide predicted to bind <500nM) of: 
    a) frameshifts 
    b) nonsynonymous substitutions (aka missense_variant in SnpEff) 
    c) in-frame insertion/deletion 
    d) frameshift truncation (aka stop_gained in SnpEff)
"""

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



def get_neoantigens(patient, hdir, kd=500):
    variant_to_neos = defaultdict(list)
    noe_files =[os.path.join(hdir, 'Neoantigens', 'pan41', f'neoantigens_{patient}.txt'), os.path.join(hdir, 'Neoantigens', 'pan41', f'neoantigens_other_{patient}.txt')]
    for file in noe_files:
        f = open(file, 'r')
        f.readline()
        for line in f.readlines():
            variant = line.split('\t')[1]
            neo = line.split('\t')[3]
            score = line.split('\t')[6]
            if float(score) >= float(kd):
                variant_to_neos[variant].append(neo)
    return variant_to_neos


def get_lesion_neo_loads(lesion_to_effectvariants, variant_to_neos):
    result = {'total': 0, 'fs': 0, 'nonsyn': 0, 'inframe_indel': 0, 'fs_trunc': 0, 'pre_stop': 0}
    for effect, variants in lesion_to_effectvariants.items():
        for var in variants[1]:
            if len(variant_to_neos[var]) >= 1:
                result[effect] += 1
    return result

def check_annotation(variant, snpeff_ann, varcode_ann, outfile):
    f = open(outfile, 'a')
    annotation_mapping = {'Substitution': 'missense_variant', 'FrameShiftTruncation': 'stop_gained', 'FrameShift': 'frameshift', 'PrematureStop': 'stop_gained'}
    snpeff = snpeff_ann.split(',')
    varcode = varcode_ann.split('--')
    snpeff_effects = [ann.split('|')[1] for ann in snpeff if ann.split('|')[1] in annotation_mapping.values()]
    varcode_effects = [ann.split('(')[0] for ann in varcode[1:] if ann.split('(')[1] in annotation_mapping.keys()]
    snpeff_effects_str = (', ').join(list(set(snpeff_effects)))
    varcode_effects_str = (', ').join(list(set(varcode_effects)))
    #print(varcode, snpeff)
    if (set([annotation_mapping[effect] for effect in varcode_effects]) != set(snpeff_effects) or
        len(list(set(snpeff_effects))) > 1 or len(list(set(varcode_effects))) > 1):
        f.write(f'{variant}\t{snpeff_effects_str}\t{varcode_effects_str}\n')


def get_lesion_variants(lesions, patients, args, outdir):
    lesion_to_effectvariants = dict()

    for lesion, patient in zip(lesions, patients):
        lesion_to_effectvariants[lesion] = {'total':[0,[]], 'fs':[0,[]], 'nonsyn':[0,[]], 'inframe_indel':[0,[]], 'fs_trunc':[0,[]], 'pre_stop':[0,[]]}
        raw_variants = get_raw_variants(lesion, patient, args.hdir) if args.check_raw else None
        total, fs, nonsyn, inframe_indel, fs_trunc, pre_stop = 0, 0, 0, 0, 0, 0
        with open(os.path.join(args.hdir, 'VCF', patient, lesion+'_ann.vcf'), 'r') as snpeff_file, open(os.path.join(args.hdir, 'VCF', patient, lesion+'_varcode.vcf'), 'r') as varcode_file:
            for snpeff_line, varcode_line in zip(snpeff_file.readlines(), varcode_file.readlines()):
                assert(snpeff_line.split('\t')[i] == varcode_line.split('\t')[i] for i in range(len(snpeff_line.split('\t'))) if i != 7)
                if snpeff_line.startswith('#'): continue
                variant = snpeff_line.split('\t')[2]
                count = snpeff_line.split('\t')[10].strip()
                
                if not args.check_raw and (count == '0:0' or count.split(':')[-1].strip() == '0'):
                    #print(f"Warning: variant {variant} has count {count} in {lesion}")
                    continue
                
                elif args.check_raw and not variant in raw_variants:
                    #print(f"Warning: variant {variant} not present with 'PASS' filter in raw (strelka/mutect) vcf file for {lesion}")
                    continue

                lesion_to_effectvariants[lesion]['total'][1].append(variant)
                lesion_to_effectvariants[lesion]['total'][0] += 1

                snpeff_ann, varcode_ann = snpeff_line.split('\t')[7], varcode_line.split('\t')[7]

                if args.annotation_check:
                    check_annotation(variant, snpeff_ann, varcode_ann, os.path.join(outdir, 'annotation_checks.txt'))
                
                if args.fs_annotation:
                    if 'frameshift' in snpeff_ann or 'FrameShift' in varcode_ann:
                        lesion_to_effectvariants[lesion]['fs'][1].append(variant)
                        lesion_to_effectvariants[lesion]['fs'][0] += 1

                else:
                    ref, alt = variant.split('_')[-2], variant.split('_')[-1]
                    if (len(ref) == 1 and len(alt) != 1) or (len(alt) == 1 and len(ref) != 1) :
                        lesion_to_effectvariants[lesion]['fs'][1].append(variant)
                        lesion_to_effectvariants[lesion]['fs'][0] += 1

                if 'missense_variant' in snpeff_ann or 'Substitution' in varcode_ann:
                    lesion_to_effectvariants[lesion]['nonsyn'][1].append(variant)
                    lesion_to_effectvariants[lesion]['nonsyn'][0] += 1
                
                if abs(len(ref) - len(alt)) % 3 == 0: 
                    lesion_to_effectvariants[lesion]['inframe_indel'][1].append(variant)
                    lesion_to_effectvariants[lesion]['inframe_indel'][0] += 1
                
                if ('stop_gained' in snpeff_ann and len(ref) != len(alt)) or 'FrameShiftTruncation' in varcode_ann:
                    lesion_to_effectvariants[lesion]['fs_trunc'][1].append(variant)
                    lesion_to_effectvariants[lesion]['fs_trunc'][0] += 1

                if ('stop_gained' in snpeff_ann and len(ref) == len(alt)) or 'PrematureStop' in varcode_ann:
                    lesion_to_effectvariants[lesion]['pre_stop'][1].append(variant)
                    lesion_to_effectvariants[lesion]['pre_stop'][0] += 1

    return lesion_to_effectvariants

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="path to home directory with /VCF folder")
    parser.add_argument('-lesions', required=True, help="comma-separated list of lesion ids")
    parser.add_argument('-patients', required=True, help="comma-separated list of patient ids, corresponding to lesion ids list")
    parser.add_argument('-fs_annotation', default=False, action='store_true')
    parser.add_argument('-check_raw', default=False, action='store_true')
    parser.add_argument('-annotation_check', default=False, action='store_true')

    args = parser.parse_args()

    lesions = args.lesions.split(',')
    patients = args.patients.split(',')

    if len(lesions) != len(patients):
        print("Error: length of lesions and patients list do not match.")
        return 

    outdir = os.path.join(args.hdir, 'LesionVariantsComparisons')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    out_file = os.path.join(outdir, 'lesion_loads.xlsx')
    out_df = pd.DataFrame(columns=['total', 'frameshift', 'nonsynonymous_substitution', 'inframe_indel', 'frameshift_truncation', 'premature_stop'])
    if os.path.exists(out_file):
        out_df = pd.read_excel(out_file, index_col=000) 

    lesion_to_effectvariants = get_lesion_variants(lesions, patients, args, outdir)
    patient_to_neos = {patient: get_neoantigens(patient, args.hdir) for patient in list(set(patients))}

    
    for patient,lesion in zip(patients,lesions):
        effects = ['total', 'fs', 'nonsyn', 'inframe_indel', 'fs_trunc', 'pre_stop' ]
        neo_loads = get_lesion_neo_loads(lesion_to_effectvariants[lesion], patient_to_neos[patient])
        data = [f'{lesion_to_effectvariants[lesion][effect][0]}, {neo_loads[effect]}' for effect in effects]
        out_df.loc[lesion] = pd.Series({'total': data[0],
                                    'frameshift': data[1], 
                                    'nonsynonymous_substitution': data[2], 
                                    'inframe_indel': data[3], 
                                    'frameshift_truncation': data[4],
                                    'premature_stop': data[5]})
    
    out_df = df[::-1]
    out_df = out_df.iloc[ np.unique( out_df.index.values, return_index = True )[1] ]

    with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
       out_df.to_excel(writer)


if __name__ == "__main__":
    main()