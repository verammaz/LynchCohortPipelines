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



"""
1) Script to input a list of lesions and then compare shared frameshift variant at loci across either VCF or neopipe (4_4636272626_GA_G)
2) Script to input a list of lesions and output which of the 46 variants are present in any lesion in that list (see excel sheet)
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


def get_fs_variants(lesion, patient, variant_to_nmd, hdir, fs_annotation, check_raw, v=False):
    fs_variants = defaultdict(list)

    raw_variants = get_raw_variants(lesion, patient, hdir) if check_raw else None
    
    with open(os.path.join(hdir, 'VCF', patient, lesion+'_ann.vcf'), 'r') as snpeff_file, open(os.path.join(hdir, 'VCF', patient, lesion+'_varcode.vcf'), 'r') as varcode_file:
        snpeff_lines = snpeff_file.readlines()
        varcode_lines = varcode_file.readlines()

        for snpeff_line, varcode_line in zip(snpeff_lines, varcode_lines):
        
            assert(snpeff_line.split('\t')[i] == varcode_line.split('\t')[i] for i in range(len(snpeff_line.split('\t'))) if i != 7)
            
            if snpeff_line.startswith('#'): continue
            
            variant = snpeff_line.split('\t')[2]
            count = snpeff_line.split('\t')[10].strip()
            
            if not check_raw and (count == '0:0' or count.split(':')[-1].strip() == count.split(':')[0].strip()):
                if v:
                    print(f"Warning: variant {variant} has count {count} in {lesion}")
                continue
            
            elif check_raw and not variant in raw_variants:
                if v:
                    print(f"Warning: variant {variant} not present with 'PASS' filter in raw (strelka/mutect) vcf file for {lesion}")
                continue

            variant_to_nmd[('_').join(variant.split('_')[0:2])] = 1 if 'nonsense_mediated_decay' in snpeff_line.split('\t')[7] else 0

            if fs_annotation:
                snpeff_ann, varcode_ann = snpeff_line.split('\t')[7], varcode_line.split('\t')[7]
                if 'frameshift' in snpeff_ann or 'FrameShift' in varcode_ann:
                    fs_variants[('_').join(variant.split('_')[0:2])].append(variant)

            else:
                ref, alt = variant.split('_')[-2], variant.split('_')[-1]
                if (len(ref) == 1 and len(alt) != 1) or (len(alt) == 1 and len(ref) != 1) :
                    fs_variants[('_').join(variant.split('_')[0:2])].append(variant)

        return fs_variants
    

def get_xl_outfile_mode(out_file):
    if os.path.exists(out_file):
        if os.path.getsize(out_file) == 0:
            return 'w'
        else:
            return 'a'
    else:
        return 'w'
    


def write_frameshifts(writer, lesions, lesion_to_fsvariants, variant_to_nmd):
     for n in range(2, len(lesions) + 1):
            combinations = itertools.combinations(lesions, n)
            for subset in combinations:
                variant_sets = (set(lesion_to_fsvariants[s].keys()) for s in subset)
                shared_variants = list(set.intersection(*variant_sets))
                df = pd.DataFrame()
                if len(shared_variants):
                    df = pd.DataFrame(columns=shared_variants, index=[s for s in subset])
                    for s in subset:
                        df.loc[s] = pd.Series({var:(',').join(lesion_to_fsvariants[s][var]) for var in shared_variants})
                    df.loc['NMD'] = pd.Series({var: variant_to_nmd[var] for var in shared_variants})
                df.to_excel(writer, index_label=f'Total: {len(shared_variants)}', index=True, sheet_name=(',').join(list(subset)))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="path to home directory with /VCF folder")
    parser.add_argument('-lesions', required=True, help="comma-separated list of lesion ids")
    parser.add_argument('-patients', required=True, help="comma-separated list of patient ids, corresponding to lesion ids list")
    parser.add_argument('-fs_annotation', default=False, action='store_true', help='count frameshift only through annotation (i.e. not just looking if indel is multiple of three)')
    parser.add_argument('-specific_fs_xl', help="excel file with specific frameshift variants and their coordinates in a dataframe")
    parser.add_argument('-check_raw', default=False, action='store_true', help="check for variant in lesion raw (mutect/strelka) vcf file before counting it in a lesion")
    parser.add_argument('-nmd', default=False, action='store_true', help='add nonsense-mediated-decay row to output table')
    parser.add_argument('-v', default=False, action='store_true', help='enable warning messages')
    args = parser.parse_args()

    lesions = args.lesions.split(',')
    patients = args.patients.split(',')

    if len(lesions) != len(patients):
        print("Error: length of lesions and patients list do not match.")
        return 

    if len(lesions) < 2:
        print("Error: please specify at least two lesion ids.")
        return

    lesion_to_fsvariants = dict()
    variant_to_nmd = dict()

    print(lesion_to_fsvariants, variant_to_nmd)

    for lesion, patient in zip(lesions, patients):
        lesion_to_fsvariants[lesion] = get_fs_variants(lesion, patient, variant_to_nmd, args.hdir, args.fs_annotation, args.check_raw, v=args.v)

    outdir = os.path.join(args.hdir, 'LesionVariantsComparisons')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    out_file = os.path.join(outdir, 'shared_frameshifts.xlsx')
   
    mode = get_xl_outfile_mode(out_file)

    
    # Write or append
    if mode == 'a':
        with pd.ExcelWriter(out_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            write_frameshifts(writer, lesions, lesion_to_fsvariants, variant_to_nmd)
    else:
        with pd.ExcelWriter(out_file, engine='openpyxl', mode='w') as writer:
            write_frameshifts(writer, lesions, lesion_to_fsvariants, variant_to_nmd)

    


    if args.specific_fs_xl is not None:
        out_file = os.path.join(outdir, 'specific_fs_variants.xlsx')

        mode = get_xl_outfile_mode(out_file)

        df = pd.read_excel(args.specific_fs_xl)
        chroms = [chrom[3:] for chrom in list(df['HG19_id'])]
        coords = list(df['COORD_HG19'])
        genes = list(df['GENE'])
        out_df = pd.DataFrame()

        if mode == 'a':
            out_df = pd.read_excel(out_file, sheet_name='specific fs variants', index_col=0) 
        else:
            out_df['chrom_pos'] = [f'{chrom}_{coord}' for chrom, coord in zip(chroms, coords)]
            out_df['gene_id'] = genes
        
        for lesion in lesions:
            if f'lesion_{lesion}' in out_df.columns: continue
            fs_presence = ['+' if var in lesion_to_fsvariants[lesion].keys() else '-' for var in out_df['chrom_pos']]
            out_df[f'lesion_{lesion}'] = fs_presence
        
        plus_df = out_df.applymap(lambda x: 1 if x == '+' else 0)
        out_df.loc['total_variants'] = plus_df.sum(axis=0)
        out_df.loc[:,'total_lesions'] = plus_df.sum(axis=1)


        if mode == 'a':
            with pd.ExcelWriter(out_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                out_df.to_excel(writer, index=True, sheet_name='specific fs variants')
        else:
            with pd.ExcelWriter(out_file, engine='openpyxl', mode='w') as writer:
                out_df.to_excel(writer, index=True, sheet_name='specific fs variants')



if __name__ == "__main__":
    main()