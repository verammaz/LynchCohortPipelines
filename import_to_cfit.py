import argparse
import os
import json
from collections import defaultdict
import pickle

from union_variants_pre import Variant
from union_variants_post import read_bamcounts

from cfit.util.Analysis import Analysis
from cfit.plot.PlotTreeAll import PlotTreeAll

def fix_vcf_format(hdir, patient_id, mapping):
    samples = []
    for pat_mapping in mapping:
        print(pat_mapping)
        for s in pat_mapping['samples']:
            samples.append(s[1])

    vcf_dir = os.path.join(hdir, 'VCF', patient_id)

    all_variants = defaultdict(list)

    for sample in samples:
        vcf_file = os.path.join(vcf_dir, f'{sample}.vcf')
        with open(vcf_file, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'): continue
                else:
                    chrom, pos, id, ref, alt = line.split('\t')[:5]
                    variant = Variant(chrom, pos, ref, alt, None, None, None, None)
                    all_variants[f"{chrom}_{pos}"].append(variant)

    
    sample_to_variants = dict()

    for sample in samples:
        bamcounts_file = os.path.join(args.data_dir, sample, f'{sample}_bamcounts.txt')
        variant_to_counts = read_bamcounts(bamcounts_file, all_variants, sample)
        sample_to_variants[sample] = variant_to_counts
    

    
    out_dir = os.path.join(args.data_dir, "..", "..","VCF", args.patient_id)

    for file in os.listdir(out_dir):
        if file.endswith('.vcf'):
            sample_name = file.split('_')[0]
            with open(file, 'r') as f_in:
                file_lines = f_in.readlines()
            f_in.close()
            with open(file, 'w') as f_out:
                for line in file_lines:
                    if line.startswith('#'):
                        f_out.write(line)
                    line_components = line.split('\t')
                    variant_id = line_components[2]
                    line_components[-1] = sample_to_variants[sample_name][variant_id] + '\n'
                    f_out.write(('\t').join(line_components))
            f_out.close()

    
    



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="Analysis data home directory")
    parser.add_argument('-patient_id', required=True, help='Patient identifier')
    parser.add_argument('-drivers', nargs='+', help='List of driver mutations to mark on clone tree plot')
    parser.add_argument('-fix_vcf_format', default=False, action='store_true',
                        help='there WAS a bug in union_variants_post.py for format of VCF counts, need to re-run CFIT and re-write VCF files')


    args = parser.parse_args()

    # Read files with driver gene names and recurring antigens

    drivers = args.drivers

    if drivers is None:
        drivers = []

        if os.path.exists(os.path.join(args.hdir, 'adenoma_drivers.txt')):
            with open(os.path.join(args.hdir, 'adenoma_drivers.txt'), 'r') as f:
                drivers_adenoma = [line.strip() for line in f.readlines()]
                drivers.extend(drivers_adenoma)

        if os.path.exists(os.path.join(args.hdir, 'carcinoma_drivers.txt')):
            with open(os.path.join(args.hdir, 'carcinoma_drivers.txt'), 'r') as f:
                drivers_carcinoma = [line.strip() for line in f.readlines()]
                drivers.extend(drivers_carcinoma)

        if os.path.exists(os.path.join(args.hdir, 'adenoma_shared_fs_interpatient.txt')):
            with open(os.path.join(args.hdir, 'adenoma_shared_fs_interpatient.txt'), 'r') as f:
                fs_recur_interpatient = [line.strip() for line in f.readlines()]
                drivers.extend(fs_recur_interpatient)

        if os.path.exists(os.path.join(args.hdir, 'adenoma_shared_fs_intrapatient.txt')):
            with open(os.path.join(args.hdir, 'adenoma_shared_fs_intrapatient.txt'), 'r') as f:
                fs_recur_intrapatient = [line.strip() for line in f.readlines()]
                drivers.extend(fs_recur_intrapatient)


    
    anl = Analysis()
    anl.set_MHC_version("pan41")
   
    with open(os.path.join(args.hdir, "cfit_config.json")) as f:
        config = json.load(f)

    with open(os.path.join(args.hdir, "cfit_mapping.json")) as f1:
        mappingjs = json.load(f1)


    mapping = []

    # get patient mapping
    for data_dict in mappingjs:
        if data_dict['name'] == args.patient_id or data_dict['pname'] == args.patient_id:
            mapping.append(data_dict)
    
    if args.fix_vcf_format:
        fix_vcf_format(args.hdir, args.patient_id, mapping)


    # initialize cfit patient(s):
    anl.initialize_config(config, mapping, dir=args.hdir, kd_thr=500, ns=[9], tree_format="pairtree")

    # to generate all graphs and save into single html file:
    for patient in anl.patients.values():
        pat_dir = os.path.join(args.hdir, config["tree_dir"], patient.name)
        plots = PlotTreeAll(patient, drivers=drivers, outdir=pat_dir, tree_format='pairtree')

