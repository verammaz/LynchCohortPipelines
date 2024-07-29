import argparse 
import os
import json
import pandas as pd


def infer_sex(variant_chroms):
    if len([chrom for chrom in variant_chroms if chrom == 'Y']) > 0:
        return 'male'
    else:
        return 'female'
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-patient_id", required=True)
    parser.add_argument("-vcf_dir", required=True)
    parser.add_argument("-patient_sex_info_file")
    parser.add_argument("-single_sample_trees", default=False, type=lambda x: (str(x).lower() == '1'))

    args = parser.parse_args()

    tree_dir = os.path.join(args.vcf_dir, "..", "..", "Pairtree", args.patient_id)

    if not os.path.exists(tree_dir):
        os.makedirs(tree_dir)
    
    vcf_files = []
    sample_names = []

    for file_name in os.listdir(args.vcf_dir):
        if file_name.startswith('.'): continue
        if file_name.endswith('.vcf'):
            file_path = os.path.join(args.vcf_dir, file_name)
            vcf_files.append(file_path)
            sample_names.append(os.path.splitext(os.path.basename(file_path))[0])

    ssm_file = os.path.join(tree_dir, f"{args.patient_id}.ssm")
    params_file = os.path.join(tree_dir, f"{args.patient_id}.params.json")
    
    sex = ''

    if args.patient_sex_info_file is not None:
        df = pd.read_csv(args.file_info, sep='\t')
        sex = df[df['patient'] == args.patient_id].loc['sex']
    
    # Create combined tree input files
    with open(ssm_file, 'w') as ssm:
        ssm.write("id\tname\tvar_reads\ttotal_reads\tvar_read_prob\n")
        file_handles = [open(file, 'r') for file in vcf_files]
        assert all(next(file_handle).strip().startswith('#') for file_handle in file_handles)
        i = 0
        if sex == '':
            lines = file_handles[0].readlines()
            sex = infer_sex(line.split('\t')[0] for line in lines if not line.startswith('#'))
            file_handles[0].seek(0)
            file_handles[0].readline()
        for lines in zip(*file_handles):
            variants = [line.split('\t')[2].strip() for line in lines]
            assert all(variant == variants[0] for variant in variants)
            var_reads = [line.split('\t')[10].split(':')[1].strip() for line in lines]
            total_reads = [line.split('\t')[10].split(':')[0].strip() for line in lines]
            var_prob = '0.001' if variants[0].split('_')[0] in ['X', 'Y'] and sex=='male' else '0.499'
            var_probs = [var_prob for line in lines]
            ssm.write(f"s{i}\t{variants[0]}\t{(', ').join(var_reads)}\t{(', ').join(total_reads)}\t{(', ').join(var_probs)}\n")
            i += 1
        for file in file_handles:
            file.close()

    D = {"clusters" : [], "garbage": [], "samples": [s for s in args.samples]}

    with open(params_file, 'w') as f:
        json.dump(D, f)


    # Create single sample tree input files
    if args.single_sample_trees:

        for vcf, sample in zip(vcf_files, sample_names):

            sample_tree_dir = os.path.join(tree_dir, sample)
            os.mkdir(sample_tree_dir)
            ssm_file = os.path.join(sample_tree_dir, f"{sample}.ssm")
            params_file = os.path.join(sample_tree_dir, f"{sample}.params.json")

            with open(vcf, 'r') as infile, open(ssm_file, 'w') as outfile:
                ssm.write("id\tname\tvar_reads\ttotal_reads\tvar_read_prob\n")
                for line in infile.readlines():
                    variant = line.split('\t')[2].strip()
                    var_reads = line.split('\t')[10].split(':')[1].strip()
                    total_reads = line.split('\t')[10].split(':')[0].strip()
                    var_prob = '0.001' if variants[0].split('_')[0] in ['X', 'Y'] and sex=='male' else '0.499'
                    outfile.write(f"s{i}\t{variant}\t{var_reads}\t{total_reads}\t{var_prob}\n")
                    i += 1
                
            D = {"clusters" : [], "garbage": [], "samples": [sample]}

            with open(params_file, 'w') as f:
                json.dump(D, f)


if __name__ == "__main__":
    main()
   

