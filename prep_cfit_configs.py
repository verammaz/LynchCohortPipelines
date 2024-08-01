import argparse
import os
import json
import pandas as pd

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True)
    parser.add_argument('-patient_id', required=True, type=str)
    parser.add_argument('-sample_info', required=True)
    parser.add_argument("-single_sample_trees", default=False, type=lambda x: (str(x).lower() == '1'))
    parser.add_argument("-hla")

    args = parser.parse_args()

    sample_info = pd.read_csv(args.sample_info, sep='\t')
    sample_info = sample_info[sample_info['Patient'] == args.patient_id]

    config_file = os.path.join(args.hdir, 'cfit_config.json')

    if not os.path.exists(config_file):
        with open(config_file, 'w') as f:
            json.dump({
                        "vcf_dir": "VCF",
                        "tree_dir": "Pairtree",
                        "aln_dir": "",
                        "neo_dir": "Neoantigens"
                        }, f)
    
    mapping_file = os.path.join(args.hdir, 'cfit_mapping.json')

    mappingjs = []

    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            mappingjs = json.load(f)
    f.close()

    new_mappingjs = []

    for data_dict in mappingjs:
        if (data_dict['pname'] != args.patient_id and 
            data_dict['name'] != args.patient_id and 
            data_dict['name'] not in [f"{s}" for s in sample_info['Sample'].tolist()]):
            new_mappingjs.append(data_dict)
    
    samples = [
            [args.patient_id, str(s), f"Sample_{s}", f"{s}_ann.vcf", 
            sample_info.loc[sample_info.Sample == s, 'TimePoint'].values[0], 
            sample_info.loc[sample_info.Sample == s, 'Tissue'].values[0], 0, 0]
            for s in sample_info['Sample'].tolist()
            ]
    
    patientjs = {
        "name": args.patient_id,
        "pname": args.patient_id,
        "cohort": "Lynch",
        "response": 0,
        "dead": 0,
        "type": "-",
        "samples": samples
    }

    new_mappingjs.append(patientjs)

    if args.single_sample_trees:
        for elem in samples:
            patient, sample, sample_name, ann_vcf, timepoint, tissue, _, _ = elem
            samplejs = {
                "name": f"{sample}",
                "pname": args.patient_id,
                "cohort": "Lynch",
                "response": 0,
                "dead": 0,
                "type": "-",
                "samples": [[patient, sample, sample_name, ann_vcf, timepoint, tissue, 0, 0]]
            }
            new_mappingjs.append(samplejs)
        
        hla_df = pd.read_csv(args.hla, sep="\t", header=None)
        hlas = hla_df[hla_df[0] == args.patient_id][1].to_list()[0]
        print(hla_df)
        print(hlas)
        for s in sample_info['Sample'].tolist():
            new_row = pd.DataFrame([[str(s), hlas]])
            hla_df = pd.concat([hla_df, new_row], ignore_index=True)
        
        hla_new_df = hla_df.drop_duplicates(keep='last')
        hla_new_df.to_csv(args.hla, sep='\t', header=None, index=False)
    
    with open(mapping_file, 'w') as f:
        json.dump(new_mappingjs, f, indent=4)




