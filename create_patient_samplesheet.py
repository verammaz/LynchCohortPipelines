import argparse
import pandas as pd
from collections import defaultdict
import os


# assumes all fastq files for patient are locatred in same directory
# if this is not the case, change samplesheet.csv manually after running this script

# assumes fastq files are named {sample}_R1_001.fastq.gz, {sample}_R2_001.fastq.gz
# can overwrite by providing -r1_suffix and -r2_suffix options (e.g. '_R2_001.fastq.gz')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-patient_dir", required=True, help="path to patient subdirectory in Raw/ folder")
    parser.add_argument("-patient_id", required=True, help="patient identifier")
    parser.add_argument("-samples", required=True, help="list of sample ids (NOT including the normal sample)")
    parser.add_argument("-normal_sample", required=True, help="normal sample id")
    parser.add_argument("-fastq_dir", required=True, help="root directory where all fastq files are located")
    parser.add_argument("-r1_suffix", default=None, help="R1 fastq file suffix")
    parser.add_argument("-r2_suffix", default=None, help="R2 fastq file suffix")

    args = parser.parse_args()

    patient_dir = os.path.abspath(args.patient_dir)
    samples = args.samples.split(',')

    if args.r1_suffix is not None:
        if args.r2_suffix is None:
            print("Error: please provide BOTH -r1_suffix and -r2_suffix options")
            return 
    
    if args.r2_suffix is not None:
        if args.r1_suffix is None:
            print("Error: please provide BOTH -r1_suffix and -r2_suffix options")
            return 
    
    d = defaultdict(list)

    r1_suffix = args.r1_suffix if args.r1_suffix is not None else "_R1_001.fastq.gz"
    r2_suffix = args.r2_suffix if args.r2_suffix is not None else "_R2_001.fastq.gz"

    for sample in samples:
        d['patient'].append(args.patient_id)
        d['sample'].append(sample)
        d['fastq_1'].append(os.path.join(args.fastq_dir, f"{sample}{r1_suffix}"))
        d['fastq_2'].append(os.path.join(args.fastq_dir, f"{sample}{r2_suffix}"))
        #d['bam'].append(os.path.join(patient_dir, sample, f"{sample}.bam"))
        #d['bai'].append(os.path.join(patient_dir, sample, f"{sample}.bai"))
        d['status'].append(1)

    d['patient'].append(args.patient_id)
    d['sample'].append('Normal')
    d['fastq_1'].append(os.path.join(args.fastq_dir, f"{args.normal_sample}{r1_suffix}"))
    d['fastq_2'].append(os.path.join(args.fastq_dir, f"{args.normal_sample}{r2_suffix}"))
    d['status'].append(0)
    
    df = pd.DataFrame(d)

    df.to_csv(os.path.join(patient_dir, 'samplesheet.csv'), index=False, columns=['patient', 'sample','fastq_1', 'fastq_2', 'status'])  



if __name__ == "__main__":
    main()