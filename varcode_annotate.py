import os
import argparse
from pyensembl import ensembl_grch37

from varcode import load_vcf
from varcode.effects import Insertion, Deletion

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-patient', required=True)
    parser.add_argument('-hdir', required=True)

    args = parser.parse_args()

    for root, dirs, files in os.walk(os.path.join(args.hdir, 'VCF', args.patient)):
        print(root, dirs)
        for file in files:
            print(file)
            if file.startswith('.'): continue
            if 'ann' not in file and 'template' not in file and file.endswith('.vcf'):
                print(file)
                #lesion = os.path.basename(file).split('.')[0]

                vcfVariants = load_vcf(os.path.join(root,file), genome=ensembl_grch37)
                vcfEffects = vcfVariants.effects()
                nonSilentMutations = vcfEffects.drop_silent_and_noncoding()
                nonSilentMutations.top_priority_effect_per_gene_id()

                nonSilentMutations.filter_by_effect_priority(Deletion).top_priority_effect_per_gene_id()

                print(nonSilentMutations)