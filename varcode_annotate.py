import os
import argparse
from pyensembl import ensembl_grch37

from varcode import Variant
from varcode import load_vcf
from varcode.effects import Insertion, Deletion

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-patient', required=True, help='patient id')
    parser.add_argument('-hdir', required=True, help='directory with VCF subdirectory')

    args = parser.parse_args()

    for root, dirs, files in os.walk(os.path.join(args.hdir, 'VCF', args.patient)):
       
        for file in files:

            if file.startswith('.'): continue 
            if file.endswith('.vcf'):
                if 'varcode' in file or 'ann' in file or 'template' in file: continue
                
                lesion = os.path.basename(file).split('.')[0]

                vcfVariants = load_vcf(os.path.join(root,file), genome=ensembl_grch37)
                vcfEffects = vcfVariants.effects()
                nonSilentMutations = vcfEffects.drop_silent_and_noncoding()
                nonSilentMutations.top_priority_effect_per_gene_id()
                nonSilentMutations.filter_by_effect_priority(Insertion)
                variant_to_effect = nonSilentMutations.groupby_variant()
            
                with open(os.path.join(root, file), 'r') as infile, open(os.path.join(root, lesion+"_varcode.vcf"), 'w') as outfile:
                    for line in infile.readlines():
                        if line.startswith('#'):
                            outfile.write(line)
                            continue
                        vcf_fields = line.split('\t')
                        variant = Variant(contig=vcf_fields[0], start=vcf_fields[1], ref=vcf_fields[3], alt=vcf_fields[4], ensembl=ensembl_grch37)
                        try:
                            effect = variant_to_effect[variant].detailed_string()
                        except KeyError:
                            effect = "silent or noncoding"
                        vcf_fields[7] = ('|').join(effect.split('\n'))
                        outfile.write(('\t').join(vcf_fields))
                
