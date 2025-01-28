# 1. Shared Frameshifts

## Goal:
- Get shared frameshift variants across lesions.
- Check for presence of specific frameshift variants in lesions.

## The Script:
Run the script as
```bash
python shared_frameshifts.py [OPTIONS] -hdir <home_dir> -lesions <lesions_list> -patients <patients_list>
```

### Required arguments:
```
-hdir                  Path to home directory with /VCF subfolder.
-lesions               Comma-separated list of lesion ids (no spaces).
-patients              Comma-separated list of patient ids, corresponding to lesion ids in <lesions> list (no spaces). 
```

> Note: The script assumes that patients[i] = lesion[i] patient. Ordering matters!


### Optional arguments:

| Parameter                 | Description   |	
| :----------------------------------------: | :------: |
| `-h` | Display usage message. |
| `-fs_annotation`| Flag to count frameshifts only through annotation (i.e. not just looking if indel is multiple of three).
| `-specific_fs_xl` | Excel file with specific frameshift variants and their coordinates in a dataframe. |
| `-specific_fs_cols` | Column names in `-specific_fs_xl` dataframe corresponding to RefChromId, RefCoord, GeneId (comma-separated, no spaces).  |
| `-check_raw` | Flag to check for variants in the lesion raw (mutect/strelka) vcf file before counting them.
| `--single_sample_trees` | Flag to add single sample tree config infos to CFIT config files. |

> Note: The XL file specified for the `specific_fs_xl` option must have at least the following columns: RefChromId (chromosome id in reference genome), RefCoord (start coordinate of variant), GeneId (name of gene). The exact names of these columns can be specified with the `-specific_fs_cols` option, in the described order (otherwise, default is 'HG19_id', 'COORD_HG19', 'GENE'). 

The script will output the following files:

*shared_frameshifts.xlsx* : This will be a dynamic file, with new sheets within the XL workbook created when the script is called with different combinations of lesions. Each sheet will correspond to a set of lesions and will contain a table of frameshift variants found in each of the lesions. The total number of shared variants and the variant NMD status is also reported.

*specific_fs.xlsx* : If the `-specific_fs_xl` option is specified, this file will also be dynamically updated on new calls of the script. It will contain a dataframe with rows corresponding variants in specified `-specific_fs_xl` and  columns 
[chrom_pos, gene_id, lesion1, lesion2, ...,	lesionN, total_lesions].

