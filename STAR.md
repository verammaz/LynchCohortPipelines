# STAR (Spliced Transcripts Alignment to a Reference)

> Source: https://github.com/alexdobin/STAR 

## Running on Minerva:

### 1: GTF annotation file

> GTF annotation files can be found here (check consistent genome build): https://www.gencodegenes.org/human/releases.html

Before running STAR genome indexing or alignment, ensure that the chromosome names in your reference FASTA and GTF file match. A mismatch (e.g., chr1 in GTF vs. 1 in FASTA) will cause STAR to exit with a fatal error.

#### Automated compatibility check (included in script)
The pipeline includes a compatibility check that compares chromosome naming conventions in both files. It will:

Extract the first chromosome name from the GTF and FASTA

Warn you if one uses chr prefixes while the other does not

Exit early with a helpful message if a mismatch is detected

#### Manual check (optional)
To manually inspect the chromosome naming formats:

```bash
# First chromosome from the FASTA
zgrep '^>' /path/to/reference.fasta | head -n 1

# First chromosome from the GTF (should be an exon line)
zgrep -v '^#' /path/to/annotation.gtf | awk '$3 == "exon" {print $1}' | head -n 1
```

#### Fixing common mismatches
If your GTF uses chr1 but your FASTA uses 1, strip the chr prefix from the GTF:

```bash
sed 's/^chr//' annotation.gtf > annotation.nochr.gtf
```
Then use annotation.nochr.gtf for your STAR command.



### 2: Set up config

Set the following file and directory paths in your `config.sh` file:

```bash
# specify file path for COMPATIBLE gtf annotations (required for STAR indexing)
GTF_FILE=""

# specify directory for STAR indexing files
STAR_DIR=""
```

> Important: .gft annotation file must be unzipped


### 3: Submit job to LSF

```bash
cd LynchCohortPipelines
./submit_star.sh -r <reads1.fastq.gz, reads2.fastq.gz> -o <output_prefix>
```

