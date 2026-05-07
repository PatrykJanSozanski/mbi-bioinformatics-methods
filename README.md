# mbi-bioinformatics-methods

Small command-line bioinformatics utilities.

## Requirements

- Python 3
- Biopython (for `gc_content.py`)
- pandas and pyranges (for `count_variants_by_gene.py`)

## Tools

### 1) `gc_content.py`

Calculates GC-content from sequences in a FASTA file.

#### Usage

```bash
python gc_content.py path/to/sequences.fasta
```

#### Output

For each sequence:
- Sequence ID
- Sequence length
- GC-content (%)

Final summary:
- Number of sequences processed
- Total sequence length
- Weighted GC-content (%)

If the FASTA file contains no sequences, the script exits with:
`No sequences found in the provided FASTA file.`

### 2) `count_variants_by_gene.py`

Counts VCF variants that overlap gene intervals from a refFlat annotation file.

#### Usage

```bash
python count_variants_by_gene.py path/to/variants.vcf path/to/annotation.refflat
```

Optional arguments:
- `-o, --output PATH` — write tab-separated output to file
- `--include-zero` — include genes with zero overlapping variants

#### Output

Tab-separated table with columns:
- `Gene`
- `VariantCount`
