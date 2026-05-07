# gc-content

Small command-line tool to calculate GC-content from sequences in a FASTA file.

## Requirements

- Python 3
- Biopython

## Usage

```bash
python gc_content.py path/to/sequences.fasta
```

## What it prints

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
