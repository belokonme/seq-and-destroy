# seq-and-destroy

Bioinformatic utilities for FASTQ filtering, BLAST output parsing and GenBank neighbouring gene extraction.

Lovingly named after Metallica's 1983 song Seek and Destroy.
The name is also a reference to the author's destroyed sleep schedule.

# Requirements

Requires Python 3.10+ for proper functioning.

# Installation

Simply `git clone` this repository to use it on your machine:

- SSH:

`git clone git@github.com:belokonme/seq-and-destroy.git`

- HTTPS:

`git clone https://github.com/belokonme/seq-and-destroy.git`

# Repository structure

```python
-/
|- README.md
|- seq_and_destroy.py   # main utilities
|- example_data/        # folder with example files to try out
    |- example_blast_results.txt
    |- example_gbk.gbk
    |- example_fastq.fastq
|- helpers/             # folder with modules imported to main script
    |- fastq.py
    |- bio_files.py
|- imgs/
    |- my_honest_reaction.png 
```

# Usage

Import the functions for use:

```python
from seq_and_destroy import (filter_fastq, parse_blast_output, select_genes_from_gbk_to_fasta)
```

Launch the utilities in the root of this repository.

## filter_fastq()

Function `filter_fastq()` filters reads in FASTQ-file by GC content, length, and quality of read. It accepts 5 arguments:

1. **input_fastq** (_str_) 

Name of your FASTQ-file you want to filter.

2. **output_fastq** (_str_)

Name of your output FASTQ-file with filtered reads.

3. **gc_bounds** (_int | float | tuple_)

A single upper limit of GC content percentage *or* a precise interval you want to set. Default value for this argument is `gc_bounds = (0, 100)`.

4. **length_bounds** (_int | float | tuple_)

Length of a read (upper limit or interval). Default value is `length_bounds = (0, 2 ** 32)`.

5. **quality_threshold** (_int | float_)

Minimal quality of a read. Default value is `quality_threshold = 0`.

The results are stored in "filtered/" folder.

### Examples

1. Default settings:

```python
filter_fastq("example_data/example_fastq.fastq", "all_reads.fastq")
```

Writes all reads in the file (in this case, 89).

2. With parameteres:

GC content: from 40 to 60%, maximum length: 50, mean quality: 30.

```python
filter_fastq("example_data/example_fastq.fastq", "filtered_reads.fastq",
               gc_bounds = (40, 60), length_bounds = 50, quality_threshold = 30)
```

Writes only reads that satisfy the parameteres (in this case, 27).

## parse_blast_output()

Function parse_blast_output() gathers BLAST best hits for genes of interest and writes them to file. It accepts 2 arguments:

1. **input_file** (_str_)

Path to your input file (output after BLAST search).

2. **output_file** (_str_)

Name of your output file with best hits.

### Examples

```
parse_blast_output("example_data/example_blast_results.txt", "best_hits.txt")
```

Writes 37 best hits in alphabetical order to file `best_hits.txt`

## select_genes_from_gbk_to_fasta()

Function select_genes_from_gbk_to_fasta() selects neighbouring genes from .gbk file and writes their respective protein sequences to file in FASTA format. Genes of interest are not written to the output; shared neighbours are written only once. 

Function accepts 5 arguments:

1. **input_gbk** (_str_)

Path to your .gbk input file.

2. **genes** (_str | list | tuple_)

Your gene or genes of interest that you want to find neighbours of.

3. **n_before** (_int_)

Number of genes that come before the gene(s) of interest in the file. Default value is `1`.

4. **n_after** (_int_)

Number of genes that come after the gene(s) of interest in the file. Default value is `1`.

5. **output_fasta** (_str_)

Name of your output file with selected neighbouring genes and their sequences. Default is `"output.fasta"`.

### Examples

1. A single gene with two neighbours (one before and one after):

```python
select_genes_from_gbk_to_fasta("example_data/example_gbk.gbk", "phrB",
                                output_fasta = "phrB_neighbours.fasta")
```

Will write _two_ neighbours of the phrB gene to file.

2. Several genes and two neighbours before and two after:

```python
select_genes_from_gbk_to_fasta("example_data/example_gbk.gbk",
                              ["phrB", "dtpD"], n_before = 2, n_after = 2,
                                output_fasta = "neighbours.fasta")
```

Will write _three_ neighbours in total for genes phrB and dtpD to file. 

# Personal note

![](imgs/my_honest_reaction.png)
