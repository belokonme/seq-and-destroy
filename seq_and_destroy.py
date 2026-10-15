from pathlib import Path

from helpers.fastq import (gc_content, mean_quality, create_bounds, is_within_bounds, read_fastq, write_read)

from helpers.bio_files import (read_gbk, get_neighbouring_genes, write_to_fasta)


def filter_fastq(
                 input_fastq: str, 
                 output_fastq: str, 
                 gc_bounds: int | float | tuple = (0, 100), 
                 length_bounds: int | float | tuple = (0, 2 ** 32), 
                 quality_threshold: int | float = 0
                 ) -> None:
    """
    Filters reads in .fastq file by GC content, length, quality of read.

    Arguments:
    input_fastq: str (name of input .fastq file)
    output_fastq: str (name of output .fastq file, will be saved to "filtered/")
    gc_bounds: int / float / tuple (interval for GC content; a single number is the upper bound)
    length_bounds: int / float / tuple (interval for read length)
    quality_threshold: int / float (minimal quality of read)

    If read in .fastq passes filtration, writes it to file in "filtered/" folder.
    """
    gc_bounds = create_bounds(gc_bounds)
    length_bounds = create_bounds(length_bounds)
    Path("filtered").mkdir(exist_ok=True)
    with open(Path("filtered")/output_fastq, "w") as out:
        for name, seq, plus, quality in read_fastq(input_fastq):
            if is_within_bounds(gc_content(seq), gc_bounds) and is_within_bounds(len(seq), length_bounds) and mean_quality(quality) >= quality_threshold:
                write_read(out,(name, seq, plus, quality))


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Gathers BLAST best hits for genes of interest and writes them to file.

    Arguments:
    input_file: str (name of input file with BLAST results)
    output_file: str (name of output file with best hits)

    Written file is sorted by alphabet.
    """
    best_hits = []
    query = "Sequences producing significant alignments:"
    with open(input_file) as blast_output:
        for line in blast_output:
            if line.strip() == query:
                blast_output.readline()
                blast_output.readline()
                hit = blast_output.readline()[:66].strip()
                best_hits.append(hit)
    with open(output_file, "w") as out:
        best_hits.sort(key = str.lower)
        for hit in best_hits:
            out.write(hit + "\n")


def select_genes_from_gbk_to_fasta(
                                   input_gbk: str, 
                                   genes: str | list | tuple, 
                                   n_before: int = 1, 
                                   n_after: int = 1, 
                                   output_fasta: str = "output.fasta"
                                   ) -> None:
    """
    Selects neighbouring genes from .gbk file and writes their protein sequences to file in .fasta format.

    Arguments:
    input_gbk: str (path to .gbk input file)
    genes: str / list / tuple (genes of interest)
    n_before: int (number of genes before gene of interest)
    n_after: int (number of genes after gene of interest)
    output_fasta: str (name of output file in .fasta format)
    """
    if isinstance(genes, str):
        genes = [genes]
    cds_list = read_gbk(input_gbk)
    neighbours = get_neighbouring_genes(cds_list, genes, n_before, n_after)
    write_to_fasta(neighbours, output_fasta)