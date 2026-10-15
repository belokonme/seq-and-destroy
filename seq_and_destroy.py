from pathlib import Path

from helpers.fastq import (gc_content, mean_quality, create_bounds, is_within_bounds, read_fastq, write_read)

from helpers.bio_files import (read_gbk, get_neighbouring_genes, write_to_fasta)

def filter_fastq(input_fastq, output_fastq, gc_bounds = (0, 100), length_bounds = (0, 2 ** 32), quality_threshold = 0):
    gc_bounds = create_bounds(gc_bounds)
    length_bounds = create_bounds(length_bounds)
    Path("filtered").mkdir(exist_ok=True)
    with open(Path("filtered")/output_fastq, "w") as out:
        for name, seq, plus, quality in read_fastq(input_fastq):
            if is_within_bounds(gc_content(seq), gc_bounds) and is_within_bounds(len(seq), length_bounds) and mean_quality(quality) >= quality_threshold:
                write_read(out,(name, seq, plus, quality))

def parse_blast_output(input_file, output_file):
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

def select_genes_from_gbk_to_fasta(input_gbk, genes, n_before = 1, n_after = 1, output_fasta = "output.fasta"):
    if isinstance(genes, str):
        genes = [genes]
    cds_list = read_gbk(input_gbk)
    neighbours = get_neighbouring_genes(cds_list, genes, n_before, n_after)
    write_to_fasta(neighbours, output_fasta)