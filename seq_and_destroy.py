from helpers.fastq import (gc_content, mean_quality, create_bounds, is_within_bounds, read_fastq, write_fastq)

def filter_fastq(input_fastq, output_fastq, gc_bounds = (0, 100), length_bounds = (0, 2 ** 32), quality_threshold = 0):
    filtered_reads = []
    gc_bounds = create_bounds(gc_bounds)
    length_bounds = create_bounds(length_bounds)
    for name, seq, plus, quality in read_fastq(input_fastq):
        if is_within_bounds(gc_content(seq), gc_bounds) and is_within_bounds(len(seq), length_bounds) and mean_quality(quality) >= quality_threshold:
            filtered_reads.append((name, seq, plus, quality))
    write_fastq(filtered_reads, output_fastq)