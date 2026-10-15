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