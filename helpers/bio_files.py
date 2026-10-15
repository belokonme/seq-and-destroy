def read_gbk(path):
    cds_list = []
    cds_info = None
    with open(path) as gbk:
        for line in gbk:
            if line.startswith("     CDS"):
                cds_info = {"gene":None, "locus_tag":None, "translation": ""}
                cds_list.append(cds_info)
            elif line.startswith("                     /gene") and cds_info is not None:
                cds_info["gene"] = line.strip().split('"')[1]
            elif line.startswith("                     /locus_tag") and cds_info is not None:
                cds_info["locus_tag"] = line.strip().split('"')[1]
            elif line.startswith("     ") and line[5] != " ":
                cds_info = None
            elif line.startswith("                     /translation") and cds_info is not None:
                seq = line.strip()
                while not seq.endswith('"'):
                    seq += gbk.readline().strip()
                cds_info["translation"] = seq.replace("/translation=", '').replace('"', '')
        return cds_list


def get_neighbouring_genes(cds_list, genes, n_before, n_after):
    neighbouring_genes = []
    genes_of_interest_indices = []
    for i, cds in enumerate(cds_list):
        if cds['gene'] in genes or cds['locus_tag'] in genes:
            genes_of_interest_indices.append(i)
    for index in genes_of_interest_indices:
        for j in range(index - n_before, index + n_after + 1):
            if j >= 0 and j < len(cds_list) and j not in genes_of_interest_indices and cds_list[j] not in neighbouring_genes:
                neighbouring_genes.append(cds_list[j])
    return neighbouring_genes


def write_to_fasta(neighbouring_genes, output_fasta):
    with open(output_fasta, "w") as out:
        for cds in neighbouring_genes:
            if cds['gene'] is not None:
                name = cds['gene']
            else:
                name = cds['locus_tag']
            out.write('>' + name + '\n')
            out.write(cds['translation'] + '\n')