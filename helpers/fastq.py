def gc_content(seq):
    if not seq:
        return "Nothing there..."
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100
