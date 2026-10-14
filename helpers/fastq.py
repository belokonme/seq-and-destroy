def gc_content(seq):
    if not seq:
        return 0
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def mean_quality(quality):
    if not quality:
        return 0
    total_quality = 0
    for symbol in quality:
        total_quality += ord(symbol) - 33
    return total_quality / len(quality)