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

def create_bounds(bounds):
    if isinstance(bounds, (float, int)):
        return (0, bounds)
    else:
        return bounds

def is_within_bounds(value, bounds):
    lower, upper = bounds
    if value <= upper and value >= lower:
        return True
    else: 
        return False

def read_fastq(path):
    reads = []
    with open(path) as file:
        while True:
            name = file.readline().strip()
            if name == "":
                break
            seq = file.readline().strip()
            plus = file.readline().strip()
            quality = file.readline().strip()
            reads.append((name, seq, plus, quality))
    return reads