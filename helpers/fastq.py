def gc_content(seq: str) -> float:
    """
    Computes GC content in a sequence.

    Arguments:
    seq: str, nucleotide sequence (case-insensitive)

    Returns float (percentage of GC content)
    If sequence is empty, returns 0.
    """
    if not seq:
        return 0
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100


def mean_quality(quality: str) -> float:
    """
    Computes mean quality of a sequence.
    
    Arguments:
    quality: str, quality of a sequence coded in phred33
    
    Returns float (mean phred33 quality score)
    If quality is empty, returns 0
    """
    if not quality:
        return 0
    total_quality = 0
    for symbol in quality:
        total_quality += ord(symbol) - 33
    return total_quality / len(quality)


def create_bounds(bounds: int | float | tuple) -> tuple:
    """
    Converts bounds inputted by user to tuple interval (lower, upper).

    Arguments:
    bounds: int / float / tuple

    Returns bounds
    If bounds are a single number, returns an interval (0, bounds) as tuple
    If bounds are a tuple, does nothing and returns the same tuple
    """
    if isinstance(bounds, (float, int)):
        return (0, bounds)
    else:
        return bounds


def is_within_bounds(value: int | float, bounds: tuple) -> bool:
    """
    Checks whether GC content of a sequence is within bounds set by user.
    Can also check if length of a sequence is within set bounds.

    Arguments:
    value: int / float
    bounds: tuple

    Returns True if value is in bounds
    Returns False if value is not in bounds
    """
    lower, upper = bounds
    if value <= upper and value >= lower:
        return True
    else:
        return False


def read_fastq(path: str):
    """
    Reads sequences in .fastq file one by one.

    Arguments:
    path: str (path to file)

    Gathers a single read with its info from file
    Yields read as tuple (name, seq, plus, quality)
    """
    with open(path) as file:
        while True:
            name = file.readline().strip()
            if name == "":
                break
            seq = file.readline().strip()
            plus = file.readline().strip()
            quality = file.readline().strip()
            yield (name, seq, plus, quality)


def write_read(file, read: tuple) -> None:
    """
    Writes reads to file one by one.

    Arguments:
    file: opened output file
    read: tuple[str, str, str, str]

    Writes a single read to file in .fastq format
    """
    name, seq, plus, quality = read
    file.write(name + "\n" + seq + "\n" + plus + "\n" + quality + "\n")