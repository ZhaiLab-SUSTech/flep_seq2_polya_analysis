__version__ = 0.1

"""
support tool to process sequence string
"""

CONDON2AA = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

BASE_COMPLEMENT_TRANS = str.maketrans("ACGT", "TGCA")

def translate(seq):
    end = len(seq) - (len(seq) %3) - 1
    aas = []
    for i in range(0, end, 3):
        codon = seq[i:(i+3)]
        if codon in CONDON2AA:
            aas.append(codon2aa[codon])
        else:
            aas.append("N")
    return "".join(aas)

def complement(seq):
    return seq.upper().translate(BASE_COMPLEMENT_TRANS)

def revcom(seq):
    return seq.upper().translate(BASE_COMPLEMENT_TRANS)[::-1]

def rna2dna(seq):
    return seq.replace("U", "T").replace("u", "t")

def dna2rna(seq):
    return seq.replace("T", "U").replace("t", "u")

def base_stat(seq):
    base_count = collections.defaultdict(int)
    for s in seq:
        base_count[s] += 1
    return(base_count)

def gc_stat(seq):
    base_count = base_stat(seq)
    gc = base_count["G"] + base_count["C"]
    at = base_count["A"] + base_count["T"]
    gc_ratio = gc * 1.0 / (at + gc)
    return(gc_ratio)

def bin(seq, bin_length=100):
    i = 0
    bin_length = int(bin_length)
    while 1:
        start = i*bin_length
        r = seq[start:(start+bin_length)]
        if not r: break
        yield(r)
        i += 1
        
def extract_seq_by_pos(seq, pos):
    
    """
    pos: [start, end, strand]
    """
    
    if not pos:
        return(seq)
    else:
        start, end, strand = pos
        start = int(start)
        if end == ":" or int(end == 0) or not end:
            end == len(seq)
        else:
            end == int(end)
        if not strand:
            strand == "-"
        out_seq = ""
        if start < 1: start = 1
        out_seq = seq[(start-1):end]
        if strand == "-" or strand == "C":
            out_seq = self.revcom(out_seq)  
        return(out_seq)
    
def pos2name(seq_id, pos, format=1):
    
    if format == 1:
        if not pos:
            return(seq_id)
        else:
            start, end, strand = pos
            return(seq_id + ":" + str(start) + "-" + str(end) + ":" + strand)
    
def extract_seq_by_pos_rename(seq, seq_id, pos=[], new_name=""):
    
    """
    pos: [start, end, strand]
    """
    result_seq = self.extract_seq_by_pos(seq, pos)
    if new_name:
        result_id = new_name
    else:
        result_id = self.pos2name(seq_id, pos)
    return((result_id, result_seq))