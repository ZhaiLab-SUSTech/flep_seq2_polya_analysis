import pandas as pd
import numpy as np
import pysam
from bugv.readers.feature import FeatureArray, Feature
import time #delete after test

"""
version 0.2.1 by jiajinbu 2021.01.04
"""

def _generate_empty_depth():
    return pd.DataFrame([], columns=["pos", "depth"])
    
def _generate_empty_junc():
    return pd.DataFrame([], columns=["start", "end", "count", "left_depth", "right_depth", "ratio"])
    
def _generate_empty_depth_junc():
    return [_generate_empty_depth(), _generate_empty_junc()]
    
def _generate_empty_depth_juncs():
    return [_generate_empty_depth_junc(), _generate_empty_depth_junc(), _generate_empty_depth_junc()]

def get_read_strand(read, mRNA_strand=True, is_FR=False):
    
    """
    is_FR is strand-specific. True: read1 is same strand to mRNA.
    For most of the strand-specific RNA-seq, is_FR is False.
    
    Return:
    read_is_plus #relative to mRNA
    You can set mRNA_strand if you only want to read strand information.
    """
    read_is_plus = not read.is_reverse
    #if read2, convert to read1
    if not read.is_read1:
        read_is_plus = not read_is_plus
    #if RF, convert
    if not is_FR:
        read_is_plus = not read_is_plus
    #if minus, convert
    if not mRNA_strand:
        read_is_plus = not read_is_plus
    return read_is_plus

def read_coverage(bam_obj, chr_name, start, end, 
                  mRNA_strand=True, 
                  is_FR=False, 
                  filter_unique=False, 
                  filter_func=None,
                  read_specific_strand=True):
    
    def _read_coverage(filter_strand=False, return_sense_strand=True):
        
        def _generate_read_func():
            #need: filter_strand return_sense_strand
            def _read_call_back_func(read):
                result = True
                if filter_unique and read.get_tag("NH") > 1:
                    result = False
                if result and filter_strand:
                    read_is_sense = get_read_strand(read, mRNA_strand, is_FR)
                    if read_is_sense != return_sense_strand:
                        result = False
                if result and filter_func is not None:
                    result = filter_func(read)
                return result
                
            return(_read_call_back_func)
        
        depth = bam_obj.count_coverage(chr_name, start-1, end, 
                                read_callback=_generate_read_func())
        
        depth = np.array(depth).max(axis=0).astype(int)
        pos = np.arange(start, end+1) #1-based
        return pd.DataFrame({"depth": depth, "pos": pos}, columns=["pos", "depth"])
             
    all_depth = _read_coverage(filter_strand=False)
    if read_specific_strand:
        sense_depth = _read_coverage(filter_strand=True, return_sense_strand=True)
        anti_depth = _read_coverage(filter_strand=True, return_sense_strand=False)
    else:
        sense_depth = _generate_empty_depth()
        anti_depth = _generate_empty_depth()
    return [all_depth, sense_depth, anti_depth]

def get_read_blocks(read):
    
    """
    1-based
    
    don't use read.get_blocks
    insert also seprate block
    
    bam_obj.find_introns !!!
    
    输入pysam的read比对对象
    输出该read比对由于splicing junction拆分成的block。
    过滤掉小的删除、插入和错配等。
    按染色质位置从小到大排列。
    是一个列表，每一个元素是一个起始位置和终止位置（均是1-based）组成的列表。
    """
    
    blocks = []
    start = read.reference_start + 1
    end = start - 1
    for (flag, length) in read.cigartuples:
        if flag == 4 or flag == 5 or flag == 1: continue
        if flag == 0 or flag == 2:
            end += length
        if flag == 3:
            blocks.append([start, end])
            start = end + length + 1
            end = start - 1
    blocks.append([start, end])
    return blocks

def get_introns(read, min_overlap=0, get_junc=False, return_exon_flag=False):
    
    """
    1-based 
    [[start, end]]
    the start and end position of intron
    """
    
    if read.is_unmapped:
        if return_exon_flag:
            return [[[]], [[]]]
        else:
            return [[]]

    blocks = get_read_blocks(read)
    introns = []
    for i, (start, end) in enumerate(blocks[:-1]):
        if min_overlap:
            if i == 0 and end - start + 1 < min_overlap:
                continue
            if i == len(blocks) - 2:
                s, e = blocks[-1]
                if e - s + 1 < min_overlap:
                    continue
        if get_junc:
            introns.append([end, blocks[i+1][0]])
        else:
            introns.append([end+1, blocks[i+1][0]-1])
    if return_exon_flag:
        return [blocks, introns]
    else:
        return introns

def get_exon_introns(read):
    
    exon_introns = []
    exons, introns = get_introns(read, min_overlap=0, get_junc=False, return_exon_flag=True)
    for exon, intron in zip(exons[:-1], introns):
        exon_introns.append([exon, "exon"])
        exon_introns.append([intron, "intron"])
    exon_introns.append([exons[-1], "exon"])
    return exon_introns

def get_reads_introns(reads, min_overlap=0, get_junc=False):
    """
    return pd.DataFrame(introns, columns=["start", "end", "count"])
    if get_junc, start = start - 1, end = end + 1
    position: 1-based
    
    #tuple can be hashed
    #bam_obj.find_introns([read for read in bam_obj.fetch("chr1", 3852, 4100) if not read.is_unmapped])
    #some read are unmapped. Must filter.
    """
    from collections import Counter
    
    introns = Counter()
    for read in reads:
        if read.is_unmapped:
            continue
        a = get_introns(read, min_overlap, get_junc)
        for start, end in a:
            introns[(start, end)] += 1
    introns = sorted(introns.items(), key=lambda x: x[0][0])
    introns = [[s, e, v]for (s, e), v in introns]
    introns = pd.DataFrame(introns, columns=["start", "end", "count"])
    return introns
    
def split_read_by_strand(bam_obj, chr_name, start, end, 
                            mRNA_strand=True, 
                            is_FR=False, 
                            filter_unique=False, 
                            filter_func=None):
    
    """
    return [reads, sense_reads, anti_reads]
    reads list.
    sense_reads and anti_reads share elements with reads.
    """
    
    reads = list(bam_obj.fetch(chr_name, start-1, end))
    sense_reads = []
    anti_reads = []
    for read in reads:
        if filter_unique and read.get_tag("NH") > 1:
            continue
        if filter_func and filter_func(read):
            continue
        if get_read_strand(read, mRNA_strand, is_FR):
            sense_reads.append(read)
        else:
            anti_reads.append(read)
    return([reads, sense_reads, anti_reads])
                
def get_depth_junc(bam_obj, chr_name, start, end, 
                  mRNA_strand=True, 
                  is_FR=False, 
                  filter_unique=False, 
                  filter_func=None,
                  read_specific_strand=True,
                  min_overlap=0):
    
    def add_depth_to_junc(junctions, depths):
        junctions["left_depth"] = pd.merge(junctions[["start"]], depths, left_on="start", right_on="pos", how="left")["depth"]
        junctions["right_depth"] = pd.merge(junctions[["end"]], depths, left_on="end", right_on="pos", how="left")["depth"]
        junctions["ratio"] = junctions["count"]/(junctions["left_depth"] + junctions["right_depth"])*2
        return junctions
    
    if mRNA_strand == "+" or mRNA_strand == "-":
        mRNA_strand = mRNA_strand == "+"
    
    reads, sense_reads, anti_reads = split_read_by_strand(bam_obj,
        chr_name, start, end, mRNA_strand, is_FR, filter_unique, filter_func)
    all_junctions = get_reads_introns(reads, min_overlap, get_junc=True)
    
    x_left = start
    x_right = end
    s = all_junctions["start"]
    if len(s) > 0:
        x_left = min(x_left, s.min())
        x_right = max(x_right, s.max())
        
    if read_specific_strand:
        if read_specific_strand:
            sense_junctions = get_reads_introns(sense_reads, min_overlap, get_junc=True)
            anti_junctions = get_reads_introns(anti_reads, min_overlap, get_junc=True)
        s = sense_junctions["start"]
        if len(s) > 0:
            x_left = min(x_left, s.min())
            x_right = max(x_right, s.max())
        s = anti_junctions["start"]
        if len(s) > 0:
            x_left = min(x_left, s.min())
            x_right = max(x_right, s.max())
    
    depth_data = read_coverage(bam_obj, chr_name, x_left, x_right, 
              mRNA_strand, 
              is_FR, 
              filter_unique, 
              filter_func,
              read_specific_strand)
              
    all_junctions = add_depth_to_junc(all_junctions, depth_data[0])
    if read_specific_strand:
        sense_junctions = add_depth_to_junc(sense_junctions, depth_data[1])
        anti_junctions = add_depth_to_junc(anti_junctions, depth_data[2])
    else:
        sense_junctions = _generate_empty_junc()
        anti_junctions = _generate_empty_junc()
    
    depth_data = [d[(d["pos"]>=start) & (d["pos"]<=end)] for d in depth_data]
    junc_data = [all_junctions, sense_junctions, anti_junctions]
    return list(zip(depth_data, junc_data))

def convert_read_to_feature(read, region_type_name="mRNA", get_strand_func=None):
    
    def default_strand_func(read):
         strand = "-" if read.is_reverse else "+"
         return strand
         
    if not get_strand_func:
        get_strand_func = default_strand_func
    
    region_start = read.reference_start + 1
    region_end = read.reference_end
    chr_name = read.reference_name
    region_strand = get_strand_func(read)
    region_id = ",".join([read.query_name,
                          chr_name, 
                          str(region_start), 
                          str(region_end),
                          region_strand,
                          str(read.query_alignment_start+1)])
    children = []
    exon_introns = get_exon_introns(read)
    for (block_start, block_end), block_type in  exon_introns:
        children.append(Feature(block_start, block_end, chr_name, region_strand, "", 1, block_type))
    feature = Feature(region_start, region_end, chr_name, region_strand, region_id, 
                        1, region_type_name, parm={"read": read}, children=children)
    return feature

def convert_reads_to_features(reads, region_type_name="mRNA", get_strand_func=None, features_type=None):
    
    features = []
    for read in reads:
        if not read.is_unmapped:
            features.append(convert_read_to_feature(read, region_type_name, get_strand_func))
    if not features_type:
        features_type = FeatureArray
    return features_type(features)    
        
class RNABam:
    
    def __init__(self, file_name="", is_FR=False, 
                  filter_unique=False, 
                  filter_func=None,
                  read_specific_strand=True,
                  min_overlap=0):
        
        """
        is_FR to set strand-specific:
        For most of strand-specific, is_FR is True.
        You can use True for non-strand-specific library. 
        """
        self.file_name = file_name
        self.is_FR = is_FR
        self.filter_unique = filter_unique
        self.filter_func = filter_func
        self.read_specific_strand = read_specific_strand
        self.min_overlap = min_overlap
        self.region_data = _generate_empty_depth_juncs()
    
    def _get_strand_func(self, read):
        if get_read_strand(read, True, self.is_FR):
            strand = "+"
        else:
            strand = "-"
        return strand
    
    def read_coverage(self, bam_obj, chr_name, start, end, 
                  mRNA_strand=True, 
                  filter_unique=False, 
                  filter_func=None,
                  read_specific_strand=True):
        return read_coverage(bam_obj, chr_name, start, end, 
                  mRNA_strand, self.is_FR, filter_unique, filter_func,
                  read_specific_strand)
                      
    def get_depth_junc(self, bam_obj, chr_name, start, end, 
                  mRNA_strand=True, 
                  is_FR=False, 
                  filter_unique=False, 
                  filter_func=None,
                  read_specific_strand=True,
                  min_overlap=0):
        return get_depth_junc(bam_obj, chr_name, start, end, 
                  mRNA_strand, 
                  is_FR, 
                  filter_unique, 
                  filter_func,
                  read_specific_strand,
                  min_overlap)
    
    def _load_region_data_func(self):
        return self.get_depth_junc(self.bam_obj, 
                            self.chr_name,
                            self.region_start,
                            self.region_end,
                            self.region_strand,
                            self.is_FR,
                            self.filter_unique, 
                            self.filter_func,
                            self.read_specific_strand,
                            self.min_overlap)
    
    def load_region_data(self, region, select_features=None):
        
        """
        region: [chr_name, region_start, region_end, region_strand]
        region_strand is like mRNA strand.
        
        return ([[[total_depth], [total_junction]],
                 [[sense_depth], [sense_junction]],
                 [[anti_depth],  [anti_junction]]])
        sense or anti_sense is relative to mRNA (region_strand).
        
        if read_specific_strand == False:
            [[sense_depth], [sense_junction]] is [[], []]
            [[anti_depth],  [anti_junction]]] is [[], []]
        """
        #print(time.time(), "open", self.file_name)
        
        chr_name, region_start, region_end, region_strand = region
        self.chr_name = chr_name
        self.region_start = region_start
        self.region_end = region_end
        self.region_strand = region_strand
        
        self.open()
        try:
            region_data = self._load_region_data_func()
            self.region_data = region_data
        finally:
            self.close()
        #print(time.time, "close", self.file_name)
        return(region_data)
        
    def open(self):
        self.close()
        self.bam_obj = pysam.AlignmentFile(self.file_name, "rb")
        
    def close(self):
        try:
            self.bam_obj.close()
        except:
            pass
            
class _RNABamOnlyReads(RNABam):
        
    def _load_region_data_func(self):
        
        reads = list(self.bam_obj.fetch(self.chr_name,
                       self.region_start - 1,
                       self.region_end))
        #after read all reads, you can close bam_obj
        return reads
