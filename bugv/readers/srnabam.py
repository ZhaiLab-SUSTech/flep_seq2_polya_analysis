from bugv.baseseq_tool import revcom
import pandas as pd
import pysam

class sRNABam():
        
    def __init__(self, file_name, fileformat="shortstack"):
        self.file_name = file_name
        self.region_data = pd.DataFrame([], columns=["seq", "start", "end", "strand", "length", "is_unique", "map_time", "counts", "ir_rel_plus"])
        self.fileformat = fileformat 
        #fileformat: shortstack, seq_count, sim_seq_count
        #暂时不区分shortstack和seq_count
        self.open()
    
    def load_region_data(self, region, select_features=None):
        #select_features无用
        fileformat = self.fileformat
        
        chr_name, region_start, region_end, region_strand = region
        self.chr_name = chr_name
        self.region_start = region_start
        self.region_end = region_end
        self.region_strand = region_strand
        origin_reads = self.bam_obj.fetch(contig=chr_name, start=region_start-1, end=region_end)
        
        reads = []
        for read in origin_reads:
            chr_ = read.reference_name
            start = read.reference_start + 1
            strand = "-" if read.is_reverse else "+"
            end = read.reference_end
            length = end - start + 1
            map_time = read.get_tag("XM") - 1 #I'm not sure it's right for all possible bowtie option.
            is_unique = map_time == 1
            seq = read.query_sequence
            if read.is_reverse: seq = revcom(seq)
            if fileformat == "sim_seq_count":
                count = int(read.query_name.split("_")[1])
            else:
                count = 1
            reads.append([seq, start, end, strand, length, is_unique, map_time, count])
        if reads:
            reads = pd.DataFrame(reads, columns=["seq", "start", "end", "strand", "length", "is_unique", "map_time", "counts"])
            reads = reads.groupby(reads.columns.tolist()[:-1])["counts"].sum().reset_index(name="counts")
            #reads = reads.sort_values(["start", "end", "strand"])
            reads["is_rel_plus"] = reads["strand"] == region_strand
        else:
            reads = pd.DataFrame(reads, columns=["seq", "start", "end", "strand", "length", "is_unique", "map_time", "counts", "is_rel_plus"])
        self.region_data = reads
        return(reads)
    
    def open(self):
        self.close()
        self.bam_obj = pysam.AlignmentFile(self.file_name, "rb")
        
    def close(self):
        try:
            self.bam_obj.close()
        except:
            pass