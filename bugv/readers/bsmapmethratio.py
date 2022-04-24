import pysam
import pandas as pd
import gzip
import numpy as np

class BsmapMethratioTbi():
    
    """
    使用完后，需用close方法关闭文件句柄
    
    load_region_data调用stat_bin_region:
    返回的是一个列表，每个元素是一个DataFrame,
    分别是CG, CHG, CHH的统计结果。
    """
    
    def __init__(self, file_name):
        self.filein = file_name
        self.is_open = False
    
    def load_region_data(self, region, bin_length=1, min_ct_count=1):
        self.open()
        try:
            self.region_data = self.stat_bin_region(region, bin_length, min_ct_count)
        finally:
            self.close()
    
    def check_format(self):
        """"
        检测是否为原始格式（列数大于等于8）或简化的格式（列数为6）
        设置属性is_origin_format
        原始格式：
        chr_name, pos, strand, context, ratio, eff_CT_count, C_count, CT_count
        unmethlyated C was converted to T
        not return eff_CT_count
        the eff_CT_count is new and is calculated by adjusting the
         methylation ratio for the C/T SNP, using the reverse strand
         mapping information.
        简化的格式：
        chr_name, pos, strand, context, C_count, CT_count
        """
        with gzip.open(self.filein, 'rt') as f:
            l = next(f)
            d = l.rstrip("\n").split("\t")
            if len(d) >= 8:
                self.is_origin_format = True
            else:
                self.is_origin_format = False
    
    def open(self):
        self.check_format()
        self.tbx = pysam.TabixFile(self.filein)
        self.is_open = True
    
    def close(self):
        try:
            self.tbx.close()
        except:
            pass
        self.is_open = False
                        
    def read_region(self, region, min_ct_count=1, to_df=False, based0=False):
        
        def _parse_tbx_lines(lines):
            is_origin_format = self.is_origin_format
            for l in lines:
                d = l.rstrip("\n").split("\t")
                if self.is_origin_format:
                    chr_name, pos, strand, context, ratio, eff_CT_count, C_count, CT_count = d[:8]
                else:
                    chr_name, pos, strand, context, C_count, CT_count = d
                pos = int(pos)
                C_count, CT_count = int(C_count), int(CT_count)
                yield([chr_name, pos, strand, context, C_count, CT_count])
        
        #only for tbi method
        chr_name, region_start, region_end = region[:3]
        if not based0:
            region_start -= 1
        if region_start < 0:
            region_start = 0
        records = self.tbx.fetch(chr_name, region_start, region_end)
        records = _parse_tbx_lines(records)
        if to_df:
            header = ["chr_name", "pos", "strand", "context", "c_count", "ct_count"]
            records = pd.DataFrame(list(records), columns=header)
        return(records)
    
    def stat_bin(self, df, region_start=None, bin_length=1):
        if region_start is None:
            region_start = df["pos"].min()
        if bin_length != 1:
            df["pos"] = np.floor_divide(df["pos"].values - region_start, bin_length)*bin_length + region_start
        df = df.groupby(["pos", "context"])[["c_count", "ct_count"]].sum().reset_index()
        df["ratio"] = df["c_count"]/df["ct_count"]
        df["bin_length"] = bin_length
        return df
    
    def stat_bin_region(self, region, bin_length=1, min_ct_count=1):
        df = self.read_region(region, min_ct_count=min_ct_count, to_df=True)
        region_start = region[1]
        df = self.stat_bin(df, region_start, bin_length)
        return [d for i, d in df.groupby("context")]
        
        
        
