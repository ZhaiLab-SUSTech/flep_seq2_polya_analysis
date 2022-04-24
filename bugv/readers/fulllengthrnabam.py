import pandas as pd
import numpy as np
import pysam
from bugv.readers.rnabam import _RNABamOnlyReads, convert_reads_to_features
from bugv.readers.feature import Feature, FeatureArray

"""
version 0.2.1 by jiajinbu 2021.01.04
"""

class RNABamFeature(_RNABamOnlyReads):
        
    def _pre_process_feature(self, feature_array):
        pass
        
    def _load_region_data_func(self):
        reads = self.bam_obj.fetch(self.chr_name,
                       self.region_start - 1,
                       self.region_end)
        region_type_name = "mRNA"
        feature_array = convert_reads_to_features(reads, region_type_name, self._get_strand_func)
        self._pre_process_feature(feature_array)
        return feature_array

class FulllengthBam(RNABamFeature):
    
    """
    The tags extract from bam are stored in feature_array.df, but not in each read feature.
    """
    
    DEFAULT_CONFIG = {"min_polya_length": 15}
    
    def __init__(self, file_name="", add_polyA_feature_flag=True, min_polya_length=-1):
        self.file_name = file_name
        self.add_polyA_feature_flag = add_polyA_feature_flag
        if min_polya_length == -1:
            self.min_polya_length = self.DEFAULT_CONFIG["min_polya_length"]
    
    def _get_strand_func(self, read):
        read_strand = ""
        try:
            read_strand = read.get_tag("rs")
        except:
            pass
        if not read_strand:
            read_strand = "-" if read.is_reverse else "+"
        return read_strand
    
    def _pre_process_feature(self, feature_array):
        
        def get_and_set_tags(feature_array,
                     keys=["mi", "ty", "rs", "ir", "pa", "pn"],
                     column_names = ["mRNA_id", "read_type", "rna_strand", "intron_retention", "polya_length", "pos_not_consistent"],
                     default_values=["", "", "", "", 0, 1]):
        
            def _get_tag(read, key, default_value):
                try:
                    return read.get_tag(key)
                except:
                    return default_value
            
            list_of_column_data = []
            #don't use list_of_column_data = [[]] *  len(column_names)
            for i in range(len(column_names)):
                list_of_column_data.append([])
            for feature in feature_array.features:
                read = feature.parm["read"]
                for key, default_value, column_name, column_data in zip(keys, default_values, column_names, list_of_column_data):
                    value = _get_tag(read, key, default_value)
                    column_data.append(value)
                    feature.__setattr__(column_name, value)                        
                
            for column_name, column_data in zip(column_names, list_of_column_data):
                feature_array.df[column_name] = column_data
                
        def add_polyA_feature(region_data, min_polya_length=15):
            for feature in region_data:
                if feature.polya_length >= min_polya_length:
                    if feature.strand == "+":
                        feature.children.append(Feature(feature.end+1, 
                                                     feature.end + feature.polya_length, 
                                                     feature.chr_name,
                                                     feature.strand, "", 1, "polyA"))
                        feature.end = feature.end + feature.polya_length                 
                    else:
                        feature.children.insert(0, Feature(feature.start - feature.polya_length, 
                                                     feature.start - 1, 
                                                     feature.chr_name,
                                                     feature.strand, "", 1, "polyA"))
                        feature.start = feature.start - feature.polya_length            
            return region_data
            
        get_and_set_tags(feature_array)
        if self.add_polyA_feature_flag:
            add_polyA_feature(feature_array, self.min_polya_length)
        