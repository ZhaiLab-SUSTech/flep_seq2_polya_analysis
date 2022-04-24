import pandas as pd
import numpy as np
import sys

"""
version 0.2.1 by jiajinbu 2021.01.04
"""

DEFAULT_FEATURE_HEADER = ["chr_name", "start", "end", "strand", "height", "id", "type"]

def exon2bed(filein, fileout, input_0_based=False):
    
    """
    It will split the first column with "_" and only include the first element.
    For example: subtilisin_exon1 to subtilisin
    
    input:
    subtilisin	chr8	8531232	8531468	+
    subtilisin	chr8	8531556	8531665	+
    subtilisin	chr8	8531754	8531851	+
    subtilisin	chr8	8531956	8532120	+
    subtilisin	chr8	8541200	8541535	+
    CHS3	chr8	8538936	8539113	+
    CHS3	chr8	8539236	8540224	+
    CHS1	chr8	8535381	8536369	-
    CHS1	chr8	8536492	8536669	-
    """
    SPLIT_MRNA = "_"
    
    def convert_mRNA_to_one_line(mRNA_df):
        mRNA = mRNA_df["mRNA"].values[0]
        chr_name = mRNA_df["chr"].values[0]
        strand = mRNA_df["strand"].values[0]
        if input_0_based:
            block_starts = mRNA_df["start"].values 
        else:
            block_starts = mRNA_df["start"].values - 1
        block_ends = mRNA_df["end"].values
        mRNA_start = block_starts[0]
        mRNA_end = block_ends[-1]
        block_sizes = block_ends - block_starts
        block_starts = block_starts - mRNA_start
    
        return [chr_name, mRNA_start, mRNA_end, mRNA, 0, strand, 
                    mRNA_start, mRNA_end, 0, len(block_starts),
                    ",".join(block_sizes.astype(str)) + ",",
                    ",".join(block_starts.astype(str)) + ","]
                    
    exon_df = pd.read_table(filein, 
                     names=["mRNA", "chr", "start", "end", "strand"])
    exon_df["mRNA"] = [s.split(SPLIT_MRNA)[0] for s in exon_df.mRNA.values]
    mRNA_df = exon_df.groupby("mRNA").apply(convert_mRNA_to_one_line)
    mRNA_df = pd.DataFrame(mRNA_df.values.tolist())
    mRNA_df.to_csv(fileout, sep="\t", header=False, index=False)

class Feature:
    
    """
    Feature attribution:
        chr_name
        start
        end
        strand
        height
        id
        type
        parm (a dict)
        children (a list of Feature objects)
    Method:
        to_list:
            return a list of features. The first element is self. If it has children, the remain elements
            are childrens. The element is also a list, default is :
            [chr_name, start, end, strand, height, id, type, gene_id, mRNA_id]
            gene_id, mRNA_id default is "."
        to_df:
            the result is dataframe of the `to_list` result.
    Potenital type name:
        gene
        mRNA
        3UTR
        CDS
        5UTR
    """
    
    def __init__(self, start, end, chr_name=None, strand=None, feature_id=None, height=1, type_name=None, parm={}, children=[]):
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.strand = strand
        self.height = height
        self.id = feature_id 
        self.type = type_name
        if parm:
            self.parm = parm
        else:
            self.parm = {}
        if children:
            self.children = children
        else:
            self.children = []
        
    def _to_list(self, header=DEFAULT_FEATURE_HEADER):
        return [self.__getattribute__(h) for h in header]
        
    def to_list(self, gene_id=".", mRNA_id=".", header=DEFAULT_FEATURE_HEADER):        
        if self.type == "gene":
            gene_id = self.id
            mRNA_id = "."            
            d = self._to_list(header)
            d.extend([gene_id, mRNA_id])
            data = [d]
            for child_feature in self.children:
                data.extend(child_feature.to_list(gene_id, mRNA_id, header))
        elif self.type == "mRNA":
            mRNA_id = self.id
            d = self._to_list(header)
            d.extend([gene_id, mRNA_id])
            data = [d]
            for child_feature in self.children:
                data.extend(child_feature.to_list(gene_id, mRNA_id, header))
        else:
            d = self._to_list(header)
            d.extend([gene_id, mRNA_id])
            data = [d]
        return data
    
    def to_df(self, header=DEFAULT_FEATURE_HEADER):
        data = self.to_list(header)
        df = pd.DataFrame(data, columns=header + ["gene_id", "mRNA_id"])        

class FeatureArray:
    
    """
    You can provide a list of Feature objects or a file_name to use for `read` method to read
    features. So you can create a new class by rewrite the `read` method.
    
    Attributes:
        header: default: ["chr_name", "start", "end", "strand", "height", "id", "type"]
        file_name
        features: a list of feature objects
        df  : records the attribution values in header of each feature elements, didn't 
              record the children feature, and didn't record gene_id and mRNA_id value.
              But some method can direct modify df value. For example, you can add tow columns:
              "y_min" and "y_max", then you can use `trans_value` method to add these attributes 
              to feature objects.
              If you change df directly, you may need to use `sort_by_index_order` method to update
              features.
        index_data: a dict, key is feature.id, value is feature.
        index_order: default is np.arrange(len(self.df))
    Methods:
        to_list/to_df: like the to_list/to_df method of Feature
        trans_value: see Attirbutes `df`.
    
    Note: The order of df and features may be not consistent. The index_order is the order of df.
        So if you sort the features, you must change the index_order of df.
    """
    
    def __init__(self, features=None, file_name=None, 
                header=DEFAULT_FEATURE_HEADER, 
                parent_index=False, children_index=False):
        self.load_data(features, file_name, header, parent_index, children_index)
        self._pre_prosess_data_after_load()
    
    def is_empty(self):
        return not len(self.df)
    
    def load_data(self, features=None, file_name=None, 
                header=DEFAULT_FEATURE_HEADER, 
                parent_index=False, children_index=False):
        self.header = header
        self.file_name = file_name
        if features is None and file_name is not None:
            features = self.read(file_name)
        self.features = features
        df = pd.DataFrame([[f.__getattribute__(h) for h in header] for f in features],
                                columns=header)
        df["index_order"] = np.arange(len(df))
        self.df = df
        
        index_data = {}
        if parent_index:
            for f in features:
                index_data[f.id] = f
                if children_index:
                    for c in f.children:
                        index_data[c.id] = c
        self.index_data = index_data
    
    def _pre_prosess_data_after_load(self):
        pass
    
    def __iter__(self):
        return(iter(self.features))
    
    def to_list(self, header=DEFAULT_FEATURE_HEADER):
        data = []
        for feature in self.features:
            data.extend(feature.to_list(header=header))
        return data
        
    def to_df(self, header=DEFAULT_FEATURE_HEADER):
        data = self.to_list(header)
        return pd.DataFrame(data, columns=header + ["gene_id", "mRNA_id"])
    
    def trans_value(self, columns):
        features = self.features
        for column in columns:
            for i, v in zip(self.df["index_order"], self.df[column]):
                features[i].__setattr__(column, v)       
            
    def sort_by_index_order(self):
        new_features = []
        features = self.features
        for i in self.df["index_order"]:
            new_features.append(features[i])
        self.features = new_features
        self.df["index_order"] = np.arange(len(self.df))
    
    def get_feature(self, id):
        return(self.index_data[id])
        
    def get_features(self, ids):
        return([self.get_feature(i) for i in ids])
        
    def select_by_index(self, indexs=None):
        if indexs is None:
            indexs = self.df["index_order"]
        return([self.features[i] for i in indexs])
    
    def _pre_prosess_data(self):
        pass
    
    def load_region_data(self, region=None, select_features=None):
        chr_name, region_start, region_end, region_strand = region
        self.chr_name = chr_name
        self.region_strt = region_start
        self.region_end = region_end
        self.region_strand = region_strand
        if select_features is not None:
            self.region_data = type(self)(self.get_features(select_features), self.header)
        else:
            f1 = self.df["chr_name"] == chr_name
            f2 = (self.df["start"] >= region_start) & (self.df["start"] <= region_end)
            f3 = (self.df["end"] >= region_start) & (self.df["end"] <= region_end)
            f4 = (self.df["start"] <= region_start) & (self.df["end"] >= region_start)
            f = f1 & (f2 | f3 | f4)
            self.region_data = type(self)(self.select_by_index(self.df["index_order"][f]), self.header)
        self._pre_prosess_data()
        
    def open(self):
        pass
        
    def close(self):
        pass

class GeneFeatures(FeatureArray):
        
    """
    features的元素是基因或者mRNA,不建议混合
    """
    
    def read(self, file_name):
        pass
    
    def generate_intron(self, features):
        
        def _generate_intron(features=None):
            #直接添加到features里面
            if features:
                chr_name = features[0].chr_name
                strand = features[0].strand
                height = features[0].height
                if len(features) == 1:
                    features[0].parm["is_single"] = True
                else:
                    if strand == "+":
                        features[0].parm["is_first"] = True
                        features[-1].parm["is_last"] = True
                    else:
                        features[-1].parm["is_first"] = True
                        features[0].parm["is_last"] = True
                    for i, feature in enumerate(features[:-1]):
                        e = feature.end
                        next_s = features[i+1].start
                        if next_s > e + 1:
                            intron_start = e + 1
                            intron_end = next_s - 1
                            features.append(Feature(intron_start, intron_end, chr_name, strand, "", height, "intron", {}, []))
        
        def _pre_process(feature):
            if feature.type == "gene":
                for child_feature in feature.children:
                    _pre_process(child_feature)
            elif feature.type == "mRNA":
                _generate_intron(feature.children)
        
        for feature in features:
            _pre_process(feature)
        return features
                    
class BedFeatures(GeneFeatures):
    
    """
    Haven't record 5 score, 9 itemRgb and track line
    https://genome.ucsc.edu/FAQ/FAQformat.html
    """
    
    def read(self, file_name, based0=True, region_type_name="mRNA"):
        
        features = []
        i = 0
        for l in open(file_name):
            if l.startswith("browser"):
                continue
            elif l.startswith("track"):
                continue
            d = l.rstrip("\n").split()
            chr_name, region_start, region_end = d[:3]
            region_start, region_end = int(region_start), int(region_end) 
            if len(d) == 3:
                region_id = f"{chr_name}:{region_start}:{region_end}"
            else:
                region_id = d[3]
            if len(d) < 6:
                strand = "+"
            else:
                strand = d[5]
            if len(d) == 12:
                block_starts = np.array(d[11].split(",")[:-1]).astype("int") + region_start
                block_ends = block_starts + np.array(d[10].split(",")[:-1]).astype("int")
                left_cds = int(d[6])
                right_cds = int(d[7])
                UTR5s = []
                UTR3s = []
                CDSs = []
                for i, (block_start, block_end) in enumerate(zip(block_starts, block_ends)):
                    if left_cds >= block_end: # >= is > for 0-based
                        UTR5s.append([block_start, block_end])
                    elif left_cds >= block_start:
                        if left_cds > block_start:
                            UTR5s.append([block_start, left_cds])
                        if right_cds < block_end:
                            CDSs.append([left_cds, right_cds])
                            if right_cds != block_end - 1: # is != for 0-based
                                UTR3s.append([right_cds, block_end])
                            UTR3s.extend(list(zip(block_starts[i+1:], block_ends[i+1:])))
                            break
                        else:
                            CDSs.append([left_cds, block_end])
                    elif right_cds >= block_end: # >= is > for 0-based
                        CDSs.append([block_start, block_end])
                    elif right_cds >= block_start:
                        if right_cds > block_start:
                            CDSs.append([block_start, right_cds])
                        if right_cds != block_end - 1:
                            UTR3s.append([right_cds, block_end])
                        UTR3s.extend(list(zip(block_starts[i+1:], block_ends[i+1:])))
                        break
                if strand == "-":
                    type_names = ["3UTR", "CDS", "5UTR"]
                else:
                    type_names = ["5UTR", "CDS", "3UTR"]
                children = []
                for xys, type_name in zip([UTR5s, CDSs, UTR3s], type_names):
                    for start, end in xys:
                        if based0: start += 1
                        children.append(Feature(start, end, chr_name, strand, "", 1, type_name))
            else:
                if based0: 
                    region_start += 1
                children = [Feature(region_start, region_end, chr_name, strand, "", 1, "CDS")]
            if based0: 
                region_start += 1
            features.append(Feature(region_start, region_end, chr_name, strand, region_id, 
                        1, region_type_name, children=children))
        features = self.generate_intron(features)
        return(features) 
        
if __name__ == "__main__":
    filein, fileout = sys.argv[1:3]
    exon2bed(filein, fileout, input_0_based=False)
            
