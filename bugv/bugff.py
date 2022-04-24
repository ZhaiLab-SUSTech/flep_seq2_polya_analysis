import re
import bugv.buseq as buseq
import sys
import os
from collections import defaultdict
import pandas as pd
import numpy as np

'''
使用gff_obj = GffRead(gtf_or_gff_file)生成gff_obj，文件后缀名是gtf，将会按照gtf格式读取，否则按照gff格式读取。

不同的读取方式依赖的是_GtfRead和_GffRead对象，它们只是读取方式不一样，但都属于_GtfGffRead对象，

因此共同的属性和方法是在_GtfGffRead对象里写的。

对象__init__时会调用_read_row_data方法。因此_GtfRead和_GffRead对象的核心不同在_read_row_data方法的不同。
核心数据存在_data属性。

1. 读取每行的基本方法_filter_row_data 以及 feature_data格式
读取每行时，它们都会用到_GtfGffRead的_filter_row_data方法，会将带有#号的行去掉，然后用"\t"分割，如果字段数不等于9则pass，
然后返回一个feature data的列表，分别对应gtf/gff文件九列：
[chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data]
注意，start和end此时是整数格式，annot_data是调用_parse_annot方法将原来的第九列数据转化为字典数据，注意gtf和gff格式parse annot时
用的正则表达式不同，分别是对应类的_re_parse_annot属性。

2. gtf读取方法
feature类型必须属于["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR"]。
对于每行，获得annot_data中的gene_id和transcript_id。
_data是一个字典，key是gene_id，value也是一个字典，肯定含有transcript属性。
transcript属性的值也是一个字典，key是transcript_id，值也是一个字典，
字典包括的基本key有["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR"], 每个key的值是一个列表，列表每个元素是一个feature data。
读取完所有行后会调用_reorganize方法，该方法将依次调用_reorganize_gtf_transcript和_reorganize_gtf_gene方法。
首先会将transcript中的各个基本key的features按起始位置由小到大排序。根据exon信息提取该transcript的chr，start，end，strand信息，
并添加到transcript这个字典上。如果没有exon，则会根据UTR和CDS信息生成exon，注意这时必须有CDS信息，而且exon的feature data会共用UTR和CDS的feature data。

_reorganize_gene方法将会根据transcript的chr,start,end,strand信息生成gene的chr,start,end,strand信息

3. gff读取方法
gff会首先按层级关系读，将信息存入_row_gff_data属性，该属性是个字典，有三个键，data（值为字典，其值为feature的ID），first_id（值为列表），sorted_ids（值为列表）。
首先获取annot_data的ID，Name和Parent信息（可以是个列表），如果没有Parent信息，则将其ID存入first_id里。每行的record_num累加
_row_gff_data["data"][ID] = [feature, parents, child, name, record_num, feature_data]
然后会统计结构。
最后再将gff文件转为gtf文件：
feature必须是这里面的：["five_prime_UTR", "exon", "start_codon", "CDS", "stop_codon", "three_prime_UTR"]
如果某feature不在这里面，那么它和它的子feature均不进行处理。
按照原始读取顺序遍历每个ID，处理时把"five_prime_UTR"和"three_prime_UTR"转化为"5UTR"和"3UTR"。同时transcript需要添加两个键，"start_codon"和"stop_codon"
然后调用_reorganize_gtf_transcript和_reorganize_gtf_gene方法。
'''

class Feature():
    
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

class FeatureArray():
        
    def __init__(self, features, header=["chr_name", "start", "end", "strand", "height", "id", "type"], parent_index=False, children_index=False):
        self.header = header
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
    
    def __iter__(self):
        return(iter(self.features))
    
    def trans_value(self, columns):
        features = self.features
        for column in columns:
            for i, v in zip(self.df["index_order"], self.df[column]):
                features[i].__setattr__(column, v)       
    
    def sort_like_igv(self, xaxis_is_invert=False, y_invert=False):
        #pandas:
        #feature_id, chr, start, end, strand, height, parm, child_feature
        #只要保证ls_start和ls_end都大于0，且ls_start小于ls_end即可。按照ls_start和ls_end进行排序
        #只要保证ls_start和ls_end都大于0，且ls_start小于ls_end即可。按照ls_start和ls_end进行排序
        df = self.df
        if not xaxis_is_invert:
            df["ls_start"] = df["start"]
            df["ls_end"] = df["end"]
        else:
            max_x = df["end"].max()
            df["ls_start"] = max_x - df["end"] + 1
            df["ls_end"] = max_x - df["start"] + 1
            df = df.sort_values(["ls_start", "ls_end"])
        
        y_mins = []
        have_overlap_regions = []
        now_max_x = 0
        for index, row in df.iterrows():
            start = row["ls_start"]
            end = row["ls_end"]
            height = row["height"]
            if start > now_max_x:
                y_min = 1
                y_max = height
                now_max_x = end
                have_overlap_regions = [[1, y_max, now_max_x]]
            else:
                for d in have_overlap_regions:
                    if d[2] < start:
                        d[2] = start - 1
                if len(have_overlap_regions) > 1:
                    new_have_overlap_regions = [have_overlap_regions[0]]
                    for d in have_overlap_regions[1:]:
                        if d[2] == new_have_overlap_regions[-1][2]:
                            new_have_overlap_regions[-1][1] = d[1]
                        else:
                            new_have_overlap_regions.append(d)
                    have_overlap_regions = new_have_overlap_regions
                have_insert = False
                for (i, d) in enumerate(have_overlap_regions):
                    x1, x2, x3 = d
                    if x3 < start and (x2 - x1 + 1) >= height:
                        have_insert = True
                        y_min = x1
                        y_max = y_min + height - 1
                        if y_max != x3:
                            d[0] = y_max + 1
                            have_overlap_regions.insert(i, [x1, y_max, end]) #可以尝试与上一个元素合并，不过意义不大，以为下一轮还要合并
                        else:
                            d[2] = end
                        break
                if not have_insert:
                    y_min = have_overlap_regions[-1][1] + 1 if have_overlap_regions else 1
                    y_max = y_min + height - 1
                    have_overlap_regions.append([y_min, y_max, end])
                if end > now_max_x:
                    now_max_x = end
            y_mins.append(y_min)

        df["y_min"] = y_mins
        df["y_min"] = df["y_min"] - 1
        df["y_max"] = df["y_min"] + df["height"]
        if y_invert:
            df["y_max"] = -df["y_max"]
            df["y_min"] = -df["y_min"]
        self.df = df
        self.trans_value(["y_min", "y_max"])
    
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
        if indexs == None:
            indexs = self.df["index_order"]
        return([self.features[i] for i in indexs])
                         
    def load_region_data(self, region=None, select_features=None):
        chr_name, region_start, region_end, region_strand = region
        self.chr_name = chr_name
        self.region_strt = region_start
        self.region_end = region_end
        self.region_strand = region_strand
        if select_features is not None:
            self.region_data = type(self)(self.get_features(select_features, self.header))
        else:
            f1 = self.df["chr_name"] == chr_name
            f2 = (self.df["start"] >= region_start) & (self.df["start"] <= region_end)
            f3 = (self.df["end"] >= region_start) & (self.df["end"] <= region_end)
            f = f1 & (f2 | f3)
            self.region_data = type(self)(self.select_by_index(self.df["index_order"][f]), self.header)
    
    def open(self):
        pass
        
    def close(self):
        pass

class GeneFeatures(FeatureArray):
        
    """
    features的元素是基因或者mRNA,不建议混合
    """
    
    def generate_intron(self, features):
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
            
    def _pre_process(self, feature):
        if feature.type == "gene":
            for child_feature in feature.children:
                self._pre_process(child_feature)
        elif feature.type == "mRNA":
            self.generate_intron(feature.children)
    
    def pre_process(self, track_obj):
        self.sort_like_igv(track_obj.xaxis_is_invert, track_obj.config["feature_plot"]["sort_reverse"])
        for feature in self.features:
            self._pre_process(feature)

class ListOrder():
    
    def __init__(self):
        self.orders = []
        self.datas = {}
    
    def add(self, s):
        if s not in self.datas:
            self.datas[s] = 1
            self.orders.append(s)
        else:
            self.datas[s] += 1

class _GtfGffRead():
    
    _re_parse_annot = re.compile("^(.+?) (.+)$")
    def __init__(self, file=""):
        if file:
            self.load(file)
        
    def load(self, file):
        self._file = file
        if file:
            self._read_row_data()
            
    def _parse_annot(self, s):
        re_parse_annot = self._re_parse_annot
        data = {}
        for l in s.split(";"):
            if not l:
                continue
            m = re_parse_annot.match(l)
            if m:
                key = m.group(1).strip()
                value = m.group(2).strip()
                value = value.strip('"')
                data[key] = value
            else:
                raise ValueError("Cann't prase gff/gtf annot %s\n" % (s))
        return data
        
    def _filter_row_data(self):
        
        file = self._file
        data = {}
        parser = self._parse_annot
        for l in open(file):
            if l.startswith("#"):
                continue
            d = l.rstrip().split("\t")
            if len(d) != 9:
                print("#filter row")
                continue
            chr_, annot_method, feature, start, end, score, strand, coding_frame, annot = d
            start, end = int(start), int(end)
            annot_data = parser(annot)
            yield([chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data])
    
    def _reorganize(self):
        self._reorganize_gtf_transcript()
        self._reorganize_gtf_gene()
    
    def cal_utr(self, transcript_data):
        #结果会和exons和cds等共享数据
        exons = transcript_data["exon"]
        cdss = transcript_data["CDS"]
        strand = transcript_data["strand"]
        if cdss:
            utr5 = []
            utr3 = []
            first_cdss_start, first_cdss_end  = cdss[0][3:5]
            last_cdss_start, last_cdss_end  = cdss[-1][3:5]
            find_f = False
            for e in exons:
                es, ee = e[3:5]
                if not find_f:
                    if ee < first_cdss_start:
                        utr5.append(e)
                    elif es < first_cdss_start:
                        utr5.append(e[0:3] + [es, first_cdss_start-1] + e[5:])
                        find_f = True
                    else:
                        find_f = True
                if find_f:
                    if es > last_cdss_end:
                        utr3.append(e)
                    elif ee > last_cdss_end:
                        utr3.append(e[0:3] + [last_cdss_end+1, ee] + e[5:])
            if strand == "-":
                utr5, utr3 = utr3, utr5
            for u in utr5:
                u[2] = "5UTR"
            for u in utr3:
                u[2] = "3UTR"
            return(utr5, cdss, utr3)            
        else:
            return(exons, [], [])

    def _reorganize_gtf_transcript(self):
        '''
        Haven't finished if exon not exist!
        '''
        data = self._data
        mRNA2dict = {}
        for gene, gene_data in data.items():
            for transcript, transcript_data in gene_data["transcript"].items():
                transcript_data["exon"].sort(key=lambda x: x[3])
                transcript_data["CDS"].sort(key=lambda x: x[3])
                transcript_data["5UTR"].sort(key=lambda x: x[3])
                transcript_data["3UTR"].sort(key=lambda x: x[3])
                exon = transcript_data["exon"]
                if not exon:
                    first_cds = transcript_data["CDS"][0]
                    exon = []
                    strand = first_cds[6]
                    if strand == "+":
                        utr1 = transcript_data["5UTR"]
                        utr2 = transcript_data["3UTR"]
                    else:
                        utr1 = transcript_data["3UTR"]
                        utr2 = transcript_data["5UTR"]
                    for chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data in utr1:
                        exon.append([chr_, annot_method, "exon", start, end, score, strand, coding_frame, annot_data])
                    chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data = first_cds
                    if exon and start == exon[-1][4] + 1:
                        exon[-1][4] = end
                        exon[-1][8] = annot_data
                    else:
                        exon.append([chr_, annot_method, "exon", start, end, score, strand, coding_frame, annot_data])
                    for chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data in transcript_data["CDS"][1:]:
                        exon.append([chr_, annot_method, "exon", start, end, score, strand, coding_frame, annot_data])
                    if utr2:
                        chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data = utr2[0]
                        if start == exon[-1][4] + 1:
                            exon[-1][4] = end
                        else:
                            exon.append([chr_, annot_method, "exon", start, end, score, strand, coding_frame, annot_data])
                        for chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data in utr2[1:]:
                            exon.append([chr_, annot_method, "exon", start, end, score, strand, coding_frame, annot_data])
                    transcript_data["exon"] = exon
                try:
                    start = exon[0][3]
                except:
                    print("good")
                    print(transcript)
                    print(transcript_data)
                end = exon[-1][4]
                strand = exon[0][6]
                chr_ = exon[0][0]
                transcript_data["chr"] = chr_
                transcript_data["start"] = start
                transcript_data["end"] = end
                transcript_data["strand"] = strand
                cal_utr5, cal_cds, cal_utr3 = self.cal_utr(transcript_data)
                transcript_data["Cal_5UTR"] = cal_utr5
                transcript_data["Cal_3UTR"] = cal_utr3
                transcript_data["gene"] = gene
                mRNA2dict[transcript] = transcript_data
        self._mRNA_data = mRNA2dict
    
    def del_genes(self, genes):
        genes = set(genes)
        data = self._data
        new_gene_orders = []
        for gene in genes:
            del data[gene]
        for gene in self._gene_orders:
            if gene not in genes: new_gene_orders.append(gene)
        self._gene_orders = new_gene_orders
                
    def _reorganize_gtf_gene(self):
        data = self._data
        no_mRNA_gene = []
        for gene, gene_data in data.items():
            transcript_poss = [[v["chr"], v["start"], v["end"], v["strand"]] for k, v in gene_data["transcript"].items()]
            if len(transcript_poss)==0:
                no_mRNA_gene.append(gene)
                continue
            start = min( [m[1] for m in transcript_poss] )
            end = max( [m[2] for m in transcript_poss] )
            strand = transcript_poss[0][3]
            chr_ = transcript_poss[0][0]
            gene_data["start"] = start
            gene_data["end"] = end
            gene_data["strand"] = strand  
            gene_data["chr"] = chr_
        self.del_genes(no_mRNA_gene)
        self.set_index()
            
    def _iter_mRNA(self):
        geneids = self._gene_orders
        data = self._data
        for gene in geneids:
            for mRNA, mRNA_data in data[gene]["transcript"].items():
               yield([gene, mRNA, mRNA_data])
    
    def get_mRNA(self, mRNA):
        return(self._mRNA_data[mRNA])
    
    def get_gene(self, gene):
        return(self._data[gene])
        
    def get_feature(self, feature):
        try:
            return((feature, self._data[feature], "gene"))
        except:
            return((feature, self._mRNA_data[feature], "mRNA"))
        
    def get_features(self, region=None, select_features=None):
        if select_features is not None:
            return([self.get_feature(feature) for feature in select_features])
        else:
            chr_name, region_start, region_end, region_strand = region
            self.chr_name = chr_name
            self.region_strt = region_start
            self.region_end = region_end
            self.region_strand = region_strand
            df = self.df
            f1 = df["chr_name"] == chr_name
            f2 = (df["start"] >= region_start) & (df["start"] <= region_end)
            f3 = (df["end"] >= region_start) & (df["end"] <= region_end)
            f = f1 & (f2 | f3)
            return(self.get_features(select_features=df["id"][f]))            
    
    def select2data(self, feature_datas):
        
        #和get_features搭配着用
                
        def get_mRNA_children(mRNA_data):
            data = []
            chr_name = mRNA_data["chr"]
            strand = mRNA_data["strand"]
            if strand == "+":
                type1s = ["Cal_5UTR", "CDS", "Cal_3UTR"]
                type2s = ["5UTR", "CDS", "3UTR"]
            else:
                type1s = ["Cal_3UTR", "CDS", "Cal_5UTR"]
                type2s = ["3UTR", "CDS", "5UTR"]
            
            for type1, type2 in zip(type1s, type2s):
                for d in mRNA_data[type1]:
                    start, end = d[3:5]
                    data.append(Feature(start, end, chr_name, strand, "", 1, type2))
            return(data)
        
        def mRNA2feature(mRNA, mRNA_data):
            return(Feature(mRNA_data["start"], mRNA_data["end"], 
                                    mRNA_data["chr"], mRNA_data["strand"], 
                                    mRNA, 1, "mRNA", children=get_mRNA_children(mRNA_data)))
        
        def get_gene_children(gene_data):            
            data = []
            for mRNA, mRNA_data in gene_data["transcript"].items():
                data.append(mRNA2feature(mRNA, mRNA_data))
            return(data)
        
        def gene2feature(gene, gene_data):
            children = get_gene_children(gene_data)
            return(Feature(gene_data["start"], gene_data["end"], 
                            gene_data["chr"], gene_data["strand"], 
                            gene, len(children), "gene", children=children))
                    
        data = []
        #id, chr_name, start, end, strand, type, height, parm, info, children
                
        for feature_id, feature_data, feature_type in feature_datas:
            if feature_type == "gene":
                feature_obj = gene2feature(feature_id, feature_data)
            elif feature_type == "mRNA":
                feature_obj = mRNA2feature(feature_id, feature_data)
            data.append(feature_obj)
            
        data = GeneFeatures(data, parent_index=True, children_index=True)
        self.region_data = data
        return(data)
    
    def load_region_data(self, region=None, select_features=None):

        feature_datas = self.get_features(region, select_features)
        self.region_data = self.select2data(feature_datas)
        return(self.region_data)
    
    def set_index(self):

        chr_names = []
        starts = []
        ends = []
        strands = []
        ids = []
                
        for gene, gene_data in self._iter_gene():
            chr_names.append(gene_data["chr"])
            starts.append(gene_data["start"])
            ends.append(gene_data["end"])
            strands.append(gene_data["strand"])
            ids.append(gene)
            
        df = pd.DataFrame({"chr_name": chr_names,
                           "start": starts,
                           "end": ends,
                           "strand": strands,
                           "id": ids})
        self.df = df

    
    def to_plot_format(self):
        
        #def __init__(self, start, end, chr_name=None, strand=None, feature_id=None, height=1, type_name=None, parm={}, children={}):
        
        def mRNA2df(mRNA_data):
            data = []
            chr_name = mRNA_data["chr"]
            strand = mRNA_data["strand"]
            if strand == "+":
                type1s = ["Cal_5UTR", "CDS", "Cal_3UTR"]
                type2s = ["5UTR", "CDS", "3UTR"]
            else:
                type1s = ["Cal_3UTR", "CDS", "Cal_5UTR"]
                type2s = ["3UTR", "CDS", "5UTR"]
            
            for type1, type2 in zip(type1s, type2s):
                for d in mRNA_data[type1]:
                    start, end = d[3:5]
                    data.append(Feature(start, end, chr_name, strand, "", 1, type2))
            return(data)
        
        def gene_mRNA2df(gene_data):            
            data = []
            for mRNA, mRNA_data in gene_data["transcript"].items():
                data.append(Feature(mRNA_data["start"], mRNA_data["end"], 
                                    mRNA_data["chr"], mRNA_data["strand"], 
                                    mRNA, 1, "mRNA", children=mRNA2df(mRNA_data)))
            return(data)
                    
        data = []
        #id, chr_name, start, end, strand, type, height, parm, info, children
        
        for gene, gene_data in self._iter_gene():
            children = gene_mRNA2df(gene_data)
            data.append(Feature(gene_data["start"], gene_data["end"], 
                                                gene_data["chr"], gene_data["strand"], 
                                                gene, len(children), "gene", children=children))
        data = GeneFeatures(data, parent_index=True, children_index=True)
        self.list_data = data
        return(data)
        
    def to_plot_format_old(self):
        
        #暂不用了
        def mRNA2df(mRNA_data):
            data = []
            chr_name = mRNA_data["chr"]
            strand = mRNA_data["strand"]
            for type1, type2 in zip(["Cal_5UTR", "CDS", "Cal_3UTR"], ["5UTR", "CDS", "3UTR"]):
                for d in mRNA_data[type1]:
                    start, end = d[3:5]
                    data.append([
                        "",
                        chr_name,
                        start,
                        end,
                        strand,
                        type2,
                        1,
                        {},
                        {},
                        {}
                    ])
            return(data)
        
        def gene_mRNA2df(gene_data):            
            data = []
            for mRNA, mRNA_data in gene_data["transcript"].items():
                data.append([
                    mRNA,
                    mRNA_data["chr"],
                    mRNA_data["start"],
                    mRNA_data["end"],
                    mRNA_data["strand"],
                    "mRNA",
                    1,
                    {},
                    {},
                    mRNA2df(mRNA_data)
                ])
            return(data)
                    
        data = []
        #id, chr_name, start, end, strand, type, height, parm, info, children
        
        for gene, gene_data in self._iter_gene():
            children = gene_mRNA2df(gene_data)
            data.append([gene, 
                gene_data["chr"], 
                gene_data["start"], 
                gene_data["end"], 
                gene_data["strand"],
                "gene",
                len(children),
                {},
                {},
                children
            ])
        self.list_data = data
        return(data)
    
    def to_df_old(self):
        
        #too slow
                
        def mRNA2df(mRNA_data):
            starts = []
            ends = []
            parms = []
            childrens = []
            types = []
            
            for type1, type2 in zip(["Cal_5UTR", "CDS", "Cal_3UTR"], ["5UTR", "CDS", "3UTR"]):
                for d in mRNA_data[type1]:
                    start, end = d[3:5]
                    starts.append(start)
                    ends.append(end)
                    types.append(type2)     
                    parms.append({})
                    childrens.append({})           
            
            df = pd.DataFrame({"start": starts,
                               "end": ends,
                               "parm": parms,
                               "type": types,
                               "children": childrens})
            df["chr_name"] = mRNA_data["chr"]
            df["strand"] = mRNA_data["strand"]
            df["height"] = 1
            return(df)
        
        def gene_mRNA2df(gene_data):
            starts = []
            ends = []
            parms = []
            childrens = []
            ids = []
            
            for mRNA, mRNA_data in gene_data["transcript"].items():
                starts.append(mRNA_data["start"])
                ends.append(mRNA_data["end"])
                parms.append({})
                ids.append(mRNA)
                childrens.append(mRNA2df(mRNA_data))
            
            df = pd.DataFrame({"start": starts,
                               "end": ends,
                               "parm": parms,
                               "id": ids,
                               "children": childrens})
            df["chr_name"] = gene_data["chr"]
            df["strand"] = gene_data["strand"]
            df["height"] = 1
            df["type"] = "mRNA"
            return(df)
        
        chr_names = []
        starts = []
        ends = []
        strands = []
        heights = []
        parms = []
        childrens = []
        ids = []
        
        i = 0
        for gene, gene_data in self._iter_gene():
            i += 1
            if i > 100:
                print(gene)
                i = 0
            chr_names.append(gene_data["chr"])
            starts.append(gene_data["start"])
            ends.append(gene_data["end"])
            strands.append(gene_data["strand"])
            parms.append({})
            ids.append(gene)
            children = gene_mRNA2df(gene_data)
            childrens.append(children)
            heights.append(len(childrens))
            
        df = pd.DataFrame({"chr_name": chr_names,
                           "start": starts,
                           "end": ends,
                           "strand": strands,
                           "parm": parms,
                           "id": ids,
                           "height": heights,
                           "children": childrens})
        df["type"] = "gene"
        self.df = df
        return(df)
    
    def _iter_mRNA_length(self):
        for gene, mRNA, mRNA_data in self._iter_mRNA():
            exons = mRNA_data["exon"]
            mRNA_length = 0
            for exon in exons:
                mRNA_length += exon[4] - exon[3] + 1
            yield([gene, mRNA, mRNA_length])
    
    def _iter_gene(self):
        geneids = self._gene_orders
        data = self._data
        for gene in geneids:
            yield([gene, data[gene]])
    
    def _iter_gene_pos(self):
        for gene, gene_data in self._iter_gene():
            chr_ = gene_data["chr"]
            start = gene_data["start"]
            end = gene_data["end"]
            strand = gene_data["strand"]
            yield([gene, [chr_, start, end, strand]])
    
    def write_mRNA_intron(self, fileout, fileout_no_intron):
        with open(fileout, 'w') as o, open(fileout_no_intron, 'w') as o_no:
            for gene, mRNA, mRNA_data in self._iter_mRNA():
                exons =  mRNA_data["exon"]
                if len(exons) == 1:
                    o_no.write(gene + "\t" + mRNA + "\n")
                    continue
                intron_start = exons[0][4] + 1
                if exons[0][6] == "+":
                    intron_ids = ["intron" + str(i+1) for i in range(len(exons)-1)]
                else:
                    intron_ids = ["intron" + str(i+1) for i in list(range(len(exons)-1))[::-1]]                    
                for (chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data), intron_id in zip(mRNA_data["exon"][1:], intron_ids):
                    intron_end = start - 1
                    o.write("\t".join([chr_, str(intron_start), str(intron_end), strand, mRNA, gene, intron_id])+"\n")
                    intron_start = end + 1
                    
    def write_mRNA_exon(self, fileout):
        with open(fileout, 'w') as o:
            for gene, mRNA, mRNA_data in self._iter_mRNA():
                exons = mRNA_data["exon"]
                if exons[0][6] == "+":
                    exon_ids = ["exon" + str(i+1) for i in range(len(exons))]
                else:
                    exon_ids = ["exon" + str(i+1) for i in list(range(len(exons)))[::-1]]
                for (chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data), exon_id in zip(exons, exon_ids):
                    o.write("\t".join([chr_, str(start), str(end), strand, mRNA, gene, exon_id])+"\n")
        
    def get_promter_pos(self, mRNA_data, feature_type="CDS", length=2000):
        #feature_type also can be "exon"
        cdss = mRNA_data[feature_type]
        strand = cdss[0][6]
        chr_ = cdss[0][0]
        start = cdss[0][3]
        end = cdss[-1][4]
        if strand == "+":
            promter_start = start - length
            promter_end = start - 1
        else:
            promter_start = end + 1
            promter_end = end + length
        
        return [chr_, promter_start, promter_end, strand]
    
    def iter_promoter(self, feature_type="CDS", length=2000):
        for gene, mRNA, mRNA_data in self._iter_mRNA():
            chr_, promter_start, promter_end, strand = self.get_promter_pos(mRNA_data, feature_type, length)
            promter_name = mRNA + ":promoter:" + chr_ + ":" + str(promter_start) + "-" + str(promter_end) + ":" + strand
            yield([chr_, promter_start, promter_end, strand, promter_name, mRNA, gene])
        
    def write_promter(self, fileout, feature_type="CDS", length=2000):
        for d in self.iter_promoter(feature_type, length):
            o.write("\t".join([str(i) for i in d])+"\n")
    
    def extract_promoter_seq(self, feature_type="CDS", length=2000, genome_fasta_file = "", genome_fasta_seqs = None):
        
        if not genome_fasta_seqs:
            genome_fasta_seqs = buseq.FastaSeqs(genome_fasta_file)
            genome_fasta_seqs.load_seq()
        
        _extract_seq_by_pos = genome_fasta_seqs.extract_seq_by_one_pos
        
        for d in self.iter_promoter(feature_type, length):
            chr_, promter_start, promter_end, strand, promter_name, mRNA, gene = d
            seq = _extract_seq_by_pos(chr_, [promter_start, promter_end, strand])
            yield(mRNA + "\t" + promter_name, seq)
        
    def write_promoter_seq(self, fileout, feature_type="CDS", length=2000, genome_fasta_file = "", genome_fasta_seqs = None, line_max_word=70):
        buseq.FastaSeqsTool().write_fasta(self.extract_promoter_seq(feature_type, length, genome_fasta_file, genome_fasta_seqs), fileout, line_max_word)
    
    def write_mRNA_length(self, fileout):
        with open(fileout, 'w') as o:
            for gene, mRNA, mRNA_length in self._iter_mRNA_length():
                o.write('\t'.join([gene, mRNA, str(mRNA_length)]) + "\n")
                    
    def write_gene_pos(self, fileout):
        gene_pos_data = self._iter_gene_pos()
        with open(fileout, 'w') as o:
            for gene, pos in gene_pos_data:
                chr_, start, end, strand = pos
                start = str(start)
                end = str(end)
                o.write('\t'.join([gene, chr_, start, end, strand]) + "\n")
                
    def extract_gene_seq(self, genome_fasta_file = "", genome_fasta_seqs = None):
        
        '''
        return generator: [gene_name, seq]
        '''
        
        if not genome_fasta_seqs:
            genome_fasta_seqs = buseq.FastaSeqs(genome_fasta_file)
            genome_fasta_seqs.load_seq()
        
        _extract_seq_by_pos = genome_fasta_seqs.extract_seq_by_one_pos
        
        for gene, pos in self._iter_gene_pos():
            chr_, start, end, strand = pos
            yield(gene, _extract_seq_by_pos(chr_, [start, end, strand]))
        
        #2019.01.06 replace
        #def _iter_pos():
        #    for gene, pos in self._iter_gene_pos():
        #        chr_, start, end, strand = pos
        #        yield([chr_, gene, [start, end, strand]])
        #
        #return genome_fasta_seqs.extract_seq_by_pos(pos_data = _iter_pos())
    
    def extract_cds_start_end(self):
        for gene, mRNA, mRNA_data in self._iter_mRNA():
            cds_data = mRNA_data["CDS"]
            try:
                chr_ = cds_data[0][0]
            except:
                print(mRNA + " extract cds wrong!")
                continue
            strand = cds_data[0][6]
            start = cds_data[0][3]
            end = cds_data[-1][4]
            if strand == "-":
                start, end = end, start
            yield([gene, mRNA, strand, start, end])
    
    def extract_cds_or_mRNA(self, is_cds=0, genome_fasta_file = "", genome_fasta_seqs = None):
        
        feature_type = "CDS" if is_cds else "exon"    
        
        if not genome_fasta_seqs:
            genome_fasta_seqs = buseq.FastaSeqs(genome_fasta_file)
            genome_fasta_seqs.load_seq()
        
        _extract_seq_by_pos = genome_fasta_seqs.extract_seq_by_one_pos
        
        for gene, mRNA, mRNA_data in self._iter_mRNA():
            cds_data = mRNA_data[feature_type]
            #chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data
            try:
                chr_ = cds_data[0][0]
            except:
                print(mRNA + " extract cds wrong!")
                continue
            strand = cds_data[0][6]
            seqs = []
            for chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data in cds_data:
                seq = _extract_seq_by_pos(chr_, [start, end, strand])
                seqs.append(seq)
            if strand == "-":
                seqs = seqs[::-1]
            cds_seq = ''.join(seqs)
            yield([gene, mRNA, cds_seq])
    
    def extract_cds_seq(self, genome_fasta_file = "", genome_fasta_seqs = None):
        return self.extract_cds_or_mRNA(1, genome_fasta_file, genome_fasta_seqs)
    
    def extract_mRNA_seq(self, genome_fasta_file = "", genome_fasta_seqs = None):
        return self.extract_cds_or_mRNA(0, genome_fasta_file, genome_fasta_seqs)
    
    def extract_protein_seq(self, genome_fasta_file = "", genome_fasta_seqs = None):
        translate_seq = buseq.FastaSeqTool().translate
        for gene, mRNA, cds_seq in self.extract_cds_seq(genome_fasta_file, genome_fasta_seqs):
            yield([gene, mRNA, translate_seq(cds_seq)])
    
    def write_cds_end(self, fileout):
        with open(fileout, 'w') as o:
            o.write("gene\tmRNA\tstrand\tcds_start\tcds_end\n")
            for gene, mRNA, strand, start, end in self.extract_cds_start_end():
                o.write("\t".join([gene, mRNA, strand, str(start), str(end)]) + "\n")
    
    def write_cds_seq(self, fileout, genome_fasta_file = "", genome_fasta_seqs = None, line_max_word=70):
        seqs = ([mRNA, cds_seq] for gene, mRNA, cds_seq in self.extract_cds_seq(genome_fasta_file, genome_fasta_seqs))
        buseq.FastaSeqsTool().write_fasta(seqs, fileout, line_max_word)        
    
    def write_mRNA_seq(self, fileout, genome_fasta_file = "", genome_fasta_seqs = None, line_max_word=70):
        seqs = ([mRNA, mRNA_seq] for gene, mRNA, mRNA_seq in self.extract_mRNA_seq(genome_fasta_file, genome_fasta_seqs))
        buseq.FastaSeqsTool().write_fasta(seqs, fileout, line_max_word)
    
    def write_protein_seq(self, fileout, genome_fasta_file = "", genome_fasta_seqs = None, line_max_word=70):
        seqs = ([mRNA, protein_seq] for gene, mRNA, protein_seq in self.extract_protein_seq(genome_fasta_file, genome_fasta_seqs))
        buseq.FastaSeqsTool().write_fasta(seqs, fileout, line_max_word)        

    def write_gene_seq(self, fileout, genome_fasta_file = "", genome_fasta_seqs = None, line_max_word=70):
        seqs = self.extract_gene_seq(genome_fasta_file, genome_fasta_seqs)
        buseq.FastaSeqsTool().write_fasta(seqs, fileout, line_max_word)
        
    def write_gtf(self, fileout):
        with open(fileout, 'w') as o:
            for gene, mRNA, mRNA_data in self._iter_mRNA():
                for feature in ["exon", "CDS", "5UTR", "3UTR"]:
                    for chr_, annot_method, feature1, start, end, score, strand, coding_frame, annot_data in mRNA_data[feature]:
                        annot = 'transcript_id "%s"; gene_id "%s";' % (mRNA, gene)
                        o.write("\t".join([chr_, annot_method, feature, str(start), str(end), score, strand, coding_frame, annot])+"\n")
                        
    def open(self):
        pass
        
    def close(self):
        pass

class _GtfRead(_GtfGffRead):
    _re_parse_annot = re.compile("^(.+?) (.+)$")
    
    def __init__(self, file=""):
        super().__init__(file)
        
    def _anotdata_to_anot(self):
        '''
        Noets: transcript_id "AT1G01010.1" --> data["transcript_id"] = AT1G01010.1
        lost "" information.
        '''
        pass
    
    def _write_line(self, data):
        pass
    
    def _read_row_data(self):
        
        '''
        Not support other feature!!!
        feature_data = [chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data]
        annot_data = { gene_id: "",
                       transcript_id: ""
                     }
        {"gene": transcript:  {"transcript_id": 5UTR:[feature_data],
                                                exon:[feature_data],
                                                start_codon:[feature_data],
                                                CDS:[feature_data],
                                                stop_codon:[feature_data],
                                                3UTR: [feature_data]
                              }
        }
        '''
        other_feature = set(["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR"])
        filted_row_data = self._filter_row_data()
        data = {}
        chr_orders = ListOrder()
        gene_orders = ListOrder()
        transcript_orders = ListOrder()
        for feature_data in filted_row_data:
            chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data = feature_data
            if feature == "five_prime_utr":
                feature = "5UTR"
            elif feature == "three_prime_utr":
                feature = "3UTR"
            if feature not in other_feature:
                continue
            gene_id = annot_data["gene_id"]
            transcript_id = annot_data["transcript_id"]
            
            if gene_id not in data:
                data[gene_id] = {"transcript":{}}
            gene_transcript_data = data[gene_id]["transcript"]
            if transcript_id not in gene_transcript_data:
                transcript_data = gene_transcript_data[transcript_id] = {
                                                                            "5UTR": [],
                                                                            "exon": [],
                                                                            "start_codon": [],
                                                                            "CDS": [],
                                                                            "stop_codon":[],
                                                                            "3UTR": []
                                                                        }
            gene_transcript_data[transcript_id][feature].append(feature_data)
            chr_orders.add(chr_)
            gene_orders.add(gene_id)
            transcript_orders.add(transcript_id)
        self._data = data
        #for transcript, transcript_data in self._data["Glyma01G000600"]["transcript"].items():
        #    print(transcript, transcript_data)
        self._chr_orders = chr_orders.orders
        self._gene_orders = gene_orders.orders
        self._transcript_orders = transcript_orders.orders
        self._reorganize()
     
class _GffRead(_GtfGffRead):
    
    '''
    默认用_read_row_data()读取数据到_row_gff_data，
    然后用_row_gff_data_2_gtf_data()将原始数据转化为gtf数据，存储到_data。
    
    gff_data含三个键，data,first_id,sorted_ids。gff文件中一般含有gene,
    mRNA,CDS等特征。除了基因外，其他特征有一个Parent属性，指示其父节点是什么，
    这里first_id记录没有父节点的行，也就是基因。
    data里是以ID为键，存储了[feature, parent, child, name, record_num, feature_data]。
    其中feature_data前几列和原gff文件一样，只是start和end已经转化为int，同时最后的注释列变为字典。
    sorted_ids按读取顺序存储了id值。
    
    '''
    
    _re_parse_annot = re.compile("^(.+?)=(.*)$")

    def __init__(self, file=""):
        super().__init__(file)
        
    def _read_row_data(self):        
        filted_row_data = self._filter_row_data()
        gff_data = {"data":{}, "first_id":[], "sorted_ids":[]}
        gff_data_data = gff_data["data"]
        gff_data_first_id = gff_data["first_id"]
        gff_data_sorted_ids = gff_data["sorted_ids"]
        record_num = 0
        rename_id = 0
        for feature_data in filted_row_data:
            feature = feature_data[2]
            annot_data = feature_data[-1]
            try:
                id_ = annot_data["ID"]
            except:
                id_ = str(rename_id)
                rename_id += 1
            name = annot_data.get("Name")
            parents = annot_data.get("Parent")
            if parents:
                parents = parents.split(",")
            else:
                parents = []
            child = []
            record_num += 1
            gff_data_data[id_] = [feature, parents, child, name, record_num, feature_data]
            gff_data_sorted_ids.append(id_)
            if parents:
                for parent in parents:
                    try:
                        gff_data_data[parent][2].append(id_)
                    except:
                        gff_data_data[parent][2].append(id_)
            else:
                gff_data_first_id.append(id_)        
        self._row_gff_data = gff_data
        self._stat_structure()
        self._row_gff_data_2_gtf_data()
    
    def _stat_structure(self):
        gff_data = self._row_gff_data
        gff_data_data = gff_data["data"]
        gff_data_first_id = gff_data["first_id"]
        gff_data_sorted_ids = gff_data["sorted_ids"]
        
        feature_structure = defaultdict(dict)
        for id_ in gff_data_sorted_ids:
            feature, parents, child, name, record_num, feature_data = gff_data_data[id_]
            feature_structure[feature][id_] = {}
            for parent in parents:
                if parent:
                    parent_feature = gff_data_data[parent][0]
                    try:
                        feature_structure[parent_feature][parent][feature] += 1
                    except:
                        feature_structure[parent_feature][parent][feature] = 1
        summary_feature_structure = {}
        for feature, data in feature_structure.items():
            total_nums = len(data.keys())
            summary_feature_structure[feature] = {"total_nums": total_nums, "child_feature": {}}
            for id_, child_features in data.items():
                for child_feature, num in child_features.items():
                    try:
                        summary_feature_structure[feature]["child_feature"][child_feature][num] += 1
                    except:
                        try:
                            summary_feature_structure[feature]["child_feature"][child_feature][num] = 1
                        except:
                            summary_feature_structure[feature]["child_feature"][child_feature] = {num:1}
            for child_feature, child_feature_data in summary_feature_structure[feature]["child_feature"].items():
                child_feature_data[0] = total_nums - sum(child_feature_data.values())
        for feature, feature_data in summary_feature_structure.items():
            print("\n#######################################")
            print(feature + "\t\t" + str(feature_data["total_nums"]))
            for child_feature, child_feature_data in feature_data["child_feature"].items():
                print("\t\t" + child_feature, child_feature_data)
            print("#######################################")
            
    def _row_gff_data_2_gtf_data(self):
        
        '''
        gff3: five_prime_UTR, three_prime_UTR
        gtf: 5UTR, 3UTR
        '''
        
        row_gff_data = self._row_gff_data
        gff_data = row_gff_data["data"]
        gff_ids = row_gff_data["sorted_ids"]
        gtf_data = {}
        chr_orders = ListOrder()
        gene_orders = ListOrder()
        transcript_orders = ListOrder()
        other_feature = set(["five_prime_UTR", "exon", "start_codon", "CDS", "stop_codon", "three_prime_UTR"])
        mRNA_indicator = {}
        pass_ids = set()
        for gff_id in gff_ids:
            feature, parents, child, name, record_num, feature_data = gff_data[gff_id]
            if feature == "gene":
                chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data = feature_data
                chr_orders.add(chr_)
                gene_orders.add(gff_id)            
                gtf_data[gff_id] = {"transcript": {},
                                     "start": start,
                                     "end": end,
                                     "strand": strand,
                                     "chr": chr_,
                                     "feature_data": feature_data }
            #<bug> Araport11_GFF3中有些gene没有mRNA，而是标记为lnc_RNA或antisense_lncRNA，
            #这里直接将lnc_RNA当成mRNA了。没有做任何标记。以后可以改一改
            #<origin>
            elif feature == "mRNA":
            #</orgin>
            #elif feature in ["mRNA", "lnc_RNA", "antisense_lncRNA"]:
            #</bug>
                chr_, annot_method, feature, start, end, score, strand, coding_frame, annot_data = feature_data
                transcript_orders.add(gff_id)
                for parent in parents:
                    if parent in pass_ids:
                        pass_ids.add(gff_id)
                        continue
                    mRNA_indicator[gff_id] = gtf_data[parent]["transcript"][gff_id] = { "5UTR": [],
                                                                                        "exon": [],
                                                                                        "start_codon": [],
                                                                                        "CDS": [],
                                                                                        "stop_codon":[],
                                                                                        "3UTR": [],
                                                                                        "start": start,
                                                                                        "end": end,
                                                                                        "strand": strand,
                                                                                        "chr": chr_
                                                                                      }
            elif feature in other_feature:
                if feature == "five_prime_UTR":
                    feature = "5UTR"
                elif feature == "three_prime_UTR":
                    feature = "3UTR"
                try:
                    for parent in parents:
                        if parent in pass_ids:
                            pass_ids.add(gff_id)
                            continue
                        mRNA_indicator[parent][feature].append(feature_data)
                    #</bug>
                except:
                    #<bug> 以前是报错退出，现在让报错不退出
                    #<origin>
                    raise Warning("%s not belong mRNA" % gff_id)
                    #</origin>
                    #print("%s not belong mRNA\n" % gff_id)
                    #</bug>
            else:
                #<bug> 以前是报错退出，现在让报错不退出
                #<origin>
                #raise Warning("Specific feature not be processed: %s " % (feature))
                #</origin>
                print("Specific feature not be processed: %s %s" % (feature, gff_id))
                pass_ids.add(gff_id)
                #</bug>
        self._data = gtf_data
        self._chr_orders = chr_orders.orders
        self._gene_orders = gene_orders.orders
        self._transcript_orders = transcript_orders.orders
        self._reorganize()
    
    def _display_structure(self):
        row_gff_data = self._row_gff_data["data"]
        parent_child = {}
        for id_, value in row_gff_data.items():
            feature = value[0]
            children = value[2]
            for child in children:
                child_feature = row_gff_data[child][0]
                key = feature + "\t" + child_feature
                try:
                    parent_child[key] += 1
                except:
                    parent_child[key] = 1
        
        for k, v in parent_child.items():
            print(k + "\t" + str(v))  

"""
            old
            class GffRead():
        
                def __init__(self, file):
                    if file:
                        self.load(file)
    
                def _check_gff_or_gtf(self, file=""):
                    if not file:
                        file = self._file
                    if file:
                        file_short_name, file_name_extension = os.path.splitext(file)
                        file_name_extension = file_name_extension.lower()
                        if file_name_extension == ".gtf":
                            return "gtf"
                        elif file_name_extension in [".gff", ".gff3"]:
                            return "gff"
                        else:
                            raise IOError("File :'%s' is not gtf or gff file. Please use correct file extension name." % (file))
                    else:
                        raise IOError("File :'%s' is not gtf or gff file. Please use correct file extension name." % (file))
        
                def load(self, filein):
                    self._file = filein
                    self.gtf_or_gff = gtf_or_gff = self._check_gff_or_gtf()
                    if gtf_or_gff == "gtf":
                        self._object = _GtfRead(self._file)
                    else:
                        self._object = _GffRead(self._file)
    
                def __getattr__(self, attribute):
                    return getattr(self._object, attribute)
            
"""

def GffRead(file_name):
    
    def _check_gff_or_gtf(file_name):
        if file_name:
            file_short_name, file_name_extension = os.path.splitext(file_name)
            file_name_extension = file_name_extension.lower()
            if file_name_extension == ".gtf":
                return "gtf"
            elif file_name_extension in [".gff", ".gff3"]:
                return "gff"
            else:
                raise IOError("File :'%s' is not gtf or gff file. Please use correct file extension name." % (file_name))
        else:
            raise IOError("File :'%s' is not gtf or gff file. Please use correct file extension name." % (file_name))
    
    gtf_or_gff = _check_gff_or_gtf(file_name)
    if gtf_or_gff == "gtf":
        obj = _GtfRead(file_name)
    else:
        obj = _GffRead(file_name)
    obj._file = file_name
    obj.gtf_or_gff = gtf_or_gff
    return(obj)

def blocks_in_blocks(block_xs, block_ys, up=0, down=0):
    
    """
    Input:
    Both block_xs and block_ys are dict:
        key:   chromatin_name
        value: a list of all read_exons in this chromatin, 
               each element is :
               [block_name, start, end, read_alignment_strand]
    
    Output:
    A iterator. Each element is a list represent a overlap relationship: 
    [x_id, y_id, to_block_dir, type5, pos5, type3, pos3, 
        chr_name, x_start, x_end, x_dir, y_start, y_end, y_dir]
    x_id is block_name in block_xs
    y_id is block_name in block_ys
    to_block_dir "+":  the strand of x_id and x_id is same, "-": not same
    
    0-based
    
    type5, pos5: the 5' indicate the 5' positon of x_id relative to the 5' of y_id
    type3, pos3: the 3' indicate the 3' positon of x_id relative to the 5' of y_id
    
    
    y_id (length: l)  ----->>>>>>>>>>>>>>>>>>>>>>---------      >>>> repesent y_id
                      |    |   |                |    |
                      |  on5(0)|            on3(l-1) |
                    up(-4)  in(4)                    |
                                                down (l+4)
    
    x_id                    
    type5: on5 type3: on3  >>>>>>>>>>>>>>>>>>>>>>
    type5: on3 type3: on5  <<<<<<<<<<<<<<<<<<<<<<
    """
    
    def pos_2_block(pos, block_start, block_end, block_dir="+", onisin = False):
        # give a block 
        type_, relative_position = "", 0
        #type list : "down","on5","in","on3","up"
        if block_dir == "+":
            relative_position = pos - block_start  # 0 for the same position
            if pos < block_start:
                type_ = "up"
            elif pos == block_start:
                type_ = "on5"
            elif pos < block_end:
                type_ = "in"
            elif pos == block_end:
                type_ = "on3"
            else:
                type_ = "down"
        else:
            relative_position = block_end - pos
            if pos > block_end:
                type_ = "up"
            elif pos == block_end:
                type_ = "on5"
            elif pos > block_start:
                type_ = "in"
            elif pos == block_start:
                type_ = "on3"
            else:
                type_ = "down"
        if onisin:
            if type_ in ["on5", "on3"]:
                type_ = "in"
        return((type_, relative_position))
        
    def blocks_in_blocks_one_chr(x, y, up=0, down=0, chr_name=""):
        '''
        Find overlapped genome features (blocks) when give two
        kind of genome features. 
        x, y are lists of two kind of genome features sorted by start position.
        Blocks are also list: [id_, start, end, dir_].
        The start position is smaller than the end position no matther 
        what the direction is.
        '''
        
        results = []
        for data in y:
            start,end = data[1:3]
            up_pos = start - up
            if up_pos < 1: up_pos = 1
            data.append(up_pos)
            data.append(end+down)
        
        for x_id,x_start,x_end,x_dir in x:
            del_pos = []
            for i, (y_id,y_start,y_end,y_dir,y_up,y_down) in enumerate(y):

                if x_end < y_up: 
                    break
                elif x_start > y_down:
                    del_pos.append(i)
                    continue

                if x_dir == "+": 
                    x_5_terminal, x_3_terminal = x_start, x_end
                else:
                    x_5_terminal, x_3_terminal = x_end, x_start
            
                if x_dir == y_dir:
                    to_block_dir = "+"
                else:
                    to_block_dir = "-"

                type5, pos5 = pos_2_block(x_5_terminal, y_start, y_end, y_dir)
                type3, pos3 = pos_2_block(x_3_terminal, y_start, y_end, y_dir)

                results.append([x_id, y_id, to_block_dir, 
                                type5, pos5, type3, pos3, 
                                chr_name, x_start, x_end, x_dir,
                                y_start, y_end, y_dir])
            del_pos.reverse()
            for i in del_pos:
                del y[i]
        return results
        
    for chr_, this_chr_block_xs in block_xs.items():
        if chr_ not in block_ys:
            continue
        this_chr_block_ys = block_ys[chr_]
        results = blocks_in_blocks_one_chr(this_chr_block_xs, this_chr_block_ys, 
                                  up, down, chr_)
        for r in results:
            yield(r)
               
def gff_write_gene_seq(file_geneome_fasta, file_gff, fileout_gene_fasta):
    gffread = GffRead(file_gff)
    gffread.write_gene_seq(fileout_gene_fasta, file_geneome_fasta)

def gff_write_mRNA_seq(file_geneome_fasta, file_gff, fileout_gene_fasta):
    gffread = GffRead(file_gff)
    gffread.write_mRNA_seq(fileout_gene_fasta, file_geneome_fasta)

def gff_write_cds_seq(file_geneome_fasta, file_gff, fileout_gene_fasta):
    gffread = GffRead(file_gff)
    gffread.write_cds_seq(fileout_gene_fasta, file_geneome_fasta)

def gff_write_protein_seq(file_geneome_fasta, file_gff, fileout_protein_fasta):
    gffread = GffRead(file_gff)
    gffread.write_protein_seq(fileout_protein_fasta, file_geneome_fasta)

def gff_write_cds_end(file_gff, fileout):
    gffread = GffRead(file_gff)
    gffread.write_cds_end(fileout)

def gff_write_gene_pos(file_gff, fileout):
    gffread = GffRead(file_gff)
    gffread.write_gene_pos(fileout)

def gff_write_intron_pos(file_gff, fileout, fileout_no_intron):
    gffread = GffRead(file_gff)
    gffread.write_mRNA_intron(fileout, fileout_no_intron)

def gff_write_mRNA_pos(file_gff, fileout):
    gffread = GffRead(file_gff)
    gffread.write_mRNA_exon(fileout)
    
def gff_write_promoter_pos(file_gff, fileout, feature_type="CDS", length=2000):
    gffread = GffRead(file_gff)
    gffread.write_promter(fileout, feature_type, length)
        
def gff_write_promoter_seq(file_geneome_fasta, file_gff, fileout, feature_type="CDS", length=2000):
    gffread = GffRead(file_gff)
    gffread.write_promoter_seq(fileout, feature_type, length, file_geneome_fasta)    

def gff2gtf(file_gff, fileout):
    gffread = GffRead(file_gff)
    gffread.write_gtf(fileout)
    
def main():
    file_geneome_fasta, file_gff, fileout_gene_fasta = sys.argv[1:]
    gff_write_gene_seq(file_geneome_fasta, file_gff, fileout_gene_fasta)
    
if __name__ == '__main__':
    main()