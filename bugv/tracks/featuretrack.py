from . track import BuGVTrack
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
from matplotlib import patches
import mpl_toolkits.axisartist as axisartist
import matplotlib.lines as lines
from bugv.bugff import GffRead
from bugv.readers.feature import BedFeatures
from bugv.readers.fulllengthrnabam import FulllengthBam, RNABamFeature
import numpy as np
import pandas as pd
import time
"""
version 0.3 by jiajinbu 2021.09.08
simple plot

version 0.2.3 by jiajinbu 2021.03.14

version 0.2.2 by jiajinbu 2021.01.25
add feature_plot.plot_pos_not_consistent paramter.

version 0.2.1 by jiajinbu 2021.01.04
"""

def add_rect_to_rects(rect, facecolor, lw, edgecolor, zorder, rects):
    parms = (facecolor, lw, edgecolor, zorder)
    try:
        rects[parms].append(rect)
    except:
        rects[parms] = [rect]

class FeaturesTrack(BuGVTrack):
    
    """
    feature是一个字典，必须含有"chr_name", "start", "end", "strand", "children", "height",
    children是一个列表，可以为空。为空表示不含有子feature。每个元素是一个feature。
    height必须指定。一般为整数。

    when plot, it first use `pre_process_data` method to handle region_data.
    and you then try to use .config["feature_plot"]["feature_process_func"] func to handle region_data.
    in `pre_process_data` it first run `filter_data` method to filter some data, then sort the data 
    based on self.config["feature_plot"]["sort_method"]. the sort method will add order, y_max and y_min
    column, and trans_value 'y_max' and 'y_min' to each feature obj. the plot will iter each feature and
    plot.
    """
    
    DEFAULT_CONFIG = {"file_type": GffRead,
                      "file_name": "",
                      "ylim_style" : None,
                      "max_cache_obj_number" : 5,
                      "y_zero_line": { "plot": False },
                      "x_base_pad": 0.5,
                      "read_obj_kwargs": {"parent_index": True},
                      "keep_select_features": True,
                      "plot_show": True,
                      ##for feature plot                      
                      "feature_plot" : {
                          "plot_gene_type" : "up_arrow", #"up_arrow",
                          "up_arrow": {
                          "cds_ratio": 0.5,
                          "arrow_length": 100,
                          "arrow_header_length": 50,
                          "utr_ratio": 0.25,
                          "arrow_height_ratio": 0.7,
                          "arrow_header_ratio": 0.3,
                          "arrow_color": "black",
                          "arrow_line_width": 0.4,
                          "arrow_header_linewidth": 0.2,
                          },
                          "sort_method": "igv",
                          "feature_process_func": None,
                          "linewidth" : 0.5,
                          "exon" : {
                              "is_first" : False,
                              "is_last" : False,
                              "is_single" : False,
                              "facecolor" : "black",
                              "edgecolor" : "black",
                              "zorder" : 10,
                              "first_exon_arrow_ratio" : 0.2,
                              "first_exon_arrow_length": 100,
                              "last_exon_arrow_ratio"  : 0.2,
                              "last_exon_arrow_length" : 100
                          },
                          "polyA" : {
                              "plot_show": True,
                              "is_first" : False,
                              "is_last" : False,
                              "is_single" : False,
                              "facecolor" : "red",
                              "edgecolor" : "red",
                              "zorder" : 10,
                              "first_exon_arrow_ratio" : 0.2,
                              "first_exon_arrow_length": 100,
                              "last_exon_arrow_ratio"  : 0.2,
                              "last_exon_arrow_length" : 100
                          },
                          "5UTR" : {
                              "is_first" : False,
                              "is_last" : False,
                              "is_single" : False,
                              "facecolor" : "black",
                              "edgecolor" : "black",
                              "zorder" : 10,
                              "first_exon_arrow_ratio" : 0.2,
                              "first_exon_arrow_length": 100,
                              "last_exon_arrow_ratio"  : 0.2,
                              "last_exon_arrow_length" : 100
                          },
                          "3UTR" : {
                              "is_first" : False,
                              "is_last" : False,
                              "is_single" : False,
                              "facecolor" : "black",
                              "edgecolor" : "black",
                              "zorder" : 10,
                              "first_exon_arrow_ratio" : 0.2,
                              "first_exon_arrow_length": 100,
                              "last_exon_arrow_ratio"  : 0.2,
                              "last_exon_arrow_length" : 100
                          },
                          "CDS" : {
                              "is_first" : False,
                              "is_last" : False,
                              "is_single" : False,
                              "facecolor" : "black",
                              "edgecolor" : "black",
                              "zorder" : 10,
                              "first_exon_arrow_ratio" : 0.2,
                              "first_exon_arrow_length": 100,
                              "last_exon_arrow_ratio"  : 0.2,
                              "last_exon_arrow_length" : 100
                          },
                          "intron" : {
                              "intron_height_ratio": 0.01,
                              "edgecolor" : "black",
                              "facecolor" : "black",
                              "zorder" : 9,
                              "linewidth": 0.5
                          },
                          "gene" : {
                              
                          },
                          "mRNA" : {
                              
                          },
                          "feature_space" : 0.2, #比例
                          "sort_reverse" : True,
                          "plot_show": True,
                      }
                      } #注意所有对象共享改属性
    
    FIRST_OR_LAST_EXON_CODES = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]  
    
    SINGLE_EXON_CODES = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]
    
    ARROW_HEADER_CODES = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]
    
    ARROW_CODES = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]
    
    def set_up_arrow_gene_style(self):
        self.load_config(
            {
                "feature_plot": {
                    "plot_gene_type": "up_arrow",
                    "up_arrow": {
                        "cds_ratio": 0.5,
                        "arrow_length": 100,
                        "arrow_header_length": 50,
                        "utr_ratio": 0.25,
                        "arrow_height_ratio": 0.7,
                        "arrow_header_ratio": 0.3,
                        "arrow_color": "black",
                        "arrow_line_width": 0.4,
                        "arrow_header_linewidth": 0.2,
                    },
                    "5UTR": {"facecolor": "black", "edgecolor": "black"},
                    "3UTR": {"facecolor": "black", "edgecolor": "black"},
                    "intron": {"intron_height_ratio": 0.01},
                }
            }
        )
    
    def set_white_black_gene_style(self):
        self.load_config(
            {
                "feature_plot": {
                    "plot_gene_type": "white_black_arrow",
                    "5UTR": {"facecolor": "white"},
                    "3UTR": {"facecolor": "white"},
                    "intron": {"intron_height_ratio": 0.25},
                }
            }
        )
    
    def plot_ax(self):
        #feature_data is a list.
        #print(time.time(), "pre_process_data")
        self.pre_process_data()
        feature_data = self.region_data
        
        if self.config["feature_plot"]["feature_process_func"] is not None:
            feature_data = self.config["feature_plot"]["feature_process_func"](feature_data)
        
        ax = self.ax
        
        if feature_data.is_empty():
            ylim_min = 0
            ylim_max = 1
        else:
            ylim_min = min([feature_data.df["y_min"].min(), feature_data.df["y_max"].min()])
            ylim_max = max([feature_data.df["y_min"].max(), feature_data.df["y_max"].max()])
        self.set_ax_ylim(ylim_min, ylim_max)
        rects = {}
        #print(time.time(), "plot")
        for feature in feature_data:
            self.plot_feature(feature, ax, y_is_plus=True, rects=rects)
        #print(time.time(), "add_rect")
        for (facecolor, lw, edgecolor, zorder), this_rects in rects.items():
            rect_collections = PatchCollection(this_rects, facecolor=facecolor, lw=lw, edgecolor=edgecolor, zorder=zorder)
            ax.add_collection(rect_collections)
        #print(time.time(), "add_rect_end")            
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.set_xticks([])
        ax.set_yticks([])

    def pre_process_data(self):
        if not self.region_data.is_empty():
            #print(time.time(), "pre_process_filter_data")
            self.filter_data()
            sort_method = self.config["feature_plot"]["sort_method"]
            #print(time.time(), "pre_process_sort_data")
            if not self.region_data.is_empty():
                if sort_method == "igv":
                    self.sort_features_like_igv()
                else:
                    self.sort_features_by_order()
            #print(time.time(), "pre_process_sort_data_end")

    def filter_data(self):
        #feataure_data = self.region_data
        pass        
        
    def sort_features_by_order(self):

        featurearray_obj = self.region_data
        xaxis_is_invert = self.xaxis_is_invert
        y_invert = self.config["feature_plot"]["sort_reverse"]

        df = featurearray_obj.df
        
        if "order" not in df:
            self.add_feature_order()
            df = featurearray_obj.df     
        
        tmp_order_index = df["order"].values.argsort()
        tmp_order_rank = tmp_order_index.argsort()
        df["y_max"] = df["height"].values[tmp_order_index].cumsum()[tmp_order_rank]
        df["y_min"] = df["y_max"] - df["height"]
        
        if y_invert:
            df["y_max"] = -df["y_max"]
            df["y_min"] = -df["y_min"]
        featurearray_obj.df = df
        featurearray_obj.trans_value(["y_min", "y_max"])
    
    def add_feature_order(self):
        
        featurearray_obj = self.region_data
        xaxis_is_invert = self.xaxis_is_invert
        
        df = featurearray_obj.df
        if not xaxis_is_invert:
            df = df.sort_values(["start", "end"])
        else:
            df = df.sort_values(["start", "end"], ascending=False)
        df["order"] = np.arange(len(df))
        featurearray_obj.df = df

    def sort_features_like_igv(self):
        #pandas:
        #feature_id, chr, start, end, strand, height, parm, child_feature
        #只要保证ls_start和ls_end都大于0，且ls_start小于ls_end即可。按照ls_start和ls_end进行排序
        #只要保证ls_start和ls_end都大于0，且ls_start小于ls_end即可。按照ls_start和ls_end进行排序

        featurearray_obj = self.region_data
        xaxis_is_invert = self.xaxis_is_invert
        y_invert = self.config["feature_plot"]["sort_reverse"]

        df = featurearray_obj.df
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
        featurearray_obj.df = df
        featurearray_obj.trans_value(["y_min", "y_max"])
        
    def plot_feature(self, feature, ax, y_is_plus=True, rects=None):
        #除了第一级可以为None，后续均不能为None
        y_min = feature.y_min
        y_max = feature.y_max
        children = feature.children
        self.plot_base_feature(feature, ax, y_min, y_max, y_is_plus, rects)
        #for child_feature in children:
        #    self.plot_base_feature(child_feature, ax, y_min, y_max, y_is_plus)
        
    def plot_base_feature(self, feature, ax, feature_y_min, feature_y_max, y_is_plus=True, rects=None):
        x_base_pad = self.config["x_base_pad"]
        plot_gene_type = self.config["feature_plot"]["plot_gene_type"]
        plot_up_arrow = self.config["feature_plot"]["up_arrow"]
        if feature.chr_name:
            if feature.chr_name != self.chr_name:
                return None
        elif (feature.end < self.region_start) or (feature.start > self.region_end):
            return None

        feature_type = feature.type
        get_parm = self.get_parm_func([feature.parm, 
                                       self.config["feature_plot"][feature_type],
                                       self.config["feature_plot"],
                                       self.config])

        origin_feature_y_min = feature_y_min
        origin_feature_y_max = feature_y_max
        ls_feature_y_space = feature_y_max - feature_y_min
        ls_feature_space = self.config["feature_plot"]["feature_space"] * ls_feature_y_space
        feature_y_max = feature_y_max - ls_feature_space
        feature_y_min = feature_y_min + ls_feature_space

        x_start = feature.start
        x_end = feature.end
        x_start -= x_base_pad
        x_end += x_base_pad

        if feature.strand == "-":
            x_start, x_end = x_end, x_start
        y_lower = feature_y_min
        y_up = feature_y_max
        if y_up < y_lower:
            y_up, y_lower = y_lower, y_up
        y_middle = (y_lower + y_up)/2.0        
        
        if not get_parm("plot_show"):
            return None
        
        if feature_type == "exon" or feature_type == "5UTR" or feature_type == "3UTR" or feature_type == "CDS" or feature_type == "polyA":
            if plot_gene_type == "white_black_arrow":
                if get_parm("is_single"):
                    middle_point_x11 = x_start + (x_end - x_start) * get_parm("first_exon_arrow_ratio")
                    if feature.strand == "+":
                        middle_point_x1 = min(middle_point_x11, x_start + get_parm("first_exon_arrow_length"))
                    else:
                        middle_point_x1 = max(middle_point_x11, x_start - get_parm("first_exon_arrow_length"))
                    middle_point_x21 = x_end - (x_end - x_start) * get_parm("last_exon_arrow_ratio")
                    if feature.strand == "+":
                        middle_point_x2 = max(middle_point_x21, x_end - get_parm("last_exon_arrow_length"))
                    else:
                        middle_point_x2 = min(middle_point_x21, x_end + get_parm("last_exon_arrow_length"))
                    verts = [
                        (x_start, y_lower),  # left, bottom
                        (middle_point_x1, y_middle),
                        (x_start, y_up),  # left, top
                        (middle_point_x2, y_up),  # right, top
                        (x_end, y_middle),
                        (middle_point_x2, y_lower),  # right, bottom
                        (0., 0.)  # ignored
                    ]
                    path = Path(verts, self.SINGLE_EXON_CODES)
                    #try plt.fill_between, maybe faster
                    patch = patches.PathPatch(path, facecolor=get_parm("facecolor"), 
                                                 lw=get_parm("linewidth"),
                                                 edgecolor=get_parm("edgecolor"), 
                                                 zorder=get_parm("zorder"))
                    ax.add_patch(patch)
                elif get_parm("is_first"):
                    middle_point_x11 = x_start + (x_end - x_start) * get_parm("first_exon_arrow_ratio")
                    if feature.strand == "+":
                        middle_point_x = min(middle_point_x11, x_start + get_parm("first_exon_arrow_length"))
                    else:
                        middle_point_x = max(middle_point_x11, x_start - get_parm("first_exon_arrow_length"))
                    verts = [
                        (x_start, y_lower),  # left, bottom
                        (middle_point_x, y_middle),
                        (x_start, y_up),  # left, top
                        (x_end, y_up),  # right, top
                        (x_end, y_lower),  # right, bottom
                        (0., 0.)  # ignored
                    ]
                    path = Path(verts, self.FIRST_OR_LAST_EXON_CODES)
                    patch = patches.PathPatch(path, facecolor=get_parm("facecolor"), 
                                                 lw=get_parm("linewidth"),
                                                 edgecolor=get_parm("edgecolor"), 
                                                 zorder=get_parm("zorder"))
                    ax.add_patch(patch)
                elif get_parm("is_last"):                
                    middle_point_x21 = x_end - (x_end - x_start) * get_parm("last_exon_arrow_ratio")
                    if feature.strand == "+":
                        middle_point_x = max(middle_point_x21, x_end - get_parm("last_exon_arrow_length"))
                    else:
                        middle_point_x = min(middle_point_x21, x_end + get_parm("last_exon_arrow_length"))
                    verts = [
                        (x_start, y_lower),  # left, bottom
                        (x_start, y_up),  # left, top
                        (middle_point_x, y_up),  # right, top
                        (x_end, y_middle),
                        (middle_point_x, y_lower),  # right, bottom
                        (0., 0.)  # ignored
                    ]
                    path = Path(verts, self.FIRST_OR_LAST_EXON_CODES)
                    patch = patches.PathPatch(path, facecolor=get_parm("facecolor"), 
                                                 lw=get_parm("linewidth"),
                                                 edgecolor=get_parm("edgecolor"), 
                                                 zorder=get_parm("zorder"))
                    ax.add_patch(patch)
                else:
                    rect = patches.Rectangle((x_start, y_lower), width=x_end - x_start, height=y_up - y_lower)
                    facecolor=get_parm("facecolor")
                    lw=get_parm("linewidth")
                    edgecolor=get_parm("edgecolor")
                    zorder=get_parm("zorder")
                    if rects is None:
                        ax.add_patch(rect)
                    else:
                        add_rect_to_rects(rect, facecolor, lw, edgecolor, zorder, rects)
            elif plot_gene_type == "up_arrow":
                origin_y_height = y_up - y_lower
                if feature_type == "5UTR" or feature_type == "3UTR":
                    finall_y_lower = y_lower + origin_y_height * (plot_up_arrow["cds_ratio"] - plot_up_arrow["utr_ratio"])/2
                    finall_height =  origin_y_height * plot_up_arrow["utr_ratio"]
                else:
                    finall_y_lower = y_lower
                    finall_height = origin_y_height * plot_up_arrow["cds_ratio"]
                rect = patches.Rectangle((x_start, finall_y_lower), width=x_end - x_start, height=finall_height)
                facecolor=get_parm("facecolor")
                lw=get_parm("linewidth")
                edgecolor=get_parm("edgecolor")
                zorder=get_parm("zorder")
                if rects is None:
                    ax.add_patch(rect)
                else:
                    add_rect_to_rects(rect, facecolor, lw, edgecolor, zorder, rects)

                if get_parm("is_single") or get_parm("is_first"):
                    #arrow_line_end_x
                    #arrow_line_width 
                    #arrow_point_x
                    #arrow_point_y
                    #arrow_width_x
                    #arrow_height_y
                
                
                    arrow_x_start = x_start
                    arrow_line_length = plot_up_arrow["arrow_length"] - plot_up_arrow["arrow_header_length"]
                    if feature.strand == "+":
                        arrow_header_x_start = arrow_x_end = x_start + arrow_line_length
                        arrow_header_x_end = arrow_x_end + plot_up_arrow["arrow_header_length"]
                    else:
                        arrow_header_x_start = arrow_x_end = x_start - arrow_line_length
                        arrow_header_x_end = arrow_x_end - plot_up_arrow["arrow_header_length"]
                    arrow_y_start = finall_y_lower + finall_height
                    arrow_y_end = y_lower + plot_up_arrow["arrow_height_ratio"] * origin_y_height
                    arrow_header_y_lower = arrow_y_end - plot_up_arrow["arrow_header_ratio"] * origin_y_height / 2
                    arrow_header_y_up = arrow_y_end + plot_up_arrow["arrow_header_ratio"] * origin_y_height / 2
                
                
                
                    print(arrow_x_start, arrow_x_end, arrow_y_start, arrow_y_end)
                    ax.plot((arrow_x_start, arrow_x_start), (arrow_y_start, arrow_y_end), lw=plot_up_arrow["arrow_line_width"], 
                                    color=plot_up_arrow["arrow_color"])
                    ax.plot((arrow_x_start, arrow_x_end), (arrow_y_end, arrow_y_end), lw=plot_up_arrow["arrow_line_width"], 
                                    color=plot_up_arrow["arrow_color"])
                
                    verts = [
                        (arrow_header_x_start, arrow_header_y_lower),
                        (arrow_header_x_start, arrow_header_y_up),
                        (arrow_header_x_end, arrow_y_end), 
                        (0., 0.)
                    ]
                    path = Path(verts, self.ARROW_HEADER_CODES)
                    #try plt.fill_between, maybe faster
                    patch = patches.PathPatch(path, facecolor=plot_up_arrow["arrow_color"], 
                                                 lw=plot_up_arrow["arrow_header_linewidth"],
                                                 edgecolor=plot_up_arrow["arrow_color"])
                    ax.add_patch(patch)
                
                    """
                    print("arrow", (arrow_x_end, arrow_y_end), (arrow_x_start, arrow_y_start))
                    ax.annotate(
                        "",
                        xy=(arrow_x_end, arrow_y_end), # arrow point position
                        xytext=(arrow_x_start, arrow_y_start), # line end position
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="angle", color=plot_up_arrow["arrow_color"], linewidth=0.2),
                    )
                    """
            else:
                pass              
        elif feature_type == "intron":
            if plot_gene_type == "white_black_arrow":
                exon_ratio = 1
            elif plot_gene_type == "up_arrow":
                exon_ratio = plot_up_arrow["cds_ratio"]
            else:
                exon_ratio = 1
            origin_y_height = y_up - y_lower
            intron_height_ratio = get_parm("intron_height_ratio")
            finall_y_lower = y_lower + origin_y_height * (exon_ratio - intron_height_ratio)/2
            finall_height =  origin_y_height * intron_height_ratio
            rect = patches.Rectangle((x_start, finall_y_lower), width=x_end - x_start, height=finall_height)
            facecolor=get_parm("facecolor")
            lw=get_parm("linewidth")
            edgecolor=get_parm("edgecolor")
            zorder=get_parm("zorder")
            if rects is None:
                ax.add_patch(rect)
            else:
                add_rect_to_rects(rect, facecolor, lw, edgecolor, zorder, rects)
            #line = lines.Line2D((x_start, x_end), (y_middle, y_middle),
            #        lw=get_parm("linewidth"), 
            #        color=get_parm("edgecolor"), axes=ax, zorder=get_parm("zorder"))
            #ax.add_line(line)
        elif feature_type == "gene":
            y_space = (origin_feature_y_max - origin_feature_y_min)/feature.height
            for i, mRNA in enumerate(feature.children):
                y_min = origin_feature_y_min + y_space * i
                y_max = origin_feature_y_min + y_space * (i + 1)
                self.plot_base_feature(mRNA, ax, y_min, y_max, y_is_plus, rects)
        elif feature_type == "mRNA":
            for child_feature in feature.children:
                self.plot_base_feature(child_feature, ax, origin_feature_y_min, origin_feature_y_max, y_is_plus, rects)
        
class BedTrack(FeaturesTrack):
    
    def init(self):
        self.load_config({"file_type": BedFeatures})
        self.set_up_arrow_gene_style()
        
        
class RNABamAlignTrack(FeaturesTrack):
    
    """
    The input is a bam file generated by 
    """
    
    def init(self):
        self.set_white_black_gene_style()
        self.load_config({"file_type": RNABamFeature, 
                          "read_obj_kwargs": {},
                          "ylim_style" : "not_change",
                          "feature_plot" : {
                              "feature_process_func": None,
                              "sort_method": "igv",
                              "linewidth" : 0,
                              "exon" : {
                                  "is_first" : False,
                                  "is_last" : False,
                                  "is_single" : False,
                                  "facecolor" : '#5D93C4',
                                  "edgecolor" : '#5D93C4',
                                  "zorder" : 10,
                                  "first_exon_arrow_ratio" : 0.2,
                                  "first_exon_arrow_length": 100,
                                  "last_exon_arrow_ratio"  : 0.2,
                                  "last_exon_arrow_length" : 100
                              },
                              "intron" : {
                                  "edgecolor" : '#A6A6A6', #change in 2021.01
                                  "zorder" : 9,
                                  "linewidth": 0,
                                  "facecolor" : "#A6A6A6",
                              },
                              "feature_space" : 0.2, #比例
                              "sort_reverse" : True
                          }
                          })
    
class FulllengthBamTrack(FeaturesTrack):
    
    """
    The input is a bam file generated by add_read_info_tag_to_bam.py in FLEPseq pipeline.
    
    You can choice feature_plot.select_read_types to determine which read types will be ploted. 
    See detail in feature_plot.select_read_types_dict. You also can directly set a list to feature_plot.select_read_types.
    You can set feature_plot.read_type_order to determine the plot order of each kind of read type.
    See detain in feature_plot.read_type_order_dict.
    You can set read_type_sort_method to determine which method will be used to sort the read in a specific read type.
    The method name includes "start_end", "start", "end", "ir".
    """
    
    def init(self):
        self.set_white_black_gene_style()
        self.load_config({"file_type": FulllengthBam, 
                          "read_obj_kwargs": {},
                          "ylim_style" : "not_change",
                          "plot_show": True,
                          "feature_plot" : {
                              "select_read_types": "high_quality",  #or None, add  = 
                              "plot_pos_not_consistent": True,
                              "select_read_types_dict": {
                                  "all": None,
                                  "high_quality": ["elongating", "polya"],
                                  "elongating": ["elongating"],
                                  "polya": ["polya"],
                                  "polya_contain_5lost": ["polya_3_not_in_last_exon", "polya_5lost", "polya"],
                                  "apa": ["polya", "polya_3_not_in_last_exon"],
                                  "mRNA": ["elongating", "splicing_intermediate", "elongating_3_mapping_low_accuracy",
                                          "elongating_5lost", "polya_3_not_in_last_exon", "polya_5lost", "polya"]
                              },
                              "read_type_order": "combine_elongating", #you can provide str to select different type, or a list of list, see read_type_order_dict config and add_feature_order method
                              "read_type_order_dict": {
                                  "seperate": [[["elongating"], "elongating"],
                                         [["splicing_intermediate"], "splicing_intermediate"],
                                         [["elongating_3_mapping_low_accuracy"], "elongating_3_mapping_low_accuracy"],
                                         [["elongating_5lost"], "elongating_5lost"],
                                         [["polya_3_not_in_last_exon"], "polya_3_not_in_last_exon"],
                                         [["polya_5lost"], "polya_5lost"],
                                         [["polya"], "polya"],
                                         [["antisense"], "antisense"],
                                         [["no_adapter"], "no_adapter"],
                                         [["duplicated"], "duplicated"],
                                         [["integenic"], "integenic"],
                                         [["other"], "other"]
                                        ],
                                  "combine_elongating": [[["elongating", "splicing_intermediate", "elongating_3_mapping_low_accuracy", ], "elongating"],
                                         [["elongating_5lost"], "elongating_5lost"],
                                         [["polya_3_not_in_last_exon"], "polya_3_not_in_last_exon"],
                                         [["polya_5lost"], "polya_5lost"],
                                         [["polya"], "polya"],
                                         [["antisense"], "antisense"],
                                         [["no_adapter", "duplicated", "integenic", "other"], "other"],
                                        ]
                              },
                              "read_type_sort_method_default": {"elongating": "end",
                                      "splicing_intermediate": "end",
                                      "elongating_3_mapping_low_accuracy": "end",
                                      "elongating_5lost": "end",
                                      "polya_3_not_in_last_exon": "end",
                                      "polya_5lost": "end",
                                      "polya": "ir",
                                      "antisense": "start_end",
                                      "no_adapter": "start_end",
                                      "duplciated": "start_end",
                                      "integenic": "start_end",
                                      "other": "start_end"
                              },
                              "read_type_sort_method_modfiy_by_read_type_order": {
                                  "apa": {"polya": "end"}
                              },
                              "read_type_sort_method": {},
                              "feature_process_func": None,
                              "sort_method": "not_order",
                              "linewidth" : 0,
                              "exon" : {
                                  "is_first" : False,
                                  "is_last" : False,
                                  "is_single" : False,
                                  "facecolor" : '#5D93C4',
                                  "edgecolor" : '#5D93C4',
                                  "zorder" : 10,
                                  "first_exon_arrow_ratio" : 0.2,
                                  "first_exon_arrow_length": 100,
                                  "last_exon_arrow_ratio"  : 0.2,
                                  "last_exon_arrow_length" : 100
                              },
                              "polyA" : {
                                  "plot_show": True,
                                  "is_first" : False,
                                  "is_last" : False,
                                  "is_single" : False,
                                  "facecolor" : 'lightcoral',
                                  "edgecolor" : 'lightcoral',
                                  "zorder" : 10,
                                  "first_exon_arrow_ratio" : 0.2,
                                  "first_exon_arrow_length": 100,
                                  "last_exon_arrow_ratio"  : 0.2,
                                  "last_exon_arrow_length" : 100
                              },
                              "intron" : {
                                  "intron_height_ratio": 0.25,
                                  "edgecolor" : '#A6A6A6', #change in 2021.01
                                  "zorder" : 9,
                                  "linewidth": 0,
                                  "facecolor" : "#A6A6A6",
                              },
                              "feature_space" : 0.2, #比例
                              "sort_reverse" : True,
                              "plot_show": True,
                          }
                          })
    
    def filter_data(self):
        #it is better to move this method to FulllengthBam class
        region_data = self.region_data

        df = region_data.df
        
        select_features = self.config["select_features"]
        if select_features:
            df = df.query('mRNA_id in @select_features')
        select_read_types = self.config["feature_plot"]["select_read_types"]
        if isinstance(select_read_types, str):
            select_read_types = self.config["feature_plot"]["select_read_types_dict"][select_read_types]
        if select_read_types:
            df = df.query('read_type in @select_read_types')
        if not self.config["feature_plot"]["plot_pos_not_consistent"] and "pos_not_consistent" in df:
            df = df.query('pos_not_consistent == 0')
        region_data.df = df
        region_data.sort_by_index_order()
        
    def add_feature_order(self):
        
        def add_type_inner_order(df):
            #you can provide a dict, key is read type, value is start, end, start_end, ir.
            def add_ir_order(tmp_df):
    
                def tmp_cal_ir_order(ir_nums, max_intron_num):
                    s = 0
                    for ir_num in ir_nums:
                        if ir_num:
                            s += 2**(max_intron_num - ir_num)
                    return s

                irs_unique_str = tmp_df.intron_retention.unique()
                irs_unique_int = [[int(i) if i else 0 for i in s.split(":")] for s in irs_unique_str]
                max_intron_num = max(max(i) for i in irs_unique_int)
                irs_order = [tmp_cal_ir_order(i, max_intron_num) for i in irs_unique_int]
                df_irs2order = pd.DataFrame({"intron_retention": irs_unique_str, "irs_order": irs_order}).sort_values("irs_order")
                tmp_df = tmp_df.merge(df_irs2order, on="intron_retention")
                return tmp_df
            
            def get_read_type_order_method_dict():
                feature_config = self.config["feature_plot"]
                read_type_order_method_dict = feature_config["read_type_sort_method_default"].copy()
                
                select_read_types = feature_config["select_read_types"]
                if isinstance(select_read_types, str):
                    if select_read_types in feature_config["read_type_sort_method_modfiy_by_read_type_order"]:
                        for read_type, sort_method in feature_config["read_type_sort_method_modfiy_by_read_type_order"][select_read_types].items():
                            read_type_order_method_dict[read_type] = sort_method
                        
                read_type_sort_method = feature_config["read_type_sort_method"]
                for read_type, sort_method in read_type_sort_method.items():
                    read_type_order_method_dict[read_type] = sort_method
                
                return read_type_order_method_dict
            
            def get_sort_func(sort_method="end"):
                
                def start_end_func(rna_strand, tmp_df):
                    if rna_strand == "+": 
                        tmp_df = tmp_df.sort_values(["start", "end"])
                    else:
                        tmp_df = tmp_df.sort_values(["end", "start"], ascending=False)
                    return tmp_df
                
                def start_func(rna_strand, tmp_df):
                    if rna_strand == "+": 
                        tmp_df = tmp_df.sort_values(["start"])
                    else:
                        tmp_df = tmp_df.sort_values(["end"], ascending=False)
                    return tmp_df
                    
                def end_func(rna_strand, tmp_df):
                    if rna_strand == "+": 
                        tmp_df = tmp_df.sort_values(["end"])
                    else:
                        tmp_df = tmp_df.sort_values(["start"], ascending=False)
                    return tmp_df
                    
                def polya_func(rna_strand, tmp_df):
                    return tmp_df.sort_values(["polya_length"])
                    
                def ir_func(rna_strand, tmp_df):
                    sort_by = "rentention_intron" #sort by retention_introns
                    tmp_df = add_ir_order(tmp_df)
                    tmp_df = tmp_df.sort_values("irs_order", ascending=True)
                    return tmp_df
                    
                sort_method_func_dict = {
                    "start_end": start_end_func,
                    "end": end_func, 
                    "start": start_func, 
                    "ir": ir_func,
                    "polya": polya_func
                }
                
                return sort_method_func_dict[sort_method]
                
            read_ids = []
            orders = []
            
            sort_by = ""
            sort_by_read_type_column_name = "new_read_type" if "new_read_type" in df else "read_type"
            
            for (mRNA_id, read_type, rna_strand), tmp_df in df.groupby(["mRNA_id", sort_by_read_type_column_name, "rna_strand"]):
                
                sort_method = get_read_type_order_method_dict()[read_type] 
                sort_func = get_sort_func(sort_method)
                tmp_df = sort_func(rna_strand, tmp_df)
                read_ids.extend(tmp_df.id.values)
                orders.extend(np.arange(len(tmp_df)))
            
            df = df.merge(pd.DataFrame({"id": read_ids, 
                                        "type_inner_order": orders}),
                          on="id")
            return df
                
        def add_read_type_order(df):
            read_type_ordered = self.config["feature_plot"]["read_type_order"]
            if isinstance(read_type_ordered, str):
                read_type_ordered = self.config["feature_plot"]["read_type_order_dict"][read_type_ordered]
            tmp_origin_read_type = []
            tmp_new_read_type = []
            tmp_type_orders = []
            for type_order, (need_to_combine_types, new_type_name) in enumerate(read_type_ordered):
                for origin_read_type in need_to_combine_types:
                    tmp_origin_read_type.append(origin_read_type)
                    tmp_new_read_type.append(new_type_name)
                    tmp_type_orders.append(type_order)
            read_type_orderd_df = pd.DataFrame({"read_type": tmp_origin_read_type, "new_read_type": tmp_new_read_type, "read_type_order": tmp_type_orders})
            df = df.merge(read_type_orderd_df, on="read_type")
            return df
        
        def add_mRNA_order(df):
            mRNA_id_orders = df.mRNA_id.values
            mRNA_id_orders[mRNA_id_orders == ""] = "ZZZZZ"
            df["mRNA_id_order"] = mRNA_id_orders
            return df
        
        featurearray_obj = self.region_data

        df = featurearray_obj.df
        df = add_read_type_order(df)
        df = add_type_inner_order(df)
        df = add_mRNA_order(df)
        
        df = df.sort_values(["mRNA_id_order", "read_type_order", "type_inner_order"])        
        df["order"] = np.arange(len(df))
        featurearray_obj.df = df        
        