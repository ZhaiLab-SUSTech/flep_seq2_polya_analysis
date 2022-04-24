from . track import BuGVTrack
import numpy as np
import matplotlib
import pandas as pd
from matplotlib.path import Path
import matplotlib.patches as patches

class InvertReaptReader:
    
    """
    Input (No header):
    chr_name, start1, end1, start2, end2
    start1 < end1 
    If invert repeat, start2 > end2
    """
    
    def __init__(self, file_name=""):
        
        self.file_name = file_name
        self.df = self.read_file(file_name)
        
    def read_file(self, file_name):
        return pd.read_table(file_name, names=["chr_name", "start1", "end1", "start2", "end2"])
        
    def load_region_data(self, region):
        chr_name, region_start, region_end, strand = region
        df = self.df
        f1 = df["chr_name"].values == chr_name
        start1 = df["start1"].values
        end1 = df["end1"].values
        start2 = df["start2"].values
        end2 = df["end2"].values
        f2 = (start1 >= region_start) & (start1 <= region_end)
        f3 = (end1 >= region_start) & (end1 <= region_end)
        f4 = (start2 >= region_start) & (start2 <= region_end)
        f5 = (end2 >= region_start) & (end2 <= region_end)
        f6 = (region_start >= start1) & (region_end <= end1)
        f7 = (region_end >= start2) & (region_end <= end2)
        f = f1 & (f2 | f3 | f4 | f5 | f6 | f7)
        df = df[f].copy()
        self.region_data = df
        
    def open(self):
        pass
    
    def close(self):
        pass


class InvertRepeatTrack(BuGVTrack):
    
    DEFAULT_CONFIG = {"file_type": InvertReaptReader,
                      "file_name": "",
                      "y_zero_line": {"plot": True},
                      "ylim_style" : None,
                      "y_zero_line": {
                          "plot": False
                      },
                      "invertrepeattrack": {
                          "color": "blue",
                          "edgecolor": None,
                          "linewidth": 2,
                          "alpha": 0.5,
                          "max_height_point_ratio": 3/4,
                          "compress_ratio": 1/2,
                          "block_ratio": 0,
                          "first_exon_arrow_ratio" : 0.2,
                          "first_exon_arrow_length": 100,
                          "last_exon_arrow_ratio"  : 0.2,
                          "last_exon_arrow_length" : 500,
                          "block_color": "black",
                          "block_alpha": 1,
                          "block_line_width": 0,
                          "block_edgecolor": None
                      }
                      } #注意所有对象共享改属性
    
    SINGLE_EXON_CODES = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]              
    
    def plot_ax(self):
        
        ls_config = self.config["invertrepeattrack"]
        color = ls_config["color"]
        edgecolor = ls_config["edgecolor"]
        linewidth = ls_config["linewidth"]
        alpha = ls_config["alpha"]
        max_height_point_ratio = ls_config["max_height_point_ratio"]
        compress_ratio = ls_config["compress_ratio"]
        block_ratio = ls_config["block_ratio"]
        first_exon_arrow_ratio = ls_config["first_exon_arrow_ratio"]
        first_exon_arrow_length = ls_config["first_exon_arrow_length"]
        last_exon_arrow_ratio = ls_config["last_exon_arrow_ratio"]
        last_exon_arrow_length = ls_config["last_exon_arrow_length"]
        block_color = ls_config["block_color"]
        block_alpha = ls_config["block_alpha"]
        block_line_width = ls_config["block_line_width"]
        block_edgecolor = ls_config["block_edgecolor"]

        ax = self.ax
        
        df = self.region_data
        df["r1"] = (df["end2"].values - df["end1"].values)/2
        df["r2"] = (df["start2"].values - df["start1"].values)/2
        y = 0
        
        max_y = max(max(df["r1"].max(), df["r2"].max())+y, 1)
        arc2 = Path.arc(0, 180)        
        for (start1, end1, start2, end2) in df[["start1", "end1", "start2", "end2"]].values:
            if block_ratio:
                y_lower = 0
                y_up = -block_ratio
                y_middle = (y_up-y_lower)/2 + y_lower
                for index, x_start, x_end in [[0, start1, end1], [1, start2, end2]]:
                    middle_point_x11 = x_start + (x_end - x_start) * first_exon_arrow_ratio
                    if x_start < x_end:
                        middle_point_x1 = min(middle_point_x11, x_start + first_exon_arrow_length)
                    else:
                        middle_point_x1 = max(middle_point_x11, x_start - first_exon_arrow_length)
                    middle_point_x21 = x_end - (x_end - x_start) * last_exon_arrow_ratio
                    if x_start < x_end:
                        middle_point_x2 = max(middle_point_x21, x_end - last_exon_arrow_length)
                    else:
                        middle_point_x2 = min(middle_point_x21, x_end + last_exon_arrow_length)
                    if index == 0:
                        end1 = middle_point_x2
                    else:
                        end2 = middle_point_x2
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
                    patch = patches.PathPatch(path, facecolor=block_color, 
                                                 lw=block_line_width,
                                                 edgecolor=block_edgecolor, 
                                                 alpha=block_alpha)
                    ax.add_patch(patch)
            r1 = (end2 - end1)/2
            r2 = (start2 - start1)/2
            x1 = end1 + r1
            x2 = start1 + r2
            if r1 < r2/2:
                r1 = r1/compress_ratio
                temp_angle1 = np.arctan(r2 * max_height_point_ratio /r1)
                temp_angle2 = 2*temp_angle1 - 1/2*np.pi
                r1 = r2 * max_height_point_ratio - r1 * np.tan(temp_angle2)
                angle1 = 360 - np.degrees(temp_angle2)
                angle2 = 180 + np.degrees(temp_angle2)
                arc1 = Path.arc(angle1, angle2)
                v1 = arc1.vertices * r1
                v1 = v1 * [compress_ratio,1]
                v1 = v1 + [x1, y-v1[0, 1]]
                c1 = arc1.codes
            else:
                v1 = arc2.vertices * r1 + [x1, y]
                c1 = arc2.codes
            v2 = arc2.vertices[::-1] * r2 + [x2, y]
            c2 = arc2.codes
            v = np.vstack([v1, v2[0, :], v2, v1[0, :], v1[0, :], (0, 0)])
            v = v * [1, 1/max_y]
            c = np.hstack([c1, Path.LINETO, c2, Path.LINETO, Path.MOVETO, Path.CLOSEPOLY])
            path=Path(v, c)
            patch = patches.PathPatch(path, facecolor=color, edgecolor=edgecolor, lw=linewidth, alpha=alpha)
            ax.add_patch(patch)
        
        self.set_ax_ylim(0-block_ratio, 1)
                
                
        
        
        
            