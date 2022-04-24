from . track import BuGVTrack
from bugv.readers.srnabam import sRNABam

class sRNAPointTrack(BuGVTrack):
    
    DEFAULT_CONFIG = {"plot_set_func" : None,
                      "file_type": sRNABam,
                      "file_name": "",
                      "read_obj_kwargs": {"fileformat": "shortstack"},
                      "is_unique": False,
                      "split_strand": True,
                      "y_zero_line": {"plot": True},
                      "ylim_style" : "change_less",
                      "srnapointtrack": {
                          "max_count": None, #150
                          "point_color": ['#71C9F8', '#11870D', '#FF8205', '#717171'], #21nt, 22nt, 24nt, 其他小RNA的颜色
                          #['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3']
                      },
                      "library_size": 1000000,
                      "plot_tpm": False,
                      "ylog": False
                      } #注意所有对象共享改属性
    
    def plot_ax(self):
        #使用point_color, point_size, point_marker, iner_x_zero_line_width和plot_set_func参数。
        
        _cols = self.config["srnapointtrack"]["point_color"]
        def length2color(x):
            if x == 21:
                return(_cols[0])
            elif x == 22:
                return(_cols[1])
            elif x == 24:
                return(_cols[2])
            else:
                return(_cols[3])
        
        ax = self.ax
        reads = self.region_data
        
        max_count = self.config["srnapointtrack"]["max_count"]
        
        #如果没有序列
        if len(reads) == 0:
            self.set_ax_ylim(0, 1)
        else:
            if reads is not None:
                reads["length_color"] = reads["length"].map(length2color)
                reads["x"] = reads["start"]
                if self.config["plot_tpm"]:
                    reads["tpm"] = reads["counts"] * 1000000 / self.config["library_size"]
                    reads["y"] = reads["tpm"]
                else:
                    reads["y"] = reads["counts"]
                if self.config["split_strand"]:
                    reads["y"] = reads["y"].values * reads["is_rel_plus"].map(lambda x: 1 if x else -1)
                reads["size"] = self.config["point_size"]
                reads["marker"] = self.config["point_marker"]
                reads["edge_width"] = self.config["point_marker_edge_width"]
                point_have_edge_center = self.config["point_have_edge_center"]
                if self.config["plot_set_func"] is not None:
                    reads = self.config["plot_set_func"](reads)
                if self.config["is_unique"]:
                    reads = reads[reads["is_unique"]]
                reads = reads.sample(frac=1).reset_index(drop=True)            
                if max_count is not None:
                    l = reads["y"].values
                    l[l > max_count] = max_count
                    l[l < -max_count] = -max_count
                    reads["y"] = l
                if self.config["ylog"]:
                    if self.config["plot_tpm"]:
                        min1 = 1000000 / self.config["library_size"]
                    else:
                        min1 = 1
                    self.y_axis_min1 = min1
                    reads["y"] = self.to_log(reads["y"], min1)
        
                for (size, marker, edge_width), d in reads.groupby(["size", "marker", "edge_width"]):
                    if point_have_edge_center == 1:
                        facecolors = d["length_color"]
                        edgecolors = "none"
                    elif point_have_edge_center == 2:
                        facecolors = "none" 
                        edgecolors = d["length_color"]
                    else:
                        facecolors = d["length_color"]
                        edgecolors = d["length_color"]
                    ax.scatter(d["x"], d["y"], facecolors=facecolors, 
                               s=size, marker=marker, edgecolors=edgecolors,
                               linewidths=edge_width) #linewidth=0.3
                self.region_data = reads
            if reads is not None:
                min_y = reads["y"].min()
                max_y = reads["y"].max()
            else:
                min_y = 0
                max_y = 0
                
            if self.config["split_strand"]:
                self.set_ax_ylim(min_y, max_y)
            else:
                self.set_ax_ylim(0, max_y)