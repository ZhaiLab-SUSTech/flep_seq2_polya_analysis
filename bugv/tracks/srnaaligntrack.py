from . track import BuGVTrack
from bugv.readers.srnabam import sRNABam

"""
unfinished
"""

class sRNAAlignTrack(BuGVTrack):
    
    DEFAULT_CONFIG = {"file_type": sRNABam,
                      "file_name": "",
                      "y_zero_line": {"plot": True},
                      "read_obj_kwargs": {"fileformat": "shortstack"},
                      "sRNA_length_color": ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3'], #21nt, 22nt, 24nt, 其他小RNA的颜色
                      "strand_color": ['#E41A1C', '#377EB8'],
                      "read_order_method": "insert" , #insert: IGV-like #stacked 简单堆叠 和 most_down 最高的放最下面，暂不用，未开发完成
                      "color_mode": "sRNA_lenght", #sRNA_length或者plus_minus
                      "ylim_style": "not_change",
                      } #注意所有对象共享该属性
            
    def plot_ax(self):
        #使用参数。
        
        def length2color(x):
            cols = self.config["sRNA_length_color"]
            if x == 21:
                return(cols[0])
            elif x == 22:
                return(cols[1])
            elif x == 24:
                return(cols[2])
            else:
                return(cols[3])
        
        ax = self.ax
        reads = self.region_data
        
        #如果没有序列
        if len(reads) == 0:
            self.set_ax_ylim(0, 1)
        else:
            if self.config["read_order_method"] == "insert":
                reads = self.read_order_method_igv_like(reads, self.region_strand)
            elif self.config["read_order_method"] == "stacked":
                reads = self.read_order_method_stacked(reads)

            if self.config["color_mode"] == "sRNA_lenght":
                reads["color"] = reads["length"].map(length2color)
            elif self.config["color_mode"] == "plus_minus":
                reads["color"] = reads["is_rel_plus"].map(lambda x: self.config["strand_color"][0] if x else self.config["strand_color"][1])
            reads["x_left"] = reads["start"]
            reads["x_right"] = reads["end"]

            codes = [
                Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
            ]

            for index, row in reads.iterrows():
                x_start = row["x_left"]
                x_end = row["x_right"]
                y_up = row["y_up"]
                y_lower = row["y_down"]
                col = row["color"]
                verts = [
                    (x_start, y_lower),  # left, bottom
                    (x_start, y_up),  # left, top
                    (x_end, y_up),  # right, top
                    (x_end, y_lower),  # right, bottom
                    (0., 0.)  # ignored
                ]
                path = Path(verts, codes)
                patch = patches.PathPatch(path, facecolor=col, lw=0)
                ax.add_patch(patch)
            self.region_data = reads
            self.set_ax_ylim(0, reads.y_up.max() + 1)
                
    def read_order_method_stacked(self, reads):
        #y可能需要调整+1或-1
        reads["y_up"] = reads.counts.cumsum()
        reads["y_lower"] = [0] + reads["y_up"][:-1].to_list()
        return(reads)
        
    def read_order_method_most_down(self, reads):
        #y可能需要调整+1或-1
        sort_by_position = True
        if sort_by_position:
            reads = reads.sort_values(["start", "end", "strand"])
        else:
            reads = reads.sort_values("counts", ascending=False)

        min_x = reads["start"].min()
        max_x = reads["end"].max()
        have_regions = np.zeros(max_x - min_x + 1)

        y_downs = []
        for index, row in reads.iterrows():
            start = row["start"]
            end = row["end"]
            count = row["counts"]
            y_down = have_regions[(start-min_x):(end-min_x+1)].max()
            y_up = y_down + count
            y_downs.append(y_down)
            have_regions[(start-min_x):(end-min_x+1)] = np.maximum(have_regions[(start-min_x):(end-min_x+1)], y_up)

        reads["y_down"] = y_downs
        reads["y_up"] = reads["y_down"] + reads["counts"]
        return(reads)
    
    def read_order_method_igv_like(self, reads, strand="+"):
        #method3
        #只要保证ls_start和ls_end都大于0，且ls_start小于ls_end即可。按照ls_start和ls_end进行排序
        if strand == "+":
            reads["ls_start"] = reads["start"]
            reads["ls_end"] = reads["end"]
        else:
            max_x = reads["end"].max()
            reads["ls_start"] = max_x - reads["end"] + 1
            reads["ls_end"] = max_x - reads["start"] + 1
            reads = reads.sort_values(["ls_start", "ls_end"])
            
        y_downs = []
        have_overlap_regions = []
        now_max_x = 0
        for index, row in reads.iterrows():
            start = row["ls_start"]
            end = row["ls_end"]
            count = row["counts"]
            if start > now_max_x:
                y_down = 1
                y_up = count
                now_max_x = end
                have_overlap_regions = [[1, y_up, now_max_x]]
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
                    if x3 < start and (x2 - x1 + 1) >= count:
                        have_insert = True
                        y_down = x1
                        y_up = y_down + count - 1
                        if y_up != x3:
                            d[0] = y_up + 1
                            have_overlap_regions.insert(i, [x1, y_up, end]) #可以尝试与上一个元素合并，不过意义不大，以为下一轮还要合并
                        else:
                            d[2] = end
                        break
                if not have_insert:
                    y_down = have_overlap_regions[-1][1] + 1 if have_overlap_regions else 1
                    y_up = y_down + count - 1
                    have_overlap_regions.append([y_down, y_up, end])
                if end > now_max_x:
                    now_max_x = end
            y_downs.append(y_down)

        reads["y_down"] = y_downs
        reads["y_down"] = reads["y_down"] - 1
        reads["y_up"] = reads["y_down"] + reads["counts"]  
        return(reads)
