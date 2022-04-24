from . track import BuGVTrack
from bugv.readers.bsmapmethratio import BsmapMethratioTbi
import numpy as np
import matplotlib

class MethylationTrack(BuGVTrack):
    
    DEFAULT_CONFIG = {"file_type": BsmapMethratioTbi,
                      "file_name": "",
                      "read_obj_kwargs": {},
                      "y_zero_line": {"plot": True},
                      "ylim_style" : "not_change",
                      "methylationtrack": {
                          "context": "CG", #"CG", "CHG", "CHH"
                          "colors": ['#E41A1C', '#377EB8', '#4DAF4A'],
                          "color": None,
                          "zero_add": 0.02,
                          "linewidth": 0.1,
                          "bin_length": 1,
                          "min_ct_count": 1
                      }
                      } #注意所有对象共享改属性
    
    def _modify_load_region_kwargs(self, kwargs):
        kwargs["bin_length"] = self.config["methylationtrack"]["bin_length"]
        kwargs["min_ct_count"] = self.config["methylationtrack"]["min_ct_count"]
        return kwargs
    
    def plot_ax(self):
        
        self.set_ax_ylim(0, 1)
        ls_config = self.config["methylationtrack"]
        context = ls_config["context"]
        context_index = {"CG":0, "CHG":1, "CHH": 2}[context]
        df = self.region_data[context_index]
        color = ls_config["color"]
        linewidth = ls_config["linewidth"]
        if color is None:
            color = ls_config["colors"][context_index]
        zero_add = ls_config["zero_add"]
        ax = self.ax

        x = df["pos"].values
        y = df["ratio"].values
        y[y<zero_add] = zero_add
        lines = [(i[:2], i[2:]) for i in np.array([x, np.zeros(len(x)), x, y]).T]
        ln_coll = matplotlib.collections.LineCollection(lines, color=color, lw=linewidth)
        ax.add_collection(ln_coll)