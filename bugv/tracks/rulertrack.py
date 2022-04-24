from . track import BuGVTrack
import numpy as np


class RulerTrack(BuGVTrack):
    
    #config:
    DEFAULT_CONFIG = {
        "xaxis_y_bottom_pos": 0.2,
        "xaxis_y_top_pos" : 1,
        "where" : "bottom",
        "ticks_number": 3,
        "chr_name_y_pos": 0.3,
        "ylim_style": None,
        "ruler" : {
            "right_space_ratio": 0.1,
            "ruler_length": 1000,
            "ruler_y_pos": 0.3,
            "tick_height": 0.1, # set minus to invert
            "label_height": 0.3, # set minus to invert
            "color": "black",
            "linewidth": 1
        }
    }
    
    DEFAULT_RULERS = np.array([1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000] )
    
    def plot_ax(self):            
        
        #most from pygenometracks: class XAxisTrack
        ax = self.ax
        
        ls_config = self.config["ruler"]
        right_space_ratio = ls_config["right_space_ratio"]
        ruler_length = ls_config["ruler_length"]
        ruler_y_pos = ls_config["ruler_y_pos"]
        tick_height = ls_config["tick_height"]
        label_height = ls_config["label_height"]
        color = ls_config["color"]
        linewidth = ls_config["linewidth"]
        size=self.config["y_axis_font"]
        tick_y_pos = ruler_y_pos + tick_height
        label_y_pos = ruler_y_pos + label_height
        
        xlim_min, xlim_max = ax.get_xlim()
        xlim_space = np.abs(xlim_max - xlim_min)
        ruler_max_length = xlim_space * (1-right_space_ratio)
        
        if ruler_max_length < ruler_length:
            ruler_length = self.DEFAULT_RULERS[(self.DEFAULT_RULERS < ruler_max_length)][-1]
        
        if xlim_min < xlim_max:
            ruler_x_max = xlim_max - right_space_ratio * xlim_space
            ruler_x_min = ruler_x_max - ruler_length
        else: # axis is inverted
            ruler_x_min = xlim_max + right_space_ratio * xlim_space
            ruler_x_max = ruler_x_min + ruler_length
        
        x_pos = [ruler_x_min, ruler_x_min, ruler_x_max, ruler_x_max]
        y_pos = [tick_y_pos, ruler_y_pos, ruler_y_pos, tick_y_pos]
        
        ax.plot(x_pos, y_pos, color=color, linewidth=linewidth)
        ax.text((ruler_x_min + ruler_x_max)/2, 
                label_y_pos,
                ruler_length,
                horizontalalignment='center',
                #verticalalignment='top',
                size=size)
        self.set_ax_ylim(0,1)