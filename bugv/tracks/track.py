import pysam
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib import patches
import mpl_toolkits.axisartist as axisartist
import matplotlib.lines as lines
import bugv.utilities 
import os
import pickle

"""
version 0.2.1 by jiajinbu 2021.01.04
"""

class BuGVTrack:
    
    """
    请查看BuGV文档。
    1. 绘图参数
    每个track的绘图参数都在该track对象的config属性里。
    
    可使用load_config(config)方法设置参数。注意这里面如果config元素值是字典的话，会先向下遍历替换。而不是直接替换为新的字典。
    由于是直接复制。config的值是列表的话，将会进行浅拷贝再传值。但是仍存在对对象，以及深拷贝才能拷贝的值有可能共享内存的问题。
    
    在新建一个track时，必须调用父类BuGVTrack的__init__方法。
    如`super(XaxisTrack, self).__init__(*args, **kwargs)`。
    这时会首先加载BuGVTrack类的ALL_DEFAULT_CONFIG默认参数。
    然后再加载DEFAULT_CONFIG参数。
    DEFAULT_CONFIG可以在子类中定义为类属性。
    注意所有对象共用类属性，因此禁止通过对象更改ALL_DEFAULT_CONFIG
    和DEFAULT_CONFIG。
    
    后续可以添加locked参数。将一些参数locked。防止后续更改。
    
    2. 每个BuGVTrack对象有一个bgv属性。存储了其所属BuGV。
    而每个BuGV对象又含有个track_objs对象，存储了其所有的BuGVTrack对象。
    可以通过BuGV的load_config(config)方法设置其所属的所有对象的参数。
    
    3. y坐标
    每个track的plot_ax用于画图。该函数内部通常要在合适位置（一般在比较靠后位置）
    调用set_ax_ylim方法。该方法分两步：
    第一步获得想用的ylim，即rel_ylim_min和rel_ylim_max，
    将它们存入同名属性。这可以通过config["fix_ylim"]来设置。如果config["fix_ylim"]为None,
    则会调用set_rel_ylim方法。该方法默认是获得当前ax的ylim。如果想更改该方法，可以再plot函数里直接
    先产生rel_ylim_min和rel_ylim_max，再调用set_ax_ylim(rel_ylim_min, rel_ylim_max)。注意更改后
    仍不是实际的的rel_ylim_min和rel_ylim_max，还根据config["ylim_style"]设定的方法进行修改。
    
    第二步调用_set_ax_ylim方法设置真正的ylim_min和ylim_max。该方法会使用config["y_margin_ratio"]的
    参数，它指轴上下两部分分别添加多大margin比例。会根据该值，将rel_ylim_min和rel_ylim_max缩放，
    并获得程序实际使用的ylim。
    因此这里有两套坐标。一套是ax实际的ylim，一套是我们关心的rel_ylim。它是去掉上下y_margin_ratio的ylim。
    
    4. 绘制y轴
    调用plot_y_axis在left_ax上绘制y轴。实际上是在rel_ylim除画线和添加标记。
    """
    
    ALL_DEFAULT_CONFIG = {"y_margin_ratio" : (0.05, 0.05),
                          "y_axis_line_width" : 0.5,
                          "y_axis_line_color" : "black",
                          "max_cache_obj_number" : 0,
                          "min_y_max": 1,
                          "fix_ylim" : None,
                          "ylim_style": "change_less", # None not plot, change_down, change_up, not_change
                          "ylim_plot_style": None, #None, "not_minus", "not_plus"
                          "ylim_plot_ratio": 1,
                          "ylim_plot_lim": None,
                          "ylim_up_down_ratio": 1, #None
                          "y_axis_label_valign_center" : False,
                          "y_axis_ticks": None,
                          "ylog": False,
                          "x_base_pad": 0.5,
                          "point_size" : 1,
                          "point_marker" : "o",
                          "point_marker_edge_width" : 1,
                          "point_have_edge_center" : 2, #1: no edge, 2: no center, 3: have edge and center 
                          "x_axis_font" : 6,
                          "y_axis_font" : 6,
                          "x_axis_line_color": "grey",
                          "y_tick_round_ndigits": None, #round(x, y_tick_round_ndigits)
                          "track_height" : 1,
                          "is_axisartist" : False,
                          
                          "null_border": True,
                          
                          "y_zero_line": { #每个track如果画y=0的x轴时用的参数
                              "plot": False, #是否绘制
                              "linewidth": 0.5, 
                              "zorder": 0,
                              "color": "grey",
                              "linestyle": '-'
                          },
                          
                          #for obj
                          "file_type" : None,
                          "file_name" : "",
                          "read_obj_kwargs": {},
                          "load_obj_region_kwargs": {},
                          "obj" : None,
                          "data" : None,
                          "fix_data" : None,
                          "region_data" : None,
                          "select_features" : None,
                          "keep_select_features": False,
                         } # (ylim_min, ylim_max)
    DEFAULT_CONFIG = {}
    
    
    def __init__(self, config={}, track_name=""):
        #设置默认属性
        self.rel_ylim_max = None
        self.rel_ylim_min = None
        self.bgv = None
        self.xaxis_is_invert = False
        self.select_features = None
        self.track_name = track_name
        
        #加载绘图参数
        self.config = {}
        self.load_config(self.ALL_DEFAULT_CONFIG)
        self.load_config(self.DEFAULT_CONFIG)
        self.init()
        self.load_config(config)
        self.load_obj()
        self.open()
        
    def init(self):
        pass
        
    def close(self):
        if self.obj is not None:
            try:
                self.obj.close()
            except:
                pass
        
    def copy(self):
        pass
    
    def get_data(self):
        pass
    
    def get_parm_func(self, parm_dict_list):
        #先找第一个元素所指向的字典里是否含有或是否为None，按顺序寻找配置，找不到返回None
        def _get_parm(parm):
            for parm_dict in parm_dict_list:
                value = parm_dict.get(parm)
                if value is not None:
                    return value
            return None
            
        return(_get_parm)
    
    def load_config(self, config):
        bugv.utilities.update_dict(self.config, config)
    
    def _modify_load_region_kwargs(self, kwargs):
        return kwargs
    
    def load_region_data(self):
        if self.config["fix_data"] is not None:
            self.region_data = self.config["fix_data"]
        elif self.obj is not None:
            if self.config["keep_select_features"]:
                ls_kwargs = self.config["load_obj_region_kwargs"].copy()
                ls_kwargs["select_features"] = self.config["select_features"]
            else:
                ls_kwargs = self.config["load_obj_region_kwargs"]
            ls_kwargs["region"] = self.region
            ls_kwargs = self._modify_load_region_kwargs(ls_kwargs)
            self.obj.load_region_data(**ls_kwargs)
            self.region_data = self.obj.region_data
    
    def _modify_load_obj_kwargs(self, kwargs):
        return kwargs
    
    def clear_cache(self):
        self._obj_cache = {}
    
    def load_obj(self):
        max_cache_obj_number = self.config["max_cache_obj_number"]
        if not hasattr(type(self), "_obj_cache"):
            _obj_cache = type(self)._obj_cache = {}
        else:
            _obj_cache = self._obj_cache 

        need_creat_obj = self.config["file_name"] or self.config["data"] is not None
        if need_creat_obj:
            ls_kwargs = self.config["read_obj_kwargs"].copy()
            if self.config["file_name"]:
                file_name = ls_kwargs["file_name"] = self.config["file_name"]
            if self.config["data"] is not None:
                ls_kwargs["data"] = self.config["data"]
            ls_kwargs = self._modify_load_obj_kwargs(ls_kwargs)
            if self.config["file_name"]:
                file_path = os.path.abspath(self.config["file_name"])
                if file_path in _obj_cache:
                    self.obj = _obj_cache[file_path]
                else:
                    if self.config["file_name"].endswith(".pickle"):
                        with open(self.config["file_name"], 'rb') as f:
                            self.obj = pickle.load(f)
                    else:
                        self.obj = self.config["file_type"](**ls_kwargs)
                    if max_cache_obj_number:
                        if len(_obj_cache) >= max_cache_obj_number:
                            _obj_cache.pop(list(_obj_cache.keys())[0])
                        _obj_cache[file_path] = self.obj
            else:
                self.obj = self.config["file_type"](**ls_kwargs)                
        elif self.config["obj"] is not None:
            self.obj = self.config["obj"]
        else:
            self.obj = None
    
    def to_log(self, n, min1=None):
        x = np.abs(n)
        if min1 is None:
            min1 = min(x[x>0])
        x = x/min1
        x[x>0] = np.log2(x[x>0])
        x = x * np.sign(n)
        return x
    
    def move(self):
        pass
            
    def open(self):
        if self.obj is not None:
            self.obj.close()
            self.obj.open()
    
    def plot_y_zero_line(self):
        if self.config["y_zero_line"]["plot"]:
            ls_config = self.config["y_zero_line"].copy()
            del ls_config["plot"]
            self.ax.axhline(y=0, **ls_config)
    
    def plot_border(self):
        ax = self.ax
        if self.config["is_axisartist"]:
            if self.config["null_border"]:
                ax.axis[:].set_visible(False)
        else:
            if self.config["null_border"]:
                ax.spines['right'].set_color('none')
                ax.spines['top'].set_color('none')
                ax.spines['bottom'].set_color('none')
                ax.spines['left'].set_color('none')
                ax.set_xticks([])
                ax.set_yticks([])
            
    def plot(self):
        self.load_region_data()
        self.get_data()
        self.plot_border()
        self.plot_ax()
        self.plot_y_zero_line()
        
        self.plot_left_ax() #don't before plot_ax()
        self.plot_right_ax()
        
    def plot_ax(self):
        pass
    
    def plot_left_ax(self):
        self.plot_y_axis()
    
    def plot_right_ax(self):
        pass
    
    def plot_y_axis(self):
        ylim_style = self.config["ylim_style"]
        ylim_plot_style = self.config["ylim_plot_style"] #None, "not_minus", "not_plus"
        ylim_plot_ratio = self.config["ylim_plot_ratio"]
        ylim_plot_lim = self.config["ylim_plot_lim"]
        
        if ylim_style is not None:
            ax = self.left_ax
            
            #get tick postion (tick_min, tick_max)
            if ylim_plot_lim is not None:
                tick_min, tick_max = ylim_plot_lim
            else:
                tick_min = self.rel_ylim_min
                tick_max = self.rel_ylim_max
                if ylim_plot_style == "not_minus":
                    tick_min = 0
                elif ylim_plot_style == "not_plus":
                    tick_max = 0
                tick_min = tick_min * ylim_plot_ratio
                tick_max = tick_max * ylim_plot_ratio
                y_tick_round_ndigits = self.config["y_tick_round_ndigits"]
                tick_min = round(tick_min, y_tick_round_ndigits)
                if not y_tick_round_ndigits: tick_min = int(tick_min)
                tick_max = round(tick_max, self.config["y_tick_round_ndigits"])
                if not self.config["y_tick_round_ndigits"]: tick_max = int(tick_max)
            
            if self.config["ylog"]:
                tick_min_pos, tick_max_pos = self.to_log([tick_min, tick_max], self.y_axis_min1)
            else:
                tick_min_pos, tick_max_pos = tick_min, tick_max
        
            tick_right_pos = 0.9
            tick_left_pos = 0.75
            tick_label_pos = 0.65
            x_pos = [tick_left_pos, tick_right_pos, tick_right_pos, tick_left_pos]
            y_pos = [tick_min_pos, tick_min_pos, tick_max_pos, tick_max_pos]
            ax.plot(x_pos, y_pos, color=self.config["y_axis_line_color"], linewidth=self.config["y_axis_line_width"])
            
            if self.config["y_axis_label_valign_center"]:
                tick_min_valign = "center"
                tick_max_valign = "center"
            else:
                tick_min_valign = 'bottom'
                tick_max_valign = "top"
        
            ax.text(tick_label_pos, tick_min_pos, tick_min, verticalalignment=tick_min_valign, horizontalalignment='right', size=self.config["y_axis_font"])
            ax.text(tick_label_pos, tick_max_pos, tick_max, verticalalignment=tick_max_valign, horizontalalignment='right', size=self.config["y_axis_font"])
        
            if tick_min*tick_max < 0:
                ax.plot([tick_left_pos, tick_right_pos], [0, 0], color=self.config["y_axis_line_color"], linewidth=self.config["y_axis_line_width"])
                ax.text(tick_label_pos, 0, 0, verticalalignment="center", horizontalalignment="right", size=self.config["y_axis_font"])
        
            ax.set_xlim(0,1)
            ax.patch.set_visible(False)    
    
    def _set_ax_ylim(self):
        rel_ylim_space = self.rel_ylim_max - self.rel_ylim_min
        self.ylim_min = self.rel_ylim_min - rel_ylim_space * self.config["y_margin_ratio"][0]
        self.ylim_max = self.rel_ylim_max + rel_ylim_space * self.config["y_margin_ratio"][1]
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.left_ax.set_ylim(self.ylim_min, self.ylim_max)
    
    def set_ax_ylim(self, rel_ylim_min=None, rel_ylim_max=None):
        """
        """
        ylim_style = self.config["ylim_style"]      
        if self.config["fix_ylim"] is not None:
            rel_ylim_min, rel_ylim_max = self.config["fix_ylim"]
        else:
            self.set_rel_ylim(rel_ylim_min, rel_ylim_max)
            if ylim_style != "not_change":
                up = max(self.rel_ylim_max, self.rel_ylim_min)
                down = min(self.rel_ylim_max, self.rel_ylim_min)
                ylim_up_down_ratio = self.config["ylim_up_down_ratio"]
                if up > 0 and down < 0 and (ylim_up_down_ratio is not None):
                    down = -down
                    if ylim_style == "change_less":
                        if up > down:
                            down = ylim_up_down_ratio * up
                        else:
                            up = ylim_up_down_ratio * down
                    elif ylim_style == "change_up":
                        up = ylim_up_down_ratio * down
                    elif ylim_style == "change_down":
                        down = ylim_up_down_ratio * up
                    down = -down
                    rel_ylim_min, rel_ylim_max = down, up
                    self.set_rel_ylim(rel_ylim_min, rel_ylim_max)
        self.set_rel_ylim(rel_ylim_min, rel_ylim_max)
        self._set_ax_ylim()
    
    def set_ax_ylim(self, rel_ylim_min=None, rel_ylim_max=None):
        """
        """
        
        def set_rel_ylim(rel_ylim_min=None, rel_ylim_max=None):
            if rel_ylim_max is None or np.isnan(rel_ylim_max):
                if self.rel_ylim_max is None:
                    self.rel_ylim_max = self.ax.get_ylim()[1]
            else:
                self.rel_ylim_max = rel_ylim_max
            if rel_ylim_min is None or np.isnan(rel_ylim_max):
                if self.rel_ylim_min is None:
                    self.rel_ylim_min = self.ax.get_ylim()[0]
            else:
                self.rel_ylim_min = rel_ylim_min
            if self.rel_ylim_min == self.rel_ylim_max and self.rel_ylim_max == 0:
                self.rel_ylim_max = self.config["min_y_max"]

        ylim_style = self.config["ylim_style"]      
        if self.config["fix_ylim"] is not None:
            rel_ylim_min, rel_ylim_max = self.config["fix_ylim"]
        else:
            set_rel_ylim(rel_ylim_min, rel_ylim_max)
            if ylim_style != "not_change":
                up = max(self.rel_ylim_max, self.rel_ylim_min)
                down = min(self.rel_ylim_max, self.rel_ylim_min)
                ylim_up_down_ratio = self.config["ylim_up_down_ratio"]
                if up >= 0 and down <= 0 and (ylim_up_down_ratio is not None):
                    down = -down
                    if ylim_style == "change_less":
                        if up > down:
                            down = ylim_up_down_ratio * up
                        else:
                            up = ylim_up_down_ratio * down
                    elif ylim_style == "change_up":
                        up = ylim_up_down_ratio * down
                    elif ylim_style == "change_down":
                        down = ylim_up_down_ratio * up
                    down = -down
                    rel_ylim_min, rel_ylim_max = down, up
                    set_rel_ylim(rel_ylim_min, rel_ylim_max)
        set_rel_ylim(rel_ylim_min, rel_ylim_max)
        self._set_ax_ylim()
    
    def set_region(self, region, select_features=None, rel_position=None, cut_region=[None, None]):
        chr_name, region_start, region_end, region_strand = region
        self.region = region
        self.chr_name = chr_name
        self.region_start = region_start
        self.region_end = region_end
        self.region_strand = region_strand
        self.select_features = select_features
        if not self.config["keep_select_features"]:
            self.config["select_features"] = select_features
        self.config["rel_position"] = rel_position
        self.config["cut_region"] = cut_region
    
    def set_rel_ylim1(self, rel_ylim_min=None, rel_ylim_max=None):
        if rel_ylim_max is None or np.isnan(rel_ylim_max):
            if self.rel_ylim_max is None:
                self.rel_ylim_max = self.ax.get_ylim()[1]
        else:
            self.rel_ylim_max = rel_ylim_max
        if rel_ylim_min is None or np.isnan(rel_ylim_max):
            if self.rel_ylim_min is None:
                self.rel_ylim_min = self.ax.get_ylim()[0]
        else:
            self.rel_ylim_min = rel_ylim_min
        if self.rel_ylim_min == self.rel_ylim_max and self.rel_ylim_max == 0:
            self.rel_ylim_max = self.config["min_y_max"]
            
    
            
    