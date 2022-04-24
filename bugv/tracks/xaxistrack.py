from . track import BuGVTrack
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines

class XaxisGeneTrack(BuGVTrack):
    
    DEFAULT_CONFIG = {
        "ylim_style": None,
        "rel_position": None, #[start, end]
        "cut_region": [None, None], #[start, end]
        "xaxis_pos": 1,
        "xaxis_tick_pos": 0.9,
        "xtick_label_pos": 0.7,
        "xaxis_line_width": 0.5,
        "xtick_line_width": 0.5,
        "x_axis_font": 3,
        "xaxis_line_color": "black",
        "xtick_line_color": "black"
    }
    
    def plot_ax(self):
        self.plot_ax_x_axis()
    
    def _plot_x_axis(self, ax, tick_poss, tick_labels, 
                     config):
        
        tick_y_pos = config["xtick_label_pos"]
        tick_font_size = config['x_axis_font']
        xaxis_pos = config['xaxis_pos']
        xaxis_line_width = config["xaxis_line_width"]
        xaxis_line_color = config["xaxis_line_color"]
        xaxis_tick_pos = config["xaxis_tick_pos"]
        xtick_line_width = config["xtick_line_width"]
        xtick_line_color = config["xtick_line_color"]
        
        for tick_pos, tick_label in zip(tick_poss, tick_labels):
            ax.text(tick_pos, tick_y_pos, tick_label, horizontalalignment='center',
                fontsize=tick_font_size, clip_on=True)
        line = lines.Line2D((tick_poss[0], tick_poss[-1]), (xaxis_pos, xaxis_pos), lw=xaxis_line_width, color=xaxis_line_color)
        ax.add_line(line)
        for tick_pos in tick_poss:
            line = lines.Line2D((tick_pos, tick_pos), (xaxis_tick_pos, xaxis_pos), lw=xtick_line_width, color=xtick_line_color)
            ax.add_line(line)
        
    def plot_ax_x_axis(self):
        
        def generate_ticks(limit=[]):
            tmp_fig, tmp_ax = plt.subplots()
            tmp_ax.set_xlim(limit)
            plt.close()
            return tmp_ax.get_xticks()
        
        #calculate tick positions and tick labels
        rel_position = self.config["rel_position"]
        cut_region = self.config["cut_region"]
        if not rel_position:
            rel_position = self.ax.get_xlim()
        rel_start, rel_end = rel_position
        if rel_end >= rel_start:
            rel_strand = "+"
        else:
            rel_strand = "-"
        tick_x,  tick_y = cut_region
        if tick_x is None:
            tick_x = 1
        if tick_y is None:
            tick_y = abs(rel_end - rel_start) + 1
        rel_tick_poss = generate_ticks([tick_x, tick_y])
        if rel_strand == "+":
            tick_poss = rel_tick_poss + rel_start - 1
        else:
            tick_poss = rel_start - rel_tick_poss + 1
        tick_labels = ["{:.0f}".format(i) for i in rel_tick_poss]
        
        #plot axis
        self._plot_x_axis(self.ax, tick_poss, tick_labels, self.config)

#Not used now
class XaxisTrack(BuGVTrack):
    
    #config:
    DEFAULT_CONFIG = {
        "xaxis_y_bottom_pos": 0.2,
        "xaxis_y_top_pos" : 0.8,
        "where" : "bottom",
        "ticks_number": 3,
        "chr_name_y_pos": 0.3,
        "is_axisartist": True, #固定值，请勿更改
        "ylim_style": None,
        "rel_position": None, #[start, end]
    }
        
    def plot_ax(self):
        #most from pygenometracks: class XAxisTrack
        ax = self.ax
        ax.patch.set_visible(False)
        if self.config['where'] == 'top':
            ax.axis["x"] = ax.new_floating_axis(0, self.config["xaxis_y_bottom_pos"])
            ax.axis["x"].set_axis_direction("top")
            vert_align = 'top'
        else:
            ax.axis["x"] = ax.new_floating_axis(0, self.config["xaxis_y_top_pos"])
            vert_align = 'bottom'
        ax.axis["x"].major_ticks.set_tick_out(True)
        
        #self.format_tick() #暂不格式化
        ax.axis["x"].major_ticklabels.set(size=self.config['x_axis_font'])
        if self.config["rel_position"]:
            self.x_tick_to_rel_position()
        else:
            ax.text(0.5, self.config["chr_name_y_pos"], self.chr_name, horizontalalignment='center',
                    fontsize=self.config['x_axis_font'],
                    verticalalignment=vert_align, transform=ax.transAxes)
        
    def x_tick_to_rel_position(self):
        
        def generate_ticks(limit=[]):
            tmp_fig, tmp_ax = plt.subplots()
            tmp_ax.set_xlim(limit)
            plt.close()
            return tmp_ax.get_xticks()
        
        rel_position = self.config["rel_position"]
        if rel_position:
            rel_start, rel_end = rel_position
            if rel_end >= rel_start:
                rel_strand = "+"
            else:
                rel_strand = "-"
            rel_tick_poss = generate_ticks([1, abs(rel_end - rel_start) + 1])
            if rel_strand == "+":
                tick_poss = rel_tick_poss + rel_start - 1
            else:
                tick_poss = rel_start - rel_tick_poss + 1
            print(self.ax.get_xlim())
            self.ax.axis["x"].axis.set_ticks(tick_poss)
            print(self.ax.get_xlim())
            self.ax.axis["x"].axis.set_ticklabels(["{:.0f}".format(i) for i in rel_tick_poss])
        
    def format_tick(self):
        #most from pygenometracks: class XAxisTrack
        #暂不使用
        ax = self.ax
        ticks = ax.get_xticks()
        
        ticks = np.linspace(ticks[0], ticks[-1], self.config["ticks_number"])
        if ticks[-1] - ticks[0] <= 1e3:
            labels = ["{:,.0f}".format((x))
                      for x in ticks]
            labels[-2] += " b"

        elif ticks[-1] - ticks[0] <= 4e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"
        ax.axis["x"].axis.set_ticks(ticks)
        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')
        ax.axis["x"].major_ticklabels.set(size=self.config['x_axis_font'])
