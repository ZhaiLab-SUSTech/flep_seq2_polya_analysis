from . track import BuGVTrack
from bugv.readers.rnabam import RNABam
import numpy as np
from scipy import interpolate
from matplotlib.patches import Polygon

def spline(x, y, n=100):
    
    #not use for now
    tck,u     = interpolate.splprep( [x,y] ,s = 0 )
    xnew,ynew = interpolate.splev( np.linspace( 0, 1, n ), tck,der = 0)
    return(xnew, ynew)
    
def bspline(x, y, n=100, degree=3):
    
    #https://stackoverflow.com/questions/28279060/splines-with-python-using-control-knots-and-endpoints
    #https://stackoverflow.com/users/1429402/fnord
    
    cv = np.asarray([x, y]).T
    count = cv.shape[0]
    
    # Prevent degree from exceeding count-1, otherwise splev will crash
    degree = np.clip(degree,1,count-1)
    
    # Calculate knot vector
    kv = np.array([0]*degree + list(range(count-degree+1)) + [count-degree]*degree,dtype='int')
    
    # Calculate query range
    u = np.linspace(0,(count-degree),n)

    # Calculate result
    xnew, ynew = np.array(interpolate.splev(u, (kv,cv.T,degree)))
    return (xnew, ynew)

def extract_true_blocks(a):
    """
    a is bool list or np.array, for example:
    np.array([True, False, True, True, False, True, False, False, True])
    return True regions:
    
    """
    #add one False in both side of a
    paded_a = np.pad(a, (1,1))
    absdiff = np.diff(paded_a)
    ranges = np.where(absdiff)[0].reshape(-1, 2)
    return ranges
    
def extract_continous_blocks(a):
    """
    a must be np.array.
    a need to support np.diff(a). numeric or bool.
    
    return:
    (np.array([[i1, j1], [i2, j2], ...]), np.array([v1, v2, ...]))
    
    for example:
    a = np.array([1, 1, 3, 4, 4])
    extract_continous_blocks(a)
    result:
    (array([[0, 2],
            [2, 3],
            [3, 5]]), 
     array([1, 3, 4]))
    """
    
    ranges = np.append(np.where(np.insert(np.diff(a), 0, True))[0], len(a))

    if ranges[0] != 0:
        ranges = np.insert(ranges, 0, 0)
    if ranges[-1] != len(a):
        ranges = np.append(ranges, len(a))
        
    return((np.array([ranges[:-1], ranges[1:]]).T, a[ranges[:-1]]))

class RNASashimiplot(BuGVTrack):
    
    DEFAULT_CONFIG = {"plot_set_func" : None,
                      "file_type": RNABam,
                      "file_name": "",
                      "point_color": ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3'], #21nt, 22nt, 24nt, 其他小RNA的颜色
                      "ylim_style" : "not_change",#"change_less",
                      "read_specific_strand": True, #whether strand specific
                      "is_FR": False, #strand specific info
                      "filter_unique": False,
                      "filter_func": None,
                      "min_overlap": 0,
                      "plot_junc": True,
                      "y_zero_line": {"plot": True},
                      "rnasashimiplot": {
                          "all_color":  '#984EA3', #brewer.pal(9, "Set1")[5]
                          "plus_color": '#984EA3', #brewer.pal(9, "Set1")[4]
                          "minus_color": '#FF7F00',
                          "density_line_width": 0.5,
                          "junction_line_width": 1,
                          "min_junction_ratio": 0,
                          "min_junction_count": 0,
                          "line_width_based_on_count": False,
                          "max_junction_line_width": 1,
                          "min_junction_line_width": 0.1,
                          #"minus_y_pad_ratio": 0
                      }
                      } #注意所有对象共享改属性                      
                      
    def _modify_load_obj_kwargs(self, kwargs):
        kwargs["is_FR"] = self.config["is_FR"]
        kwargs["filter_unique"] = self.config["filter_unique"]
        kwargs["filter_func"] = self.config["filter_func"]
        kwargs["read_specific_strand"] = self.config["read_specific_strand"]
        kwargs["min_overlap"] = self.config["min_overlap"]
        return kwargs
        
    def plot_ax(self):
        
        MIN_YLIM = 1
        
        def plot_depth(depth_data, ax, y_is_minus=False, color="purple", x_base_pad=0.5, density_line_width=0.5):
            
            x = depth_data["pos"].values
            y = depth_data["depth"].values
            depth_more0_regions = extract_true_blocks(y>0)
            
            for i, j in depth_more0_regions:
                x0_start = x[i]
                y0 = y[i:j]
                y0ijs, y0s = extract_continous_blocks(y0)
                xs = np.pad(y0ijs.flatten()-x_base_pad, (1,1), "edge") + x0_start
                ys = np.pad(np.repeat(y0s, 2), (1,1))
                if y_is_minus:
                    ys = -ys
                xys = np.array([xs, ys]).T
                polygon = Polygon(xys, facecolor=color, edgecolor=color, lw=density_line_width)
                ax.add_patch(polygon)
        
        def _filter_cal_junc(junc_data):
            #free variable: min_junction_ratio, min_junction_count
            r = junc_data
            if min_junction_ratio:
                r = junc_data[junc_data["ratio"]>=min_junction_ratio] 
            if min_junction_count:
                r = junc_data[junc_data["count"]>=min_junction_count]
            return r
            
        def _set_parm_junc(junc_data):
            #free variable:  line_width_based_on_count
            #junction_line_width, max_junction_line_width, min_junction_line_width
            junc_data = junc_data.copy()
            if line_width_based_on_count:
                junc_data["lw"] = (min_junction_line_width + 
                                        (max_junction_line_width-min_junction_line_width) * junc_data["ratio"]
                                  )
            else:
                junc_data["lw"] = junction_line_width
            return junc_data
        
        def plot_junc_core(junc_data, ax, ymax, y_is_minus=False, color="purple", xaxis_lenght=0):
                                                    
            for i, row_data in junc_data.iterrows():
                xleft = row_data["start"]
                xright = row_data["end"]
                yleft = row_data["left_depth"]
                yright = row_data["right_depth"]
                count = row_data["count"] #not used for now
                lw = row_data["lw"]

                xmiddle = (xright + xleft)/2
                xleft_dis = xmiddle - xleft
                xright_dis = xright - xmiddle
                xleft_right_dis = xright - xleft
                
                add_height__ratio = 0.1
                if xaxis_lenght:
                    junc_length_ratio = (xright - xleft) / xaxis_lenght
                    if junc_length_ratio >= 0.3:
                        add_height__ratio = 0.3
                    elif junc_length_ratio >= 0.2:
                        add_height__ratio = 0.2
                    elif junc_length_ratio >= 0.1:
                        add_height__ratio = 0.15
                
                ymiddle = max([yleft, yright]) + add_height__ratio*ymax
                yleft_dis = ymiddle - yleft
                yright_dis = yright - ymiddle

                x0 = np.array([xleft, xleft + 0.3*xleft_dis, xleft+0.8*xleft_dis, xmiddle, 
                      xmiddle + 0.2*xright_dis, xmiddle+0.7*xleft_dis, xright])
                y0 = np.array([yleft, yleft + 0.7*yleft_dis, yleft + 0.98*yleft_dis, ymiddle, 
                      ymiddle + 0.02*yright_dis, ymiddle + 0.3*yright_dis, yright])
                if y_is_minus:
                    y0 = -y0
                x, y = bspline(x0, y0) #spline(x0, y0) not used
                ax.plot(x, y, zorder=0, color=color, lw=lw)
                
        def plot_junc(junc_data, ax, ymax, y_is_minus=False, color="purple", xaxis_lenght=0):
            junc_data = _filter_cal_junc(junc_data)
            junc_data = _set_parm_junc(junc_data)
            plot_junc_core(junc_data, ax, ymax, y_is_minus, color, xaxis_lenght)

        ax = self.ax
        depth_junc_data_all, depth_junc_data_sense, depth_junc_data_anti = self.region_data
        ymax = MIN_YLIM
        ymin = 0
        
        ls_config = self.config["rnasashimiplot"]
        all_color = ls_config["all_color"]
        plus_color = ls_config["plus_color"]
        minus_color = ls_config["minus_color"]
        junction_line_width = ls_config["junction_line_width"]
        min_junction_ratio = ls_config["min_junction_ratio"]
        min_junction_count = ls_config["min_junction_count"]
        line_width_based_on_count = ls_config["line_width_based_on_count"]
        max_junction_line_width = ls_config["max_junction_line_width"]
        min_junction_line_width = ls_config["min_junction_line_width"]
        density_line_width = ls_config["density_line_width"]
        x_base_pad = self.config["x_base_pad"]
        
        xaxis_lenght = np.abs(self.xlim[1] - self.xlim[0])        
        
        if self.config["read_specific_strand"]:
            plot_depth(depth_junc_data_sense[0], ax, 
                       y_is_minus=False, 
                       color=plus_color, 
                       x_base_pad=x_base_pad,
                       density_line_width=density_line_width)
            plus_y_max = depth_junc_data_sense[0]["depth"].max() 
            plot_depth(depth_junc_data_anti[0], ax, 
                       y_is_minus=True,
                       color=minus_color,
                       density_line_width=density_line_width)
            minus_y_max = depth_junc_data_anti[0]["depth"].max()
            if self.config["plot_junc"]:
                ymax = max(plus_y_max*1.4, ymax)
                ymin = -minus_y_max*1.4
            else:
                ymax = max(plus_y_max, ymax)
                ymin = -minus_y_max
            self.set_ax_ylim(ymin, ymax)
            junc_y_lim = max(self.rel_ylim_max, -self.rel_ylim_min)
            #另一种方法是正负分别用不同的junc_y_lim
            if self.config["plot_junc"]:
                plot_junc(depth_junc_data_sense[1], ax, junc_y_lim, False, plus_color, xaxis_lenght)
                plot_junc(depth_junc_data_anti[1], ax, junc_y_lim, True, minus_color, xaxis_lenght)
        else:
            plot_depth(depth_junc_data_all[0], ax, 
                       y_is_minus=False,
                       color=all_color,
                       x_base_pad=x_base_pad,
                       density_line_width=density_line_width)
            all_y_max = depth_junc_data_all[0]["depth"].max()
            if self.config["plot_junc"]:
                ymax = max(ymax, all_y_max*1.4)
            else:
                ymax = max(ymax, all_y_max)
            self.set_ax_ylim(ymin, ymax)
            if self.config["plot_junc"]:
                plot_junc(depth_junc_data_all[1], ax, self.rel_ylim_max, False, all_color, xaxis_lenght)