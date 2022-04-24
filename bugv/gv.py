from bugv.tracks.track import BuGVTrack
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist
import bugv.utilities 

"""

Note: 
png or jupter print: may be lost some reads when reads are too many
pdf: color will become light when reads are too many.
svg: color remain norm when reads are too many. but may be freeze when open in firefox.
you can open svg in chrome. 
You can open 
if you save the figure into png or print in jupter, some reads may be lost when reads are too many.
You can save it to pdf. 
You can open figure in chrome.

version 0.2.3 by jiajinbu 2021.03.14
speed full-length by ax.add_collection(PatchCollection(rects))

version 0.2.1 by jiajinbu 2021.01.04
"""

class BuGV:
    
    """
    learn from pygenometracks
    1. 首先创建BuGV对象。
    bgv = BuGV()
    
    2. 添加track
    bgv.add_track(track_type, config)
    track_type为BuGVTrack子类的类名（字符串格式）。
    config为绘图参数字典。
    track添加后会将track对象添加到bgv.track_objs里。
    同时每个track对象的bgv属性设为改bgv对象。
    可以通过track_obj.load_config(config)更改track_obj绘图属性。
    也可以通过bgv.load_track_config(config)更改所有track_obj绘图属性。
    
    2.1 每次添加track时，除特定track(如XaxisTrack)外，都需提供file_name属性。
    这会调用track默认的config["file_type"]值（类的名字），生成对象存入track对象的
    obj属性，该对象需具有open，close，load_region_data方法。并且load_region_data后将
    值设为region_data。每次plot时都会调用load_region_data方法，并用region_data绘图。
    
    3. set_region([chr_name, start, end, strand]) # 1-based
    设置region，chr_name，region_start，region_end，region_strand，xlim属性
    同时调用每个track_obj.set_region(region)。
    
    4. plot方法
    绘图。也可以plot(region)，相当于先set_region(region),再plot()。
    绘制时，根据track_obj数量n，利用matplotlib.gridspec.GridSpec生成n*3个网格。
    每行对应一个track，每行三列的ax分别传给track_obj的left_ax, ax, right_ax属性。
    每个track的高度由track_obj.config["track_height"]指定。所有track的三列的宽度比例是一样的。
    left_width, right_width在default_config中有设置。也可以用set_width方法设置。具体使用
    见get_track_width。
    
    5. x轴坐标为原始坐标。如果region_strand为负，则反转轴坐标。
    具体为在plot函数里ax.invert_xaxis()。
    
    6. y轴见BuGVTrack文档。
    
    2021.03.14 test the speed:
    rna-seq for a gene, max depth is 13490, load data about 6.1s, load and plot about 6.9s
    full-length read for a gene,  max depth is 8980, load data about 1.6s, load and plot about 70s.
        load data and add rect about 20s, real plot about 50s.
    iter the data structure only need 0.2s. add rect ~16s, real plot ~50s.
    full-length read for a gene, max depth is 89, load and plot about 1.7s
    Thus, you need to reduce the time of rna-seq load and full-length plot.
    use `ax.add_collection(PatchCollection(rects))` to speed, from 70s to 10s.
    open file, read read to feature_array about 1s.
    pre_process, filter_data about 0.04s, sort about 0.2-0.3s.
    about 4-5s to generate react, 5s to add rect, 1-2s to plot.
    """
    
    DEFAULT_CONFIG = {"margin": {'left': 0.04, 'right': 0.92, 'bottom': 0.03, 'top': 0.97, 'wspace': 0, 'hspace': 0},
                      "left_width": 0.3,
                      "right_width": 1}

    def __init__(self):
        
        self.__track_classes = {}
        self.update_track_classes() #generate __track_classes
        self.config = {}
        self.load_bgv_config(self.DEFAULT_CONFIG)
        self.tracks = []
        self.tracks_dict = {}
        self.fig = None
        self.xaxis_is_invert = False
        self.region = None
        self.select_features = None
    
    def update_track_classes(self, clses=None):
        """
        Usually not be used by User.
        clses is None or a list of class which is the subclass of BuGVTrack.
        If clses is None, will iter all subclass of BuGVTrack.
        Add all BuGVTrack subclasses to self.__track_classes dict.
        The name is class name (cls.__name__), the value is class.
        """
        def find_subclass(class_name):
            track_classes = {}
            for cls in class_name.__subclasses__():
                track_classes[cls.__name__] = cls
                for k, v in find_subclass(cls).items():
                    track_classes[k] = v
            return track_classes
        if clses:
            for cls in clese:
                track_class[cls.__name__] = cls
        else:
            track_classes = find_subclass(BuGVTrack)
        self.__track_classes = track_classes
        
    def get_track_class(self, track_type):
        """
        Usually not be used by User.
        """
        return self.__track_classes[track_type]
    
    def get_track(self, track_index):
        try:
            return self.tracks_dict[track_index]
        except:
            return self.tracks[track_index]
        
    def get_tracks_by_type(self, track_type):
        select_tracks = []
        if isinstance(track_type, str):
            track_type = self.get_track_class(track_type)
        for track_obj in self.tracks:
            if isinstance(track_obj, track_type):
                select_tracks.append(track_obj)
        return select_tracks
        
    def get_bed(self):
        return self.get_tracks_by_type("BedTrack")[0].obj
    
    def get_feature(self, feature):
        try:
            bed_obj = self.bed
        except:
            bed_obj = self.get_bed()
        return bed_obj.get_feature(feature)
    
    def load_track_config(self, config={}, track_class_names=[]):
        """
        Load config to track obj.
        If not provide track_class_names, load config to all track obj.
        If provide track_class_names, only load config to the track objs of specific classes.
        """
        if track_class_names:
            select_tracks = []
            for track_class_name in track_class_names:
                for track_obj in self.tracks:
                    if isinstance(track_obj, self.get_track_class(track_class_name)):
                        select_tracks.append(track_obj)
        else:
            select_tracks = self.tracks
        for track_obj in select_tracks:
            track_obj.load_config(config)
    
    def add_track(self, track_type, file_name=None, config=None, track_name=None):
        if not config:
            config = {}
        if file_name:
            config["file_name"] = file_name
        if not track_name:
            track_name = str(len(self.tracks) + 1)
        track_obj = self.get_track_class(track_type)(config, track_name)
        track_obj.bgv = self
        self.tracks.append(track_obj)
        self.tracks_dict[track_name] = track_obj
        return(track_obj)
    
    def move_track(self):
        pass
    
    def load_bgv_config(self, config={}):
        bugv.utilities.update_dict(self.config, config)
        
    def set_width(self, left_width=None, right_width=None):
        if left_width is not None:
            self.config["left_width"] = left_width
        if right_width is not None:
            self.config["right_width"] = right_width        
    
    def set_region(self, region=None, select_features=None, rel_position=None, min_x_axis_length=None, cut_region=[None, None]):
        #注意region和select_features为None，set_region等于什么也没做
        #不会覆盖以前的值
        if region is not None:
            chr_name, region_start, region_end, region_strand = region
            self.region = region
            self.chr_name = chr_name
            self.region_start = region_start
            self.region_end = region_end
            self.region_strand = region_strand
            self.xlim = [self.region_start, self.region_end]#[0, self.region_end - self.region_start + 1]
        self.select_features = select_features
        self.config["rel_position"] = rel_position
        self.config["cut_region"] = rel_position
        for track_obj in self.tracks:
            track_obj.set_region(self.region, self.select_features, self.config["rel_position"], cut_region)
            track_obj.xlim = [self.xlim[0], self.xlim[1]]
        self.min_x_axis_length = min_x_axis_length
        
    def get_track_height(self):
        track_height = [o.config["track_height"] for o in self.tracks]
        self.track_height = track_height
        return(track_height)        
    
    def get_track_width(self):
        min_x_axis_length = self.min_x_axis_length
        x_rel_length = self.region_end - self.region_start + 1 # 1 based
        left_plot_width = self.config["left_width"]
        right_plot_width = self.config["right_width"]
        plot_width = self.fig.get_figwidth() - left_plot_width - right_plot_width
        if min_x_axis_length and min_x_axis_length > x_rel_length:
            core_plot_width = plot_width * x_rel_length / min_x_axis_length
            right_plot_width = right_plot_width + (plot_width - core_plot_width)
        else:
            core_plot_width = plot_width
        return([left_plot_width, core_plot_width, right_plot_width])
    
    def add_highlight(self, tracks=None, xs=[], zorder=-1, color="grey", **key_values):
        
        if tracks is None:
            tracks = self.tracks
        
        for track in tracks:
            for x1, x2 in xs:
                track.ax.axvspan(x1, x2, zorder=zorder, color=color, **key_values)        
        
    def plot_mRNA(self, mRNA_name, upstream_length=500, downstream_length=500, cut_region=[None, None], file_pdf="", figsize=[4, 2.2], select_feature_flag=True):
        mRNA_feature = self.get_feature(mRNA_name)
        if mRNA_feature.strand == "-":
            upstream_length, downstream_length = downstream_length, upstream_length
        select_reigion_start, select_reigion_end = mRNA_feature.start - upstream_length, mRNA_feature.end + downstream_length
        if cut_region[0] is not None:
            if mRNA_feature.strand == "+":
                select_reigion_start =  mRNA_feature.start + cut_region[0] - upstream_length - 1
            else:
                select_reigion_end = mRNA_feature.end - cut_region[0] + downstream_length + 1
        if cut_region[1] is not None:
            if mRNA_feature.strand == "+":
                select_reigion_end =  mRNA_feature.start + cut_region[1] - upstream_length - 1
            else:
                select_reigion_start = mRNA_feature.end - cut_region[1] + downstream_length + 1
        select_region = [mRNA_feature.chr_name, select_reigion_start, select_reigion_end, mRNA_feature.strand]
        if select_feature_flag:
            select_features = [mRNA_name]
        else:
            select_features = None
        if mRNA_feature.strand == "+":
            rel_position = [mRNA_feature.start, mRNA_feature.end]
        else:
            rel_position = [mRNA_feature.end, mRNA_feature.start]
        self.plot(region=select_region, 
                 select_features=select_features, 
                 rel_position=rel_position,
                 cut_region=cut_region,
                 figsize=figsize, 
                 dpi=300)
        if file_pdf:
            self.fig.savefig(file_pdf)
    
    def plot(self, fig=None, region=None, select_features=None, rel_position=None, min_x_axis_length=None, cut_region=[None, None], **kwargs):
        
        if fig is None:
            fig = plt.figure(**kwargs)
        self.fig = fig
        
        #注意region和select_features为None，set_region等于什么也没做
        #不会覆盖以前的值
        self.set_region(region, select_features, rel_position, min_x_axis_length, cut_region)
        
        fig.subplots_adjust(wspace=self.config["margin"]['wspace'], 
                            hspace=self.config["margin"]['hspace'],
                            left=self.config["margin"]['left'],
                            right=self.config["margin"]['right'],
                            bottom=self.config["margin"]['bottom'],
                            top=self.config["margin"]['top'])
        
        track_height = self.get_track_height()
        width_ratios = self.get_track_width()
                
        grids = matplotlib.gridspec.GridSpec(len(self.tracks), 3,
                                             height_ratios=track_height,
                                             width_ratios=width_ratios,
                                             wspace=0.01)
        
        first_ax = None
        for track_index, track_obj in enumerate(self.tracks):
            left_ax = fig.add_subplot(grids[track_index, 0])
            left_ax.set_axis_off()
            right_ax = fig.add_subplot(grids[track_index, 2])
            right_ax.set_axis_off()
            if track_obj.config["is_axisartist"]:
                ax = axisartist.Subplot(fig, grids[track_index, 1])
                if first_ax:
                    fig.add_subplot(ax, sharex=first_ax)
                else:
                    first_ax = ax
                    fig.add_subplot(ax)
            else:
                if first_ax:
                    ax = fig.add_subplot(grids[track_index, 1], sharex=first_ax)
                else:
                    ax = fig.add_subplot(grids[track_index, 1])
                    first_ax = ax
            ax.set_xlim(self.xlim)
            if self.region_strand == "-":
                ax.invert_xaxis()
                track_obj.xaxis_is_invert = True
                self.xaxis_is_invert = True
            else:
                track_obj.xaxis_is_invert = False
                self.xaxis_is_invert = False
            track_obj.left_ax = left_ax
            track_obj.right_ax = right_ax
            track_obj.ax = ax
            track_obj.plot()                    
        
    def read_config(self):
        pass
    
    def export_config(self):
        pass
    
    def open(self):
        for track_obj in self.tracks:
            track_obj.open()
    
    def close(self):
        for track_obj in self.tracks:
            track_obj.close()
            
    def set_track_height(self, track_heights):
        for track_obj, track_height in zip(self.tracks, track_heights):
            track_obj.track_height = track_height
            
    def copy(self):
        pass