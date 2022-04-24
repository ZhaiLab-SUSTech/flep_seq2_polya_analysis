from . track import BuGVTrack

class SpaceTrack(BuGVTrack):
    
    DEFAULT_CONFIG = {
        "ylim_style": None
    }
            
    def plot_ax(self):
        ax = self.ax
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.set_xticks([])
        ax.set_yticks([])
        