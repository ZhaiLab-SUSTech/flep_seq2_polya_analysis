__version__ = "0.2"

"""
version 0.2.1 by jiajinbu 2021.01.04
version 0.2

需解决绘图完成后，重新绘制，以前的参数还在对象里保存的问题。

引用其他模块:
from baseseq_tool import revcom #baseseq_tool
bugff
"""

from .gv import BuGV
from .tracks.track import BuGVTrack

#import sys
#from . readers import * #加载所有的类
#from . tracks import * #加载所有的类

#__all__ = []
#mod = sys.modules[__name__]
#classes = [getattr(mod, x) for x in dir(mod) if isinstance(getattr(mod, x), type)]
#for cls in classes:
#    __all__.append(cls.__name__)