import os
import sys
import pkgutil
import importlib

__all__ = []

for (module_loader, name, ispkg) in pkgutil.iter_modules([os.path.dirname(__file__)]):
    mod = importlib.import_module('.' + name, __package__)
    classes = [getattr(mod, x) for x in dir(mod) if isinstance(getattr(mod, x), type)]
    #from . track import BuGVTrack
    #for cls in BuGVTrack.__subclasses__():
    for cls in classes:
        setattr(sys.modules[__name__], cls.__name__, cls)
        __all__.append(cls.__name__)