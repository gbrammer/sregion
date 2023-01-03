from .sregion import SRegion, patch_from_polygon

try:
    from .version import __version__
except ImportError:
    __version__ = '0.1'
