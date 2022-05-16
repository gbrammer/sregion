from .sregion import SRegion

try:
    from .version import __version__
except ImportError:
    return '0.1'
    