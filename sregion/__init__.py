from .sregion import SRegion, patch_from_polygon
from .utils import moc_patch_union

try:
    from .version import __version__
except ImportError:
    __version__ = "0.1"
