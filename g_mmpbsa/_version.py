try:
    from .g_mmpbsa import gmx_version
    gmx_version = gmx_version()
except:
    pass

import importlib.metadata
__version__ = importlib.metadata.version('g_mmpbsa')
