from importlib.metadata import version

import holoviews as hv

from . import pl, pp, tl
from .pl.utils import AnnDataInterface

__all__ = ["pl"]

__version__ = version("anndata-plot")

# register the AnnDataInterface with holoviews
if AnnDataInterface.datatype not in hv.core.data.datatypes:
    hv.core.data.datatypes.append(AnnDataInterface.datatype)

hv.core.Interface.register(AnnDataInterface)
