from importlib.metadata import version

try:
    import holoviews as hv
except ImportError:
    raise ImportError("holoviews is required for anndata-plot. Please install it with:\npip install holoviews")

from . import pl
from .pl.utils import AnnDataInterface

__all__ = ["pl"]

__version__ = version("anndata-plot")

# register the AnnDataInterface with holoviews
if AnnDataInterface.datatype not in hv.core.data.datatypes:
    hv.core.data.datatypes.append(AnnDataInterface.datatype)

hv.core.Interface.register(AnnDataInterface)
