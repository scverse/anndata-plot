from ._utils import ColorLike
from typing import Literal, Sequence
import numpy as np
import pandas as pd

from typing import Collection
from numpy.typing import NDArray

from .render_utils import ColorOpts
from .render_utils import LegendOpts
from .render_utils import SizeOpts
from .render_utils import AxisOpts
from .render_utils import _LegendLoc
from .render_utils import _FontWeight
from .render_utils import _FontSize

from dataclasses import asdict

from matplotlib.colors import is_color_like

from anndata import AnnData

import holoviews as hv

# copied from scanpy, pl._anndata line +- 235
# adapted to work with 1 color, that can only be in names or be a column
def _check_if_annotations(
    adata: AnnData,
    axis_name: Literal["obs", "var"],
    x: str | None = None,
    y: str | None = None,
    color: Collection[str | ColorLike] | None = None,
    use_raw: bool | None = None,
) -> bool:
    """Check if `x`, `y`, and `colors` are annotations of `adata`.

    If `axis_name` is `obs`, checks in `adata.obs.columns` and `adata.var_names`,
    if `axis_name` is `var`, checks in `adata.var.columns` and `adata.obs_names`.
    """
    annotations: pd.Index[str] = getattr(adata, axis_name).columns
    other_ax_obj = (
        adata.raw if use_raw and axis_name == "obs" else adata
    )
    names: pd.Index[str] = getattr(
        other_ax_obj, "var" if axis_name == "obs" else "obs"
    ).index

    def is_annotation(needle: pd.Index) -> NDArray[np.bool_]:
        return needle.isin({None}) | needle.isin(annotations) | needle.isin(names)

    if not is_annotation(pd.Index([x, y])).all():
        return False

    return bool(is_annotation(pd.Index([color])).all())

def scatter(
    adata: AnnData,
    x: str | None = None,
    y: str | None = None,
    basis: str | None = None,
    color_by: str | None = None,
    title: str | None = None,
    color_opts: ColorOpts | dict | None = None,
    legend_opts: LegendOpts | dict | None = None,
    size_opts: SizeOpts | dict | None = None,
    interactive: bool = False,
):

    # determine which dims to use
    if basis is not None:
        kdims = [f"obsm.X_{basis}.0", f"obsm.X_{basis}.1"]
        vdims = []
        if color_by is not None and color_by in adata.obs.columns:
            vdims = [f"obs.{color_by}"] if color_by is not None else []
        elif color_by is not None and color_by in adata.var_names:
            vdims = [color_by]
    elif _check_if_annotations(adata, "obs", x=x, y=y, color = color_by):
        kdims = [f"obs.{x}", f"obs.{y}"]
        vdims = [f"obs.{color_by}"] if color_by is not None else []
    elif _check_if_annotations(adata, "var", x=x, y=y, color = color_by):
        kdims = [f"var.{x}", f"var.{y}"]
        vdims = [f"var.{color_by}"] if color_by is not None else []
    else:
        msg = (
            "`x`, `y`, and potential `color` inputs must all "
            "come from either `.obs` or `.var`"
        )
        raise ValueError(msg)

    if title is None and color_by is not None:
        title = color_by.replace("_", " ")

    title_opts = {"title": title}

    if interactive:
        allopts = get_interactive_opts(vdims[0])
    else:
        allopts = get_static_opts(vdims[0])

    all_opts = {
        **title_opts,
        **allopts
    }

    return hv.Points(adata, kdims, vdims).opts(**all_opts)

def get_color_opts(kdims, vdims, interactive, color_by = None, cmap = None, palette = None, **kwargs):
    args = {}

    if len(vdims) != 0:
        args["color"] = vdims[0]
    if palette is not None:
        args["cmap"] = palette
    if cmap is not None:
        args["cmap"] = cmap
    return ColorOpts(**args)

def get_legend_opts(legend_loc = None, interactive = False, **kwargs):
    if legend_loc is not None:
        return LegendOpts(legend_position = legend_loc, **kwargs)
    return LegendOpts(**kwargs)


def get_interactive_opts(color_by):
    return {
        "cmap": "viridis",
        "color": color_by,
        "width": 550,
        "height": 550,
        "legend_position": "bottom_left"
    }

def get_static_opts(color_by):
    return {
        "cmap": "viridis",
        "color": color_by,
        "fig_size": 250,
    }
