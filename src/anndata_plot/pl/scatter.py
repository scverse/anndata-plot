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

# def scatter(
#     adata, #Y: np.ndarray,
#     Y: Sequence[str],
#     title: str | None = None,
#     color_by: str | None = None,
#     legend_opts: LegendOpts = LegendOpts(),
#     color_opts: ColorOpts = ColorOpts(),
#     size_opts: SizeOpts = SizeOpts(),
#     aixs_opts: AxisOpts = AxisOpts(),
#     backend_opts: dict = None,

#     # *,
#     # colors: str | Sequence[ColorLike | np.ndarray] = "blue", #should probably be a colormapping?
#     # sort_order=True,
#     # alpha=None,
#     # highlights=(),
#     # right_margin=None,
#     # left_margin=None,
#     # projection: Literal["2d", "3d"] = "2d",
#     # title=None,
#     # component_name="DC",
#     # component_indexnames=(1, 2, 3),
#     # axis_labels=None,
#     # colorbars=(False,),
#     # sizes=(1,),
#     # markers=".",
#     # color_map="viridis",
#     # show_ticks=True,
#     # ax=None,
# ):
#     # fig = hv.render(hv.Points(adata, Y), backend="matplotlib")
#     # points = hv.Points(adata, Y, colors).opts(color=colors)
#     # points.opts(color = colors)

#     # if colors is a column in obs --> hv.Points(adata, Y, colors).opts(color=colors)
#     # if colors is just one color --> hv.Points(adata, Y).opts(color=colors) --> but is this needed?
#     # if colors is a var_name --> hv.Points(adata, Y, colors).opts(color=colors)

#     if title is None:
#         title = f"Scatter plot of {Y[0]} and {Y[1]}"

#     # merge opts dicts
#     opts = {**asdict(legend_opts),
#             **asdict(color_opts),
#             **asdict(size_opts),
#             "title": title,
#             }

#     if color_by is None:
#         return (hv.Points(adata, Y).opts(**opts))

#     return (hv.Points(adata, Y, color_by).opts(**opts))

# copied from scanpy, pl._anndata line +- 235
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
    cmap: str | None = None,
    palette: Sequence[ColorLike] | None = None,
    legend_opts: LegendOpts | dict | None = None,
    legend_loc: _LegendLoc | None = "right margin",
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

    # check if color_opts is a dict
    if isinstance(color_opts, dict):
        color_opts = get_color_opts(kdims, vdims, color_by, **color_opts)
    elif color_opts is None:
        color_opts = get_color_opts(kdims, vdims, color_by, cmap, palette)

    all_opts = {
        **title_opts,
        **asdict(color_opts),
    }

    print(locals())

    print(kdims)
    print(vdims)
    print(all_opts)

    return hv.Points(adata, kdims, vdims).opts(**all_opts)

def get_color_opts(kdims, vdims, color_by = None, cmap = None, palette = None):
    args = {}

    if len(vdims) != 0:
        args["color"] = vdims[0]
    if palette is not None:
        args["cmap"] = palette
    if cmap is not None:
        args["cmap"] = cmap
    return ColorOpts(**args)
