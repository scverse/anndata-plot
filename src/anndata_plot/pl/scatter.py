from ._utils import ColorLike
from typing import Literal, Sequence
import numpy as np

from .render_utils import ColorOpts
from .render_utils import LegendOpts
from .render_utils import SizeOpts

from dataclasses import asdict

import holoviews as hv

def scatter(
    adata, #Y: np.ndarray,
    Y: Sequence[str],
    *,
    colors: str | Sequence[ColorLike | np.ndarray] = "blue", #should probably be a colormapping?
    sort_order=True,
    alpha=None,
    highlights=(),
    right_margin=None,
    left_margin=None,
    projection: Literal["2d", "3d"] = "2d",
    title=None,
    component_name="DC",
    component_indexnames=(1, 2, 3),
    axis_labels=None,
    colorbars=(False,),
    sizes=(1,),
    markers=".",
    color_map="viridis",
    show_ticks=True,
    ax=None,
):
    # fig = hv.render(hv.Points(adata, Y), backend="matplotlib")
    # points = hv.Points(adata, Y, colors).opts(color=colors)
    # points.opts(color = colors)

    # if colors is a column in obs --> hv.Points(adata, Y, colors).opts(color=colors)
    # if colors is just one color --> hv.Points(adata, Y).opts(color=colors) --> but is this needed?
    # if colors is a var_name --> hv.Points(adata, Y, colors).opts(color=colors)

    if color_map is None:
        color_map = "viridis"

    legend_opts = LegendOpts()
    color_opts = ColorOpts(color = colors, cmap = color_map)
    size_opts = SizeOpts()

    # merge opts dicts
    opts = {**asdict(legend_opts),
            **asdict(color_opts),
            **asdict(size_opts)}

    return (hv.Points(adata, Y, colors).opts(**opts))
