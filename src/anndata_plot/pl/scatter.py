from ._utils import ColorLike
from typing import Literal, Sequence
import numpy as np

import holoviews as hv

def scatter(
    adata, #Y: np.ndarray,
    Y: Sequence[str],
    *,
    colors: str | Sequence[ColorLike | np.ndarray] = "blue", # should probably be a colormapping?
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
    print("HV scatter")
    hv.Points(adata, Y)
    return None

