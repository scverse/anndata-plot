from dataclasses import dataclass

from typing import Literal

_LegendLoc = Literal[
    "none",
    "right margin",
    "on data",
    "on data export",
    "best",
    "upper right",
    "upper left",
    "lower left",
    "lower right",
    "right",
    "center left",
    "center right",
    "lower center",
    "upper center",
    "center",
]
_FontWeight = Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]
_FontSize = Literal[
    "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"
]

@dataclass
class LegendOpts:
    legend_position: str = "inner"
    legend_cols: int = 1
    show_legend: bool = True

@dataclass
class LegendOptsMpl(LegendOpts):
    legend_font_weight: _FontWeight = "normal"
    legend_font_size: _FontSize = "medium"



@dataclass
class ColorOpts:
    color: str = "blue"
    cmap: str = "viridis"

@dataclass
class SizeOpts:
    size: int = 100

@dataclass
class AxisOpts:
    xlabel: str = ""
    ylabel: str = ""
