from typing import Literal

_FontWeight = Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]
_FontSize = Literal[
    "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"
]
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
ColorLike = str | tuple[float, ...]
