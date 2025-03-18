from dataclasses import dataclass

@dataclass
class LegendOpts:
    legend_position: str = "right"

@dataclass
class ColorOpts:
    color: str = "blue"
    cmap: str = "viridis"

@dataclass
class SizeOpts:
    height: int = 400
    width: int = 600
