import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from typing import Union, Literal, Optional


def heatmap_triangle(
        matrix: Union[pd.DataFrame, np.array],
        direction: Literal["up", "down", "left", "right"] = "up",
        names=None,
        ticks_rotation: float = 0,
        ax: Optional["matplotlib.axes._axes.Axes"] = None,
        show_cbar: bool = False,
        cmap: str = "coolwarm",
        cbar_label: str = "",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        **cbar_params,
    ) -> None:
    """Plot upper diagonal of a `matrix` as a triangular heatmap
    
    Parameters
    ----------
    matrix : Union[pd.DataFrame, np.array]
        A square matrix to plot. Only the upper diagonal is plotted assumng that the matrix is symmetrical,
        e.g. containing pairwise distances
    direction : Literal["up", "down", "left", "right"]
        A direction in which the triangle will be pointed
    names : Optional
        Names of the observations to use as tick labels. If `matrix` is a DataFrame, its index will be used
    ticks_rotation : float = None
        Rotation of ticks if `names` are set
    ax : Optional["matplotlib.axes._axes.Axes"]
        Axis to plot the heatmap on
    show_cbar : bool = False
        Whether to display colorbar
    cmap : str = "coolwarm"
        Color map for the plot
    cbar_label : str = ""
        Title of the colorbar
    vmin : Optional[float]
        Minimum value for colormap
    vmax : Optional[float]
        Maximum value for colormap
    **cbar_params
        Optional colorbar parameters
    """
    N = len(matrix)

    if names is None and type(matrix) is pd.DataFrame:
        names = matrix.index

    # Define a mesh for heatmap. It will be centered at 0 and take values from -N to N
    # Each square will thus have length 2
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)]) - N / 2

    # Rotate the matrix by 45 degrees
    rotation = [[1, 1], [-1, 1]]
    A = A @ rotation

    # Additionally rotate the matrix according to the direction
    match direction:
        case "left":
            x = np.rot90(matrix, k=0)
        case "up":
            x = np.rot90(matrix, k=1)
        case "right":
            x = np.rot90(matrix, k=2)
        case "down":
            x = np.rot90(matrix, k=3)
        case _:
            raise ValueError('Wrong value of "rotation". Please set one of "up", "down", "left", "right"')

    # Create the mesh and plot it
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)
    caxes = plt.pcolormesh(X, Y, x, axes=ax, cmap=cmap, vmin=vmin, vmax=vmax)

    if ax is None:
        ax = plt.gca()

    # Remove meaningless ticks and put sample names to the correct side
    match direction:
        case "left":
            ax.set_xlim(right=0)
            ax.set_xticks([])
            ax.yaxis.tick_right()
        case "up":
            ax.set_ylim(bottom=0)
            ax.set_yticks([])
        case "right":
            ax.set_xlim(left=0)
            ax.set_xticks([])
        case "down":
            ax.set_ylim(top=0)
            ax.set_yticks([])
            ax.xaxis.tick_top()

    
    # If names are given, set them as ticks
    label_loc = np.arange(-N, N, 2) + 1
    match direction:
        case "left" | "right":
            if names is not None:
                ax.set_yticks(label_loc, names, rotation=ticks_rotation)
            else:
                ax.set_yticks([])
        case "up" | "down":
            if names is not None:
                ax.set_xticks(label_loc, names, rotation=ticks_rotation)
            else:
                ax.set_xticks([])
    
    # Add a colorbar below the heatmap triangle
    if show_cbar:
        if cbar_params is None:
            cbar_params = {}

        if "orientation" not in cbar_params:
            match direction:
                case "up" | "down":
                    cbar_params["orientation"] = "horizontal"
                case "left" | "right":
                    cbar_params["orientation"] = "vertical"

        cb = plt.colorbar(
            caxes,
            ax=ax,
            **cbar_params,
        )
        cb.set_label(cbar_label)
