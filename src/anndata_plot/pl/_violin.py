import holoviews as hv
import numpy as np
import pandas as pd


def violin(
    data: pd.DataFrame,
    var_names: list[str] | str,
    groupby: str | None = None,
    **kwargs,
):
    # If the column is categorical, use its categories mapped to continuous numbers.
    if pd.api.types.is_categorical_dtype(data[groupby]):
        categories = list(data[groupby].cat.categories)
        mapping = {cat: i for i, cat in enumerate(categories)}
        data[groupby + "_num"] = data[groupby].map(mapping)
        xdim = groupby + "_num"
        xticks = list(zip(list(mapping.values()), categories, strict=False))
        x_values = list(mapping.values())
        min_x = min(x_values)
        max_x = max(x_values)
    else:
        # If the column is numeric, use the native values.
        if pd.api.types.is_numeric_dtype(data[groupby]):
            xdim = groupby
            sorted_unique = sorted(data[groupby].unique())
            # Create xticks from unique values (converted to strings for labels)
            xticks = list(zip(sorted_unique, [str(val) for val in sorted_unique], strict=False))
            min_x = data[groupby].min()
            max_x = data[groupby].max()
        else:
            # For non-numeric, non-categorical columns, factorize the column.
            categories, codes = pd.factorize(data[groupby])
            data[groupby + "_num"] = codes
            xdim = groupby + "_num"
            xticks = list(zip(range(len(categories)), categories, strict=False))
            min_x = min(codes)
            max_x = max(codes)

    print(xticks)

    y_min = 0
    # Add a bit of padding on top (e.g. 5% extra space)
    y_max = data[var_names].max()
    y_padding = (y_max - y_min) * 0.05 if y_max > 0 else 0
    y_lim = (y_min - y_padding, y_max + y_padding)

    # Compute zoom and pan intervals based on the tick positions.
    if xticks is not None and len(xticks) > 1:
        tick_positions = sorted([t[0] for t in xticks])
        diffs = np.diff(tick_positions)
        min_interval = np.min(diffs) * 0.9
    else:
        min_interval = 0.5  # Default if only one tick.

    x_range_val = max_x - min_x
    max_interval = x_range_val + 2  # Add padding (1 unit each side)

    def _limit_x_range(plot, element):
        # Lock the x-axis bounds and zoom intervals.
        x_range_obj = plot.handles.get("x_range")
        if x_range_obj is not None:
            x_range_obj.bounds = (min_x - 0.5, max_x + 0.5)
            x_range_obj.min_interval = min_interval
            x_range_obj.max_interval = max_interval

    def _limit_y_range(plot, element):
        y_range_obj = plot.handles.get("y_range")
        if y_range_obj is not None:
            y_range_obj.bounds = y_lim
            y_range_obj.start = y_lim[0]
            y_range_obj.end = y_lim[1]

        # Set custom tick positions and labels on the Bokeh axes.
        # if xticks is not None:
        #     ticks_values = [t[0] for t in xticks]
        #     label_overrides = {t[0]: t[1] for t in xticks}
        #     xaxis = plot.handles.get("xaxis")
        #     if xaxis is not None:
        #         try:
        #             for axis in xaxis:
        #                 try:
        #                     axis.ticker = FixedTicker(ticks=ticks_values)
        #                 except Exception:
        #                     pass
        #                 axis.major_label_overrides = label_overrides
        #         except TypeError:
        #             try:
        #                 xaxis.ticker = FixedTicker(ticks=ticks_values)
        #             except Exception:
        #                 pass
        #             xaxis.major_label_overrides = label_overrides

    # Create the violin plot using the chosen dimension for the x-axis.
    vdim = hv.Dimension(var_names, range=(-y_padding / 5, y_lim[1] + y_padding / 5))

    violin_plot = hv.Violin(
        data,
        kdims=[xdim],
        vdims=[vdim],
    ).opts(
        # clip=(0, y_lim[1] + 1),
        show_legend=False,
        xlabel=groupby,
        ylabel=var_names,
        xlim=(min_x - 0.5, max_x + 0.5),
        ylim=y_lim,
        hooks=[_limit_x_range, _limit_y_range],
        # **kwargs,
    )
    return violin_plot
