import holoviews as hv
import numpy as np
from bokeh.models import HoverTool
from holoviews import dim


def volcano_plot(results_df, pval_threshold=0.05, LFC_threshold=2):
    """
    Create a volcano plot from a results dataframe.

    Parameters
    ----------
    results_df : pd.DataFrame
        The results dataframe.
    pval_threshold : float, optional
        The p-value threshold.
    LFC_threshold : float, optional
        The log2 fold change threshold.

    Returns
    -------
    p : holoviews.Scatter
        The volcano plot.
    """

    def map_DE(a):
        log2FoldChange, pval = a
        if pval < pval_threshold:
            if log2FoldChange > LFC_threshold:
                return "Up-reg"
            elif log2FoldChange < -LFC_threshold:
                return "Down-reg"
        return "Non-sig"

    # Create conditions for each category
    results_df["de_status"] = results_df[["log2FoldChange", "padj"]].apply(map_DE, axis=1)

    hover = HoverTool(
        tooltips=[
            ("gene", "@gene"),
        ]
    )

    p = (
        hv.Scatter(results_df, kdims="log2FoldChange", vdims=["padj", "de_status", "gene"], label="Volcano Plot")
        .transform(padj=-dim("padj").log10())
        .opts(
            color="de_status",
            cmap={"Up-reg": "indianred", "Down-reg": "cornflowerblue", "Non-sig": "lightgrey"},
            size=5,
            tools=[hover],
            width=600,
            height=600,
            xlabel="log2 Fold Change",
            ylabel="-log10(adjusted p-value)",
            title="Volcano Plot",
        )
        * hv.VLine(-LFC_threshold).opts(color="black", line_dash="dashed")
        * hv.VLine(LFC_threshold).opts(color="black", line_dash="dashed")
        * hv.HLine(-np.log10(pval_threshold)).opts(color="black", line_dash="dashed")
    )
    return p
