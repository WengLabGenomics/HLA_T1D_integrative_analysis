import numpy as np
import scipy as sp
import pandas as pd
import numbers
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import scanpy as sc
import os
import anndata
import matplotlib.transforms as mtrans
import matplotlib.pyplot as plt
from typing import List, Dict
from anndata import read_h5ad
import fnmatch
import matplotlib.patches as patches
import matplotlib.pyplot as plt

# ================================================================================
# Plotting utilities
# ================================================================================
def small_squares(ax, pos, size=1, linewidth=0.8):
    """
    Draw many small squares on ax, given the positions of
    these squares.
    """
    for xy in pos:
        x, y = xy
        margin = (1 - size) / 2
        rect = patches.Rectangle(
            (x + margin, y + margin),
            size,
            size,
            linewidth=linewidth,
            edgecolor="k",
            facecolor="none",
            zorder=20,
        )
        ax.add_patch(rect)


def plot_group_stats(
    dict_df_stats: Dict = None,
    df_fdr_prop: pd.DataFrame = None,
    df_assoc_fdr: pd.DataFrame = None,
    df_hetero_fdr: pd.DataFrame = None,
    plot_kws: Dict = None,
):
    """plot group-level statistics for scDRS results
    Parameters
    ----------
    dict_df_stats: Dict
        dictionary from trait -> pd.DataFrame
        group-level statistics from `scdrs perform-downstream: group-analysis`
        see https://martinjzhang.github.io/scDRS/notebooks/quickstart.html for an
        example. Alternatively, df_fdr_prop, df_assoc_fdr, df_hetero_fdr can be
        separately provided.
    df_fdr_prop : pd.DataFrame
        dataframe of proportion of cells with FDR < 0.1
    df_assoc : pd.DataFrame
        dataframe of group-level association statistics
    df_hetero : pd.DataFrame
        dataframe of group-level heterogeneity statistics
    plot_kws : Dict
        dictionary of plotting parameters (you can adjust them by scaling them with a factor to the default value), containing
        - cb_location: location of colorbar, default="top"
        - cb_fraction: fraction of colorbar, default=0.1
        - cb_pad: padding of colorbar, default=0.05
        - cb_vmin: minimum value of colorbar, default=0
        - cb_vmax: maximum value of colorbar, default=0.2
        - cb_n_bin: number of bins of colorbar, default=5
        - square_size: size of each heatmap grid, default=30
        - hetero_size: size of each cross denoting heterogeneity, default=8
        - signif_size: size of each small square denoting significant cell type, default=0.6
        - signif_width: width of each small square denoting significant cell type, default=0.5
    """
    DEFAULT_PLOT_KWS = {
        "cb_location": "top",
        "cb_fraction": 0.1,
        "cb_pad": 0.05,
        "cb_vmin": 0,
        "cb_vmax": 0.2,
        "cb_n_bin": 5,
        "square_size": 30,
        "hetero_size": 8,
        "signif_size": 0.6,
        "signif_width": 0.5,
    }
    if plot_kws is None:
        plot_kws = DEFAULT_PLOT_KWS
    else:
        # add default values
        for k, v in DEFAULT_PLOT_KWS.items():
            if k not in plot_kws:
                plot_kws[k] = v

    if dict_df_stats is not None:
        trait_list = list(dict_df_stats.keys())
        # compile df_fdr_prop, df_assoc_fdr, df_hetero_fdr from dict_df_stats
        df_fdr_prop = pd.concat(
            [
                dict_df_stats[trait]["n_fdr_0.1"] / dict_df_stats[trait]["n_cell"]
                for trait in trait_list
            ],
            axis=1,
        ).T
        df_assoc_fdr = pd.concat(
            [dict_df_stats[trait]["assoc_mcp"] for trait in trait_list], axis=1
        ).T
        df_assoc_fdr = pd.DataFrame(
            multipletests(df_assoc_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_assoc_fdr.shape
            ),
            index=df_assoc_fdr.index,
            columns=df_assoc_fdr.columns,
        )
        df_hetero_fdr = pd.concat(
            [dict_df_stats[trait]["hetero_mcp"] for trait in trait_list], axis=1
        ).T
        df_hetero_fdr = pd.DataFrame(
            multipletests(df_hetero_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_hetero_fdr.shape
            ),
            index=df_hetero_fdr.index,
            columns=df_hetero_fdr.columns,
        )
        df_fdr_prop.index = trait_list
        df_assoc_fdr.index = trait_list
        df_hetero_fdr.index = trait_list
    else:
        assert (
            (df_fdr_prop is not None)
            and (df_assoc_fdr is not None)
            and (df_hetero_fdr is not None)
        ), "If dict_df_stats is not provided, df_fdr_prop, df_assoc_fdr, df_hetero_fdr must be all provided."

    df_hetero_fdr = df_hetero_fdr.applymap(lambda x: "×" if x < 0.05 else "")
    df_hetero_fdr[df_assoc_fdr > 0.05] = ""

    fig, ax = plot_heatmap(
        df_fdr_prop,
        squaresize=plot_kws["square_size"],
        heatmap_annot=df_hetero_fdr,
        heatmap_annot_kws={"color": "black", "size": plot_kws["hetero_size"]},
        heatmap_cbar_kws=dict(
            use_gridspec=False,
            location=plot_kws["cb_location"],
            fraction=plot_kws["cb_fraction"],
            pad=plot_kws["cb_pad"],
            drawedges=True,
        ),
        heatmap_vmin=plot_kws["cb_vmin"],
        heatmap_vmax=plot_kws["cb_vmax"],
        colormap_n_bin=plot_kws["cb_n_bin"],
    )
    ax.set_xlabel(None)

    small_squares(
        ax,
        pos=[(y, x) for x, y in zip(*np.where(df_assoc_fdr < 0.05))],
        size=plot_kws["signif_size"],
        linewidth=plot_kws["signif_width"],
    )

    cb = ax.collections[0].colorbar
    cb.ax.tick_params(labelsize=8)

    cb.ax.set_title("Prop. of sig. cells (FDR < 0.1)", fontsize=8)
    cb.outline.set_edgecolor("black")
    cb.outline.set_linewidth(1)
    return fig, ax


def discrete_cmap(N, base_cmap=None, start_white=True):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    if start_white:
        color_list[0, :] = 1.0
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def plot_heatmap(
    df,
    dpi=150,
    squaresize=20,
    heatmap_annot=None,
    heatmap_annot_kws={"color": "black", "size": 4},
    heatmap_linewidths=0.5,
    heatmap_linecolor="gray",
    heatmap_xticklabels=True,
    heatmap_yticklabels=True,
    heatmap_cbar=True,
    heatmap_cbar_kws=dict(use_gridspec=False, location="top", fraction=0.03, pad=0.01),
    heatmap_vmin=0.0,
    heatmap_vmax=1.0,
    xticklabels_rotation=45,
    colormap_n_bin=5,
):
    figwidth = df.shape[1] * squaresize / float(dpi)
    figheight = df.shape[0] * squaresize / float(dpi)
    fig, ax = plt.subplots(1, figsize=(figwidth, figheight), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_facecolor("silver")
    sns.heatmap(
        df,
        annot=heatmap_annot,
        annot_kws=heatmap_annot_kws,
        fmt="",
        cmap=discrete_cmap(colormap_n_bin, "Reds"),
        linewidths=heatmap_linewidths,
        linecolor=heatmap_linecolor,
        square=True,
        ax=ax,
        xticklabels=heatmap_xticklabels,
        yticklabels=heatmap_yticklabels,
        cbar=heatmap_cbar,
        cbar_kws=heatmap_cbar_kws,
        vmin=heatmap_vmin,
        vmax=heatmap_vmax,
    )

    plt.yticks(fontsize=8)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=xticklabels_rotation,
        va="top",
        ha="right",
        fontsize=8,
    )
    ax.tick_params(left=False, bottom=False, pad=-2)
    trans = mtrans.Affine2D().translate(5, 0)
    for t in ax.get_xticklabels():
        t.set_transform(t.get_transform() + trans)
    return fig, ax


def plot_assoc_matrix(pval_dict, pval_index, meta_df, stratify_by, fdr_level=0.2):
    """
    plot (tissue / tissue-celltype) x traits association matrix
    # Arguments:
    - pval_dict: dictionary of the p-value: trait -> array, array has been ordered
    - pval_index: index of the p-values
    - meta_df: providing metadata for `stratify_by` use
    - stratify_by: the column in `meta_df` to stratify the association.
    """

    def num2str(x):
        if x > 1000:
            return "%0.1fk" % (x / 1000)
        elif x > 0:
            return "%d" % x
        else:
            return ""

    trait_list = list(pval_dict.keys())
    pval_df = pd.DataFrame(pval_dict, index=pval_index).join(meta_df, how="inner")

    assert stratify_by in ["tissue", "tissue_celltype", "celltype"]
    stratify_list = sorted(pval_df[stratify_by].unique())

    # Dataframe for plotting
    df_plot = pd.DataFrame(index=stratify_list, columns=trait_list, data=0)

    for trait in trait_list:
        pval = pval_df[trait].values
        fdr = multipletests(pval, method="fdr_bh")[1]
        temp_df = (
            pval_df.loc[fdr < fdr_level].groupby([stratify_by]).agg({stratify_by: len})
        )
        temp_df = temp_df.loc[~temp_df[stratify_by].isna()]
        df_plot.loc[temp_df.index, trait] = temp_df[stratify_by].values

    df_plot = df_plot.loc[df_plot.max(axis=1) > 10]
    df_plot = df_plot.T
    df_plot[df_plot < 10] = 0
    if df_plot.size == 0:
        print("No association")
        return
    mat_annot = np.zeros(df_plot.shape, dtype=object)
    for i_col, col in enumerate(df_plot.columns):
        mat_annot[:, i_col] = [num2str(x) for x in df_plot[col].values]
    df_plot = np.log10(df_plot + 1)

    plt.figure(figsize=[0.4 * df_plot.shape[1] + 3, 0.25 * df_plot.shape[0] + 3])
    sns.heatmap(df_plot, annot=mat_annot, fmt="s", cbar=False)
    plt.xticks(
        np.arange(df_plot.shape[1]) + 0.5, df_plot.columns, rotation=45, ha="right"
    )
    plt.show()

def plot_heatmap(
    df,
    dpi=150,
    squaresize=20,
    heatmap_annot=None,
    heatmap_annot_kws={"color": "black", "size": 4},
    heatmap_linewidths=0.8,
    heatmap_linecolor="black",
    
    heatmap_xticklabels=True,
    heatmap_yticklabels=True,
    heatmap_cbar=True,
    heatmap_cbar_kws=dict(use_gridspec=False, location="top", fraction=0.03, pad=0.01),
    heatmap_vmin=0.0,
    heatmap_vmax=1.0,
    xticklabels_rotation=90,
    colormap_n_bin=5,
):
    figwidth = df.shape[1] * squaresize / float(dpi)
    figheight = df.shape[0] * squaresize / float(dpi)
    fig, ax = plt.subplots(1, figsize=(figwidth, figheight), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_facecolor("white")
    sns.heatmap(
        df,
        annot=heatmap_annot,
        annot_kws=heatmap_annot_kws,
        fmt="",
        cmap=discrete_cmap(colormap_n_bin, "Reds"),
        linewidths=heatmap_linewidths,
        linecolor=heatmap_linecolor,
        square=True,
        ax=ax,
        xticklabels=heatmap_xticklabels,
        yticklabels=heatmap_yticklabels,
        cbar=heatmap_cbar,
        cbar_kws=heatmap_cbar_kws,
        vmin=heatmap_vmin,
        vmax=heatmap_vmax,
    )

    plt.yticks(fontsize=6)
    plt.xticks(fontsize=6)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=xticklabels_rotation,
        va="top",
        ha="right",
        fontsize=6,
    )
    ax.tick_params(left=False, bottom=False, pad=-2)
    trans = mtrans.Affine2D().translate(5, 0)
    for t in ax.get_xticklabels():
        t.set_transform(t.get_transform() + trans)
    return fig, ax


## PBMC
dict_df_stats = {
    trait: pd.read_csv(f"./out/step5-PBMC/{trait}.scdrs_group.celltype", sep="\t", index_col=0)
    for trait in ["T1D"]
}
df = dict_df_stats['T1D']
df['df_fdr_prop'] = df["n_fdr_0.1"]/df["n_cell"]
df['df_assoc_fdr'] = multipletests(df['assoc_mcp'], method="fdr_bh")[1]
df['df_hetero_fdr'] = multipletests(df["hetero_mcp"], method="fdr_bh")[1]
df = df[['df_fdr_prop','assoc_mcp','df_assoc_fdr','hetero_mcp','df_hetero_fdr']]
trait_list = list(dict_df_stats.keys())
df_fdr_prop= pd.concat(
            [
                dict_df_stats[trait]["n_fdr_0.1"] / dict_df_stats[trait]["n_cell"]
                for trait in trait_list
            ],
            axis=1,
        ).T
#df_fdr_prop = df_fdr_prop[list_order]
df_assoc_fdr = pd.concat(
            [dict_df_stats[trait]["assoc_mcp"] for trait in trait_list], axis=1
        ).T
df_assoc_fdr = pd.DataFrame(
            multipletests(df_assoc_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_assoc_fdr.shape
            ),
            index=df_assoc_fdr.index,
            columns=df_assoc_fdr.columns,
        )
#df_assoc_fdr=df_assoc_fdr[list_order]
df_hetero_fdr = pd.concat(
            [dict_df_stats[trait]["hetero_mcp"] for trait in trait_list], axis=1
        ).T
df_hetero_fdr = pd.DataFrame(
            multipletests(df_hetero_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_hetero_fdr.shape
            ),
            index=df_hetero_fdr.index,
            columns=df_hetero_fdr.columns,
        )
#df_hetero_fdr = df_hetero_fdr[list_order]
df_fdr_prop.index = trait_list
df_assoc_fdr.index = trait_list
df_hetero_fdr.index = trait_list

df_assoc_fdr = df_assoc_fdr.applymap(lambda x: "×" if x < 0.01 else "")
df_assoc_fdr[df_assoc_fdr > 0.01] = ""

fig, ax = plot_heatmap(
        df_fdr_prop,
        squaresize=100,
        heatmap_annot=df_assoc_fdr,
        heatmap_annot_kws={"color": "black", "size": 8},
        heatmap_cbar_kws=dict(
            use_gridspec=False, location="top", fraction=0.04, pad=0.7, drawedges=True
        ),
        heatmap_vmin=0,
        heatmap_vmax=0.2,
        colormap_n_bin=10,
    )
cb = ax.collections[0].colorbar
cb.ax.tick_params(labelsize=8)
cb.ax.set_title("Prop. of sig. cells (FDR < 0.1)", fontsize=6)
cb.outline.set_edgecolor("black")
cb.outline.set_linewidth(1)
plt.tight_layout()


# Pancreatic islet
dict_df_stats = {
    trait: pd.read_csv(f"./out/step5-Islet/{trait}.scdrs_group.celltype", sep="\t", index_col=0)
    for trait in ["T1D"]
}
df = dict_df_stats['T1D']
df['df_fdr_prop'] = df["n_fdr_0.1"]/df["n_cell"]
df['df_assoc_fdr'] = multipletests(df['assoc_mcp'], method="fdr_bh")[1]
df['df_hetero_fdr'] = multipletests(df["hetero_mcp"], method="fdr_bh")[1]
df = df[['df_fdr_prop','assoc_mcp','df_assoc_fdr','hetero_mcp','df_hetero_fdr']]
trait_list = list(dict_df_stats.keys())
df_fdr_prop= pd.concat(
            [
                dict_df_stats[trait]["n_fdr_0.1"] / dict_df_stats[trait]["n_cell"]
                for trait in trait_list
            ],
            axis=1,
        ).T
#df_fdr_prop = df_fdr_prop[list_order]
df_assoc_fdr = pd.concat(
            [dict_df_stats[trait]["assoc_mcp"] for trait in trait_list], axis=1
        ).T
df_assoc_fdr = pd.DataFrame(
            multipletests(df_assoc_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_assoc_fdr.shape
            ),
            index=df_assoc_fdr.index,
            columns=df_assoc_fdr.columns,
        )
#df_assoc_fdr=df_assoc_fdr[list_order]
df_hetero_fdr = pd.concat(
            [dict_df_stats[trait]["hetero_mcp"] for trait in trait_list], axis=1
        ).T
df_hetero_fdr = pd.DataFrame(
            multipletests(df_hetero_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_hetero_fdr.shape
            ),
            index=df_hetero_fdr.index,
            columns=df_hetero_fdr.columns,
        )
#df_hetero_fdr = df_hetero_fdr[list_order]
df_fdr_prop.index = trait_list
df_assoc_fdr.index = trait_list
df_hetero_fdr.index = trait_list

df_assoc_fdr = df_assoc_fdr.applymap(lambda x: "×" if x < 0.01 else "")
df_assoc_fdr[df_assoc_fdr > 0.01] = ""

fig, ax = plot_heatmap(
        df_fdr_prop,
        squaresize=100,
        heatmap_annot=df_assoc_fdr,
        heatmap_annot_kws={"color": "black", "size": 8},
        heatmap_cbar_kws=dict(
            use_gridspec=False, location="top", fraction=0.04, pad=0.7, drawedges=True
        ),
        heatmap_vmin=0,
        heatmap_vmax=0.2,
        colormap_n_bin=10,
    )
cb = ax.collections[0].colorbar
cb.ax.tick_params(labelsize=8)
cb.ax.set_title("Prop. of sig. cells (FDR < 0.1)", fontsize=6)
cb.outline.set_edgecolor("black")
cb.outline.set_linewidth(1)
plt.tight_layout()