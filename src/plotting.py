"""Visualisations pour l'analyse différentielle GSE5281.

Trois figures de référence :
    * ``volcano_plot``        — vue d'ensemble (ampleur vs significativité) ;
    * ``top_genes_heatmap``   — profil d'expression des gènes les plus marquants ;
    * ``pca_plot``            — structure globale des échantillons.

Particularité de ce projet : on annote nommément la **signature Alzheimer
canonique** (gènes connus de la pathologie), pour rendre les figures parlantes
au-delà du simple « top p-value ».
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Signature Alzheimer canonique du cortex entorhinal (gènes bien établis).
CANONICAL_AD_GENES = ["GFAP", "SST", "GABRA1", "NEFL", "CALB1", "VGF", "RTN3"]

# Palette cohérente pour tout le projet.
GROUP_PALETTE = {"Control": "#27ae60", "AD": "#c0392b"}
_UP_COLOR = "#c0392b"     # surexprimé chez AD
_DOWN_COLOR = "#2980b9"   # sous-exprimé chez AD
_NS_COLOR = "#b0b0b0"     # non significatif


def _classify(results: pd.DataFrame, alpha: float, lfc: float) -> pd.Series:
    """Étiquette chaque gène : 'up' / 'down' / 'ns'."""
    up = (results["p_adj"] < alpha) & (results["log2FC"] >= lfc)
    down = (results["p_adj"] < alpha) & (results["log2FC"] <= -lfc)
    cat = pd.Series("ns", index=results.index)
    cat[up] = "up"
    cat[down] = "down"
    return cat


def volcano_plot(
    results: pd.DataFrame,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    highlight_genes: list[str] | None = None,
    ax: plt.Axes | None = None,
):
    """Volcano plot : log2 fold-change (x) vs −log10(p-value) (y).

    Les gènes sont colorés selon leur statut (sur-/sous-exprimé chez AD vs non
    significatif), aux seuils ``FDR < alpha`` et ``|log2FC| >= lfc_threshold``.
    Les gènes de `highlight_genes` présents dans les résultats sont annotés
    nommément.
    """
    if highlight_genes is None:
        highlight_genes = CANONICAL_AD_GENES
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7))

    cat = _classify(results, alpha, lfc_threshold)
    y = -np.log10(results["p_value"].clip(lower=1e-300))
    colors = {"up": _UP_COLOR, "down": _DOWN_COLOR, "ns": _NS_COLOR}

    for label in ("ns", "down", "up"):  # ns au fond
        mask = cat == label
        ax.scatter(results.loc[mask, "log2FC"], y[mask], s=10, alpha=0.5,
                   c=colors[label], edgecolors="none",
                   label={"up": "↑ chez AD", "down": "↓ chez AD", "ns": "n.s."}[label])

    # Lignes de seuil
    p_thresh = results.loc[results["p_adj"] < alpha, "p_value"].max()
    if pd.notna(p_thresh):
        ax.axhline(-np.log10(p_thresh), color="grey", ls="--", lw=0.8)
    ax.axvline(lfc_threshold, color="grey", ls="--", lw=0.8)
    ax.axvline(-lfc_threshold, color="grey", ls="--", lw=0.8)

    # Annotation des gènes de la signature canonique
    for gene in highlight_genes:
        if gene in results.index:
            x_g = results.loc[gene, "log2FC"]
            y_g = -np.log10(max(results.loc[gene, "p_value"], 1e-300))
            ax.scatter(x_g, y_g, s=40, facecolors="none", edgecolors="black", linewidths=1.2)
            ax.annotate(gene, (x_g, y_g), xytext=(5, 4), textcoords="offset points",
                        fontsize=9, fontweight="bold")

    ax.set_xlabel("log2 fold-change (AD − contrôle)")
    ax.set_ylabel("−log10(p-value)")
    ax.set_title(f"Volcano plot — cortex entorhinal (AD vs contrôle)\n"
                 f"seuils : FDR < {alpha}, |log2FC| ≥ {lfc_threshold}")
    ax.legend(title="Statut", loc="upper right", framealpha=0.9)
    return ax


def top_genes_heatmap(
    gene_expr: pd.DataFrame,
    results: pd.DataFrame,
    meta: pd.DataFrame,
    n_top: int = 40,
    highlight_genes: list[str] | None = None,
):
    """Heatmap (z-score par gène) des gènes les plus différentiels.

    Sélectionne les `n_top` gènes les plus significatifs (FDR) puis y ajoute
    les gènes de la signature canonique (`highlight_genes`) s'ils n'y sont pas
    déjà. Les échantillons (colonnes) portent une barre de couleur par groupe ;
    les étiquettes des gènes canoniques sont mises en évidence (rouge gras).

    Returns
    -------
    seaborn.matrix.ClusterGrid
    """
    if highlight_genes is None:
        highlight_genes = CANONICAL_AD_GENES

    top = results.sort_values("p_adj").head(n_top).index.tolist()
    extra = [g for g in highlight_genes if g in gene_expr.index and g not in top]
    selected = top + extra

    data = gene_expr.loc[selected]
    # z-score par gène (ligne) : centre chaque gène pour comparer les profils.
    zdata = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

    col_colors = meta.loc[data.columns, "group"].map(GROUP_PALETTE)
    col_colors.name = "Groupe"

    g = sns.clustermap(
        zdata, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
        col_colors=col_colors, xticklabels=False, yticklabels=True,
        figsize=(11, max(8, 0.22 * len(selected))),
        cbar_kws={"label": "expression (z-score)"},
        dendrogram_ratio=(0.12, 0.08),
    )
    g.ax_heatmap.set_xlabel(f"{data.shape[1]} échantillons (colorés par groupe)")
    g.ax_heatmap.set_ylabel("")

    # Mettre en évidence les gènes canoniques dans les étiquettes
    highlight_set = set(highlight_genes)
    for tick in g.ax_heatmap.get_yticklabels():
        if tick.get_text() in highlight_set:
            tick.set_color(_UP_COLOR)
            tick.set_fontweight("bold")

    # Légende des groupes
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in GROUP_PALETTE.values()]
    g.ax_heatmap.legend(handles, GROUP_PALETTE.keys(), title="Groupe",
                        bbox_to_anchor=(1.02, 1.12), loc="upper left", frameon=False)
    g.figure.suptitle(f"Top {n_top} gènes différentiels + signature canonique\n"
                      "(z-score par gène ; gènes AD connus en rouge)", y=1.02)
    return g


def pca_plot(
    gene_expr: pd.DataFrame,
    meta: pd.DataFrame,
    n_top_var: int = 2000,
    ax: plt.Axes | None = None,
):
    """PCA des échantillons sur les `n_top_var` gènes les plus variables.

    On restreint aux gènes les plus variables (les plus informatifs), on
    standardise chaque gène (z-score), puis on décompose par SVD. On projette
    les échantillons sur les deux premières composantes, colorés par groupe.

    Returns
    -------
    (ax, scores_df, var_ratio) : axes, scores PC1/PC2 par échantillon, et part
    de variance expliquée par chaque composante.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6.5))

    top_var = gene_expr.var(axis=1).sort_values(ascending=False).head(n_top_var).index
    sub = gene_expr.loc[top_var]

    # Échantillons en lignes, gènes en colonnes ; standardisation par gène.
    X = sub.T.to_numpy()
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    scores = U[:, :2] * S[:2]
    var_ratio = (S ** 2 / (S ** 2).sum())[:2]

    scores_df = pd.DataFrame(scores, index=sub.columns, columns=["PC1", "PC2"])
    scores_df["group"] = meta.loc[scores_df.index, "group"].to_numpy()

    for grp, color in GROUP_PALETTE.items():
        m = scores_df["group"] == grp
        ax.scatter(scores_df.loc[m, "PC1"], scores_df.loc[m, "PC2"],
                   s=70, c=color, label=grp, edgecolors="white", linewidths=0.8)

    ax.set_xlabel(f"PC1 ({100*var_ratio[0]:.1f} % de variance)")
    ax.set_ylabel(f"PC2 ({100*var_ratio[1]:.1f} % de variance)")
    ax.set_title(f"PCA des échantillons (cortex entorhinal)\n{n_top_var} gènes les plus variables")
    ax.legend(title="Groupe")
    ax.axhline(0, color="grey", lw=0.5, alpha=0.5)
    ax.axvline(0, color="grey", lw=0.5, alpha=0.5)
    return ax, scores_df, var_ratio
