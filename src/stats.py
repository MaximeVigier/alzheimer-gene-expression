"""Analyse d'expression différentielle gène par gène.

Pour chaque gène, on compare l'expression (log2, normalisée) entre le groupe
Alzheimer et le groupe contrôle :

* **Test statistique** : t-test de Welch (Student à variances inégales). On ne
  suppose pas l'égalité des variances entre groupes — hypothèse rarement
  vérifiée en pratique, et les deux groupes ont des effectifs différents
  (10 vs 13). Welch est plus robuste que le t-test classique dans ce cas, sans
  coût lorsque les variances sont en réalité égales.
* **log fold-change** : différence des moyennes en échelle log2
  (``moyenne_AD − moyenne_contrôle``). Comme les données sont déjà en log2,
  cette différence EST directement le log2 fold-change. Convention : une valeur
  positive = gène **surexprimé** chez les patients Alzheimer.
* **Correction des tests multiples** : ~10 000+ gènes testés simultanément
  gonflent les faux positifs. On contrôle le **taux de fausses découvertes**
  (FDR) par la procédure de Benjamini-Hochberg.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def differential_expression(
    expr: pd.DataFrame,
    samples_group: list[str],
    samples_ref: list[str],
    group_label: str = "AD",
    ref_label: str = "Control",
) -> pd.DataFrame:
    """Analyse différentielle de chaque ligne (gène) entre deux groupes.

    Parameters
    ----------
    expr : matrice gènes × échantillons (valeurs log2, normalisées).
    samples_group : colonnes du groupe d'intérêt (ex. patients AD).
    samples_ref : colonnes du groupe de référence (ex. contrôles).
    group_label, ref_label : étiquettes, pour la documentation des colonnes.

    Returns
    -------
    pandas.DataFrame indexé par gène, trié par p-value ajustée croissante :
        ``mean_<group>, mean_<ref>, log2FC, t_stat, p_value, p_adj``.
    Le ``log2FC`` est défini comme ``mean_<group> − mean_<ref>``.
    """
    a = expr[samples_group].to_numpy()
    b = expr[samples_ref].to_numpy()

    mean_a = a.mean(axis=1)
    mean_b = b.mean(axis=1)
    log2fc = mean_a - mean_b  # données déjà en log2

    # t-test de Welch (equal_var=False), ligne par ligne.
    t_stat, p_value = stats.ttest_ind(a, b, axis=1, equal_var=False)

    # Correction Benjamini-Hochberg (contrôle du FDR).
    p_adj = stats.false_discovery_control(p_value, method="bh")

    results = pd.DataFrame(
        {
            f"mean_{group_label}": mean_a,
            f"mean_{ref_label}": mean_b,
            "log2FC": log2fc,
            "t_stat": t_stat,
            "p_value": p_value,
            "p_adj": p_adj,
        },
        index=expr.index,
    )
    return results.sort_values("p_adj")


def summarize_hits(
    results: pd.DataFrame,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.Series:
    """Compte les gènes significatifs selon le FDR et un seuil de |log2FC|.

    ``lfc_threshold = 1.0`` correspond à un fold-change de 2× (échelle linéaire).
    """
    sig = results["p_adj"] < alpha
    strong = results["log2FC"].abs() >= lfc_threshold
    return pd.Series(
        {
            "genes_testes": len(results),
            f"FDR<{alpha}": int(sig.sum()),
            f"FDR<{alpha} & |log2FC|>={lfc_threshold}": int((sig & strong).sum()),
            "surexprimes_AD": int((sig & strong & (results["log2FC"] > 0)).sum()),
            "sousexprimes_AD": int((sig & strong & (results["log2FC"] < 0)).sum()),
        }
    )
