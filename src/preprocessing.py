"""Prétraitement des données GSE5281 : métadonnées, matrice d'expression,
transformation log2 et normalisation inter-échantillons.

Le pipeline suit l'ordre standard pour des microarrays Affymetrix MAS5 :
    1. construire la matrice d'expression (sondes × échantillons) ;
    2. transformer en log2 (les intensités MAS5 sont linéaires et très
       asymétriques) ;
    3. normaliser les distributions entre échantillons (quantile
       normalization) pour rendre les puces comparables.

Le choix du filtrage des sondes (ABS_CALL) et de la collapse sonde → gène
est laissé en aval, car il relève de décisions méthodologiques.
"""
from __future__ import annotations

import re
import unicodedata

import numpy as np
import pandas as pd

# Correspondance "Disease State" GEO → étiquette de groupe analytique.
GROUP_MAP = {
    "normal": "Control",
    "alzheimer's disease": "AD",
}


def clean_text(s: str) -> str:
    """Normalise une chaîne issue des métadonnées GEO.

    Les valeurs GSE5281 traînent des caractères non-ASCII parasites en fin
    de champ (ex. ``'Human\\xa0'``). On normalise en NFKC puis on retire tout
    caractère hors ASCII imprimable, et on rogne les espaces.
    """
    s = unicodedata.normalize("NFKC", s)
    return re.sub(r"[^\x20-\x7E]", "", s).strip()


def _parse_age(value: str) -> float:
    """Extrait l'âge numérique d'un champ libre comme ``'63 years'``."""
    m = re.search(r"\d+", value)
    return float(m.group()) if m else np.nan


def build_sample_metadata(gse) -> pd.DataFrame:
    """Construit une table de métadonnées (une ligne par échantillon GSM).

    Gère les deux casses de clés présentes dans GSE5281 (Title Case et
    minuscules) en normalisant toutes les clés en minuscules.

    Returns
    -------
    pandas.DataFrame indexé par ``sample_id`` avec les colonnes
    ``title, region, disease, group, sex, age``.
    """
    rows = []
    for name, gsm in gse.gsms.items():
        fields = {}
        for ch in gsm.metadata.get("characteristics_ch1", []):
            if ":" in ch:
                key, val = ch.split(":", 1)
                fields[clean_text(key).lower()] = clean_text(val)
        disease = fields.get("disease state", "")
        rows.append(
            {
                "sample_id": name,
                "title": clean_text(gsm.metadata.get("title", [""])[0]),
                "region": fields.get("organ region", ""),
                "disease": disease,
                "group": GROUP_MAP.get(disease.lower(), "Unknown"),
                "sex": fields.get("sex", "").lower(),
                "age": _parse_age(fields.get("age", "")),
            }
        )
    return pd.DataFrame(rows).set_index("sample_id")


def filter_region(meta: pd.DataFrame, region: str = "Entorhinal Cortex") -> pd.DataFrame:
    """Filtre la table de métadonnées sur une région cérébrale (insensible à la casse)."""
    return meta[meta["region"].str.lower() == region.lower()].copy()


def build_expression_matrix(gse, sample_ids, value_col: str = "VALUE") -> pd.DataFrame:
    """Assemble la matrice d'expression (sondes × échantillons) pour `sample_ids`.

    S'appuie sur ``gse.pivot_samples`` qui aligne tous les échantillons sur
    l'index des sondes, puis restreint aux colonnes demandées.
    """
    matrix = gse.pivot_samples(value_col)[list(sample_ids)]
    matrix.index.name = "probe"
    return matrix


def log2_transform(matrix: pd.DataFrame, offset: float = 1.0) -> pd.DataFrame:
    """Transforme en log2(x + offset) pour stabiliser la variance.

    L'offset évite ``log2(0)`` ; ici les valeurs sont strictement positives
    mais l'offset reste une sécurité standard.
    """
    return np.log2(matrix + offset)


def quantile_normalize(matrix: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalization : aligne la distribution de chaque échantillon
    sur une distribution de référence commune.

    Principe : on trie les valeurs de chaque colonne, on moyenne ligne à ligne
    ces colonnes triées (= distribution de référence), puis on remappe chaque
    valeur d'origine via son rang sur cette référence. Après normalisation,
    tous les échantillons partagent exactement la même distribution.
    """
    reference = np.sort(matrix.to_numpy(), axis=0).mean(axis=1)
    positions = np.arange(1, len(reference) + 1)
    ranks = matrix.rank(method="average", axis=0)
    normed = ranks.apply(lambda col: np.interp(col, positions, reference), axis=0)
    return pd.DataFrame(normed, index=matrix.index, columns=matrix.columns)


def filter_by_abscall(
    expr: pd.DataFrame, abscall: pd.DataFrame, min_present: int
) -> pd.DataFrame:
    """Conserve les sondes « Present » (call ``P``) dans au moins `min_present`
    échantillons.

    Le seuil `min_present` est typiquement calé sur la taille du plus petit
    groupe : une sonde marqueur d'un seul groupe (ex. surexprimée uniquement
    chez les patients) reste ainsi conservée.

    Parameters
    ----------
    expr : matrice d'expression (sondes × échantillons).
    abscall : matrice des detection calls, mêmes index/colonnes que `expr`.
    min_present : seuil k (nombre minimal d'échantillons « Present »).
    """
    abscall = abscall.reindex(index=expr.index, columns=expr.columns)
    present_counts = (abscall == "P").sum(axis=1)
    kept = present_counts[present_counts >= min_present].index
    return expr.loc[kept]


def collapse_probes_to_genes(
    expr: pd.DataFrame,
    probe_to_symbol: pd.Series,
    multi_sep: str = "///",
) -> pd.DataFrame:
    """Réduit la matrice sondes × échantillons à une matrice gènes × échantillons.

    Stratégie « sonde la plus exprimée » (max mean) : pour chaque gène, on
    conserve la sonde dont l'intensité moyenne (sur tous les échantillons) est
    la plus forte — celle au meilleur rapport signal/bruit.

    Les sondes sont écartées proprement (jamais mappées au hasard) si :
      - elles n'ont **pas** de symbole de gène (NaN ou vide) ;
      - elles portent **plusieurs** symboles (ex. ``'DDR1 /// MIR4640'``),
        car l'attribution à un gène unique serait arbitraire.

    Parameters
    ----------
    expr : matrice d'expression (sondes × échantillons).
    probe_to_symbol : Series indexée par ID de sonde → symbole de gène.
    multi_sep : séparateur des symboles multiples dans l'annotation.

    Returns
    -------
    pandas.DataFrame indexé par symbole de gène (un gène = une ligne).
    """
    symbols = probe_to_symbol.reindex(expr.index).astype("string").str.strip()

    # Sondes valides : symbole présent, non vide, et non multiple.
    valid = symbols.notna() & (symbols != "") & (~symbols.str.contains(multi_sep, regex=False))
    expr_valid = expr.loc[valid]
    symbols_valid = symbols.loc[valid]

    # Pour chaque gène, sélectionner la sonde d'intensité moyenne maximale.
    mean_expr = expr_valid.mean(axis=1)
    best_probe = mean_expr.groupby(symbols_valid.to_numpy()).idxmax()  # gène -> ID sonde

    gene_matrix = expr_valid.loc[best_probe.to_numpy()].copy()
    gene_matrix.index = best_probe.index
    gene_matrix.index.name = "gene"
    return gene_matrix.sort_index()
