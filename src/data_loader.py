"""Téléchargement et mise en cache du dataset GEO GSE5281.

GEOparse récupère le fichier SOFT compressé (``.soft.gz``) depuis les
serveurs du NCBI GEO. Ce fichier est mis en cache dans ``data/raw/`` :
tout appel ultérieur réutilise la copie locale au lieu de re-télécharger
(le fichier pèse plusieurs dizaines de Mo).
"""
from __future__ import annotations

from pathlib import Path

import GEOparse

# Racine du projet = dossier parent de src/ ; permet d'appeler ce module
# depuis n'importe où (notebooks/, scripts, tests) sans casser les chemins.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"

GSE_ID = "GSE5281"


def load_gse(geo_id: str = GSE_ID, destdir: Path | str = RAW_DIR, silent: bool = True):
    """Charge un GSE depuis le cache local ; le télécharge s'il est absent.

    Parameters
    ----------
    geo_id : str
        Identifiant GEO du dataset (par défaut ``GSE5281``).
    destdir : Path | str
        Dossier de cache des fichiers bruts.
    silent : bool
        Réduit la verbosité de GEOparse.

    Returns
    -------
    GEOparse.GSE
        L'objet GSE contenant métadonnées, plateforme(s) (GPL) et
        échantillons (GSM).
    """
    destdir = Path(destdir)
    destdir.mkdir(parents=True, exist_ok=True)
    return GEOparse.get_GEO(geo=geo_id, destdir=str(destdir), silent=silent)


if __name__ == "__main__":
    gse = load_gse()
    print(f"{GSE_ID} chargé : {len(gse.gsms)} échantillons, "
          f"{len(gse.gpls)} plateforme(s).")
