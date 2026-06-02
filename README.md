# Analyse d'expression génique différentielle — Maladie d'Alzheimer (GSE5281)

> 🚧 Projet en cours de construction. Le README complet (contexte biologique,
> méthode, résultats et visuels) sera rédigé à l'issue de l'analyse.

Analyse d'expression différentielle entre tissu cérébral de patients atteints
de la maladie d'Alzheimer et témoins sains âgés, à partir du dataset public
[GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)
(microarray Affymetrix HG-U133 Plus 2.0). L'analyse est restreinte au
**cortex entorhinal**, première structure cérébrale touchée dans Alzheimer.

## Stack
Python · pandas · numpy · scipy · matplotlib · seaborn · GEOparse

## Reproduire l'analyse
```bash
pip install -r requirements.txt
```
Puis exécuter les notebooks de `notebooks/` dans l'ordre (01 → 04). Les données
GEO sont téléchargées et mises en cache automatiquement dans `data/raw/`
(non versionné).

## Structure
```
notebooks/   Analyse narrée, étape par étape (acquisition → interprétation)
src/         Fonctions réutilisables (chargement, prétraitement, stats, plots)
data/        Données brutes et transformées (non versionnées)
results/     Figures et tables de gènes différentiels
```
