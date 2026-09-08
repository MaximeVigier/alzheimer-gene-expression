import streamlit as st
import pandas as pd
import numpy as np
from PIL import Image
import os

# Set page config
st.set_page_config(
    page_title="Alzheimer Gene Expression Analysis",
    page_icon="🧠",
    layout="wide"
)

# Title and introduction
st.title("🧠 Alzheimer's Disease Gene Expression Analysis")
st.markdown("""
This interactive dashboard presents the results of a differential gene expression analysis in the entorhinal cortex 
of Alzheimer's disease (AD) patients compared to healthy aged controls. The analysis uses data from GSE5281 (Liang et al., 2008).

The study focuses on the earliest region affected by AD pathology - the entorhinal cortex - to capture a clear biological signal.
""")

# Load figures
figures_dir = "results/figures"

# Display volcano plot
st.subheader("Volcano Plot")
volcano_img_path = os.path.join(figures_dir, "volcano_entorhinal.png")
if os.path.exists(volcano_img_path):
    volcano_img = Image.open(volcano_img_path)
    st.image(volcano_img, caption="Volcano plot of differential gene expression", use_column_width=True)
else:
    st.warning("Volcano plot image not found")

# Display heatmap
st.subheader("Heatmap - Top Differential Genes")
heatmap_img_path = os.path.join(figures_dir, "heatmap_top_genes.png")
if os.path.exists(heatmap_img_path):
    heatmap_img = Image.open(heatmap_img_path)
    st.image(heatmap_img, caption="Heatmap of top differentially expressed genes", use_column_width=True)
else:
    st.warning("Heatmap image not found")

# Display PCA
st.subheader("PCA of Samples")
pca_img_path = os.path.join(figures_dir, "pca_samples.png")
if os.path.exists(pca_img_path):
    pca_img = Image.open(pca_img_path)
    st.image(pca_img, caption="PCA plot showing separation of AD vs control samples", use_column_width=True)
else:
    st.warning("PCA image not found")

# Interactive table
st.subheader("Interactive Differential Expression Table")

# Load the CSV file
csv_path = "results/tables/de_results_entorhinal.csv"
if os.path.exists(csv_path):
    df = pd.read_csv(csv_path)
    
    # Show basic info about the data
    st.markdown(f"Total genes in analysis: {len(df)}")
    
    # Add filters for FDR and log2FC
    col1, col2 = st.columns(2)
    
    with col1:
        fdr_threshold = st.slider("FDR Threshold", 0.0, 0.1, 0.05, 0.01)
        
    with col2:
        log2fc_threshold = st.slider("log2FoldChange Threshold", 0.0, 3.0, 1.0, 0.1)
    
    # Apply filters
    filtered_df = df[(df['FDR'] < fdr_threshold) & (abs(df['log2FoldChange']) > log2fc_threshold)]
    
    # Show filtered table
    st.markdown(f"Showing {len(filtered_df)} genes with FDR < {fdr_threshold} and |log2FC| > {log2fc_threshold}")
    st.dataframe(filtered_df.reset_index(drop=True))
else:
    st.warning("Differential expression table not found")

# Canonical AD genes table
st.subheader("Canonical AD Genes")

canonical_genes_data = {
    'Gene': ['GFAP', 'CD44', 'SPP1', 'GABRA1', 'SST', 'NEFL', 'CALB1', 'RTN3', 'VGF'],
    'log2FC (AD − ctrl)': [1.78, 2.91, 2.91, -2.76, -2.11, -1.41, -1.38, -1.34, -1.18],
    'FDR': ['0.0014', '0.001', '<0.001', '0.0009', '0.0078', '0.0059', '0.0125', '0.0069', '0.0012'],
    'Expected role': ['Reactive astrogliosis', 'Glial / inflammatory', 'Glial / inflammatory', 'GABAergic synaptic', 'Interneuron marker', 'Neurofilament (neuronal)', 'Calcium-binding neuronal', 'Neuronal, APP processing', 'Neuropeptide / synaptic'],
    'Direction': ['↑ as expected', '↑ as expected', '↑ as expected', '↓ as expected', '↓ as expected', '↓ as expected', '↓ as expected', '↓ as expected', '↓ as expected']
}

canonical_genes_df = pd.DataFrame(canonical_genes_data)
st.dataframe(canonical_genes_df)

# Post-translational genes table
st.subheader("Post-Translational Genes (Not Differentially Expressed)")

post_translational_data = {
    'Gene': ['APP', 'MAPT', 'PSEN1', 'APOE'],
    'log2FC': [0.34, 0.00, -0.08, -0.41],
    'FDR': [0.080, 0.993, 0.841, 0.396],
    'Significant?': ['no', 'no', 'no', 'no']
}

post_translational_df = pd.DataFrame(post_translational_data)
st.dataframe(post_translational_df)

# Biological interpretation
st.subheader("Biological Interpretation")

st.markdown("""
### The Glial/Neuronal Switch

The transcriptomic signature observed in the entorhinal cortex reveals a dual, mirror-image dynamic:

1. **Neuronal collapse**: Neurons die, synapses are lost, and with them go their characteristic markers:
   - SST (somatostatin, inhibitory interneurons)
   - GABRA1 (GABA-A receptor)
   - NEFL (neurofilament light, axonal integrity)
   - CALB1 (calbindin, specific neuronal populations)
   - VGF (a neuropeptide involved in synaptic plasticity)

2. **Glial reaction**: Under assault, astrocytes activate and enter a state of reactive astrogliosis:
   - They proliferate, hypertrophy,
   - Overexpress GFAP (Glial Fibrillary Acidic Protein)
   
This is not a sign of repair — it is an inflammatory response to ongoing damage.
""")

# Methodological caveats
st.subheader("Methodological Caveats")

st.markdown("""
- Modest sample size (10 vs 13) 
- Microarray data (2007), less sensitive than current RNA-seq 
- Age/sex verified comparable between groups but not explicitly modelled  
- The high proportion of significant genes (~39 %) reflects a genuinely strong biological effect and the absence of a multivariate model
- For publication-grade analysis, `limma` (linear model + empirical-Bayes moderated variance) would be used

These trade-offs are deliberate for a pipeline-focused portfolio project.
""")

st.markdown("---")
st.caption("Data from NCBI GEO accession [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281). Code and dataset available on [GitHub](https://github.com/MaximeVigier/alzheimer-gene-expression)")