# Multi-Omics Integration App
import os
import tempfile
import requests
import numpy as np
import pandas as pd
import streamlit as st
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px

from pyvis.network import Network
import networkx as nx
from gseapy import enrichr
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap
from matplotlib.colors import Normalize

# -----------------------------
# App Configuration
# -----------------------------
st.set_page_config(page_title="Multi-Omics App", layout="wide")
st.image("https://raw.githubusercontent.com/priyadarshinikp1/Multiomics-Integrator-app/main/logo.png", width=200)
st.title("🧬 Multi-Omics Integration Vizzhy App")

with st.sidebar:
    st.markdown("---")
    st.markdown("**👨‍💻 Created by: PRIYADARSHINI**")
    st.markdown("[LinkedIn](https://www.linkedin.com/in/priyadarshini24) | [GitHub](https://github.com/priyadarshinikp1)")


# -----------------------------
# File Upload Section
# -----------------------------
st.header("📁 Upload Omics Data")

genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")


if genomics:
    gdf = pd.read_csv(genomics)

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)

if proteomics:
    pdf = pd.read_csv(proteomics)

# -----------------------------
# Sidebar Filters
# -----------------------------
st.sidebar.header("⚙️ Settings")

cadd_thresh = float(st.sidebar.text_input("Min CADD Score (Genomics)", value="20"))
logfc_thresh = float(st.sidebar.text_input("Min |logFC| (Transcriptomics)", value="1"))
t_pval_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="1000"))

run_enrichment = st.sidebar.checkbox("Run Enrichment Analyses", value=True)
show_network = st.sidebar.checkbox("Show Network Visualization", value=True)
show_association_table = st.sidebar.checkbox("Show Association Table", value=True)
num_pathways_to_show = st.sidebar.slider("Number of Pathways to Display in Network", min_value=1, max_value=100, value=10)

# -----------------------------
# Preview Filtered Data: Top N Rows
# -----------------------------
preview_n = st.sidebar.slider("Preview Top N Filtered Rows", 5, 50, 10)

# Display Filtered Data Previews
st.subheader("🔍 Filtered Data Preview")

if genomics and transcriptomics and proteomics:
    try:
        gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
        tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh)]
        
    # Filter Transcriptomics data based on logFC threshold
        if logfc_thresh > 0:
    # For positive logFC threshold, filter for values >= logFC_thresh
              tdf_filtered = tdf[tdf['logFC'] >= logfc_thresh]
        elif logfc_thresh < 0:
    # For negative logFC threshold, filter for values < logFC_thresh
              tdf_filtered = tdf[tdf['logFC'] < logfc_thresh]
        else:
    # If logFC_thresh is 0, no filtering needed
              tdf_filtered = tdf
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
        
        # Display filtered data for all three omics
        st.markdown("**Genomics**")
        st.dataframe(gdf_filtered.head(preview_n))  # Preview the top N filtered rows for genomics
        st.markdown("**Transcriptomics**")
        st.dataframe(tdf_filtered.head(preview_n))  # Preview the top N filtered rows for transcriptomics
        st.markdown("**Proteomics**")
        st.dataframe(pdf_filtered.head(preview_n))  # Preview the top N filtered rows for proteomics

    except Exception as e:
        st.error(f"Integration error: {e}")

# -----------------------------
# Filtering and Integration
# -----------------------------
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    try:
        gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
        tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh) & (tdf['logFC'].abs() >= logfc_thresh)]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

        union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])

        def extract_uniprot_ids(protein_series):
            ids = set()
            for entry in protein_series.dropna():
                for pid in str(entry).split(";"):
                    if pid.strip():
                        ids.add(pid.strip())
            return ids

        def map_uniprot_to_gene(uniprot_ids):
            mapping = {}
            ids = list(uniprot_ids)
            for i in range(0, len(ids), 100):
                chunk = ids[i:i+100]
                query = " OR ".join([f"accession:{id_}" for id_ in chunk])
                url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
                try:
                    r = requests.get(url)
                    if r.status_code == 200:
                        lines = r.text.strip().split('\n')[1:]
                        for line in lines:
                            acc, genes = line.split('\t')
                            mapping[acc] = genes.split()[0] if genes else acc
                except Exception as e:
                    st.warning(f"UniProt API error: {e}")
            return mapping

        unique_uniprot_ids = extract_uniprot_ids(pdf_filtered['Protein'])
        uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)

        expanded_rows = []
        for _, row in pdf_filtered.iterrows():
            for pid in str(row['Protein']).split(';'):
                pid = pid.strip()
                gene = uniprot_gene_map.get(pid)
                if gene:
                    expanded_rows.append({'Protein': pid, 'GeneName': gene})

        expanded_protein_df = pd.DataFrame(expanded_rows)
        protein_gene_map = dict(zip(expanded_protein_df['Protein'], expanded_protein_df['GeneName']))

        all_entities = union_genes | set(protein_gene_map.values())

        results = {}
        raw_assoc_data = []

        if run_enrichment:
            st.header("📊 Enrichment Analyses")
            libraries = {
                "Reactome Pathways": "Reactome_2016",
                "Disease Associations": "DisGeNET",
                "HMDB Metabolites": "HMDB_Metabolites"
            }

            for name, lib in libraries.items():
                try:
                    gene_list_clean = [str(g).strip() for g in union_genes if pd.notna(g)]
                    enr = enrichr(gene_list=gene_list_clean, gene_sets=lib, outdir=None)

                    if enr.results.empty:
                        continue

                    df = enr.results.copy()
                    df['-log10(pval)'] = -np.log10(df['P-value'])
                    df = df.rename(columns={"Term": "Pathway", "Genes": "Genes_Involved"})
                    results[name] = df

                    fig = px.bar(df.head(10), x="Pathway", y="-log10(pval)", title=f"Top 10 {name}")
                    st.plotly_chart(fig)
                except Exception as e:
                    st.error(f"Error in {name}: {e}")

        if show_network and results:
            st.subheader("🧠 Interactive Omics Network")
            net = Network(height='800px', width='100%', directed=False)
            net.force_atlas_2based()

            legend_items = {
                "Gene": 'gray', "Protein": 'gold',
                "Pathway": 'skyblue', "Metabolite": 'lightgreen', "Disease": 'lightcoral'
            }

            for i, (label, color) in enumerate(legend_items.items()):
                net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=-i*50, physics=False, fixed=True)

            color_map = {
                "Reactome Pathways": "skyblue",
                "Disease Associations": "lightcoral",
                "HMDB Metabolites": "lightgreen"
            }

            for name, df in results.items():
                color = color_map.get(name, "gray")
                for _, row in df.head(num_pathways_to_show).iterrows():
                    term = row['Pathway']
                    net.add_node(term, label=term, color=color)
                    for gene in row['Genes_Involved'].split(';'):
                        gene = gene.strip()
                        if not gene:
                            continue
                        net.add_node(gene, label=gene, color='gray')
                        net.add_edge(gene, term)
                        matched_proteins = [prot for prot, gname in protein_gene_map.items() if gname == gene]
                        for prot in matched_proteins:
                            net.add_node(prot, label=prot, color='gold')
                            net.add_edge(gene, prot)
                        raw_assoc_data.append({
                            'Gene': gene,
                            'Protein': ';'.join(matched_proteins),
                            'Pathway': term if name == 'Reactome Pathways' else '',
                            'Metabolite': term if name == 'HMDB Metabolites' else '',
                            'Disease': term if name == 'Disease Associations' else ''
                        })

            with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
                net.save_graph(tmp_file.name)
                st.components.v1.html(open(tmp_file.name, 'r', encoding='utf-8').read(), height=800)

        if show_association_table and raw_assoc_data:
            st.subheader("📄 Gene-Protein-Term Association Summary")
            df = pd.DataFrame(raw_assoc_data)
            assoc_df = df.groupby('Gene').agg({
                'Protein': lambda x: ';'.join(set(filter(None, x))),
                'Pathway': lambda x: ';'.join(set(filter(None, x))),
                'Disease': lambda x: ';'.join(set(filter(None, x))),
                'Metabolite': lambda x: ';'.join(set(filter(None, x)))
            }).reset_index()

            assoc_df['non_nulls'] = assoc_df.notnull().sum(axis=1)
            assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')
            st.dataframe(assoc_df)

    except Exception as e:
        st.error(f"Integration error: {e}")

# -----------------------------
# UMAP + KMeans Clustering
# -----------------------------
st.header("🧬 UMAP + KMeans Clustering (Multi-Omics)")

try:
    merged_df = pd.merge(gdf_filtered, tdf_filtered, on='Gene', how='inner')
    merged_df = pd.merge(merged_df, pdf_filtered, on='Gene', how='inner')

    genomics_data = merged_df[['CADD']].values
    transcriptomics_data = merged_df[['logFC', 'AveExpr', 'B']].values
    proteomics_data = merged_df[['Intensity']].values

    combined_data = np.concatenate([genomics_data, transcriptomics_data, proteomics_data], axis=1)
    scaler = StandardScaler()
    normalized_data = scaler.fit_transform(combined_data)

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.3, n_components=2, random_state=42)
    umap_results = reducer.fit_transform(normalized_data)

    n_clusters = st.sidebar.slider("Number of KMeans Clusters", min_value=2, max_value=10, value=5)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(umap_results)

    merged_df['Cluster'] = clusters

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(umap_results[:, 0], umap_results[:, 1], c=merged_df['CADD'], cmap='Spectral', alpha=0.7, s=40)
    ax.set_title('UMAP Projection Colored by CADD Score')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    fig.colorbar(scatter, label='CADD Score')
    st.pyplot(fig)

    fig2, ax2 = plt.subplots(figsize=(10, 8))
    sns.scatterplot(x=umap_results[:, 0], y=umap_results[:, 1], hue=clusters, palette='viridis', ax=ax2, s=60)
    ax2.set_title('UMAP Projection with KMeans Clustering')
    st.pyplot(fig2)

    st.subheader("🔍 Clustered Multi-Omics Data")
    st.dataframe(merged_df[['Gene', 'Cluster'] + [col for col in merged_df.columns if col not in ['Gene', 'Cluster']]].head(20))

except Exception as e:
    st.error(f"UMAP clustering error: {e}")
