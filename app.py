import pandas as pd
import streamlit as st
import networkx as nx
from gseapy import enrichr
from pyvis.network import Network
import tempfile
import streamlit.components.v1 as components
from bs4 import BeautifulSoup

st.set_page_config(page_title="Int_Enridb App", layout="wide")
st.image("https://raw.githubusercontent.com/priyadarshinikp1/Multiomics-Integrator-app/main/logo.png", width=200)
st.title("🧬 Omics Integration & Over-Representation Explorer - Vizzhy")
with st.sidebar:
    st.markdown("---")
    st.markdown("**👨‍💻 Created by:PRIYADARSHINI")
    st.markdown("[LinkedIn](www.linkedin.com/in/priyadarshini24) | [GitHub](https://github.com/priyadarshinikp1)")

# --- File Upload Section ---
st.sidebar.header("📁 Upload Your Omics Data")
genomics_file = st.sidebar.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics_file = st.sidebar.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics_file = st.sidebar.file_uploader("Upload Proteomics CSV", type="csv")

if not all([genomics_file, transcriptomics_file, proteomics_file]):
    st.warning("Please upload all three omics files to proceed.")
    st.stop()

# --- Load Datasets ---
genomics_data = pd.read_csv(genomics_file)
transcriptomics_data = pd.read_csv(transcriptomics_file)
proteomics_data = pd.read_csv(proteomics_file)

# --- Extract and Process Genes ---
genomics_genes = set(genomics_data['Gene'].str.upper())
transcriptomics_genes = set(transcriptomics_data['Gene'].str.upper())
proteomics_genes = set(proteomics_data['Gene'].str.upper())

common_genes = list(genomics_genes & transcriptomics_genes & proteomics_genes)
st.success(f"✅ Found {len(common_genes)} common genes across all omics layers.")

# --- Over-Representation Libraries ---
libraries = [
    "HMDB_Metabolites", "DisGeNET", "OMIM_Disease",
    "KEGG_2021_Human", "Reactome_2016",
    "TRANSFAC_and_JASPAR_PWMs", "GO_Biological_Process_2021", "KEA_2015"
]

lib_to_type = {
    "HMDB_Metabolites": "metabolite",
    "DisGeNET": "disease",
    "OMIM_Disease": "disease",
    "KEGG_2021_Human": "pathway",
    "Reactome_2016": "pathway",
    "TRANSFAC_and_JASPAR_PWMs": "regulator",
    "GO_Biological_Process_2021": "process",
    "KEA_2015": "enzyme"
}

# --- Perform Over-Representation Analysis ---
st.header("🧠 Over-Representation Analysis Results")
results = {}
for lib in libraries:
    try:
        ora = enrichr(gene_list=common_genes, gene_sets=lib, outdir=None)
        results[lib] = ora.results
        st.success(f"✓ Over-Representation completed: {lib}")
    except Exception as e:
        st.warning(f"[ERROR] {lib}: {e}")

# --- Display Top ORA Terms ---
for lib, df in results.items():
    st.subheader(f"{lib} - Top Results")
    st.dataframe(df[['Term', 'P-value', 'Adjusted P-value', 'Genes']].head(10))

# === NETWORK ===
top_n = st.slider("Select number of top terms to visualize:", 5, 50, 5)
degree_threshold = st.slider("Minimum node degree to include:", 1, 10, 1)

filtered_results = {
    lib: df.sort_values("P-value").head(top_n)
    for lib, df in results.items()
}

legend_items = {
    "Gene": 'rgb(169,169,169)',
    "Protein": 'rgb(138,43,226)',
    "Enzyme": 'rgb(255,165,0)',
    "Metabolite": 'rgb(152,251,152)',
    "Pathway": 'rgb(135,206,250)',
    "Process": 'rgb(255,182,193)',
    "Disease": 'rgb(255,99,71)',
    "Regulator": 'rgb(205,133,63)',
}

type_color_map = {
    "gene": legend_items["Gene"],
    "protein": legend_items["Protein"],
    "enzyme": legend_items["Enzyme"],
    "metabolite": legend_items["Metabolite"],
    "pathway": legend_items["Pathway"],
    "process": legend_items["Process"],
    "disease": legend_items["Disease"],
    "regulator": legend_items["Regulator"],
    "term": "rgb(200,200,200)"
}

st.subheader("🧩 Interactive Omics Network")

try:
    net = Network(height='800px', width='100%', notebook=False, directed=False)
    net.force_atlas_2based()
    net.show_buttons(filter_=['physics'])

    y_pos = 0
    for label, color in legend_items.items():
        net.add_node(f"legend_{label}", label=label, shape='box', color=color,
                     size=20, x=-1000, y=y_pos, physics=False, fixed=True)
        y_pos -= 50

    temp_graph = nx.Graph()
    for gene in common_genes:
        temp_graph.add_node(gene, type="gene")

    for lib, df in filtered_results.items():
        node_type = lib_to_type.get(lib, "term")
        for _, row in df.iterrows():
            term = row['Term']
            genes = [g.strip().upper() for g in row['Genes'].split(';')]
            temp_graph.add_node(term, type=node_type)
            for gene in genes:
                if gene in common_genes:
                    temp_graph.add_edge(gene, term)

    for _, row in proteomics_data.iterrows():
        gene = row['Gene'].strip().upper()
        protein = row['Protein'].strip()
        if gene in temp_graph.nodes:
            temp_graph.add_node(protein, type="protein")
            temp_graph.add_edge(gene, protein)

    nodes_to_keep = [n for n, d in temp_graph.degree() if d >= degree_threshold]
    temp_graph = temp_graph.subgraph(nodes_to_keep)

    for node in temp_graph.nodes:
        n_type = temp_graph.nodes[node].get("type", "term")
        net.add_node(
            node,
            label=node,
            color=type_color_map.get(n_type, "gray"),
            size=15 if n_type == "gene" else 25
        )

    for s, t in temp_graph.edges:
        net.add_edge(s, t)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
        net.save_graph(tmp_file.name)

        with open(tmp_file.name, 'r', encoding='utf-8') as f:
            html = f.read()

        soup = BeautifulSoup(html, 'html.parser')
        script_tag = soup.find_all("script")[-1]

        highlight_js = """
        <script type=\"text/javascript\">
        network.on(\"click\", function (params) {
            if (params.nodes.length === 0) return;
            var clickedNodeId = params.nodes[0];
            nodes.update(
                nodes.get().map(function (node) {
                    return {
                        id: node.id,
                        color: node.originalColor || node.color
                    };
                })
            );
            var connectedNodes = network.getConnectedNodes(clickedNodeId);
            connectedNodes.push(clickedNodeId);
            nodes.update(
                connectedNodes.map(function (id) {
                    var node = nodes.get(id);
                    node.originalColor = node.color;
                    return {
                        id: id,
                        color: '#FFFF00'
                    };
                })
            );
        });
        </script>
        """
        soup.body.append(BeautifulSoup(highlight_js, 'html.parser'))

        updated_html_path = tmp_file.name.replace(".html", "_highlight.html")
        with open(updated_html_path, "w", encoding='utf-8') as f:
            f.write(str(soup))

        components.html(open(updated_html_path, 'r', encoding='utf-8').read(), height=800)

except Exception as e:
    st.error(f"Network rendering failed: {e}")



# === SUMMARY TABLE OF ALL ASSOCIATIONS ===
st.subheader("🧾 Summary Table of Gene Associations")

# Initialize summary dictionary
summary_dict = {gene: {
    "Transcription Factor": [],
    "Protein": [],
    "Enzyme": [],
    "Metabolite": [],
    "Pathway": [],
    "Process": [],
    "Disease": []
} for gene in common_genes}

# Fill in the enrichment-based associations
# Fill in the enrichment-based associations
for lib, df in results.items():
    assoc_type = lib_to_type.get(lib, None)
    if not assoc_type:
        continue
    for _, row in df.iterrows():
        term = row["Term"]
        # Filter out mouse or non-human terms
        if any(x in term.lower() for x in ["mouse", "mus musculus", "mmu", "murine"]):
            continue
        genes = [g.strip().upper() for g in row["Genes"].split(";")]
        for gene in genes:
            if gene in summary_dict:
                if assoc_type == "regulator":
                    summary_dict[gene]["Transcription Factor"].append(term)
                elif assoc_type == "enzyme":
                    summary_dict[gene]["Enzyme"].append(term)
                elif assoc_type == "metabolite":
                    summary_dict[gene]["Metabolite"].append(term)
                elif assoc_type == "pathway":
                    summary_dict[gene]["Pathway"].append(term)
                elif assoc_type == "process":
                    summary_dict[gene]["Process"].append(term)
                elif assoc_type == "disease":
                    summary_dict[gene]["Disease"].append(term)


# Fill in proteomics-based associations
for _, row in proteomics_data.iterrows():
    gene = row['Gene'].strip().upper()
    protein = row['Protein'].strip()
    if gene in summary_dict:
        summary_dict[gene]["Protein"].append(protein)

# Convert to DataFrame
summary_df = pd.DataFrame.from_dict(summary_dict, orient='index').reset_index()
summary_df.rename(columns={'index': 'Gene'}, inplace=True)

# Optional: Join list values into semicolon-separated strings
for col in summary_df.columns[1:]:
    summary_df[col] = summary_df[col].apply(lambda x: '; '.join(set(x)) if isinstance(x, list) else '')

st.dataframe(summary_df)

st.markdown("### 🔍 Filter & Customize Summary Table Display")

# Search bar for disease terms
search_term = st.text_input("Search for disease or term (case-insensitive):", "")

# Option to show top N or all genes
max_rows = st.slider("Select number of genes to display (use slider to show top N):", 
                     min_value=5, max_value=len(summary_df), value=10, step=1)

# Optional sorting: by number of associations
summary_df["Total Associations"] = summary_df.iloc[:, 1:].apply(lambda x: sum(bool(i) for i in x), axis=1)

# Filter by search term (e.g., disease name)
if search_term:
    filtered_df = summary_df[summary_df.apply(
        lambda row: any(search_term.lower() in str(row[col]).lower() for col in summary_df.columns[1:]),
        axis=1)]
else:
    filtered_df = summary_df

# Sort and slice
filtered_df = filtered_df.sort_values("Total Associations", ascending=False).drop(columns=["Total Associations"])
filtered_df = filtered_df.head(max_rows)

# Display the updated summary table
st.dataframe(filtered_df)
