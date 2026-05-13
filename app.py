"""
Project: Betalain Pathway Explorer
Description: Interactive Streamlit App for Bioinformatics Sequence Analysis of DOD Gene in Mirabilis jalapa
Author: Zahidul Hasan
Department: Botany, University of Chittagong, Bangladesh
Email: zahidulhasan.botany.cu@gmail.com
LinkedIn: linkedin.com/in/zahidulhasan-botany-cu
Status: Professional Modular Version for GitHub/Streamlit
Date: May 2026
"""

import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
from stmol import showmol
import py3Dmol

# --- Page configuration ---
st.set_page_config(
    page_title="Betalain Pathway Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Custom Style ---
st.markdown("""
    <style>
    .main {
        background-color: #f5f7f9;
    }
    .stMetric {
        background-color: #ffffff;
        padding: 15px;
        border-radius: 10px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    h1 {
        color: #2c3e50;
        font-family: 'Helvetica Neue', sans-serif;
    }
    </style>
    """, unsafe_allow_html=True)

# --- Path Setup ---
BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# Ensure file paths are correct
GBK_FILE = os.path.join(DATA_DIR, "mirabilis_dod.gbk")
BLAST_CSV = os.path.join(RESULTS_DIR, "blast_results.csv")
FLOWER_IMG = os.path.join(DATA_DIR, "mirabilis_jalapa_flower.png")
LOGO_IMG = os.path.join(DATA_DIR, "logo.png")
PHYLO_TREE_IMG = os.path.join(RESULTS_DIR, "phylogenetic_tree_dod.png") # Image obtained from MEGA12 analysis
GENOMIC_ANALYSIS_IMG = os.path.join(RESULTS_DIR, "genomic_analysis.png")
SEQ_LENGTH_IMG = os.path.join(RESULTS_DIR, "sequence_length_comparison.png") 


# --- Sidebar ---
if os.path.exists(LOGO_IMG):
    st.sidebar.image(LOGO_IMG, width=250)
else:
    # If the image is not found, display an emoji or text
    st.sidebar.title("🧬 DOD Explorer") 

st.sidebar.title("Analysis Control Panel")
app_mode = st.sidebar.selectbox("Select Analysis Section", 
    ["Dashboard Overview", "Sequence Analysis", "BLAST Homology", "Phylogenetics", "Domain Architecture"])

st.sidebar.info("Project: Mirabilis jalapa DOD Gene Analysis (Automated Pipeline)")
st.sidebar.write("Lead Researcher & Developer: **Zahidul Hasan**")
st.sidebar.write("Dept. of Botany, University of Chittagong, Bangladesh")


# --- Data loading function ---
@st.cache_data
def load_gbk_data():
    if os.path.exists(GBK_FILE):
        return SeqIO.read(GBK_FILE, "genbank")
    return None

@st.cache_data
def load_blast_data():
    if os.path.exists(BLAST_CSV):
        return pd.read_csv(BLAST_CSV)
    return None

# --- 1. Dashboard Overview ---
if app_mode == "Dashboard Overview":
    st.title("🧬 Betalain Pathway Explorer")
    st.subheader("In silico Identification and Sequence Analysis of DOD Gene in Mirabilis jalapa")
    
    col1, col2, col3 = st.columns(3)
    
    record = load_gbk_data()
    if record:
        dna_len = len(record.seq)
        gc_val = gc_fraction(record.seq) * 100
        # Logic Update: Now displays length according to biological translation
        protein_seq = record.seq.translate(to_stop=True)
        protein_len = len(protein_seq)
        
        with col1:
            st.metric("Gene Length", f"{dna_len} bp", delta="Target Gene")
        with col2:
            st.metric("GC Content", f"{gc_val:.2f}%", delta="Normal Range")
        with col3:
            st.metric("Protein Length", f"{protein_len} aa", delta="Translated (to stop)")
    else:
        st.warning("GenBank file not found in 'data/mirabilis_dod.gbk'.")


    st.markdown("---")
    st.markdown("""
    ### The Significance and Objective of this Project: 
    This dashboard was developed to analyze the molecular characteristics of the **DOD (DOPA 4,5-dioxygenase)** gene, a key enzyme involved in the **Betalain Biosynthesis Pathway** responsible for flower pigmentation in the plant *Mirabilis jalapa*.
    
    **Relationship with the Betalain Pathway:** The DOD enzyme catalyzes the conversion of DOPA into Betalamic Acid, a crucial precursor required for the biosynthesis of red and yellow betalain pigments (betacyanins = red-violet and betaxanthins = yellow-orange). This in silico analysis provides insights into the gene’s structural features and evolutionary relationships.
    """)
     
    
    if os.path.exists(FLOWER_IMG):
        st.image(FLOWER_IMG, caption="Mirabilis jalapa - Flower Pigmentation Study", width=500)
    else:
        st.warning("Flower image not found in the data directory.")

# --- 2. Sequence analysis ---
elif app_mode == "Sequence Analysis":
    st.title("🔍 Nucleotide and Protein Analysis")
    
    
    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Genomic Summary", 
        "📏 Length Comparison", 
        "🧬 DNA Composition", 
        "🧪 Amino Acid Frequency"
    ])
    
    # Tab 1: Genomic Summary (Static Image)
    with tab1:
        if os.path.exists(GENOMIC_ANALYSIS_IMG):
            st.subheader("Genomic Composition and BLAST Identity Overview")
            st.image(GENOMIC_ANALYSIS_IMG, use_container_width=True)
            st.caption("Visual Summary of DNA Composition and Homology Search")
        else:
            st.error("Genomic Analysis image not found.")
            
    # Tab 2: Sequence Length Comparison Image
    with tab2:
        if os.path.exists(SEQ_LENGTH_IMG):
            st.subheader("Protein Sequence Length Comparison")
            st.image(SEQ_LENGTH_IMG, use_container_width=True)
            st.info("Observation: Here we can see the protein length (267 aa) of Mirabilis jalapa is similar to other homologous sequences.")
        else:
            st.error("Sequence Length Comparison image not found.")

    # Data loading (for subsequent tabs)
    record = load_gbk_data()
    if record:
        # Tab 3: Interactive DNA Pie Chart
        with tab3:
            st.write("#### Interactive DNA Base Distribution")
            bases = dict(Counter(record.seq))
            fig_dna = px.pie(values=list(bases.values()), names=list(bases.keys()), 
                           title="Base Composition (A, T, G, C)", 
                           color_discrete_sequence=px.colors.qualitative.Pastel)
            st.plotly_chart(fig_dna)
            
        # Tab 4: Interactive Protein Bar Chart
        with tab4:
            st.write("#### Protein Sequence Statistics")
            protein_seq = record.seq.translate(to_stop=True)
            aa_counts = dict(Counter(protein_seq))
            aa_df = pd.DataFrame(list(aa_counts.items()), columns=['Amino Acid', 'Frequency']).sort_values(by='Frequency', ascending=False)
            
            fig_aa = px.bar(aa_df, x='Amino Acid', y='Frequency', color='Frequency', 
                          title="Top Amino Acids in DOD Protein")
            st.plotly_chart(fig_aa)
    else:
        st.error("GenBank file not found in 'data/mirabilis_dod.gbk'.")


# --- 3. BLAST Homology Analysis ---
elif app_mode == "BLAST Homology":
    st.title("💥 BLASTp Search Results")
    df = load_blast_data()
    if df is not None:
        st.dataframe(df.style.highlight_max(axis=0, subset=['Identity (%)'], color='lightgreen'))
        fig_blast = px.scatter(df, x="Identity (%)", y="E-value", size="Identity (%)", color="Organism Name",
                             hover_name="Organism Name", log_y=True, title="Blast Hits: Identity vs E-value")
        st.plotly_chart(fig_blast, use_container_width=True)
    else:
        st.warning("BLAST data file not found in 'results/blast_results.csv'.")

# --- 4. Phylogenetics ---
elif app_mode == "Phylogenetics":
    st.title("🌳 Evolutionary Relationships")
    st.markdown("#### Phylogenetic Tree Analysis of DOD Gene")

    # First image: General/Pipeline Tree
    tree_img_1 = os.path.join(RESULTS_DIR, "phylogenetic_tree.png")
    if os.path.exists(tree_img_1):
        st.subheader("Pipeline Generated Tree")
        st.image(tree_img_1, caption="Phylogenetic Tree (Auto-generated)", use_container_width=True)
    
    st.markdown("---") 

    # Second image: MEGA specific Tree
    # Ensure that PHYLO_TREE_IMG = "phylogenetic_tree_dod.png" is included in the path setup.
    if os.path.exists(PHYLO_TREE_IMG):
        st.subheader("Custom MEGA Analysis")
        st.image(PHYLO_TREE_IMG, caption="Phylogenetic Tree constructed using MEGA 12", use_container_width=True)
        
        # Detailed analysis description
        st.info("""
        **Methodology Details (MEGA 12):**
        - **Method:** Neighbor-Joining (NJ)
        - **Bootstrap:** 1000 replicates
        - **Model:** Poisson correction
        """)
    
    # If no image is found
    if not os.path.exists(tree_img_1) and not os.path.exists(PHYLO_TREE_IMG):
        st.error("Sorry, no phylogenetic tree image found.")
        


# --- 5. Protein Structure Analysis ---
elif app_mode == "Domain Architecture":
    st.title("🏗️ Protein Structure & Domain Architecture")
    
    # Set up two tabs — 3D view and domain image display
    tab_3d, tab_domain = st.tabs(["🧬 3D Interactive Model", "🖼️ Conserved Domains"])
    
    with tab_3d:
        st.subheader("Interactive 3D Structure (SWISS-MODEL Prediction)")
        
        # PDB file path (ensure model_01.pdb file is in the results folder)
        pdb_file = os.path.join(RESULTS_DIR, "model_01.pdb") 
        
        if os.path.exists(pdb_file):
            with open(pdb_file, "r") as f:
                pdb_data = f.read()
            
            # 3D Viewer Setup
            view = py3Dmol.view(width=None, height=400)
            view.addModel(pdb_data, 'pdb')
            view.setStyle({'cartoon': {'color': 'spectrum'}}) # Colorful cartoon style
            view.zoomTo()
            view.spin(True) # The model will rotate slowly.
            
            showmol(view, height=400, width='100%')
            
            st.info("The model can be interactively zoomed and rotated using the mouse.")
            
            # Download button
            st.download_button("Download PDB File", pdb_data, "DOD_model.pdb")
        else:
            st.error("PDB file not found in 'results/model_01.pdb'.")
            st.info("Please download the PDB file from SWISS-MODEL and place it in the results folder.")
    with tab_domain:
        domain_img = os.path.join(RESULTS_DIR, "conserved_domain.png")
        if os.path.exists(domain_img):
            st.image(domain_img, caption="Conserved Functional Domains", use_container_width=True)
        else:
            st.error("Domain image not found.")

st.markdown("---")
st.caption("© 2026 Zahidul Hasan | Department of Botany, University of Chittagong")