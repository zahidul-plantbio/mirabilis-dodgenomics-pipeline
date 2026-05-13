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

# --- পেজ কনফিগারেশন ---
st.set_page_config(
    page_title="Betalain Pathway Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- কাস্টম স্টাইল ---
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

# --- পাথ সেটআপ ---
BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# ফাইল পাথ নিশ্চিত করা
GBK_FILE = os.path.join(DATA_DIR, "mirabilis_dod.gbk")
BLAST_CSV = os.path.join(RESULTS_DIR, "blast_results.csv")
FLOWER_IMG = os.path.join(DATA_DIR, "mirabilis_jalapa_flower.png") # এক্সটেনশন png নিশ্চিত করুন
LOGO_IMG = os.path.join(DATA_DIR, "logo.png")
PHYLO_TREE_IMG = os.path.join(RESULTS_DIR, "phylogenetic_tree_dod.png") # MEGA থেকে প্রাপ্ত ইমেজ
GENOMIC_ANALYSIS_IMG = os.path.join(RESULTS_DIR, "genomic_analysis.png")
SEQ_LENGTH_IMG = os.path.join(RESULTS_DIR, "sequence_length_comparison.png") # নাম ঠিক করা হয়েছে


# --- সাইডবার ---
if os.path.exists(LOGO_IMG):
    st.sidebar.image(LOGO_IMG, width=250)
else:
    # যদি ইমেজ না পাওয়া যায় তবে একটি ইমোজি বা টেক্সট দেখাবে
    st.sidebar.title("🧬 DOD Explorer") 

st.sidebar.title("Analysis Control Panel")
app_mode = st.sidebar.selectbox("Select Analysis Section", 
    ["Dashboard Overview", "Sequence Analysis", "BLAST Homology", "Phylogenetics", "Domain Architecture"])
  # ["ড্যাশবোর্ড ওভারভিউ", "সিকোয়েন্স এনালাইসিস", "BLAST হোমোলজি", "ফাইলোজেনেটিক্স", "ডোমেইন আর্কিটেকচার"])
st.sidebar.info("Project: Mirabilis jalapa DOD Gene Analysis (Automated Pipeline)")
st.sidebar.write("Built by: **Zahidul Hasan**")
st.sidebar.write("Department of Botany, University of Chittagong, Bangladesh")


# --- ডাটা লোড করার ফাংশন ---
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

# --- ১. ড্যাশবোর্ড ওভারভিউ ---
if app_mode == "Dashboard Overview":
    st.title("🧬 Betalain Pathway Explorer")
    st.subheader("In silico Identification and Sequence Analysis of DOD Gene in Mirabilis jalapa")
    
    col1, col2, col3 = st.columns(3)
    
    record = load_gbk_data()
    if record:
        dna_len = len(record.seq)
        gc_val = gc_fraction(record.seq) * 100
        # লজিক আপডেট: এখন বায়োলজিক্যাল ট্রান্সলেশন অনুযায়ী দৈর্ঘ্য দেখাবে
        protein_seq = record.seq.translate(to_stop=True)
        protein_len = len(protein_seq)
        
        with col1:
            st.metric("Gene Length", f"{dna_len} bp", delta="Target Gene")
        with col2:
            st.metric("GC Content", f"{gc_val:.2f}%", delta="Normal Range")
        with col3:
            st.metric("Protein Length", f"{protein_len} aa", delta="Translated (to stop)")
    else:
        st.warning("GenBank ফাইলটি 'data/mirabilis_dod.gbk' পাথে পাওয়া যায়নি।")

# সঠিক ইন্ডেন্টেশন (লাইন ১১০ থেকে শুরু)
#এই ড্যাশবোর্ডটি *Mirabilis jalapa* গাছের ফুলের রঙের জন্য দায়ী **Betalain Biosynthesis Pathway**-এর একটি মূল এনজাইম **DOD (DOPA 4,5-dioxygenase)** জিনের আণবিক বৈশিষ্ট্য বিশ্লেষণের জন্য তৈরি। 
#**Betalain Pathway-র সাথে সম্পর্ক:** DOD এনজাইমটি DOPA-কে Betalamic Acid-এ রূপান্তর করে, যা লাল ও হলুদ পিগমেন্ট তৈরির জন্য অপরিহার্য। এই ইন-সিলিকো বিশ্লেষণ জিনের গঠন ও বিবর্তনীয় সম্পর্ক বুঝতে সাহায্য করে।
    st.markdown("---")
    st.markdown("""
    ### The Significance and Objective of this Project: 
    This dashboard was developed to analyze the molecular characteristics of the **DOD (DOPA 4,5-dioxygenase)** gene, a key enzyme involved in the **Betalain Biosynthesis Pathway** responsible for flower pigmentation in the plant *Mirabilis jalapa*.
    
    **Relationship with the Betalain Pathway:** The DOD enzyme catalyzes the conversion of DOPA into Betalamic Acid, a crucial precursor required for the biosynthesis of red and yellow betalain pigments. This in silico analysis provides insights into the gene’s structural features and evolutionary relationships.
    """)
     
    # লক্ষ্য করুন: এই অংশটি এখন 'if app_mode == "ড্যাশবোর্ড ওভারভিউ":' এর ভেতরে আছে
    if os.path.exists(FLOWER_IMG):
        st.image(FLOWER_IMG, caption="Mirabilis jalapa - Flower Pigmentation Study", width=500)
    else:
        st.warning("ফুলের ইমেজটি data ফোল্ডারে পাওয়া যায়নি।")

# --- ২. সিকোয়েন্স এনালাইসিস ---
elif app_mode == "Sequence Analysis":
    st.title("🔍 Nucleotide and Protein Analysis")
    
    # সব এনালাইসিসকে একটি মাত্র ট্যাব সেটের অধীনে আনা (আরও প্রফেশনাল)
    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Genomic Summary", 
        "📏 Length Comparison", 
        "🧬 DNA Composition", 
        "🧪 Amino Acid Frequency"
    ])
    
    # ট্যাব ১: স্ট্যাটিক জেনোমিক সামারি ইমেজ
    with tab1:
        if os.path.exists(GENOMIC_ANALYSIS_IMG):
            st.subheader("Genomic Composition and BLAST Identity Overview")
            st.image(GENOMIC_ANALYSIS_IMG, use_container_width=True)
            st.caption("Visual Summary of DNA Composition and Homology Search")
        else:
            st.error("Genomic Analysis ইমেজটি পাওয়া যায়নি।")
            
    # ট্যাব ২: সিকোয়েন্স লেন্থ কম্পারিজন ইমেজ
    with tab2:
        if os.path.exists(SEQ_LENGTH_IMG):
            st.subheader("Protein Sequence Length Comparison")
            st.image(SEQ_LENGTH_IMG, use_container_width=True)
            st.info("পর্যবেক্ষণ: এখানে দেখা যাচ্ছে Mirabilis jalapa এর প্রোটিন দৈর্ঘ্য (২৬৭ aa) অন্যান্য হোমোলোগাস সিকোয়েন্সের কাছাকাছি।")
        else:
            st.error("Sequence Length Comparison ইমেজটি পাওয়া যায়নি।")

    # ডাটা লোড করা (পরের ট্যাবগুলোর জন্য)
    record = load_gbk_data()
    if record:
        # ট্যাব ৩: ইন্টারঅ্যাক্টিভ DNA পাই চার্ট
        with tab3:
            st.write("#### Interactive DNA Base Distribution")
            bases = dict(Counter(record.seq))
            fig_dna = px.pie(values=list(bases.values()), names=list(bases.keys()), 
                           title="Base Composition (A, T, G, C)", 
                           color_discrete_sequence=px.colors.qualitative.Pastel)
            st.plotly_chart(fig_dna)
            
        # ট্যাব ৪: ইন্টারঅ্যাক্টিভ প্রোটিন বার চার্ট
        with tab4:
            st.write("#### প্রোটিন সিকোয়েন্স স্ট্যাটিস্টিকস")
            protein_seq = record.seq.translate(to_stop=True)
            aa_counts = dict(Counter(protein_seq))
            aa_df = pd.DataFrame(list(aa_counts.items()), columns=['Amino Acid', 'Frequency']).sort_values(by='Frequency', ascending=False)
            
            fig_aa = px.bar(aa_df, x='Amino Acid', y='Frequency', color='Frequency', 
                          title="Top Amino Acids in DOD Protein")
            st.plotly_chart(fig_aa)
    else:
        st.error("GenBank ফাইল লোড করা সম্ভব হয়নি।")

# বাকি সেকশনগুলো অপরিবর্তিত রাখা হলো (BLAST, ফাইলোজেনেটিক্স, ইত্যাদি)
# --- ৩. BLAST হোমোলজি ---
elif app_mode == "BLAST Homology":
    st.title("💥 BLASTp Search Results")
    df = load_blast_data()
    if df is not None:
        st.dataframe(df.style.highlight_max(axis=0, subset=['Identity (%)'], color='lightgreen'))
        fig_blast = px.scatter(df, x="Identity (%)", y="E-value", size="Identity (%)", color="Organism Name",
                             hover_name="Organism Name", log_y=True, title="Blast Hits: Identity vs E-value")
        st.plotly_chart(fig_blast, use_container_width=True)
    else:
        st.warning("BLAST ডাটা ফাইলটি 'results/blast_results.csv' পাথে পাওয়া যায়নি।")

# --- ৪. ফাইলোজেনেটিক্স ---
elif app_mode == "Phylogenetics":
    st.title("🌳 Evolutionary Relationships")
    st.markdown("#### Phylogenetic Tree Analysis of DOD Gene")

    # প্রথম ইমেজ: General/Pipeline Tree
    tree_img_1 = os.path.join(RESULTS_DIR, "phylogenetic_tree.png")
    if os.path.exists(tree_img_1):
        st.subheader("Pipeline Generated Tree")
        st.image(tree_img_1, caption="Phylogenetic Tree (Auto-generated)", use_container_width=True)
    
    st.markdown("---") # একটি ডিভাইডার লাইন

    # দ্বিতীয় ইমেজ: MEGA specific Tree
    # নিশ্চিত করুন পাথ সেটআপে PHYLO_TREE_IMG = "phylogenetic_tree_dod.png" আছে
    if os.path.exists(PHYLO_TREE_IMG):
        st.subheader("Custom MEGA Analysis")
        st.image(PHYLO_TREE_IMG, caption="Phylogenetic Tree constructed using MEGA 12", use_container_width=True)
        
        # এনালাইসিসের বিস্তারিত বর্ণনা
        st.info("""
        **Methodology Details (MEGA 12):**
        - **Method:** Neighbor-Joining (NJ)
        - **Bootstrap:** 1000 replicates
        - **Model:** Poisson correction
        """)
    
    # যদি কোনো ইমেজই না পাওয়া যায়
    if not os.path.exists(tree_img_1) and not os.path.exists(PHYLO_TREE_IMG):
        st.error("দুঃখিত, কোনো ফাইলোজেনেটিক ট্রি ইমেজ পাওয়া যায়নি।")
        

# --- ৫. ডোমেইন আর্কিটেকচার ---
# elif app_mode == "ডোমেইন আর্কিটেকচার":
#     st.title("🏗️ Protein Domain Architecture")
#     domain_img = os.path.join(RESULTS_DIR, "conserved_domain.png")
#     if os.path.exists(domain_img):
#         st.image(domain_img, caption="Conserved Functional Domains", use_container_width=True)
#     else:
#         st.error("ডোমেইন ইমেজটি পাওয়া যায়নি।")

# --- ৫. প্রোটিন স্ট্রাকচার এনালাইসিস ---
elif app_mode == "Domain Architecture":
    st.title("🏗️ Protein Structure & Domain Architecture")
    
    # দুটি ট্যাব তৈরি করা (একটি ৩ডি ভিউয়ের জন্য, একটি ডোমেইন ইমেজের জন্য)
    tab_3d, tab_domain = st.tabs(["🧬 3D Interactive Model", "🖼️ Conserved Domains"])
    
    with tab_3d:
        st.subheader("Interactive 3D Structure (SWISS-MODEL Prediction)")
        
        # PDB ফাইল পাথ (নিশ্চিত করুন model_01.pdb ফাইলটি results ফোল্ডারে আছে)
        pdb_file = os.path.join(RESULTS_DIR, "model_01.pdb") 
        
        if os.path.exists(pdb_file):
            with open(pdb_file, "r") as f:
                pdb_data = f.read()
            
            # ৩ডি ভিউ সেটআপ
            view = py3Dmol.view(width=None, height=400)
            view.addModel(pdb_data, 'pdb')
            view.setStyle({'cartoon': {'color': 'spectrum'}}) # রঙিন কার্টুন স্টাইল
            view.zoomTo()
            view.spin(True) # মডেলটি ধীরে ধীরে ঘুরবে
            
            showmol(view, height=400, width='100%')
            
            st.info("The model can be interactively zoomed and rotated using the mouse.")
            
            # ডাউনলোড বাটন
            st.download_button("Download PDB File", pdb_data, "DOD_model.pdb")
        else:
            st.error("PDB ফাইলটি 'results/model_01.pdb' পাথে পাওয়া যায়নি।")
            st.info("আপনার SWISS-MODEL থেকে PDB ফাইলটি ডাউনলোড করে results ফোল্ডারে রাখুন।")

    with tab_domain:
        domain_img = os.path.join(RESULTS_DIR, "conserved_domain.png")
        if os.path.exists(domain_img):
            st.image(domain_img, caption="Conserved Functional Domains", use_container_width=True)
        else:
            st.error("ডোমেইন ইমেজটি পাওয়া যায়নি।")

st.markdown("---")
st.caption("© 2026 Zahidul Hasan | Department of Botany, University of Chittagong")