"""
Project: Mirabilis jalapa DOD Gene Analysis Dashboard
Description: Interactive Streamlit App for Bioinformatics Sequence Analysis
Author: Zahidul Hasan
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

# --- পেজ কনফিগারেশন ---
st.set_page_config(
    page_title="Mirabilis DOD Analyzer",
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

# --- সাইডবার ---
st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/thumb/e/e5/Phylo_tree.svg/1200px-Phylo_tree.svg.png", width=100)
st.sidebar.title("অ্যানালাইসিস কন্ট্রোল")
app_mode = st.sidebar.selectbox("সেকশন নির্বাচন করুন", 
    ["ড্যাশবোর্ড ওভারভিউ", "সিকোয়েন্স এনালাইসিস", "BLAST হোমোলজি", "ফাইলোজেনেটিক্স", "ডোমেইন আর্কিটেকচার"])

st.sidebar.info("প্রজেক্ট: Mirabilis jalapa DOD Gene Analysis")
st.sidebar.write("নির্মাণে: **Zahidul Hasan**")

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
if app_mode == "ড্যাশবোর্ড ওভারভিউ":
    st.title("🧬 Mirabilis jalapa DOD Gene Analysis Dashboard")
    st.subheader("In silico Identification and Sequence Analysis")
    
    col1, col2, col3 = st.columns(3)
    
    record = load_gbk_data()
    if record:
        dna_len = len(record.seq)
        gc_val = gc_fraction(record.seq) * 100
        
        with col1:
            st.metric("জিনের দৈর্ঘ্য", f"{dna_len} bp", delta="Target Gene")
        with col2:
            st.metric("GC কন্টেন্ট", f"{gc_val:.2f}%", delta="Normal Range")
        with col3:
            st.metric("প্রোটিনের দৈর্ঘ্য", f"{int(dna_len/3)} aa", delta="Translated")

    st.markdown("---")
    st.markdown("""
    ### প্রজেক্টের উদ্দেশ্য:
    এই ড্যাশবোর্ডটি *Mirabilis jalapa* গাছের ফুলের রঙের জন্য দায়ী **DOD (DOPA 4,5-dioxygenase)** জিনের সিকোয়েন্স এনালাইসিস এবং বিবর্তনীয় সম্পর্ক প্রদর্শনের জন্য তৈরি করা হয়েছে।
    
    **ব্যবহৃত প্রযুক্তি:** Python, Biopython, Streamlit, Plotly, Seaborn.
    """)
    
    st.image("https://upload.wikimedia.org/wikipedia/commons/0/0d/Mirabilis_jalapa_Flower.jpg", caption="Mirabilis jalapa - Four O'Clock Flower", width=500)

# --- ২. সিকোয়েন্স এনালাইসিস ---
elif app_mode == "সিকোয়েন্স এনালাইসিস":
    st.title("🔍 নিউক্লিওটাইড এবং প্রোটিন এনালাইসিস")
    record = load_gbk_data()
    
    if record:
        tab1, tab2 = st.tabs(["DNA কম্পোজিশন", "অ্যামিনো অ্যাসিড ফ্রিকোয়েন্সি"])
        
        with tab1:
            st.write("#### DNA বেস ডিস্ট্রিবিউশন")
            bases = dict(Counter(record.seq))
            fig_dna = px.pie(values=list(bases.values()), names=list(bases.keys()), 
                           title="Base Composition (A, T, G, C)", color_discrete_sequence=px.colors.qualitative.Pastel)
            st.plotly_chart(fig_dna)
            
        with tab2:
            st.write("#### প্রোটিন সিকোয়েন্স স্ট্যাটিস্টিকস")
            protein_seq = record.seq.translate(to_stop=True)
            aa_counts = dict(Counter(protein_seq))
            aa_df = pd.DataFrame(list(aa_counts.items()), columns=['Amino Acid', 'Frequency']).sort_values(by='Frequency', ascending=False)
            
            fig_aa = px.bar(aa_df, x='Amino Acid', y='Frequency', color='Frequency', title="Top Amino Acids in DOD Protein")
            st.plotly_chart(fig_aa)

# --- ৩. BLAST হোমোলজি ---
elif app_mode == "BLAST হোমোলজি":
    st.title("💥 BLASTp Search Results")
    df = load_blast_data()
    
    if df is not None:
        st.write("NCBI SwissProt ডাটাবেস থেকে প্রাপ্ত টপ হিটসমূহ:")
        
        # ইন্টারেক্টিভ টেবিল
        st.dataframe(df.style.highlight_max(axis=0, subset=['Identity (%)'], color='lightgreen'))
        
        # ইন্টারেক্টিভ চার্ট
        fig_blast = px.scatter(df, x="Identity (%)", y="E-value", size="Identity (%)", color="Organism Name",
                             hover_name="Organism Name", log_y=True, title="Blast Hits: Identity vs E-value")
        st.plotly_chart(fig_blast, use_container_width=True)
    else:
        st.warning("BLAST ডাটা ফাইল খুঁজে পাওয়া যায়নি। দয়া করে আগে পাইপলাইন রান করুন।")

# --- ৪. ফাইলোজেনেটিক্স ---
elif app_mode == "ফাইলোজেনেটিক্স":
    st.title("🌳 Evolutionary Relationships")
    tree_img = os.path.join(RESULTS_DIR, "phylogenetic_tree.png")
    
    if os.path.exists(tree_img):
        st.image(tree_img, caption="Phylogenetic Tree (Neighbor-Joining Method)", use_container_width=True)
        st.markdown("""
        **বিশ্লেষণ:**
        - এই গাছটি দেখায় কিভাবে *Mirabilis jalapa* এর DOD প্রোটিন অন্যান্য প্রজাতির (যেমন *Beta vulgaris*, *Spinacia oleracea*) সাথে সম্পর্কিত।
        - জেনেটিক ডিস্ট্যান্স যত কম, সম্পর্ক তত কাছাকাছি।
        """)
    else:
        st.error("ট্রি ইমেজটি জেনারেট করা হয়নি।")

# --- ৫. ডোমেইন আর্কিটেকচার ---
elif app_mode == "ডোমেইন আর্কিটেকচার":
    st.title("🏗️ Protein Domain Architecture")
    domain_img = os.path.join(RESULTS_DIR, "conserved_domain.png")
    
    if os.path.exists(domain_img):
        st.image(domain_img, caption="Conserved Functional Domains", use_container_width=True)
        st.info("LigB Domain এবং বিশেষ কনসার্ভড মোটিফগুলো চিহ্নিত করা হয়েছে যা Betalain বায়োসিন্থেসিসে ভূমিকা রাখে।")
    else:
        st.error("ডোমেইন ইমেজটি পাওয়া যায়নি।")

# --- ফুটার ---
st.markdown("---")
st.caption("© 2026 Zahidul Hasan | Department of Botany, University of Chittagong | Powered by Streamlit & Biopython")