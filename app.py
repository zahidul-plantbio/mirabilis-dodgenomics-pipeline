"""
Project: Betalain Pathway Explorer
Description: Interactive Streamlit App for Bioinformatics Sequence Analysis of DOD Gene in Mirabilis jalapa
Author: Zahidul Hasan
"""

import streamlit as st
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter

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
BASE_DIR = os.path.dirname(__file__) 
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# ফাইল পাথ নিশ্চিত করা
GBK_FILE = os.path.join(DATA_DIR, "mirabilis_dod.gbk")
BLAST_CSV = os.path.join(RESULTS_DIR, "blast_results.csv")
TREE_IMG = os.path.join(RESULTS_DIR, "phylogenetic_tree.png")
DOMAIN_IMG = os.path.join(RESULTS_DIR, "conserved_domain.png")
FLOWER_IMG = os.path.join(DATA_DIR, "mirabilis_jalapa_flower.png")

# --- সাইডবার ---
LOGO_URL = "https://cdn-icons-png.flaticon.com/128/3015/3015511.png"
st.sidebar.image(LOGO_URL, width=80)
st.sidebar.title("অ্যানালাইসিস কন্ট্রোল")
app_mode = st.sidebar.selectbox("সেকশন নির্বাচন করুন", 
    ["ড্যাশবোর্ড ওভারভিউ", "সিকোয়েন্স এনালাইসিস", "BLAST হোমোলজি", "ফাইলোজেনেটিক্স", "ডোমেইন আর্কিটেকচার"])

st.sidebar.info("প্রজেক্ট: Betalain Pathway Explorer")
st.sidebar.write("বিভাগ: উদ্ভিদবিজ্ঞান বিভাগ, চট্টগ্রাম বিশ্ববিদ্যালয়")
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
    st.title("🧬 Betalain Pathway Explorer")
    st.subheader("In silico Identification and Sequence Analysis of DOD Gene in Mirabilis jalapa")
    
    record = load_gbk_data()
    if record:
        dna_len = len(record.seq)
        gc_val = gc_fraction(record.seq) * 100
        # বায়োলজিক্যাল ট্রান্সলেশন অনুযায়ী দৈর্ঘ্য (২৬৭ aa)
        protein_seq = record.seq.translate(to_stop=True)
        protein_len = len(protein_seq)
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("জিনের দৈর্ঘ্য", f"{dna_len} bp", delta="Target Gene")
        with col2:
            st.metric("GC কন্টেন্ট", f"{gc_val:.2f}%", delta="Normal Range")
        with col3:
            st.metric("প্রোটিনের দৈর্ঘ্য", f"{protein_len} aa", delta="Translated")

    st.markdown("---")
    
    # লোকাল ইমেজ প্রদর্শনের লজিক (ইন্ডেন্টেশন ঠিক করা হয়েছে)
    if os.path.exists(FLOWER_IMG):
        st.image(FLOWER_IMG, caption="Mirabilis jalapa - Flower Pigmentation Study", width=500)
    else:
        st.warning("ফুলের ইমেজটি data ফোল্ডারে পাওয়া যায়নি।")

    st.markdown("""
    ### প্রজেক্টের গুরুত্ব ও উদ্দেশ্য:
    এই ড্যাশবোর্ডটি *Mirabilis jalapa* গাছের ফুলের রঙের জন্য দায়ী **Betalain Biosynthesis Pathway**-এর একটি মূল এনজাইম **DOD (DOPA 4,5-dioxygenase)** জিনের আণবিক বৈশিষ্ট্য বিশ্লেষণের জন্য তৈরি। 
    """)

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
    else:
        st.error("সিকোয়েন্স ডাটা লোড করা সম্ভব হয়নি।")

# --- ৩. BLAST হোমোলজি ---
elif app_mode == "BLAST হোমোলজি":
    st.title("💥 BLASTp Search Results")
    df = load_blast_data()
    if df is not None:
        st.dataframe(df.style.highlight_max(axis=0, subset=['Identity (%)'], color='lightgreen'))
        fig_blast = px.scatter(df, x="Identity (%)", y="E-value", size="Identity (%)", color="Organism Name",
                             hover_name="Organism Name", log_y=True, title="Blast Hits: Identity vs E-value")
        st.plotly_chart(fig_blast, use_container_width=True)
    else:
        st.warning("BLAST ডাটা ফাইলটি 'results/blast_results.csv' পাওয়া যায়নি।")

# --- ৪. ফাইলোজেনেটিক্স ---
elif app_mode == "ফাইলোজেনেটিক্স":
    st.title("🌳 Evolutionary Relationships")
    if os.path.exists(TREE_IMG):
        st.image(TREE_IMG, caption="Phylogenetic Tree (Neighbor-Joining Method)", use_container_width=True)
    else:
        st.error(f"ইমেজ ফাইলটি পাওয়া যায়নি। পাথ: {TREE_IMG}")

# --- ৫. ডোমেইন আর্কিটেকচার ---
elif app_mode == "ডোমেইন আর্কিটেকচার":
    st.title("🏗️ Protein Domain Architecture")
    if os.path.exists(DOMAIN_IMG):
        st.image(DOMAIN_IMG, caption="Conserved Functional Domains (LigB Superfamily)", use_container_width=True)
    else:
        st.warning(f"ডোমেইন ইমেজটি পাওয়া যায়নি: {DOMAIN_IMG}")

# --- ফুটার ---
st.markdown("---")
st.caption("© 2026 Zahidul Hasan | Department of Botany, University of Chittagong")