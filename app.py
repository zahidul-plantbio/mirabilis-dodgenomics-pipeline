"""
Project: Betalain Pathway Explorer
Description: Updated Image Path Handling and UI Fix
Author: Zahidul Hasan
"""

import streamlit as st
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import plotly.express as px
from collections import Counter

# --- পেজ কনফিগারেশন ---
st.set_page_config(
    page_title="Betalain Pathway Explorer",
    page_icon="🧬",
    layout="wide"
)

# --- কাস্টম স্টাইল ---
st.markdown("""
    <style>
    .stMetric {
        background-color: #ffffff;
        padding: 15px;
        border-radius: 10px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    </style>
    """, unsafe_allow_html=True)

# --- পাথ সেটআপ (Relative Paths are key for Streamlit Cloud) ---
BASE_DIR = os.path.dirname(__file__) # ফাইলের বর্তমান ডিরেক্টরি
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# ফাইল পাথ
GBK_FILE = os.path.join(DATA_DIR, "mirabilis_dod.gbk")
BLAST_CSV = os.path.join(RESULTS_DIR, "blast_results.csv")
TREE_IMG = os.path.join(RESULTS_DIR, "phylogenetic_tree.png")
DOMAIN_IMG = os.path.join(RESULTS_DIR, "conserved_domain.png")

# --- সাইডবার ---
LOGO_URL = "https://cdn-icons-png.flaticon.com/512/3015/3015511.png"
st.sidebar.image(LOGO_URL, width=80)
st.sidebar.title("অ্যানালাইসিস কন্ট্রোল")
app_mode = st.sidebar.selectbox("সেকশন নির্বাচন করুন", 
    ["ড্যাশবোর্ড ওভারভিউ", "সিকোয়েন্স এনালাইসিস", "BLAST হোমোলজি", "ফাইলোজেনেটিক্স", "ডোমেইন আর্কিটেকচার"])

st.sidebar.markdown(f"""
**প্রজেক্ট:** Mirabilis jalapa DOD Gene Analysis  
**নির্মাণে:** Zahidul Hasan  
""")

# --- ডাটা লোডিং ---
@st.cache_data
def get_seq_record():
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
    record = get_seq_record()
    
    if record:
        dna_len = len(record.seq)
        gc_val = gc_fraction(record.seq) * 100
        protein_len = len(record.seq.translate(to_stop=True))
        
        c1, c2, c3 = st.columns(3)
        c1.metric("জিনের দৈর্ঘ্য", f"{dna_len} bp", "Target")
        c2.metric("GC কন্টেন্ট", f"{gc_val:.2f}%", "Normal")
        c3.metric("প্রোটিনের দৈর্ঘ্য", f"{protein_len} aa", "Translated")

    st.info("Betalain Pathway-র বিশ্লেষণে DOD জিনের ভূমিকা অপরিসীম।")
    
    st.image("https://upload.wikimedia.org/wikipedia/commons/0/0d/Mirabilis_jalapa_Flower.jpg", 
             caption="Mirabilis jalapa - Flower", width=400)

# --- ২. সিকোয়েন্স এনালাইসিস ---
elif app_mode == "সিকোয়েন্স এনালাইসিস":
    st.title("🔍 নিউক্লিওটাইড এবং প্রোটিন এনালাইসিস")
    record = get_seq_record()
    
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
        st.write("NCBI SwissProt ডাটাবেস থেকে প্রাপ্ত টপ হিটসমূহ:")
        st.dataframe(df.style.highlight_max(axis=0, subset=['Identity (%)'], color='lightgreen'))
        
        fig_blast = px.scatter(df, x="Identity (%)", y="E-value", size="Identity (%)", color="Organism Name",
                             hover_name="Organism Name", log_y=True, title="Blast Hits: Identity vs E-value")
        st.plotly_chart(fig_blast, use_container_width=True)
    else:
        st.warning("BLAST ডাটা ফাইলটি 'results/blast_results.csv' পাথে পাওয়া যায়নি।")

# --- ৪. ফাইলোজেনেটিক্স (ইমেজ সেকশন) ---
elif app_mode == "ফাইলোজেনেটিক্স":
    st.title("🌳 Evolutionary Relationships")
    if os.path.exists(TREE_IMG):
        st.image(TREE_IMG, caption="Phylogenetic Tree (Neighbor-Joining Method)", use_container_width=True)
    else:
        st.error(f"ইমেজ ফাইলটি পাওয়া যায়নি! নিশ্চিত করুন যে এটি এই পাথে আছে: {TREE_IMG}")
        st.info("পরামর্শ: গিটহাবে 'results/phylogenetic_tree.png' ফাইলটি আপলোড করেছেন কি না চেক করুন।")

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