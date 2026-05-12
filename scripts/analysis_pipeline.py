"""
Project: Mirabilis jalapa DOD Gene Analysis
Author: Zahidul Hasan
Department: Botany, University of Chittagong, Bangladesh
Email: zahidulhasan.botany.cu@gmail.com
LinkedIn: linkedin.com/in/zahidulhasan-botany-cu
Status: Professional Modular Version for GitHub/Streamlit
"""

import os
import logging
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

# ১. লগিং সেটআপ (টার্মিনালে সুন্দর আউটপুট দেখার জন্য)
logging.basicConfig(
    level=logging.INFO, 
    format='%(levelname)s: %(message)s'
)

# ২. ডাটা ফোল্ডার চেক করার ফাংশন
def setup_directories():
    """প্রজেক্টের জন্য প্রয়োজনীয় ফোল্ডার তৈরি করে।"""
    for folder in ['data', 'results']:
        if not os.path.exists(folder):
            os.makedirs(folder)
            logging.info(f"'{folder}' ফোল্ডার তৈরি করা হয়েছে।")

# ৩. ডাটা সংগ্রহের ফাংশন (ধাপ ১)
def fetch_ncbi_data(accession_id, email):
    """NCBI থেকে GenBank ফাইল সংগ্রহ করে।"""
    Entrez.email = email
    file_path = f"data/{accession_id}.gbk"
    
    try:
        if not os.path.exists(file_path):
            logging.info(f"NCBI থেকে {accession_id} ডাউনলোড হচ্ছে...")
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            SeqIO.write(record, file_path, "genbank")
            return record
        else:
            logging.info(f"লোকাল ফাইল '{file_path}' থেকে ডাটা লোড করা হচ্ছে।")
            return SeqIO.read(file_path, "genbank")
    except Exception as e:
        logging.error(f"ডাটা সংগ্রহে ত্রুটি: {e}")
        return None

# ৪. সিকোয়েন্স অ্যানালাইসিসের ফাংশন (ধাপ ২)
def perform_analysis(record):
    """GC Content এবং প্রোটিন ট্রান্সলেশন সম্পন্ন করে।"""
    try:
        dna_seq = record.seq
        gc_val = gc_fraction(dna_seq) * 100
        protein_seq = dna_seq.translate(to_stop=True)
        
        logging.info(f"বিশ্লেষণ সম্পন্ন: দৈর্ঘ্য {len(dna_seq)} bp, GC {gc_val:.2f}%")
        
        # প্রোটিন ফাইলটি সেভ করা
        with open("data/mirabilis_protein.fasta", "w") as f:
            f.write(f">{record.id} protein\n{protein_seq}\n")
            
        return protein_seq, gc_val
    except Exception as e:
        logging.error(f"অ্যানালাইসিসে ত্রুটি: {e}")
        return None, None

# ==========================================
# ৫. মেইন এক্সিকিউশন (কন্ট্রোল রুম)
# ==========================================
if __name__ == "__main__":
    print("\n" + "="*50)
    print("MIRABILIS JALAPA GENETICS PIPELINE v2.0")
    print("="*50 + "\n")

    # ভেরিয়েবল সেটআপ
    MY_EMAIL = "zahidulhasan.botany.cu@gmail.com"
    TARGET_ACCESSION = "AB435372.1"

    # ধাপে ধাপে কল করা
    setup_directories()
    
    genbank_record = fetch_ncbi_data(TARGET_ACCESSION, MY_EMAIL)
    
    if genbank_record:
        protein, gc = perform_analysis(genbank_record)
        
        if protein:
            print("\n" + "-"*30)
            print(f"✅ ফলাফল সংগ্রহ সম্পন্ন!")
            print(f"🧬 জিনের নাম: {genbank_record.description}")
            print(f"🧪 প্রোটিন সিকোয়েন্সের প্রথম ১০টি অ্যামিনো অ্যাসিড: {protein[:10]}...")
            print("-"*30)
            
            logging.info("প্রজেক্টের প্রথম অংশ সফলভাবে কাজ করছে।")
    else:
        logging.error("পাইপলাইন রান করা সম্ভব হয়নি।")
