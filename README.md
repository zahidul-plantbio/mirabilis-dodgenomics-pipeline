# In silico Identification and Sequence Analysis of Flower Color Genes in Mirabilis jalapa

## 📌 Project Overview
This repository contains a professional bioinformatics pipeline designed to analyze the **4,5-DOPA Dioxygenase (DOD)** gene in *Mirabilis jalapa* (Four o'clock flower). This gene plays a critical role in the betalain biosynthetic pathway, contributing to the species' unique floral pigmentation and Mendelian incomplete dominance.

## 🧬 Key Features
- **Automated Data Retrieval:** Fetches gene and protein sequences directly from NCBI Entrez.
- **Structural Bioinformatic Analysis:** Performs sequence quality checks, GC content analysis, and protein translation.
- **Homology Search & Evolutionary Mapping:** Integrated BLASTp search using the Swiss-Prot database and construction of a Neighbor-Joining phylogenetic tree.
- **Structural Modeling:** High-reliability 3D protein structure prediction with a **GMQE score of 1.00**.

## 📂 Repository Structure
```text
mirabilis-dodgenomics-pipeline/
├── data/               # Raw GenBank and FASTA sequences
├── results/            # BLAST reports, phylogenetic trees, and diagrams
├── scripts/
│   └── analysis_pipeline.py  # Main modular analysis engine
├── requirements.txt    # List of dependencies
└── README.md           # Project documentation
```
## 📊 Results Summary
The pipeline successfully identified a highly conserved LigB superfamily domain in the MjDOD protein. The phylogenetic analysis confirms its evolutionary relationship with other Caryophyllales, such as Beta vulgaris.

## 👨‍🔬 Author
**Zahidul Hasan** 
Department of Botany, University of Chittagong, Bangladesh.  
*Specialization: Plant Bioinformatics*

## 📧 Contact
If you have any questions, suggestions, or potential collaborations, feel free to reach out:

* **LinkedIn:** [Zahidul Hasan](https://www.linkedin.com/in/zahidulhasan-botany-cu)
* **Email:** [zahidulhasan.botany.cu@gmail.com](mailto:zahidul.hasan.botany@gmail.com)