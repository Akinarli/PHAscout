"""
PHAscout Konfigürasyon Modülü
=============================
Aşama 0'dan gelen kalibre edilmiş eşikler ve sistem sabitleri.
Bu değerler bilimsel olarak (ROC/F1 analizi ile) belirlenmiştir.
"""

import os

# ============================================================
# PROJE YOL YAPISI
# ============================================================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
HMM_PROFILES_DIR = os.path.join(DATA_DIR, "hmm_profiles")
PHAC_CLASS_HMM_DIR = os.path.join(HMM_PROFILES_DIR, "phac_classes")
PFAM_HMM_DIR = os.path.join(HMM_PROFILES_DIR, "pfam")
REFERENCE_SEQ_DIR = os.path.join(DATA_DIR, "reference_sequences")
POSITIVE_REF_DIR = os.path.join(REFERENCE_SEQ_DIR, "positive")
NEGATIVE_REF_DIR = os.path.join(REFERENCE_SEQ_DIR, "negative")

# ============================================================
# PFAM PROFİL EŞLEMELERİ (Gen → PFAM Accession)
# ============================================================
PFAM_PROFILES = {
    "phaC": ["PF07167", "PF00561"],   # Alpha/beta hydrolase + Abhydrolase
    "phaA": ["PF00108", "PF02803"],   # Thiolase N-term + C-term
    "phaB": ["PF13561"],               # adh_short_C2 (SDR ailesi)
    "phaJ": ["PF00767"],               # Enoyl-CoA hydratase
    "phaG": ["PF00121"],               # Acyltransferase
    "phaP": ["PF09361"],               # Phasin
    "phaR": ["PF07879"],               # Regulator (transkripsiyon düzenleyici)
    "phaE": ["PF09712"],               # PhaE (Class III alt birim)
}

# ============================================================
# BLOSUM62 NORMALİZE SKOR EŞİKLERİ (Aşama 0 Kalibrasyonu)
# Bu değerler ROC/F1 analizi ile hesaplanmıştır.
# ============================================================
BLOSUM62_THRESHOLDS = {
    "phaB": 0.4778,  # F1 = 0.73 | Negatif: FabG + SDR broad
    "phaA": 0.7435,  # F1 = 0.82 | Negatif: FadA
    "phaG": 0.35,    # Henüz kalibre edilmedi (az veri), varsayılan
    "phaJ": 0.35,    # Henüz kalibre edilmedi (az veri), varsayılan
}

# ============================================================
# UZUNLUK FİLTRELERİ (Amino Asit Sayısı Aralıkları)
# Literatürden ve UniProt reviewed kayıtlardan belirlenmiştir.
# ============================================================
LENGTH_FILTERS = {
    "phaB": (230, 280),   # Acetoacetyl-CoA reductase
    "phaJ": (130, 160),   # Enoyl-CoA hydratase
    "phaG": (280, 350),   # 3-hydroxyacyl-ACP:CoA transacylase
    "phaA": (350, 420),   # Beta-ketothiolase
}

# ============================================================
# HMM TARAMA PARAMETRELERİ
# ============================================================
HMM_EVALUE_THRESHOLD = 1e-5  # E-value üst sınır
USE_GATHERING_THRESHOLD = True  # PFAM GA eşiği kullan

# ============================================================
# PhaC SINIF PROFİL DOSYALARI
# ============================================================
PHAC_CLASS_PROFILES = {
    "Class_I": os.path.join(PHAC_CLASS_HMM_DIR, "phac_class_I.hmm"),
    "Class_II": os.path.join(PHAC_CLASS_HMM_DIR, "phac_class_II.hmm"),
    "Class_III": os.path.join(PHAC_CLASS_HMM_DIR, "phac_class_III.hmm"),
    "Class_IV": os.path.join(PHAC_CLASS_HMM_DIR, "phac_class_IV.hmm"),
}

# ============================================================
# KATALİTİK TRİAD POZİSYONLARI (HMM Hizalama Referans)
# PhaC Box motifi: G[GS].C.[GA]G
# ============================================================
CATALYTIC_TRIAD = {
    "Cys": "C",
    "Asp": "D",
    "His": "H",
}
PHAC_BOX_REGEX = r"[GSY].C.[GSA]"  # Minimal lipase box: G-x-C-x-G and variants like S-Y-C-V-G (Halomonas)

# ============================================================
# CLASS III/IV ALT BİRİM KONTROL
# ============================================================
SUBUNIT_PENALTY = {
    "full_score": 40,      # Alt birim tamam
    "missing_score": 15,   # Alt birim eksik (Polimeraz indeksi düşürülür)
}

# ============================================================
# SEZGİSEL İNDEKS SABİTLERİ
# ============================================================
HEURISTIC_INDEX = {
    "polymerase_full": 40,
    "polymerase_no_subunit": 15,
    "polymerase_none": 0,
    "monomer_full": 40,
    "monomer_partial_phaA_only": 10,
    "monomer_partial_phaB_only": 10,
    "monomer_beta_no_phaJ": 30,
    "yield_phaP": 12,
    "max_score": 92,
}

# ============================================================
# NCBI API AYARLARI
# ============================================================
NCBI_DATASETS_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")  # Kullanıcı env'den sağlayabilir
NCBI_TIMEOUT = 60  # saniye

# ============================================================
# METABOLIK YOLAK TANIMLARI (Boolean Mantığı)
# ============================================================
PATHWAYS = {
    "alpha": {
        "name": "Şekerden SCL-PHA (P3HB)",
        "required_genes": ["phaC", "phaA", "phaB"],
        "valid_phac_classes": ["Class_I", "Class_III", "Class_IV"],
        "carbon_sources": ["glucose", "fructose", "sucrose", "glycerol", "acetate"],
        "product_tendency": "P(3HB)",
    },
    "beta": {
        "name": "Yağ Asidinden MCL-PHA",
        "required_genes": ["phaC"],
        "valid_phac_classes": ["Class_II"],
        "carbon_sources": ["octanoate", "decanoate", "dodecanoate", "fatty_acids"],
        "product_tendency": "MCL-PHA (C6-C14 monomerleri)",
    },
    "gamma": {
        "name": "Şekerden MCL-PHA (de novo FAS)",
        "required_genes": ["phaC", "phaG"],
        "valid_phac_classes": ["Class_II"],
        "carbon_sources": ["glucose", "fructose", "sucrose", "glycerol"],
        "product_tendency": "MCL-PHA (3HD baskın)",
    },
    "delta": {
        "name": "Şeker + VFA Ko-substrat",
        "required_genes": ["phaC", "phaA", "phaB"],
        "optional_genes": ["phaJ"],
        "valid_phac_classes": ["Class_I", "Class_III", "Class_IV"],
        "carbon_sources": ["glucose+propionate", "glucose+valerate"],
        "product_tendency": "P(3HB-co-3HV)",
    },
}
