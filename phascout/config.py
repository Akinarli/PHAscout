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
    "phaZ": ["PF05898"],               # PHA depolymerase
}

# ============================================================
# BLOSUM62 NORMALİZE SKOR EŞİKLERİ (Aşama 0 Kalibrasyonu)
# Bu değerler min(qq, rr) normalizasyonu ile kalibre edilmiştir.
# Kontamine referanslar (FabG/MabA/FabM/antiporter/PhaC) temizlendikten
# sonra scripts/calibrate_thresholds.py ile yeniden hesaplanmıştır.
# NOT: Pozitif referans setleri küçük (n=3-9); eşikler tek başına
# kesin değil, double-layer'ın tüm gücü için ileride genişletilmeli.
# ============================================================
BLOSUM62_THRESHOLDS = {
    "phaB": 0.4899,  # F1 = 0.923 (n=7 pos / 40 neg: FabG + SDR broad). Eski kirli set: F1=0.50
    "phaA": 0.7161,  # F1 = 0.941 (n=9 pos / 96 neg: FadA)
    "phaG": 0.5823,  # F1 = 1.000 (n=3 pos / 37 neg: FabG) - az veri
    "phaJ": 0.4186,  # F1 = 1.000 (n=3 pos / 5 neg: crotonase ailesi) - cok az veri, dusuk guven
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
# ============================================================
CATALYTIC_TRIAD = {
    "Cys": "C",
    "Asp": "D",
    "His": "H",
}

# Sınıf-spesifik katalitik HMM MATCH-state kolonları.
# Bu kolonlar, deneysel olarak doğrulanmış katalitik kalıntılar
# (UniProt ACT_SITE Cys + literatür Asp/His) bir referans diziye
# hizalanarak ve scripts/derive_catalytic_columns.py ile mevcut
# HMM modellerine eşlenerek türetilmiştir (yeniden üretilebilir).
#   Class_I : C. necator P23608  C319/D480/H508
#   Class_II: P. oleovorans P26494 C296/D451/H479
#   Class_III: A. vinosum P45370  C149/D302/H331
#   Class_IV: B. albus A0A1J9SXW8 C151/D306/H335
CATALYTIC_HMM_COLUMNS = {
    "Class_I":   {"Cys": 326, "Asp": 481, "His": 510},
    "Class_II":  {"Cys": 297, "Asp": 452, "His": 480},
    "Class_III": {"Cys": 162, "Asp": 315, "His": 345},
    "Class_IV":  {"Cys": 150, "Asp": 304, "His": 333},
}

# Triad kolon eşleştirme toleransı (HMM match-state kolonu cinsinden).
# Eski deger 25 (~50 kalintilik pencere) idi ve katalitik spesifikligi
# yok ediyordu. Katalitik kalinti TEK bir hizalama kolonundadir;
# kucuk bir pay sadece HMM build varyasyonlari icindir.
TRIAD_TOLERANCE = 2

# PhaC Box (lipase box) motifi: G-x-C-x-G-G ve varyantlari.
# Box icindeki Cys, katalitik nukleofil ile AYNI kalinti olmalidir;
# bu kosul validator'da hmm-kolon hizalamasiyla zorlanir.
PHAC_BOX_REGEX = r"G.C.G[GA]"  # G-x-C-x-G-[G/A]: katalitik Cys nukleofil dirseginde

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
    "penalty_phaZ": 8,
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
