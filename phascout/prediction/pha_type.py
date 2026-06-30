"""
PHA Tipi Potansiyeli Yorumlayıcı
=================================
PhaC sınıfı + tespit edilen gen donanımı + operon kanıtını, bilinen biyokimyaya
(Rehm 2003; Steinbüchel) göre bir PHA tipi POTANSİYELİNE eşler.

ÖNEMLİ: Çıktı bir üretim hükmü DEĞİL, gen-temelli bir potansiyel ifadesidir.
Gen varlığı; ifade düzeyini, karbon kaynağını veya büyüme koşullarını garanti etmez.

Sınıf -> temel PHA tipi:
    Class I / III / IV -> SCL-PHA (kısa zincir, C3-C5)
    Class II           -> MCL-PHA (orta zincir, C6-C14)

Yardımcı genler rotayı/monomeri netleştirir:
    phaA + phaB -> asetil-CoA'dan 3HB (şekerden, SCL)
    phaG        -> de novo yağ asidi sentezinden MCL (Class II, şekerden)
    phaJ        -> beta-oksidasyondan (R)-3-hidroksiaçil-CoA besler. SENTAZA BAĞLI:
                   Class II ile MCL; SCL-class (I/III/IV) ile yağ-asidinden SCL
                   (3HHx ko-polimeri DEĞİL — sınıf substrat aralığıyla sınırlı).
"""

# Sınıf -> temel PHA tipi
_SCL_CLASSES = {"Class_I", "Class_III", "Class_IV"}
_MCL_CLASSES = {"Class_II"}

CAVEAT = (
    "Gen varlığı üretim POTANSİYELİNİ gösterir; gerçek PHA birikimi ifade düzeyine, "
    "karbon kaynağına ve büyüme koşullarına bağlıdır."
)


def _synteny_evidence(gene, gene_vector, operon_result):
    """phaA/phaB için dürüst, gen-bazlı kanıt etiketi (operon sintenisi ile).

    PhaB↔FabG (ve PhaA↔diğer thiolazlar) ayrımı dizi-tek-başına bir tavana
    sahiptir; en güvenilir ayırt edici phaC ile operon sintenisidir.
    """
    if not gene_vector.get(gene):
        return "tespit edilmedi"
    gev = (operon_result or {}).get("genes", {}).get(gene, {})
    rescued = (operon_result or {}).get("rescued", {}).get(gene)
    family = "SDR" if gene == "phaB" else "thiolaz"
    confuser = "FabG" if gene == "phaB" else "diğer thiolazlar"
    if rescued:
        return f"operon-kurtarma (phaC'ye {rescued['distance']} bp, BLOSUM {rescued['blosum']})"
    if gev.get("syntenic"):
        return f"operon-destekli (phaC'ye {gev.get('distance')} bp sintenik)"
    if not (operon_result or {}).get("available"):
        return f"BLOSUM-onaylı (operon/GFF verisi yok — {confuser} sinteni ile dışlanamadı)"
    return f"aday {family} — operon yok, {confuser}'den kesin ayrılamaz"


def classify_pha_potential(phac_result, gene_vector, operon_result=None):
    """
    PHA tipi potansiyelini belirle.

    Returns:
        dict: {
            "potential": "none" | "SCL" | "MCL" | "SCL-co-MCL",
            "products": [str, ...],
            "routes": [ {name, substrate, genes, product}, ... ],
            "confidence": "yok" | "dusuk" | "orta" | "yuksek",
            "gene_evidence": {gene: label},
            "caveat": str,
            "notes": [str, ...],
        }
    """
    operon_result = operon_result or {}
    phac_class = phac_result.get("best_class")
    functional = phac_result.get("is_functional", False)

    result = {
        "potential": "none",
        "products": [],
        "routes": [],
        "confidence": "yok",
        "gene_evidence": {},
        "caveat": CAVEAT,
        "notes": [],
    }

    # Gen kanıt etiketleri (dürüst)
    ev = {}
    for g in ["phaC", "phaG", "phaJ", "phaP", "phaR", "phaE"]:
        ev[g] = "tespit edildi" if gene_vector.get(g) else "tespit edilmedi"
    ev["phaA"] = _synteny_evidence("phaA", gene_vector, operon_result)
    ev["phaB"] = _synteny_evidence("phaB", gene_vector, operon_result)
    result["gene_evidence"] = ev

    # 1. Fonksiyonel PhaC yoksa potansiyel yok
    if not functional or phac_class not in (_SCL_CLASSES | _MCL_CLASSES):
        result["notes"].append(
            "Fonksiyonel PhaC (katalitik triad + box) doğrulanamadı; PHA potansiyeli raporlanmaz."
        )
        return result

    has_ab = gene_vector.get("phaA") and gene_vector.get("phaB")
    has_g = gene_vector.get("phaG")
    has_j = gene_vector.get("phaJ")

    routes = []

    if phac_class in _SCL_CLASSES:
        products = []
        potential = "SCL"

        if has_ab:
            products = ["P(3HB)", "P(3HB-co-3HV) (propiyonat/valerat ko-substratıyla)"]
            routes.append({
                "name": "Şekerden SCL (de novo)",
                "substrate": "şeker (glukoz, fruktoz, sükroz, gliserol)",
                "genes": ["phaA", "phaB", "phaC"],
                "product": "P(3HB)",
            })

        if has_j:
            # BIYOKIMYA: phaJ ((R)-spesifik enoyl-CoA hidrataz) beta-oksidasyondan
            # (R)-3-hidroksiaçil-CoA besler. SCL-SPESIFIK bir sentazla (Class I/III/IV,
            # substrat aralığı C4-C5) bu, yağ asidinden SCL (P3HB/P3HV) verir;
            # 3HHx (mcl) KOPOLIMERI VERMEZ. 3HHx ancak MCL-yetenekli (Class II) ya da
            # geniş-substratlı bir sentazla (ör. Aeromonas caviae PhaC) polimerleşir.
            # Bu yüzden phaJ TEK BAŞINA SCL-class organizmada SCL-co-MCL iddiasını
            # TETIKLEMEMELI (eski davranış sistematik over-claim üretiyordu). phaJ
            # yalnızca yağ-asidinden bir SCL besleme rotası olarak raporlanır.
            routes.append({
                "name": "Yağ asidinden SCL (β-oksidasyon / PhaJ)",
                "substrate": "yağ asitleri (oktanoat, dekanoat...)",
                "genes": ["phaJ", "phaC"],
                "product": "P(3HB)",
            })
            if "P(3HB)" not in products:
                products.insert(0, "P(3HB)")
            result["notes"].append(
                "phaJ tespit edildi: yağ asidi substratında β-oksidasyon besleme rotası. "
                "3HHx içeren SCL-co-MCL ko-polimeri yalnızca MCL-yetenekli (Class II) veya "
                "geniş-substratlı bir sentazla mümkündür; sınıf tek başına bunu göstermez "
                "(genuine ko-polimer üreticileri burada SCL olarak raporlanabilir — bilinen sınır)."
            )

        if routes:
            result["potential"] = potential
            result["products"] = products
            result["confidence"] = "yuksek" if has_ab else "orta"
        else:
            # Sentaz fonksiyonel ama besleme rotası YOK -> tip iddia edilmez
            result["potential"] = "belirsiz"
            result["confidence"] = "dusuk"
            result["notes"].append(
                "Fonksiyonel PhaC var; ancak monomer-sağlayan gen (phaA/phaB veya phaJ) "
                "tespit edilemedi — PHA tipi potansiyeli iddia edilemez."
            )

    else:  # Class II -> MCL
        if has_g:
            routes.append({
                "name": "Şekerden MCL (de novo yağ asidi sentezi / PhaG)",
                "substrate": "şeker (glukoz, fruktoz, gliserol)",
                "genes": ["phaG", "phaC"],
                "product": "MCL-PHA (3HD baskın)",
            })
        if has_j:
            routes.append({
                "name": "Yağ asidinden MCL (β-oksidasyon / PhaJ)",
                "substrate": "yağ asitleri (C6-C14)",
                "genes": ["phaJ", "phaC"],
                "product": "MCL-PHA (C6-C14)",
            })

        if routes:
            result["potential"] = "MCL"
            result["products"] = ["MCL-PHA (C6-C14 monomerleri)"]
            result["confidence"] = "yuksek"
        else:
            result["potential"] = "belirsiz"
            result["confidence"] = "dusuk"
            result["notes"].append(
                "Fonksiyonel Class II PhaC var; ancak monomer-sağlayan gen (phaG/phaJ) "
                "tespit edilemedi — MCL potansiyeli iddia edilemez."
            )

    result["routes"] = routes
    return result
