"""
Operon / Synteni Kanıt Modülü
==============================
Seçilen (fonksiyonel) phaC'nin genomik komşuluğunu analiz ederek phaA/phaB
gibi yardımcı genler için SINTENI kanıtı üretir.

Biyolojik gerekçe: PhaB'yi FabG'den (ya da PhaA'yı diğer thiolazlardan) tek
başına dizi benzerliğiyle ayırmak bir tavana sahiptir (her ikisi de aynı
üst-aileden). En güvenilir genomik ayırt edici, genin phaC ile aynı operonda
(birkaç kb içinde, aynı kontig) bulunmasıdır. Bu modül bunu iki amaçla kullanır:
  1. KANIT ETİKETİ: bir yardımcı gen phaC ile sintenik mi? (operon-destekli vs aday)
  2. KURTARMA: BLOSUM eşiğinin biraz altında kalan ama phaC ile SIKI sintenik
     (operon içinde) bir aday, operon kanıtıyla onaylanabilir.

Not: Operonu olmayan organizmalarda (ör. Class IV Bacillus, phaA/phaB operon-bağlı
değildir) hiçbir gen kurtarılmaz — bu doğru ve dürüst davranıştır.
"""

OPERON_MAX_BP = 3500       # operon komşuluğu eşiği (kb mertebesinde sıkı)
# Operon-kurtarma için minimum BLOSUM (gürültüyü ele; pipeline kullanır).
RESCUE_BLOSUM_FLOOR = {
    "phaA": 0.45,
    "phaB": 0.35,
}


def _interval_distance(c1, c2):
    """İki genomik özellik arasındaki minimum mesafe (bp). Farklı kontig -> inf."""
    if not c1 or not c2:
        return None
    if c1["contig"] != c2["contig"]:
        return float("inf")
    if c1["end"] < c2["start"]:
        return c2["start"] - c1["end"]
    if c2["end"] < c1["start"]:
        return c1["start"] - c2["end"]
    return 0  # örtüşen


def analyze_operon_evidence(phac_protein_id, raw_hmm_results, gff_data,
                            max_distance_bp=OPERON_MAX_BP):
    """
    Seçilen phaC'ye göre phaA/phaB için en yakın sintenik adayı bulur.

    Args:
        phac_protein_id: Seçilen (fonksiyonel) phaC protein_id'si.
        raw_hmm_results: HMMScanner.scan() ham çıktısı (gen -> hit listesi).
        gff_data: protein_id -> {contig, start, end, strand}.

    Returns:
        dict: {
          "available": bool,                # GFF + phaC koordinatı var mı
          "phac_protein_id": str,
          "is_class_i_operon": bool,        # phaA veya phaB phaC ile sintenik mi
          "genes": {
             gene: {"syntenic": bool, "distance": int|None, "protein_id": str|None}
          }
        }
    """
    result = {
        "available": False,
        "phac_protein_id": phac_protein_id,
        "is_class_i_operon": False,
        "genes": {},
    }

    if not gff_data or not phac_protein_id:
        return result
    phac_coords = gff_data.get(phac_protein_id)
    if not phac_coords:
        return result
    result["available"] = True

    for gene in ("phaA", "phaB"):
        best = {"syntenic": False, "distance": None, "protein_id": None}
        seen = set()
        for hit in raw_hmm_results.get(gene, []):
            pid = hit.get("protein_id")
            if not pid or pid in seen:
                continue
            seen.add(pid)
            dist = _interval_distance(phac_coords, gff_data.get(pid))
            if dist is None:
                continue
            if dist <= max_distance_bp and (best["distance"] is None or dist < best["distance"]):
                best["distance"] = dist
                best["protein_id"] = pid
                best["syntenic"] = True
        result["genes"][gene] = best
        if best["syntenic"]:
            result["is_class_i_operon"] = True

    return result


# Geriye uyumluluk: eski imza (yalnızca BLOSUM-onaylı genlerin sintenisi)
def analyze_operon(detected_genes: dict, gff_data: dict, max_distance_bp=OPERON_MAX_BP) -> dict:
    results = {"is_class_i_operon": False, "distances": {}}
    if not gff_data:
        return results

    def coords(name):
        gi = detected_genes.get(name, {})
        if gi.get("detected"):
            return gff_data.get(gi.get("protein_id"))
        return None

    pc = coords("phaC")
    if not pc:
        return results
    for g in ("phaA", "phaB"):
        d = _interval_distance(pc, coords(g))
        if d is not None:
            results["distances"][f"phaC-{g}"] = d
            if d <= max_distance_bp:
                results["is_class_i_operon"] = True
    return results
