"""
Halomonas Cins-İçi PHA Genomik-Potansiyel Taraması (dürüst uygulama)
====================================================================
74 Halomonas genomunu (table3_dataset.csv) PHAscout pipeline'ından geçirir ve
HER suş için PHA TİPİ POTANSİYELİNİ (SCL / MCL / belirsiz / none) + PhaC sınıfı +
tespit edilen rota genlerini içeren DÜRÜST bir tablo üretir.

ÖNEMLİ ÇERÇEVE (asla ihlal etme):
  - Çıktı GENOMİK POTANSİYELDİR, üretim DEĞİLDİR. Gen varlığı PHA birikimini
    kanıtlamaz; deneysel doğrulama gerekir.
  - Halomonas zaten bilinen bir PHA-üretici cinsidir; yüksek "potansiyel" oranı
    BEKLENEN bir sonuçtur, bir "keşif" değildir. Bu bir KARAKTERİZASYON/HARİTALAMA
    çalışmasıdır.
  - Etiketsiz set → accuracy/precision İDDİA EDİLEMEZ.

Önbellek: her rapor analysis/halomonas/cache/{accession}.json'a yazılır; tekrar
koşunca NCBI'a gidilmez. Yeniden indirmek için --force.

Kullanım:
  python -m scripts.screen_halomonas
  python -m scripts.screen_halomonas --force
"""

import os
import csv
import sys
import json
import argparse
import logging

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
if BASE not in sys.path:
    sys.path.insert(0, BASE)  # 'phascout' paketini bul (hem -m hem path ile çalışsın)
DATASET = os.path.join(BASE, "table3_dataset.csv")
OUT_DIR = os.path.join(BASE, "analysis", "halomonas")
CACHE_DIR = os.path.join(OUT_DIR, "cache")
RESULTS_CSV = os.path.join(OUT_DIR, "results.csv")

FIELDS = [
    "accession", "species", "status",
    "phac_class", "phac_functional", "phac_confidence",
    "pha_potential", "potential_confidence", "pha_products",
    "detected_genes", "ml_probability", "error",
]


def _extract(report):
    s = report.get("summary", {})
    pa = report.get("phac_analysis", {})
    genes = report.get("genes", {}).get("detected", [])
    products = report.get("pha_potential", {}).get("products", [])
    return {
        "phac_class": s.get("phac_class", "N/A"),
        "phac_functional": pa.get("functional", False),
        "phac_confidence": pa.get("confidence", 0),
        "pha_potential": s.get("pha_potential", "none"),
        "potential_confidence": s.get("potential_confidence", "yok"),
        "pha_products": " ; ".join(products),
        "detected_genes": "|".join(genes),
        "ml_probability": s.get("ml_probability", 0.0),
    }


def run(force=False):
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(DATASET, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    from phascout.pipeline import PHAscoutPipeline
    pipeline = PHAscoutPipeline()

    results = []
    n_ok = n_err = 0
    for i, row in enumerate(rows, 1):
        acc = (row.get("Accession") or "").strip()
        species = (row.get("Species") or "").strip()
        if not acc:
            continue
        rec = {"accession": acc, "species": species, "status": "ok", "error": ""}
        cache_path = os.path.join(CACHE_DIR, f"{acc}.json")

        if os.path.exists(cache_path) and not force:
            with open(cache_path, encoding="utf-8") as cf:
                report = json.load(cf)
            rec.update(_extract(report))
            rec["status"] = "cached"
            print(f"[{i}/{len(rows)}] CACHE {acc} {species} -> {rec['pha_potential']} ({rec['phac_class']})")
            results.append(rec)
            n_ok += 1
            continue

        try:
            report = pipeline.run(accession=acc)
            with open(cache_path, "w", encoding="utf-8") as cf:
                json.dump(report, cf, ensure_ascii=False, default=str)
            rec.update(_extract(report))
            print(f"[{i}/{len(rows)}] OK    {acc} {species} -> {rec['pha_potential']} "
                  f"({rec['phac_class']}, fonk={rec['phac_functional']})")
            n_ok += 1
        except Exception as e:  # noqa: BLE001 - tek genom patlarsa koşu durmaz
            rec["status"] = "error"
            rec["error"] = str(e)
            print(f"[{i}/{len(rows)}] HATA  {acc} {species}: {e}")
            n_err += 1
        results.append(rec)

    with open(RESULTS_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for rec in results:
            w.writerow({k: rec.get(k, "") for k in FIELDS})

    print(f"\n{'='*56}")
    print(f"Tarama bitti: {n_ok} başarılı, {n_err} hata.")
    print(f"Sonuç tablosu: {RESULTS_CSV}")
    print("Sonra: python -m scripts.summarize_halomonas")


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    run(force=args.force)
