"""
PHAscout Dürüst Benchmark — Kör Koşu Runner
============================================
dataset.csv'deki her accession için TAM genom pipeline'ını BİR KEZ koşturur.

Disiplin:
  - Ham rapor cache/{accession}.json'a yazılır. Tekrar çalıştırıldığında
    cache'ten okunur (NCBI'a tekrar gidilmez, "tek koşu" garanti). Yeniden
    koşmak için --force.
  - Dürüstlük doğrulaması: atıfsız veya evidence_method='annotation' satırlar
    REDDEDİLİR (benchmark/README.md demir kuralları).
  - Pipeline bir genomda patlarsa kayıt 'error' olarak işaretlenir, koşu durmaz.

Kullanım:
  python -m benchmark.run --dataset benchmark/dataset.csv
  python -m benchmark.run --dataset benchmark/dataset.csv --force
"""

import os
import csv
import json
import argparse
import logging
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(HERE, "cache")
PREDICTIONS_CSV = os.path.join(HERE, "predictions.csv")

# Dürüstlük: bu evidence_method değerleri ıslak-lab DEĞİLDİR -> reddet.
FORBIDDEN_EVIDENCE = {"", "annotation", "anotasyon", "ncbi", "genome", "in_silico", "prediction"}

PRED_FIELDS = [
    "accession", "organism_resolved", "status",
    "pred_potential", "pred_confidence", "pred_produces",
    "pred_class", "phac_functional", "phac_confidence",
    "detected_genes", "error",
]


def _validate_row(row):
    """Demir kuralları zorla. Sorun varsa hata mesajı döndür, yoksa None.

    Epistemik ayrım:
      - POZİTİF iddia (label_produces=yes): ıslak-lab kanıtı ZORUNLU. Anotasyon/
        in-silico kanıtla "üretir" demek döngüseldir (bizim aracın yaptığını
        başka bir gen-bulucuyla kıyaslamak) -> yasak.
      - NEGATİF iddia (label_produces=no): genomik-yokluk + rapor-yokluğu meşru
        bir temeldir (zayıf olsa da). Atıf yine zorunlu, ama evidence_method
        kısıtı uygulanmaz.
    """
    if not row.get("accession", "").strip():
        return "accession boş"
    if not row.get("citation", "").strip():
        return "atıf (citation) yok — atıfsız satır koşturulmaz"
    produces = row.get("label_produces", "").strip().lower()
    if produces not in ("yes", "no"):
        return "label_produces yes/no olmalı"
    if produces == "yes":
        ev = row.get("evidence_method", "").strip().lower()
        if ev in FORBIDDEN_EVIDENCE:
            return f"POZİTİF iddia için evidence_method='{ev}' ıslak-lab değil (annotation/in_silico yasak)"
    return None


def _extract_prediction(report):
    """Pipeline raporundan benchmark için karşılaştırılabilir alanları çek."""
    s = report.get("summary", {})
    pa = report.get("phac_analysis", {})
    genes = report.get("genes", {}).get("detected", [])
    return {
        "organism_resolved": report.get("organism", {}).get("organism_name", "N/A"),
        "pred_potential": s.get("pha_potential", "none"),
        "pred_confidence": s.get("potential_confidence", "yok"),
        "pred_produces": s.get("produces_pha", False),
        "pred_class": s.get("phac_class", "N/A"),
        "phac_functional": pa.get("functional", False),
        "phac_confidence": pa.get("confidence", 0),
        "detected_genes": "|".join(genes),
    }


def run(dataset_path, force=False):
    os.makedirs(CACHE_DIR, exist_ok=True)

    with open(dataset_path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print(f"[!] {dataset_path} boş. Önce held-out genomlarla doldur.")
        return

    # Pipeline'ı bir kez kur (HMM/referans yüklemesi pahalı)
    from phascout.pipeline import PHAscoutPipeline
    pipeline = PHAscoutPipeline()

    predictions = []
    n_ok = n_err = n_skip = 0

    for i, row in enumerate(rows, 1):
        acc = row.get("accession", "").strip()
        problem = _validate_row(row)
        if problem:
            print(f"[{i}/{len(rows)}] ATLANDI {acc or '(boş)'}: {problem}")
            n_skip += 1
            continue

        cache_path = os.path.join(CACHE_DIR, f"{acc}.json")
        rec = {"accession": acc, "status": "ok", "error": ""}

        if os.path.exists(cache_path) and not force:
            with open(cache_path, encoding="utf-8") as cf:
                report = json.load(cf)
            rec.update(_extract_prediction(report))
            rec["status"] = "cached"
            print(f"[{i}/{len(rows)}] CACHE {acc} -> {rec['pred_potential']}")
            predictions.append(rec)
            n_ok += 1
            continue

        try:
            report = pipeline.run(accession=acc)
            with open(cache_path, "w", encoding="utf-8") as cf:
                json.dump(report, cf, ensure_ascii=False, default=str)
            rec.update(_extract_prediction(report))
            print(f"[{i}/{len(rows)}] OK    {acc} -> {rec['pred_potential']} "
                  f"({rec['pred_class']}, fonk={rec['phac_functional']})")
            n_ok += 1
        except Exception as e:
            rec["status"] = "error"
            rec["error"] = str(e)
            print(f"[{i}/{len(rows)}] HATA  {acc}: {e}")
            n_err += 1

        predictions.append(rec)

    # predictions.csv yaz
    with open(PREDICTIONS_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=PRED_FIELDS)
        w.writeheader()
        for rec in predictions:
            w.writerow({k: rec.get(k, "") for k in PRED_FIELDS})

    print(f"\n{'='*50}")
    print(f"Koşu bitti: {n_ok} başarılı, {n_err} hata, {n_skip} atlandı (kural ihlali).")
    print(f"Tahminler: {PREDICTIONS_CSV}")
    print(f"Ham raporlar: {CACHE_DIR}/")
    print("Şimdi: python -m benchmark.metrics")


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.path.join(HERE, "dataset.csv"))
    ap.add_argument("--force", action="store_true", help="Cache'i yok say, yeniden koştur")
    args = ap.parse_args()
    run(args.dataset, force=args.force)
