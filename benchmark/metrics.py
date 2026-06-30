"""
PHAscout Dürüst Benchmark — Metrik Hesaplama
=============================================
dataset.csv (etiket) + predictions.csv (tahmin) -> dürüst metrikler.

Dürüstlük ilkeleri:
  - İKİ AYRI SORU: (1) Tespit: üretici mi değil mi? (2) Tip: üreticiler
    arasında SCL/MCL/ko-polimer doğru mu? Bunları karıştırmak yanıltır.
  - ÇEKİMSERLİK: "belirsiz" (sentaz var, besleme rotası yok) bir tahmin
    DEĞİL, bir "bilmiyorum"dur. Ayrı raporlanır; pozitif/negatif'e zorlanmaz.
  - Negatif-set zayıf: gerçek deneysel negatif azdır. Özgüllük/precision
    raporlanır AMA açık uyarıyla — bu sayılara güvenme.
  - Wilson %95 GA: küçük N'de dürüst belirsizlik aralığı.
  - Filum + sınıf kırılımı: ortalama bir sayı yalan söyler; NEREDE patladığını göster.
"""

import os
import csv
import math
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATASET = os.path.join(HERE, "dataset.csv")
PREDICTIONS = os.path.join(HERE, "predictions.csv")
# Çalışma zamanında --dataset ile override edilebilir (smoke testi için).

POSITIVE_POTENTIALS = {"SCL", "MCL", "SCL-co-MCL"}


def wilson_ci(k, n, z=1.96):
    """Wilson skor %95 güven aralığı (oran için)."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2*n)) / denom
    half = (z * math.sqrt(p*(1-p)/n + z**2/(4*n**2))) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def _norm_type(t):
    """Tip etiketini normalize et (SCL-co-MCL ve eşanlamlıları)."""
    t = (t or "").strip()
    low = t.lower().replace("_", "-").replace(" ", "")
    if low in ("scl-co-mcl", "sclcomcl", "co-polymer", "kopolimer", "scl-mcl"):
        return "SCL-co-MCL"
    if low.startswith("scl"):
        return "SCL"
    if low.startswith("mcl"):
        return "MCL"
    if low in ("none", "yok", ""):
        return "none"
    return t


def load_joined():
    if not os.path.exists(PREDICTIONS):
        raise SystemExit(f"[!] {PREDICTIONS} yok. Önce: python -m benchmark.run")
    with open(DATASET, encoding="utf-8-sig") as f:
        labels = {r["accession"].strip(): r for r in csv.DictReader(f)}
    with open(PREDICTIONS, encoding="utf-8-sig") as f:
        preds = {r["accession"].strip(): r for r in csv.DictReader(f)}

    joined = []
    for acc, lab in labels.items():
        pr = preds.get(acc)
        if not pr or pr.get("status") == "error":
            continue
        joined.append((acc, lab, pr))
    return joined, labels, preds


def report():
    joined, labels, preds = load_joined()
    if not joined:
        raise SystemExit("[!] Eşleşen (etiket+tahmin) satır yok. dataset doldurulmuş ve run koşmuş olmalı.")

    # ---- TESPIT (binary, çekimserlik ayrı) ----
    TP = TN = FP = FN = ABST_POS = ABST_NEG = 0
    detection_rows = []
    for acc, lab, pr in joined:
        truth = lab["label_produces"].strip().lower() == "yes"
        pot = pr["pred_potential"].strip()
        if pot == "belirsiz":
            if truth: ABST_POS += 1
            else: ABST_NEG += 1
            detection_rows.append((acc, truth, "ABSTAIN", pot))
            continue
        pred_pos = pot in POSITIVE_POTENTIALS
        if truth and pred_pos: TP += 1; tag = "TP"
        elif truth and not pred_pos: FN += 1; tag = "FN"
        elif not truth and pred_pos: FP += 1; tag = "FP"
        else: TN += 1; tag = "TN"
        detection_rows.append((acc, truth, tag, pot))

    n_decided = TP + TN + FP + FN
    n_abstain = ABST_POS + ABST_NEG

    print("=" * 64)
    print("PHAscout DÜRÜST BENCHMARK — SONUÇ")
    print("=" * 64)
    print(f"Toplam eşleşen genom: {len(joined)}  (karar: {n_decided}, çekimser: {n_abstain})")

    print("\n--- 1) TESPİT (üretici mi?) ---")
    print(f"  TP={TP}  FN={FN}  FP={FP}  TN={TN}")
    print(f"  Çekimser (belirsiz): pozitiflerde {ABST_POS}, negatiflerde {ABST_NEG}")

    if TP + FN > 0:
        rec = TP / (TP + FN)
        lo, hi = wilson_ci(TP, TP + FN)
        print(f"  Duyarlılık (recall)  = {rec:.3f}  [%95 GA {lo:.3f}-{hi:.3f}]  (n={TP+FN})")
    if TP + FP > 0:
        prec = TP / (TP + FP)
        print(f"  Kesinlik (precision) = {prec:.3f}   (n={TP+FP})")
    if TN + FP > 0:
        spec = TN / (TN + FP)
        lo, hi = wilson_ci(TN, TN + FP)
        print(f"  Özgüllük (specificity)= {spec:.3f}  [%95 GA {lo:.3f}-{hi:.3f}]  (n={TN+FP})")
        print(f"    [!] DİKKAT: negatif seti (n={TN+FP}) küçük/zayıf olabilir — bu sayıya az güven.")
    if n_decided > 0:
        acc_v = (TP + TN) / n_decided
        print(f"  Doğruluk (karar verilenlerde) = {acc_v:.3f}")

    # ---- TİP DOĞRULUĞU (gerçek üreticiler arasında) ----
    print("\n--- 2) TİP DOĞRULUĞU (yalnızca gerçek üreticiler) ---")
    type_correct = type_total = 0
    type_confusion = defaultdict(lambda: defaultdict(int))
    for acc, lab, pr in joined:
        if lab["label_produces"].strip().lower() != "yes":
            continue
        lt = _norm_type(lab.get("label_type"))
        pt = _norm_type(pr.get("pred_potential"))
        if pt in ("none", "belirsiz"):
            continue  # tip iddiası yok -> tip metriğine girmez (tespit'te zaten FN/çekimser)
        type_total += 1
        type_confusion[lt][pt] += 1
        if lt == pt:
            type_correct += 1
    if type_total > 0:
        ta = type_correct / type_total
        lo, hi = wilson_ci(type_correct, type_total)
        print(f"  Tip doğruluğu = {type_correct}/{type_total} = {ta:.3f}  [%95 GA {lo:.3f}-{hi:.3f}]")
        print("  Karışıklık (satır=gerçek, sütun=tahmin):")
        for lt in sorted(type_confusion):
            for pt in sorted(type_confusion[lt]):
                flag = "" if lt == pt else "  <-- HATA"
                print(f"    {lt:>12} -> {pt:<12} : {type_confusion[lt][pt]}{flag}")
    else:
        print("  (Tip iddiası olan doğru-pozitif yok)")

    # ---- KIRILIM: filum ve sınıf ----
    def breakdown(key):
        agg = defaultdict(lambda: {"correct": 0, "total": 0})
        for acc, lab, pr in joined:
            grp = (lab.get(key) or "bilinmiyor").strip() or "bilinmiyor"
            truth = lab["label_produces"].strip().lower() == "yes"
            pot = pr["pred_potential"].strip()
            if pot == "belirsiz":
                continue
            pred_pos = pot in POSITIVE_POTENTIALS
            agg[grp]["total"] += 1
            if pred_pos == truth:
                agg[grp]["correct"] += 1
        return agg

    for key, title in (("phylum", "FİLUM"), ("synthase_class", "BEKLENEN SINIF")):
        print(f"\n--- KIRILIM: {title} (tespit doğruluğu) ---")
        agg = breakdown(key)
        for grp in sorted(agg):
            c, n = agg[grp]["correct"], agg[grp]["total"]
            print(f"  {grp:>16}: {c}/{n}" + (f" ({c/n:.0%})" if n else ""))

    # ---- HATA LİSTESİ (teşhis için en değerli kısım) ----
    print("\n--- YANLIŞLAR (hata analizi) ---")
    any_err = False
    for acc, truth, tag, pot in detection_rows:
        if tag in ("FP", "FN"):
            any_err = True
            org = labels[acc].get("organism", "")
            lt = labels[acc].get("label_type", "")
            print(f"  [{tag}] {acc} {org} | gerçek={'üretir' if truth else 'üretmez'}/{lt} | tahmin={pot}")
    if not any_err:
        print("  (Karar verilenlerde FP/FN yok — ama N'e ve çekimserlere bak.)")

    print("\n" + "=" * 64)
    print("Hatırlatma: Bu sette ayar yapma. İyileştirmeyi TAZE held-out sette doğrula.")
    print("=" * 64)


if __name__ == "__main__":
    import sys
    import argparse
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=DATASET)
    args = ap.parse_args()
    DATASET = args.dataset
    report()
