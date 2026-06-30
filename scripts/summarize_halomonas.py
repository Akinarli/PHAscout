"""
Halomonas Taraması — Özet + Görselleştirme (dürüst)
====================================================
analysis/halomonas/results.csv'i okur; cins-içi PHA genomik-potansiyel
özetini ve iki portfolyo figürünü üretir:
  1) gene_heatmap.png  — 74 suş × PHA genleri varlık ısı haritası (potansiyele göre sıralı)
  2) potential_distribution.png — PHA-tipi potansiyel dağılımı

Islak-lab overlay: analysis/halomonas/wetlab_known.csv varsa (accession,label_type,
citation), o suşlar özet ve figürde '★' ile işaretlenir — accuracy iddiası DEĞİL,
yalnızca destekleyici örnekler.

ÇERÇEVE: Çıktı genomik POTANSİYELDİR, üretim değildir. Halomonas zaten bilinen bir
PHA-üretici cinsidir; yüksek potansiyel oranı beklenendir (keşif değil).

Kullanım:  python scripts/summarize_halomonas.py
"""

import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
OUT_DIR = os.path.join(BASE, "analysis", "halomonas")
RESULTS = os.path.join(OUT_DIR, "results.csv")
WETLAB = os.path.join(OUT_DIR, "wetlab_known.csv")

GENES = ["phaC", "phaA", "phaB", "phaJ", "phaG", "phaP", "phaR", "phaE"]
# Sıralama önceliği: somut tip -> belirsiz -> none
POT_ORDER = {"SCL": 0, "MCL": 1, "SCL-co-MCL": 2, "belirsiz": 3, "none": 4}


def _short_species(s):
    s = (s or "").strip()
    return (s[:38] + "…") if len(s) > 39 else s


def main():
    if not os.path.exists(RESULTS):
        sys.exit(f"[!] {RESULTS} yok. Önce: python scripts/screen_halomonas.py")
    df = pd.read_csv(RESULTS).fillna("")
    df = df[df["status"].isin(["ok", "cached"])].copy()
    if df.empty:
        sys.exit("[!] Başarılı kayıt yok.")

    # Islak-lab overlay
    wet = {}
    if os.path.exists(WETLAB):
        w = pd.read_csv(WETLAB).fillna("")
        wet = {str(r["accession"]).strip(): r for _, r in w.iterrows()}

    # Gen varlık matrisi
    for g in GENES:
        df[g] = df["detected_genes"].apply(lambda s, g=g: int(g in str(s).split("|")))

    df["pot_rank"] = df["pha_potential"].map(lambda p: POT_ORDER.get(p, 5))
    df = df.sort_values(["pot_rank", "phac_class", "species"]).reset_index(drop=True)
    df["wet"] = df["accession"].apply(lambda a: a in wet)

    # ---- DÜRÜST ÖZET (stdout) ----
    n = len(df)
    print("=" * 60)
    print(f"HALOMONAS CİNS-İÇİ PHA GENOMİK-POTANSİYEL TARAMASI (N={n})")
    print("=" * 60)
    print("UYARI: genomik potansiyel ≠ üretim. Deneysel doğrulama gerekir.")
    print("Halomonas bilinen PHA-üretici cinsidir; yüksek oran BEKLENİR.\n")

    print("--- PHA tipi potansiyeli dağılımı ---")
    for pot, c in df["pha_potential"].value_counts().items():
        print(f"  {pot:>12}: {c:>3}  ({c/n:.0%})")
    print("\n--- PhaC sınıfı (fonksiyonel olanlar) ---")
    fn = df[df["phac_functional"] == True]  # noqa: E712
    for cls, c in fn["phac_class"].value_counts().items():
        print(f"  {cls:>10}: {c:>3}")
    print(f"\nFonksiyonel PhaC: {len(fn)}/{n} ({len(fn)/n:.0%})")
    if wet:
        print(f"\nIslak-lab destekli örnek (★): {df['wet'].sum()} suş")

    # ---- FIGÜR 1: gen ısı haritası ----
    mat = df[GENES].to_numpy(dtype=float)
    fig_h = max(6.0, 0.22 * n + 1.5)
    fig, ax = plt.subplots(figsize=(8.5, fig_h))
    ax.imshow(mat, aspect="auto", cmap="Greens", vmin=0, vmax=1)
    ax.set_xticks(range(len(GENES)))
    ax.set_xticklabels(GENES, rotation=45, ha="right", fontsize=9)
    labels = []
    for _, r in df.iterrows():
        star = "★ " if r["wet"] else ""
        labels.append(f"{star}{_short_species(r['species'])}  [{r['pha_potential']}/{r['phac_class']}]")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=6.2)
    # grid
    ax.set_xticks(np.arange(-.5, len(GENES), 1), minor=True)
    ax.set_yticks(np.arange(-.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    ax.set_title("Halomonas (N=%d): PHA biyosentez geni varlığı\n"
                 "(genomik potansiyel — üretim kanıtı değil)" % n, fontsize=10)
    fig.tight_layout()
    f1 = os.path.join(OUT_DIR, "gene_heatmap.png")
    fig.savefig(f1, dpi=150)
    plt.close(fig)

    # ---- FIGÜR 2: potansiyel dağılımı ----
    counts = df["pha_potential"].value_counts()
    order = sorted(counts.index, key=lambda p: POT_ORDER.get(p, 5))
    vals = [counts[o] for o in order]
    colors = {"SCL": "#2e7d32", "MCL": "#1565c0", "SCL-co-MCL": "#6a1b9a",
              "belirsiz": "#f9a825", "none": "#9e9e9e"}
    fig2, ax2 = plt.subplots(figsize=(6.5, 4))
    bars = ax2.bar(order, vals, color=[colors.get(o, "#777") for o in order])
    ax2.bar_label(bars)
    ax2.set_ylabel("Suş sayısı")
    ax2.set_title("Halomonas PHA-tipi genomik potansiyeli dağılımı (N=%d)" % n, fontsize=10)
    ax2.set_xlabel("'belirsiz' = fonksiyonel PhaC var, monomer rotası tespit edilmedi (çekimser)")
    fig2.tight_layout()
    f2 = os.path.join(OUT_DIR, "potential_distribution.png")
    fig2.savefig(f2, dpi=150)
    plt.close(fig2)

    # Sıralı tablo (raporda kullanmak için)
    df.drop(columns=["pot_rank"]).to_csv(
        os.path.join(OUT_DIR, "results_sorted.csv"), index=False)

    print(f"\nFigürler: {f1}\n          {f2}")
    print(f"Sıralı tablo: {os.path.join(OUT_DIR, 'results_sorted.csv')}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    main()
