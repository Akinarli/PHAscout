"""
pathway_engine — PHBV/3HV + birleşik-karar (pha_type otoritesi) regresyon testleri.

İki ana koruma:
  1. PHBV: delta yolağı phaA+phaB TEK BAŞINA P(3HB-co-3HV)'yi DOĞRULANMIŞ çıktı
     olarak iddia etmemeli (CONDITIONAL + 3HV-öncül uyarısı).
  2. Birleşik karar: yolak motoru pha_type'a TABİDİR. pha_type 'belirsiz'/'none'
     derse veya PhaC fonksiyonel değilse HİÇBİR yolak aktif olmamalı; böylece
     rapor kendisiyle çelişmez.

Çalıştırma:  python -m pytest tests/test_pathway_phbv.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phascout.prediction.pathway_engine import PathwayEngine
from phascout.prediction.pha_type import classify_pha_potential


def _gv(**kw):
    base = {g: False for g in ["phaC", "phaA", "phaB", "phaJ", "phaG", "phaP", "phaR", "phaE"]}
    base.update({k: bool(v) for k, v in kw.items()})
    return base


def _phac(cls, functional=True):
    return {"best_class": cls, "is_functional": functional}


def _pathways(phac, gv):
    """Pipeline ile AYNI çağrı: pha_type otorite, sonra yolak motoru."""
    pot = classify_pha_potential(phac, gv, {})
    return PathwayEngine().determine_pathways(
        gv, phac.get("best_class"),
        functional=phac.get("is_functional", False),
        pha_potential=pot,
    ), pot


def _get(pathways, pid):
    return next(p for p in pathways if p["pathway_id"] == pid)


# ---- PHBV / delta ----

def test_delta_does_not_bare_assert_phbv():
    pw, pot = _pathways(_phac("Class_I"), _gv(phaC=1, phaA=1, phaB=1))
    assert pot["potential"] == "SCL"
    delta = _get(pw, "delta")
    assert delta["active"] is True
    prod = delta["product_tendency"]
    assert "koşullu" in prod.lower() or "yalnızca" in prod.lower() or "P(3HB)" in prod, prod
    assert delta.get("confidence") == "CONDITIONAL", delta.get("confidence")


def test_delta_note_warns_about_3hv_precursor():
    pw, _ = _pathways(_phac("Class_I"), _gv(phaC=1, phaA=1, phaB=1))
    note = _get(pw, "delta")["note"].lower()
    assert "3hv" in note and ("vfa" in note or "propiyonat" in note), note


def test_alpha_still_p3hb():
    pw, _ = _pathways(_phac("Class_I"), _gv(phaC=1, phaA=1, phaB=1))
    alpha = _get(pw, "alpha")
    assert alpha["active"] is True
    assert "P(3HB)" in alpha["product_tendency"], alpha["product_tendency"]


# ---- Birleşik karar: pha_type otoritesi ----

def test_class_ii_phac_only_no_active_pathway():
    # R. opacus deseni: Class II + sadece phaC -> pha_type 'belirsiz'.
    # ESKI HATA: beta (required=[phaC]) AKTIF -> MCL diyordu (çelişki).
    # Artık: pha_type belirsiz olduğu için HİÇBİR yolak aktif değil.
    pw, pot = _pathways(_phac("Class_II"), _gv(phaC=1))
    assert pot["potential"] == "belirsiz"
    assert all(not p["active"] for p in pw), [p["pathway_id"] for p in pw if p["active"]]
    assert any("cekimser" in p["note"].lower() or "çekimser" in p["note"].lower() for p in pw)


def test_nonfunctional_phac_no_active_pathway():
    # E. coli deseni: phaC HMM hit var ama fonksiyonel değil.
    # ESKI RİSK: yolak motoru phaC VARLIĞINA bakıp aktif diyebilirdi.
    pw, pot = _pathways(_phac("Class_I", functional=False), _gv(phaC=1, phaA=1, phaB=1))
    assert pot["potential"] == "none"
    assert all(not p["active"] for p in pw)
    assert any("fonksiyonel" in p["note"].lower() for p in pw)


def test_class_ii_mcl_activates_only_mcl():
    # Class II + phaG -> pha_type MCL. SCL yolakları (alpha/delta) aktif olmamalı.
    pw, pot = _pathways(_phac("Class_II"), _gv(phaC=1, phaG=1))
    assert pot["potential"] == "MCL"
    assert _get(pw, "gamma")["active"] is True
    assert _get(pw, "alpha")["active"] is False
    assert _get(pw, "delta")["active"] is False


def test_no_phac_no_active_pathway():
    pw, pot = _pathways(_phac(None, functional=False), _gv())
    assert pot["potential"] == "none"
    assert all(not p["active"] for p in pw)


def test_consistency_active_implies_positive_potential():
    # Genel değişmez: herhangi bir aktif yolak varsa, pha_type pozitif olmalı.
    for phac, gv in [
        (_phac("Class_I"), _gv(phaC=1, phaA=1, phaB=1)),
        (_phac("Class_II"), _gv(phaC=1, phaG=1)),
        (_phac("Class_II"), _gv(phaC=1)),             # belirsiz
        (_phac("Class_I", functional=False), _gv(phaC=1, phaA=1, phaB=1)),  # none
    ]:
        pw, pot = _pathways(phac, gv)
        if any(p["active"] for p in pw):
            assert pot["potential"] in ("SCL", "MCL", "SCL-co-MCL"), (gv, pot["potential"])


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
