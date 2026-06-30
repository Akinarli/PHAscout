"""
pathway_engine PHBV/3HV regresyon testleri.

Ana koruma: delta yolagi phaA+phaB TEK BASINA P(3HB-co-3HV)'yi DOGRULANMIS bir
cikti olarak iddia ETMEMELI. 3HV monomeri propiyonil-CoA rotasi (C5-kabul eden
tiyolaz + tek-karbon-sayili VFA) gerektirir; PHAscout 3HV-onculu sinyal aramaz.
Bu yuzden delta aktif olsa bile ciktisi KOSULLU olarak isaretlenmeli ve
intrinsik egilim P(3HB) olarak ifade edilmeli.

Calistirma:  python -m pytest tests/test_pathway_phbv.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phascout.prediction.pathway_engine import PathwayEngine


def _gv(**kw):
    base = {g: False for g in ["phaC", "phaA", "phaB", "phaJ", "phaG", "phaP", "phaR", "phaE"]}
    base.update({k: bool(v) for k, v in kw.items()})
    return base


def _delta(pathways):
    return next(p for p in pathways if p["pathway_id"] == "delta")


def test_delta_does_not_bare_assert_phbv():
    # Class I + phaA+phaB (3HV-onculu gen YOK) -> delta aktif olabilir AMA
    # ciktisi kosulsuz "P(3HB-co-3HV)" OLMAMALI.
    pw = PathwayEngine().determine_pathways(_gv(phaC=1, phaA=1, phaB=1), "Class_I")
    delta = _delta(pw)
    assert delta["active"] is True, "delta phaA+phaB ile aktif olmali"
    prod = delta["product_tendency"]
    # Kosulsuz PHBV iddiasi olmamali: ya kosul ifadesi icermeli ya da P(3HB) one cikmali
    assert "koşullu" in prod.lower() or "yalnızca" in prod.lower() or "P(3HB)" in prod, prod
    assert delta.get("confidence") == "CONDITIONAL", delta.get("confidence")


def test_delta_note_warns_about_3hv_precursor():
    pw = PathwayEngine().determine_pathways(_gv(phaC=1, phaA=1, phaB=1), "Class_I")
    note = _delta(pw)["note"].lower()
    assert "3hv" in note and ("vfa" in note or "propiyonat" in note), note


def test_alpha_still_p3hb():
    # alpha (sekerden SCL) degismemeli: P(3HB)
    pw = PathwayEngine().determine_pathways(_gv(phaC=1, phaA=1, phaB=1), "Class_I")
    alpha = next(p for p in pw if p["pathway_id"] == "alpha")
    assert alpha["active"] is True
    assert "P(3HB)" in alpha["product_tendency"], alpha["product_tendency"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"  [OK ] {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  [HATA] {fn.__name__}: {e}")
    print(f"\n{passed}/{len(fns)} gecti")
