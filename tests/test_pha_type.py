"""
pha_type.classify_pha_potential regresyon testleri.

Ana koruma: phaJ TEK BAŞINA SCL-class (Class I/III/IV) bir organizmayı
SCL-co-MCL'e YÜKSELTMEMELI. 3HHx (mcl monomer) yalnızca MCL-yetenekli
(Class II) ya da geniş-substratlı sentazla polimerleşir. Bu kural ilk dürüst
benchmark'ta (N=15) ortaya çıkan sistematik over-claim'i giderir; biyokimyadan
türetilmiştir, benchmark vakalarına fit edilmemiştir.

Çalıştırma:  python -m pytest tests/test_pha_type.py   (veya doğrudan python tests/test_pha_type.py)
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phascout.prediction.pha_type import classify_pha_potential


def _gv(**kw):
    base = {g: False for g in ["phaC", "phaA", "phaB", "phaJ", "phaG", "phaP", "phaR", "phaE"]}
    base.update({k: bool(v) for k, v in kw.items()})
    return base


def _phac(cls):
    return {"best_class": cls, "is_functional": True}


def test_scl_class_with_phaj_does_not_overclaim_copolymer():
    # P. sacchari deseni: Class I + phaA+phaB+phaJ -> SCL (SCL-co-MCL DEGIL)
    r = classify_pha_potential(_phac("Class_I"), _gv(phaC=1, phaA=1, phaB=1, phaJ=1), {})
    assert r["potential"] == "SCL", r["potential"]


def test_class_iv_phab_phaj_is_scl():
    # B. thuringiensis deseni: Class IV + phaB+phaJ (phaA yok) -> SCL (phaJ besleme rotasi)
    r = classify_pha_potential(_phac("Class_IV"), _gv(phaC=1, phaB=1, phaJ=1), {})
    assert r["potential"] == "SCL", r["potential"]


def test_class_ii_phaj_still_mcl():
    # Class II + phaG/phaJ -> MCL (degismemeli)
    r = classify_pha_potential(_phac("Class_II"), _gv(phaC=1, phaG=1, phaJ=1), {})
    assert r["potential"] == "MCL", r["potential"]


def test_class_ii_phac_only_abstains():
    # R. opacus deseni: Class II + sadece phaC -> belirsiz (FP'yi onler, degismemeli)
    r = classify_pha_potential(_phac("Class_II"), _gv(phaC=1), {})
    assert r["potential"] == "belirsiz", r["potential"]


def test_scl_full_route_unchanged():
    r = classify_pha_potential(_phac("Class_I"), _gv(phaC=1, phaA=1, phaB=1), {})
    assert r["potential"] == "SCL", r["potential"]


def test_scl_phac_only_abstains():
    r = classify_pha_potential(_phac("Class_I"), _gv(phaC=1), {})
    assert r["potential"] == "belirsiz", r["potential"]


def test_nonfunctional_phac_no_potential():
    r = classify_pha_potential({"best_class": "Class_I", "is_functional": False}, _gv(phaC=1, phaA=1, phaB=1), {})
    assert r["potential"] == "none", r["potential"]


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
