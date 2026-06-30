"""
PHAscout Alt Birim Kontrolcüsü
================================
Class III PhaC için PhaE, Class IV PhaC için PhaR alt birimlerinin
genomda bulunup bulunmadığını kontrol eder.

Alt birim eksikse polimeraz indeksi 40'tan 15'e düşürülür (config.py).

Kullanım:
    from phascout.classification.subunit_checker import SubunitChecker
    checker = SubunitChecker()
    result = checker.check(phac_class, hmm_scan_results)
"""

import logging
from phascout.config import SUBUNIT_PENALTY

logger = logging.getLogger(__name__)


class SubunitChecker:
    """
    PhaC sınıfına göre gerekli alt birimlerin varlığını kontrol eder.

    - Class III -> PhaE (PF08333) gereklidir.
    - Class IV  -> PhaR (PF07879) gereklidir.
    - Class I/II -> Alt birim gerekmez.
    """

    # Her sınıf için zorunlu alt birim
    REQUIRED_SUBUNITS = {
        "Class_III": "phaE",
        "Class_IV": "phaR",
    }

    def check(self, phac_class: str, hmm_scan_results: dict) -> dict:
        """
        Alt birim kontrolü yap.

        Args:
            phac_class: PhaC sınıfı ("Class_I", "Class_II", "Class_III", "Class_IV")
            hmm_scan_results: HMMScanner.scan() çıktısı (gen -> hit listesi)

        Returns:
            dict: {
                'subunit_required': bool,
                'subunit_name': str or None,
                'subunit_found': bool,
                'polymerase_score': int,
                'note': str,
            }
        """
        result = {
            "subunit_required": False,
            "subunit_name": None,
            "subunit_found": False,
            "polymerase_score": SUBUNIT_PENALTY["full_score"],
            "note": "",
        }

        required = self.REQUIRED_SUBUNITS.get(phac_class)

        if required is None:
            # Class I veya Class II: alt birim gerekmez
            result["note"] = f"{phac_class}: Alt birim gerektirmiyor."
            logger.info(result["note"])
            return result

        result["subunit_required"] = True
        result["subunit_name"] = required

        # HMM tarama sonuçlarında alt birimi ara
        subunit_hits = hmm_scan_results.get(required, [])

        if subunit_hits and len(subunit_hits) > 0:
            result["subunit_found"] = True
            result["polymerase_score"] = SUBUNIT_PENALTY["full_score"]
            result["note"] = (
                f"{phac_class}: {required} alt birimi BULUNDU "
                f"({subunit_hits[0]['protein_id']}). Polimeraz skoru: {SUBUNIT_PENALTY['full_score']}"
            )
        else:
            result["subunit_found"] = False
            result["polymerase_score"] = SUBUNIT_PENALTY["missing_score"]
            result["note"] = (
                f"{phac_class}: {required} alt birimi BULUNAMADI. "
                f"Polimeraz skoru dusuruldu: {SUBUNIT_PENALTY['full_score']} -> {SUBUNIT_PENALTY['missing_score']}"
            )

        logger.info(result["note"])
        return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    checker = SubunitChecker()

    # Test: Class I (alt birim gerektirmez)
    r1 = checker.check("Class_I", {})
    print(f"Class_I: {r1['note']}")

    # Test: Class III (PhaE gerekli, ama yok)
    r2 = checker.check("Class_III", {})
    print(f"Class_III (eksik): {r2['note']}")

    # Test: Class III (PhaE var)
    r3 = checker.check("Class_III", {"phaE": [{"protein_id": "WP_001234"}]})
    print(f"Class_III (var): {r3['note']}")
