"""
PHAscout Yolak Motoru (Pathway Engine)
========================================
Tespit edilen genlere ve PhaC sınıfına göre hangi metabolik yolakların
aktif olduğunu Boolean mantığıyla belirler.

4 Yolak:
  alpha: Sekerden SCL-PHA (phaC + phaA + phaB)
  beta:  Yag asidinden MCL-PHA (phaC Class II)
  gamma: Sekerden MCL-PHA (phaC Class II + phaG)
  delta: Seker + VFA ko-substrat (phaC + phaA + phaB + phaJ)

Kullanım:
    from phascout.prediction.pathway_engine import PathwayEngine
    engine = PathwayEngine()
    active = engine.determine_pathways(gene_vector, phac_class)
"""

import logging
from phascout.config import PATHWAYS

logger = logging.getLogger(__name__)


class PathwayEngine:
    """
    Boolean mantik ile metabolik yolak aktivasyonu belirler.
    """

    def determine_pathways(self, gene_vector: dict, phac_class: str) -> list:
        """
        Aktif yolakları belirle.

        Args:
            gene_vector: Gen tespit sonuclari.
                         dict: gen_adi -> bool (True = tespit edildi)
                         Ornek: {"phaC": True, "phaA": True, "phaB": True, ...}
            phac_class: PhaC sinifi ("Class_I", "Class_II", "Class_III", "Class_IV")
                        veya None (PhaC bulunamadiysa)

        Returns:
            list[dict]: Aktif yolak listesi. Her eleman:
                {
                    'pathway_id': str,
                    'name': str,
                    'active': bool,
                    'carbon_sources': list[str],
                    'product_tendency': str,
                    'missing_genes': list[str],
                    'note': str,
                }
        """
        results = []

        if not gene_vector.get("phaC", False):
            logger.info("PhaC bulunamadi, hicbir yolak aktif degil.")
            return [{
                "pathway_id": pathway_id,
                "name": pdef["name"],
                "active": False,
                "carbon_sources": [],
                "product_tendency": "Uretim beklenmez",
                "missing_genes": ["phaC"],
                "note": "PhaC bulunamadi.",
            } for pathway_id, pdef in PATHWAYS.items()]

        for pathway_id, pdef in PATHWAYS.items():
            # PhaC sinif kontrolu
            valid_classes = pdef.get("valid_phac_classes", [])
            if phac_class and phac_class not in valid_classes:
                results.append({
                    "pathway_id": pathway_id,
                    "name": pdef["name"],
                    "active": False,
                    "carbon_sources": [],
                    "product_tendency": pdef["product_tendency"],
                    "missing_genes": [],
                    "note": f"PhaC sinifi ({phac_class}) bu yolak icin uygun degil.",
                })
                continue

            # Zorunlu genlerin kontrolu (Boolean AND)
            required = pdef["required_genes"]
            missing = [g for g in required if not gene_vector.get(g, False)]

            if missing:
                results.append({
                    "pathway_id": pathway_id,
                    "name": pdef["name"],
                    "active": False,
                    "carbon_sources": [],
                    "product_tendency": pdef["product_tendency"],
                    "missing_genes": missing,
                    "note": f"Eksik genler: {', '.join(missing)}",
                })
            else:
                # Yolak aktif!
                optional = pdef.get("optional_genes", [])
                optional_found = [g for g in optional if gene_vector.get(g, False)]

                note = "Yolak AKTIF."
                confidence = "HIGH"

                if pathway_id == "beta":
                    phaj_present = gene_vector.get("phaJ", False)
                    confidence = "HIGH" if phaj_present else "MEDIUM"
                    if not phaj_present:
                        note += " phaJ bulunamadi — yag asidi substrat saglama kapasitesi belirsiz. Deneysel dogrulama onerilir."
                elif pdef.get("requires_3hv_precursor"):
                    # PHBV (3HV) over-claim'i onle: phaA+phaB TEK BASINA 3HV
                    # SAGLAMAZ. 3HV-onculu gen (or. bktB) taranmadigindan, bu
                    # yolak P(3HB-co-3HV)'yi DOGRULANMIS bir cikti olarak iddia
                    # edemez; yalnizca disaridan VFA ko-substrat altinda KOSULLU
                    # bir olasilik olarak raporlanir.
                    confidence = "CONDITIONAL"
                    note += (
                        " UYARI: P(3HB-co-3HV) potansiyeli KOSULLUDUR — 3HV monomeri "
                        "tek-karbon-sayili VFA (propiyonat/valerat) ko-substrat + C5-kabul "
                        "eden tiyolaz gerektirir. PHAscout 3HV-onculu sinyal (or. bktB) "
                        "ARAMAZ; phaA/phaB varligi PHBV'yi kanitlamaz. Bu organizmanin "
                        "intrinsik cikti egilimi P(3HB)'dir."
                    )
                    if optional_found:
                        note += f" Opsiyonel genler de mevcut: {', '.join(optional_found)}"
                else:
                    if optional_found:
                        note += f" Opsiyonel genler de mevcut: {', '.join(optional_found)}"

                results.append({
                    "pathway_id": pathway_id,
                    "name": pdef["name"],
                    "active": True,
                    "confidence": confidence,
                    "carbon_sources": pdef["carbon_sources"],
                    "product_tendency": pdef["product_tendency"],
                    "missing_genes": [],
                    "note": note,
                })

                logger.info(
                    f"Yolak {pathway_id} ({pdef['name']}): AKTIF -> {pdef['product_tendency']}"
                )

        return results

    def get_carbon_recommendations(self, active_pathways: list, user_carbon: str = None) -> list:
        """
        Aktif yolaklara gore karbon kaynagi onerileri sun.

        Args:
            active_pathways: determine_pathways() ciktisi.
            user_carbon: Kullanicinin sectigi karbon kaynagi (opsiyonel).

        Returns:
            list[dict]: Oneri listesi.
        """
        recommendations = []

        for pw in active_pathways:
            if not pw["active"]:
                continue

            if user_carbon:
                # Kullanici karbon kaynagi belirlediyse uyumluluk kontrol et
                compatible = user_carbon.lower() in [c.lower() for c in pw["carbon_sources"]]
                recommendations.append({
                    "pathway": pw["name"],
                    "carbon_source": user_carbon,
                    "compatible": compatible,
                    "product": pw["product_tendency"],
                    "note": "Uyumlu" if compatible else "Bu karbon kaynagi bu yolakla uyumsuz.",
                })
            else:
                # Tum uyumlu karbon kaynaklarini oner
                for carbon in pw["carbon_sources"]:
                    recommendations.append({
                        "pathway": pw["name"],
                        "carbon_source": carbon,
                        "compatible": True,
                        "product": pw["product_tendency"],
                        "note": "Onerilen karbon kaynagi.",
                    })

        return recommendations


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    engine = PathwayEngine()

    # Test: Cupriavidus necator benzeri (Class I, phaC + phaA + phaB)
    gene_vec = {"phaC": True, "phaA": True, "phaB": True, "phaJ": False, "phaG": False, "phaP": True}
    result = engine.determine_pathways(gene_vec, "Class_I")

    print("\n=== Yolak Analizi ===")
    for pw in result:
        status = "AKTIF" if pw["active"] else "INAKTIF"
        print(f"  {pw['pathway_id']} ({pw['name']}): {status}")
        if pw["active"]:
            print(f"    Urun: {pw['product_tendency']}")
            print(f"    Karbon: {', '.join(pw['carbon_sources'])}")
        elif pw["missing_genes"]:
            print(f"    Eksik: {', '.join(pw['missing_genes'])}")
