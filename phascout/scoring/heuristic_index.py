"""
PHAscout Sezgisel PHA Adaylik Indeksi
=======================================
Endustriyel potansiyeli siralamak icin lineer toplama indeksi.
Bu bir mutlak verim (g/L) tahmini DEGILDIR, aday suslari
kiyaslamak icin uzman degerlendiremesidir.

S_index = S_polymerase + S_monomer_supply + S_yield_booster
Aralik: [0, 92]

Kullanım:
    from phascout.scoring.heuristic_index import HeuristicIndex
    scorer = HeuristicIndex()
    score = scorer.calculate(phac_result, gene_vector, pathway_results, subunit_result)
"""

import logging
from phascout.config import HEURISTIC_INDEX

logger = logging.getLogger(__name__)


class HeuristicIndex:
    """
    Sezgisel PHA Adaylik Indeksi hesaplayicisi.
    3 katmanli puanlama: Polimeraz + Monomer Tedarik + Verim Artirici
    """

    def calculate(
        self,
        phac_result: dict,
        gene_vector: dict,
        pathway_results: list,
        subunit_result: dict,
    ) -> dict:
        """
        Tam skor hesaplama.

        Args:
            phac_result: PhaCValidator.full_analysis() ciktisi.
            gene_vector: Gen tespit sonuclari (gen_adi -> bool).
            pathway_results: PathwayEngine.determine_pathways() ciktisi.
            subunit_result: SubunitChecker.check() ciktisi.

        Returns:
            dict: {
                'total_score': int,
                'max_possible': int,
                'percentage': float,
                'tier': str,           # "Yuksek", "Orta", "Dusuk", "Potansiyelsiz"
                's1_polymerase': int,
                's2_monomer': int,
                's3_yield': int,
                'breakdown': list[str],
            }
        """
        s1 = 0  # Polimeraz skoru
        s2 = 0  # Monomer tedarik skoru
        s3 = 0  # Verim artirici
        breakdown = []

        # ============================================
        # KATMAN 1: POLIMERAZ (0-40)
        # ============================================
        if phac_result.get("phac_confirmed", False):
            s1 = subunit_result.get("polymerase_score", HEURISTIC_INDEX["polymerase_full"])
            if s1 == HEURISTIC_INDEX["polymerase_full"]:
                breakdown.append(f"S1 Polimeraz: +{s1} (PhaC dogrulanmis, alt birim tamam)")
            else:
                breakdown.append(f"S1 Polimeraz: +{s1} (PhaC dogrulanmis, alt birim EKSIK)")
        else:
            s1 = HEURISTIC_INDEX["polymerase_none"]
            breakdown.append(f"S1 Polimeraz: +{s1} (PhaC dogrulanamadi)")

        # ============================================
        # KATMAN 2: MONOMER TEDARIK (0-40)
        # ============================================
        active_pathways = [pw for pw in pathway_results if pw["active"]]

        if active_pathways:
            # En yuksek puanli yolagi sec
            best_pathway = None
            best_monomer_score = 0

            for pw in active_pathways:
                pid = pw["pathway_id"]

                if pid == "alpha":
                    # phaA + phaB -> tam skor
                    if gene_vector.get("phaA") and gene_vector.get("phaB"):
                        score = HEURISTIC_INDEX["monomer_full"]
                    elif gene_vector.get("phaA") or gene_vector.get("phaB"):
                        score = HEURISTIC_INDEX["monomer_partial_phaA_only"]
                    else:
                        score = 0

                elif pid == "beta":
                    # Beta-oksidasyon yolagi: phaJ varsa tam, yoksa 30
                    if gene_vector.get("phaJ"):
                        score = HEURISTIC_INDEX["monomer_full"]
                    else:
                        score = HEURISTIC_INDEX["monomer_beta_no_phaJ"]

                elif pid == "gamma":
                    # phaG varsa tam
                    if gene_vector.get("phaG"):
                        score = HEURISTIC_INDEX["monomer_full"]
                    else:
                        score = 0

                elif pid == "delta":
                    # phaA + phaB + (phaJ opsiyonel)
                    if gene_vector.get("phaA") and gene_vector.get("phaB"):
                        score = HEURISTIC_INDEX["monomer_full"]
                    else:
                        score = HEURISTIC_INDEX["monomer_partial_phaA_only"]
                else:
                    score = 0

                if score > best_monomer_score:
                    best_monomer_score = score
                    best_pathway = pw

            s2 = best_monomer_score
            if best_pathway:
                breakdown.append(
                    f"S2 Monomer: +{s2} (Yolak: {best_pathway['name']}, "
                    f"Urun: {best_pathway['product_tendency']})"
                )
        else:
            s2 = 0
            breakdown.append("S2 Monomer: +0 (Aktif yolak yok)")

        # ============================================
        # KATMAN 3: VERIM ARTIRICI (0-12)
        # ============================================
        if gene_vector.get("phaP", False):
            s3 = HEURISTIC_INDEX["yield_phaP"]
            breakdown.append(f"S3 Verim: +{s3} (PhaP granul stabilizatoru tespit edildi)")
        else:
            s3 = 0
            breakdown.append("S3 Verim: +0 (PhaP bulunamadi)")

        # ============================================
        # TOPLAM
        # ============================================
        total = s1 + s2 + s3
        max_score = HEURISTIC_INDEX["max_score"]
        percentage = (total / max_score * 100) if max_score > 0 else 0

        # Kademe belirleme
        if percentage >= 80:
            tier = "Yuksek Potansiyel"
        elif percentage >= 50:
            tier = "Orta Potansiyel"
        elif percentage >= 20:
            tier = "Dusuk Potansiyel"
        else:
            tier = "Potansiyelsiz"

        result = {
            "total_score": total,
            "max_possible": max_score,
            "percentage": round(percentage, 1),
            "tier": tier,
            "s1_polymerase": s1,
            "s2_monomer": s2,
            "s3_yield": s3,
            "breakdown": breakdown,
        }

        logger.info("Sezgisel Indeks: %d/%d (%.1f%%) -> %s", total, max_score, percentage, tier)
        return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    scorer = HeuristicIndex()

    # Test: Cupriavidus necator benzeri (tam donanim)
    phac = {"phac_confirmed": True}
    genes = {"phaC": True, "phaA": True, "phaB": True, "phaP": True, "phaJ": False, "phaG": False}
    pathways = [
        {"pathway_id": "alpha", "active": True, "name": "Sekerden SCL-PHA", "product_tendency": "P(3HB)"},
        {"pathway_id": "beta", "active": False, "name": "Yag asidinden MCL-PHA", "product_tendency": "MCL-PHA"},
    ]
    subunit = {"polymerase_score": 40}

    result = scorer.calculate(phac, genes, pathways, subunit)

    print(f"\nToplam: {result['total_score']}/{result['max_possible']} ({result['percentage']}%)")
    print(f"Kademe: {result['tier']}")
    for line in result["breakdown"]:
        print(f"  {line}")
