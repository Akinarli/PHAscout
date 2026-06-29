"""
PHAscout Rapor Ureticisi
=========================
Tum pipeline sonuclarini birlestirip JSON ve metin raporu ureten modul.

Kullanım:
    from phascout.reporting.report_generator import ReportGenerator
    report = ReportGenerator()
    output = report.generate(analysis_results)
"""

import json
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class ReportGenerator:
    """
    Pipeline sonuclarini birlestirip raporlayan sinif.
    """

    def generate(
        self,
        organism_info: dict,
        gene_vector: dict,
        hmm_details: dict,
        phac_result: dict,
        subunit_result: dict,
        pathway_results: list,
        heuristic_result: dict,
        carbon_recommendations: list = None,
    ) -> dict:
        """
        Tam rapor olustur.

        Returns:
            dict: JSON-uyumlu rapor.
        """
        # Tespit edilen genleri listele
        detected_genes = [g for g, found in gene_vector.items() if found]
        missing_genes = [g for g, found in gene_vector.items() if not found]

        # Aktif yolaklar
        active_pathways = [pw for pw in pathway_results if pw["active"]]

        # PHA uretim karari
        produces_pha = (
            phac_result.get("phac_confirmed", False) and len(active_pathways) > 0
        )

        report = {
            "meta": {
                "tool": "PHAscout",
                "version": "0.1.0",
                "timestamp": datetime.now().isoformat(),
            },
            "organism": organism_info,
            "summary": {
                "produces_pha": produces_pha,
                "pha_type": active_pathways[0]["product_tendency"] if active_pathways else "N/A",
                "phac_class": phac_result.get("best_class", "N/A"),
                "heuristic_score": heuristic_result.get("total_score", 0),
                "heuristic_max": heuristic_result.get("max_possible", 92),
                "heuristic_tier": heuristic_result.get("tier", "Potansiyelsiz"),
            },
            "genes": {
                "detected": detected_genes,
                "missing": missing_genes,
                "details": {},
            },
            "phac_analysis": {
                "class": phac_result.get("best_class"),
                "score": phac_result.get("best_score", 0),
                "confidence": phac_result.get("confidence", 0),
                "triad_found": phac_result.get("triad_found", False),
                "triad_residues": phac_result.get("triad_residues", {}),
                "box_found": phac_result.get("box_found", False),
                "box_match": phac_result.get("box_match"),
                "functional": phac_result.get("is_functional", False),
                "all_class_scores": phac_result.get("all_scores", {}),
                "notes": phac_result.get("notes", []),
            },
            "subunit_check": subunit_result,
            "pathways": pathway_results,
            "active_pathways": [pw["name"] for pw in active_pathways],
            "carbon_recommendations": carbon_recommendations or [],
            "heuristic_index": heuristic_result,
        }

        # Gen detaylari
        for gene, found in gene_vector.items():
            detail = {"detected": found}
            if gene in hmm_details and hmm_details[gene]:
                best = hmm_details[gene]
                if isinstance(best, list):
                    best = best[0] if best else {}
                detail["protein_id"] = best.get("protein_id", "N/A")
                detail["evalue"] = best.get("evalue", None)
                detail["score"] = best.get("score", None)
                detail["blosum_score"] = best.get("blosum_score", None)
                detail["filter_note"] = best.get("filter_note", "")
            report["genes"]["details"][gene] = detail

        return report

    def to_json(self, report: dict, filepath: str = None) -> str:
        """Raporu JSON string olarak dondur, opsiyonel olarak dosyaya kaydet."""
        json_str = json.dumps(report, indent=2, ensure_ascii=False, default=str)
        if filepath:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(json_str)
            logger.info(f"Rapor kaydedildi: {filepath}")
        return json_str

    def to_text(self, report: dict) -> str:
        """Raporu okunabilir metin formatina cevir."""
        lines = []
        lines.append("=" * 60)
        lines.append("PHAscout ANALIZ RAPORU")
        lines.append("=" * 60)

        # Organizma
        org = report.get("organism", {})
        lines.append(f"\nOrganizma: {org.get('organism_name', 'Bilinmiyor')}")
        lines.append(f"Accession: {org.get('accession', 'N/A')}")

        # Ozet
        s = report["summary"]
        lines.append(f"\n{'=' * 40}")
        lines.append("SONUC")
        lines.append(f"{'=' * 40}")

        if s["produces_pha"]:
            lines.append(f"PHA Uretimi: EVET")
            lines.append(f"PHA Tipi: {s['pha_type']}")
        else:
            lines.append(f"PHA Uretimi: HAYIR")

        lines.append(f"PhaC Sinifi: {s['phac_class']}")
        lines.append(f"Sezgisel Indeks: {s['heuristic_score']}/{s['heuristic_max']} ({s['heuristic_tier']})")

        # Genler
        g = report["genes"]
        lines.append(f"\n{'=' * 40}")
        lines.append("TESPIT EDILEN GENLER")
        lines.append(f"{'=' * 40}")
        for gene in g["detected"]:
            detail = g["details"].get(gene, {})
            pid = detail.get("protein_id", "")
            lines.append(f"  [+] {gene}: {pid}")
        for gene in g["missing"]:
            lines.append(f"  [-] {gene}: Bulunamadi")

        # PhaC
        pa = report["phac_analysis"]
        lines.append(f"\n{'=' * 40}")
        lines.append("PhaC ANALIZI")
        lines.append(f"{'=' * 40}")
        lines.append(f"  Sinif: {pa['class']} (skor: {pa['score']:.1f})")
        lines.append(f"  Guven: {pa['confidence']:.1f}")
        lines.append(f"  Triad: {'Bulundu' if pa['triad_found'] else 'Bulunamadi'}")
        lines.append(f"  Box: {'Bulundu' if pa['box_found'] else 'Bulunamadi'} ({pa['box_match']})")
        lines.append(f"  Fonksiyonel: {'Evet' if pa['functional'] else 'Hayir'}")

        # Yolaklar
        lines.append(f"\n{'=' * 40}")
        lines.append("METABOLIK YOLAKLAR")
        lines.append(f"{'=' * 40}")
        for pw in report["pathways"]:
            status = "AKTIF" if pw["active"] else "INAKTIF"
            lines.append(f"  {pw['pathway_id']} ({pw['name']}): {status}")
            if pw["active"]:
                lines.append(f"    Urun egilimi: {pw['product_tendency']}")

        # Skor detay
        hi = report["heuristic_index"]
        lines.append(f"\n{'=' * 40}")
        lines.append("SEZGISEL INDEKS DETAY")
        lines.append(f"{'=' * 40}")
        for line in hi.get("breakdown", []):
            lines.append(f"  {line}")

        lines.append(f"\n{'=' * 60}")
        lines.append(f"Rapor tarihi: {report['meta']['timestamp']}")
        lines.append(f"{'=' * 60}")

        return "\n".join(lines)


if __name__ == "__main__":
    gen = ReportGenerator()
    print("ReportGenerator basariyla yuklendi.")
