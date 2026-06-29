"""
PHAscout Ana Pipeline Orkestratoru
====================================
Tum modulleri sirasyla calistirip, sonuclari birlestiren ana modul.

Kullanım:
    from phascout.pipeline import PHAscoutPipeline
    pipeline = PHAscoutPipeline()
    report = pipeline.run(accession="GCF_000009285.1")
    # veya
    report = pipeline.run(fasta_file="proteins.faa")
"""

import logging
from phascout.input.ncbi_datasets import fetch_proteome, get_organism_info
from phascout.input.fasta_parser import FastaInput
from phascout.detection.hmm_scanner import HMMScanner
from phascout.detection.double_layer import DoubleLayerFilter
from phascout.detection.phac_validator import PhaCValidator
from phascout.classification.subunit_checker import SubunitChecker
from phascout.prediction.pathway_engine import PathwayEngine
from phascout.scoring.heuristic_index import HeuristicIndex
from phascout.reporting.report_generator import ReportGenerator

logger = logging.getLogger(__name__)


class PHAscoutPipeline:
    """
    PHAscout ana pipeline sinifi.
    Tum analiz adimlari tek bir run() cagrisiyla yururulur.
    """

    def __init__(self):
        logger.info("PHAscout Pipeline baslatiliyor...")
        self.hmm_scanner = HMMScanner()
        self.double_layer = DoubleLayerFilter()
        self.phac_validator = PhaCValidator()
        self.subunit_checker = SubunitChecker()
        self.pathway_engine = PathwayEngine()
        self.heuristic = HeuristicIndex()
        self.reporter = ReportGenerator()
        logger.info("Tum moduller yuklendi.")

    def run(
        self,
        accession: str = None,
        fasta_file: str = None,
        fasta_text: str = None,
        carbon_source: str = None,
    ) -> dict:
        """
        Tam PHAscout analizini calistir.

        Args:
            accession: NCBI Assembly accession (GCF_/GCA_).
            fasta_file: Yerel FASTA dosya yolu.
            fasta_text: FASTA metin girdisi.
            carbon_source: Kullanicinin belirledigi karbon kaynagi (opsiyonel).

        Returns:
            dict: Tam analiz raporu.
        """
        # =======================================
        # ADIM 1: GIRDI HAZIRLAMA
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 1: Girdi hazirlaniyor...")

        organism_info = {"organism_name": "Bilinmiyor", "accession": "N/A"}

        if accession:
            organism_info = get_organism_info(accession)
            logger.info(f"Organizma: {organism_info.get('organism_name')}")
            records = fetch_proteome(accession)
            fasta_input = FastaInput.from_records(records, source=accession)
        elif fasta_file:
            fasta_input = FastaInput.from_file(fasta_file)
        elif fasta_text:
            fasta_input = FastaInput.from_text(fasta_text)
        else:
            raise ValueError("accession, fasta_file veya fasta_text gereklidir.")

        logger.info(f"Girdi: {len(fasta_input)} protein ({fasta_input.mode} modu)")

        # =======================================
        # ADIM 2: HMM TARAMA (Katman 1 - Genis Ag)
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 2: HMM taramasi basliyor...")

        hmm_results = self.hmm_scanner.scan(fasta_input.records)

        # =======================================
        # ADIM 3: BLOSUM62 FILTRE (Katman 2 - Spesifiklik)
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 3: Double-Layer filtre uygulanıyor...")

        filtered_results = self.double_layer.filter_candidates(hmm_results)

        # =======================================
        # ADIM 4: GEN VEKTORU OLUSTUR
        # =======================================
        gene_vector = {}
        all_genes = ["phaC", "phaA", "phaB", "phaJ", "phaG", "phaP", "phaR", "phaE"]

        for gene in all_genes:
            hits = filtered_results.get(gene, [])
            gene_vector[gene] = len(hits) > 0

        logger.info(f"Gen vektoru: {gene_vector}")

        # =======================================
        # ADIM 5: PhaC SINIFLANDIRMA + DOGRULAMA
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 5: PhaC analizi...")

        phac_result = {
            "phac_confirmed": False,
            "best_class": None,
            "best_score": 0,
            "confidence": 0,
            "triad_found": False,
            "box_found": False,
            "is_functional": False,
            "notes": [],
            "all_scores": {},
        }

        if gene_vector["phaC"]:
            phac_hits = filtered_results.get("phaC", [])
            if phac_hits:
                best_phac = phac_hits[0]
                phac_seq = best_phac.get("sequence", "")
                if phac_seq:
                    phac_result = self.phac_validator.full_analysis(phac_seq)
                    logger.info(f"PhaC sinifi: {phac_result['best_class']}")
                    logger.info(f"PhaC fonksiyonel: {phac_result['is_functional']}")

        phac_class = phac_result.get("best_class")

        # =======================================
        # ADIM 6: ALT BIRIM KONTROLU
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 6: Alt birim kontrolu...")

        subunit_result = self.subunit_checker.check(phac_class, hmm_results)

        # =======================================
        # ADIM 7: YOLAK ANALIZI
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 7: Yolak analizi...")

        pathway_results = self.pathway_engine.determine_pathways(gene_vector, phac_class)
        carbon_recs = self.pathway_engine.get_carbon_recommendations(
            pathway_results, carbon_source
        )

        # =======================================
        # ADIM 8: SEZGISEL INDEKS
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 8: Sezgisel indeks hesaplaniyor...")

        heuristic_result = self.heuristic.calculate(
            phac_result, gene_vector, pathway_results, subunit_result
        )

        # =======================================
        # ADIM 9: RAPOR URETIMI
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 9: Rapor uretiliyor...")

        report = self.reporter.generate(
            organism_info=organism_info,
            gene_vector=gene_vector,
            hmm_details=filtered_results,
            phac_result=phac_result,
            subunit_result=subunit_result,
            pathway_results=pathway_results,
            heuristic_result=heuristic_result,
            carbon_recommendations=carbon_recs,
        )

        logger.info("=" * 50)
        logger.info("ANALIZ TAMAMLANDI.")
        return report


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    acc = sys.argv[1] if len(sys.argv) > 1 else "GCF_000009285.1"

    pipeline = PHAscoutPipeline()
    report = pipeline.run(accession=acc)

    # Metin raporu yazdir
    text = pipeline.reporter.to_text(report)
    print(text)
