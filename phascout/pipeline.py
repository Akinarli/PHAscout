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
from phascout.scoring.ml_scorer import MLScorer
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
        self.ml_scorer = MLScorer()
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
            try:
                from phascout.input.ncbi_datasets import fetch_proteome_and_gff
                records, gff_data = fetch_proteome_and_gff(accession)
                self.gff_data = gff_data
            except Exception as e:
                logger.warning(f"GFF3 alinamadi, sadece FASTA indirilecek. Hata: {e}")
                records = fetch_proteome(accession)
                self.gff_data = {}
            fasta_input = FastaInput.from_records(records, source=accession)
        elif fasta_file:
            # Ham (anote edilmemis) nukleotit genom/kontig girdisi otomatik
            # tespit edilir ve Prodigal ile gen cagrilir; protein FASTA'si ise
            # dogrudan okunur.
            from phascout.input.gene_caller import is_nucleotide_input, call_genes
            if is_nucleotide_input(path=fasta_file):
                logger.info("Nukleotit girdi tespit edildi -> Prodigal ile gen cagriliyor...")
                records, gff = call_genes(nucleotide_fasta_path=fasta_file)
                fasta_input = FastaInput.from_records(records, source=f"{fasta_file} (Prodigal)")
                self.gff_data = gff
            else:
                fasta_input = FastaInput.from_file(fasta_file)
                self.gff_data = {}
        elif fasta_text:
            from phascout.input.gene_caller import is_nucleotide_input, call_genes
            if is_nucleotide_input(text=fasta_text):
                logger.info("Nukleotit metin girdi tespit edildi -> Prodigal ile gen cagriliyor...")
                records, gff = call_genes(fasta_text=fasta_text)
                fasta_input = FastaInput.from_records(records, source="text_input (Prodigal)")
                self.gff_data = gff
            else:
                fasta_input = FastaInput.from_text(fasta_text)
                self.gff_data = {}
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

        phac_seq = ""
        best_phac_hit = None
        if gene_vector["phaC"]:
            phac_hits = filtered_results.get("phaC", [])
            # KRITIK: Katman 1 (genis PFAM agi) yuzlerce alpha/beta-hidrolazi
            # yakalar ve jenerik hidrolazlar gercek PhaC'ten DAHA IYI E-value
            # alabilir. Bu yuzden sadece phac_hits[0]'a bakmak gercek PhaC'i
            # kacirir (sessiz yanlis negatif). Tum adaylari dogrulayip,
            # ONAYLANAN (phac_confirmed) ve EN YUKSEK skorlu olani seceriz;
            # fonksiyonel (triad+box) olanlar onceliklidir.
            # Coklu-sentaz genomlarinda (or. C. necator H16'da Class I PhaC +
            # ikinci bir sentaz-homologu) ham bit-skoru farkli sinif HMM'leri
            # arasinda kiyaslanamaz. Bu yuzden secimde biyolojik baglami
            # kullaniriz: organizmanin GERCEKTEN kullanabilecegi, yani besleme
            # ROTASI tam olan sentazi tercih ederiz (Class I PhaC phaA/phaB
            # operonunda; yetim bir Class II ORF'una gore onceliklidir).
            def _route_complete(cls):
                if cls == "Class_II":
                    return gene_vector.get("phaG") or gene_vector.get("phaJ")
                return (gene_vector.get("phaA") and gene_vector.get("phaB")) or gene_vector.get("phaJ")

            best_rank = (-1, -1, -1.0)  # (is_functional, route_complete, best_score)
            n_candidates = len(phac_hits)
            for hit in phac_hits:
                seq = hit.get("sequence", "")
                if not seq:
                    continue
                cand = self.phac_validator.full_analysis(seq)
                if not cand.get("phac_confirmed"):
                    continue
                rank = (
                    1 if cand.get("is_functional") else 0,
                    1 if _route_complete(cand.get("best_class")) else 0,
                    cand.get("best_score", 0.0),
                )
                if rank > best_rank:
                    best_rank = rank
                    phac_result = cand
                    best_phac_hit = hit
                    phac_seq = seq

            if best_phac_hit is not None:
                # Secilen gercek PhaC'i listenin basina al ki rapor/operon
                # dogru proteini gostersin (jenerik hidrolazi degil).
                filtered_results["phaC"].remove(best_phac_hit)
                filtered_results["phaC"].insert(0, best_phac_hit)
                logger.info(
                    f"PhaC: {n_candidates} aday tarandi, secilen {best_phac_hit['protein_id']} "
                    f"-> sinif={phac_result['best_class']}, fonksiyonel={phac_result['is_functional']}"
                )
            else:
                logger.info(f"PhaC: {n_candidates} aday taranadi, hicbiri PhaC olarak onaylanmadi.")

        # =======================================
        # ADIM 5.5: OPERON / SINTENI KANITI (birinci-sinif kanit)
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 5.5: Operon/sinteni analizi...")

        from phascout.prediction.operon_analyzer import (
            analyze_operon_evidence, RESCUE_BLOSUM_FLOOR,
        )
        gff = self.gff_data if hasattr(self, "gff_data") else {}
        phac_pid = best_phac_hit["protein_id"] if best_phac_hit else None
        operon_result = analyze_operon_evidence(phac_pid, hmm_results, gff)

        # OPERON-KURTARMA: BLOSUM esigini gecememis ama secilen phaC ile SIKI
        # sintenik (operon icinde) bir phaA/phaB adayi, sinteni kanitiyla onaylanir.
        # Biyolojik gerekce: operon icindeki bir thiolaz/SDR neredeyse kesin
        # PhaA/PhaB'dir (FabG fatty-acid operonundadir, phaC'nin yaninda degil).
        # Operonu olmayan organizmalarda (Class IV Bacillus) hicbir sey kurtarilmaz.
        operon_result["rescued"] = {}
        for gene in ("phaA", "phaB"):
            gev = operon_result.get("genes", {}).get(gene, {})
            if gene_vector.get(gene) or not gev.get("syntenic"):
                continue
            cand_pid = gev.get("protein_id")
            cand_hit = next((h for h in hmm_results.get(gene, []) if h["protein_id"] == cand_pid), None)
            if not cand_hit:
                continue
            bl = self.double_layer._max_score_against_refs(cand_hit.get("sequence", ""), gene)
            if bl >= RESCUE_BLOSUM_FLOOR.get(gene, 0.4):
                cand_hit = dict(cand_hit)
                cand_hit["verified"] = True
                cand_hit["blosum_score"] = round(bl, 4)
                cand_hit["filter_note"] = (
                    f"OPERON-KURTARMA: BLOSUM {bl:.3f} (esik alti) ama phaC ile "
                    f"{gev['distance']} bp sintenik -> operon kaniti ile onaylandi."
                )
                filtered_results.setdefault(gene, []).insert(0, cand_hit)
                gene_vector[gene] = True
                operon_result["rescued"][gene] = {"protein_id": cand_pid, "distance": gev["distance"], "blosum": round(bl, 4)}
                logger.info(f"  {gene} OPERON-KURTARMA: {cand_pid} ({gev['distance']} bp, BLOSUM {bl:.3f})")

        # Sinteni SINIFLANDIRMAYI DEGISTIRMEZ; yalnizca kanit/etiket uretir.
        if operon_result.get("is_class_i_operon"):
            phac_result.setdefault("notes", []).append(
                "Sinteni kaniti: phaC, phaA/phaB ile operonik yakinlikta (PhaB↔FabG ayrimi icin destekleyici)."
            )

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
        # ADIM 8: MAKINE OGRENMESI (ML) SKORU
        # =======================================
        logger.info("=" * 50)
        logger.info("ADIM 8: ML Biyolojik Olasilik hesaplaniyor...")

        ml_result = self.ml_scorer.predict(
            phac_result, gene_vector, operon_result, phac_seq
        )

        # =======================================
        # ADIM 8.5: GUVEN + ALT BIRIM ENTEGRASYONU
        # =======================================
        # Karar (produces_pha) DETERMINISTIK kalir: fonksiyonel PhaC (triad+box)
        # + aktif yolak. ML ve alt birim, GUVENI ayarlar ve uyari uretir; ML
        # yapisal triad kanitini GECERSIZ KILMAZ (kara kutu, deterministik
        # katalitik kanittan daha zayif bir kanittir).
        if phac_result.get("phac_confirmed"):
            # Alt birim eksikse Class III/IV icin guven dusur + uyari
            if subunit_result.get("subunit_required") and not subunit_result.get("subunit_found"):
                phac_result["confidence"] = round(phac_result.get("confidence", 0.0) * 0.6, 1)
                warn = (
                    f"UYARI: {phac_class} sentazi {subunit_result.get('subunit_name')} "
                    f"alt birimi gerektirir ancak tespit edilmedi; aktif sentaz suphelidir."
                )
                phac_result.setdefault("notes", []).append(warn)
                phac_result["subunit_ok"] = False
            else:
                phac_result["subunit_ok"] = True

            # ML korroborasyonu (deterministik karari gecersiz kilmaz, sadece raporlanir)
            ml_p = ml_result.get("ml_probability", 0.0)
            phac_result["ml_corroborated"] = bool(ml_p >= 50.0)
            if phac_result.get("is_functional") and not phac_result["ml_corroborated"]:
                phac_result.setdefault("notes", []).append(
                    f"NOT: Fonksiyonel triad bulundu ancak ML guveni dusuk ({ml_p}%) - elle gozden gecirilebilir."
                )

        # =======================================
        # ADIM 8.7: PHA TIPI POTANSIYELI (yorum)
        # =======================================
        from phascout.prediction.pha_type import classify_pha_potential
        pha_potential = classify_pha_potential(phac_result, gene_vector, operon_result)
        logger.info(f"PHA tipi potansiyeli: {pha_potential['potential']} (guven: {pha_potential['confidence']})")

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
            heuristic_result=ml_result, # Use ml_result here to avoid rewriting reporter immediately
            carbon_recommendations=carbon_recs,
            pha_potential=pha_potential,
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
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    print(text)
