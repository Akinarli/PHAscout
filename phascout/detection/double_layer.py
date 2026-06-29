"""
PHAscout Double-Layer Filtre Modülü (Katman 2: Spesifiklik)
=============================================================
PFAM taramasından gelen adayları, BLOSUM62 normalize hizalama skoru
ile filtreler. Bu katman, FabG↔PhaB ve FadA↔PhaA gibi false positive'leri
elemeye yarar.

Modern Bio.Align.PairwiseAligner kullanır (deprecated pairwise2 değil).

Kullanım:
    from phascout.detection.double_layer import DoubleLayerFilter
    dl = DoubleLayerFilter()
    verified = dl.filter_candidates(hmm_hits)
"""

import os
import logging
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

from phascout.config import (
    BLOSUM62_THRESHOLDS,
    LENGTH_FILTERS,
    POSITIVE_REF_DIR,
)

logger = logging.getLogger(__name__)


class DoubleLayerFilter:
    """
    BLOSUM62 tabanlı spesifiklik filtresi.

    Katman 1 (HMM) adaylarını referanslara karşı hizalayarak,
    normalize skor eşiğini geçen genleri onaylar.
    """

    def __init__(self):
        self.aligner = PairwiseAligner()
        self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.mode = "local"  # Multi-domain proteinler için local alignment

        # Referans dizileri gen bazlı yükle
        self.references = {}
        self._load_references()

    def _load_references(self):
        """
        data/reference_sequences/positive/ klasöründen
        her gen için referans FASTA dizilerini yükle.
        """
        ref_files = {
            "phaB": "phab.fasta",
            "phaA": "phaa.fasta",
            "phaG": "phag.fasta",
            "phaJ": "phaj.fasta",
        }

        for gene, filename in ref_files.items():
            filepath = os.path.join(POSITIVE_REF_DIR, filename)
            if os.path.exists(filepath):
                refs = list(SeqIO.parse(filepath, "fasta"))
                self.references[gene] = refs
                logger.info(f"  {gene}: {len(refs)} referans yüklendi.")
            else:
                raise FileNotFoundError(
                    f"KRİTİK HATA: {gene} geni için Double-Layer filtreleme referans dosyası "
                    f"bulunamadı: {filepath}. PHAscout'un çalışması için 'data/reference_sequences/positive/' "
                    f"klasöründe referans FASTA dosyaları eksiksiz olmalıdır."
                )

    def _normalized_score(self, query_seq: str, ref_seq: str) -> float:
        """
        Normalize BLOSUM62 hizalama skoru hesapla.

        Formül: Score(query, ref) / Score(ref, ref)
        Sonuç 0-1 arasında olur.
        """
        try:
            score_qr = self.aligner.score(query_seq, ref_seq)
            score_rr = self.aligner.score(ref_seq, ref_seq)
            if score_rr == 0:
                return 0.0
            return score_qr / score_rr
        except Exception as e:
            logger.debug(f"Hizalama hatası: {e}")
            return 0.0

    def _max_score_against_refs(self, query_seq: str, gene: str) -> float:
        """
        Bir aday dizisinin, ilgili gen referanslarına karşı
        aldığı en yüksek normalize skoru döndür.
        """
        refs = self.references.get(gene, [])
        if not refs:
            return 0.0

        max_score = 0.0
        for ref in refs:
            score = self._normalized_score(query_seq, str(ref.seq))
            if score > max_score:
                max_score = score
        return max_score

    def _check_length(self, seq_length: int, gene: str) -> bool:
        """
        Dizi uzunluğunun beklenen aralıkta olup olmadığını kontrol et.
        Aralık dışındaysa tolerans payı (%20) ile tekrar dene.
        """
        if gene not in LENGTH_FILTERS:
            return True  # Uzunluk filtresi tanımlı değilse geç

        min_len, max_len = LENGTH_FILTERS[gene]
        tolerance = 0.20  # %20 tolerans

        lower = int(min_len * (1 - tolerance))
        upper = int(max_len * (1 + tolerance))

        return lower <= seq_length <= upper

    def filter_candidates(self, hmm_hits: dict) -> dict:
        """
        HMM taramasından gelen adayları BLOSUM62 + uzunluk filtresiyle ele.

        Bu metot sadece Double-Layer gerektiren genleri (phaA, phaB, phaG, phaJ)
        filtreler. phaC, phaP, phaR, phaE için Katman 2 gerekmez.

        Args:
            hmm_hits: HMMScanner.scan() çıktısı.
                      dict: gen_adı -> hit listesi

        Returns:
            dict: gen_adı -> onaylanmış hit listesi (+ 'blosum_score' ve 'verified' eklenir)
        """
        verified_results = {}

        # Double-Layer gerektirmeyen genler doğrudan geçer
        passthrough_genes = ["phaC", "phaP", "phaR", "phaE"]

        for gene, hits in hmm_hits.items():
            if gene in passthrough_genes:
                # Filtre yok, direkt onayla
                for hit in hits:
                    hit["verified"] = True
                    hit["blosum_score"] = None
                    hit["filter_note"] = "PFAM-only (Double-Layer gerekmez)"
                verified_results[gene] = hits
                continue

            # Double-Layer filtre uygula
            threshold = BLOSUM62_THRESHOLDS.get(gene, 0.35)
            verified = []

            for hit in hits:
                seq = hit.get("sequence", "")
                seq_len = hit.get("seq_length", len(seq))

                # Uzunluk kontrolü
                if not self._check_length(seq_len, gene):
                    hit["verified"] = False
                    hit["blosum_score"] = None
                    hit["filter_note"] = (
                        f"Uzunluk filtresi BAŞARISIZ: {seq_len} aa "
                        f"(beklenen: {LENGTH_FILTERS.get(gene, 'N/A')})"
                    )
                    logger.debug(
                        f"  {gene}/{hit['protein_id']}: Uzunluk reddedildi ({seq_len} aa)"
                    )
                    continue

                # BLOSUM62 normalize skor
                blosum_score = self._max_score_against_refs(seq, gene)
                hit["blosum_score"] = round(blosum_score, 4)

                if blosum_score >= threshold:
                    hit["verified"] = True
                    hit["filter_note"] = (
                        f"BLOSUM62 ONAYLANDI: {blosum_score:.4f} >= {threshold}"
                    )
                    verified.append(hit)
                    logger.info(
                        f"  {gene}/{hit['protein_id']}: ONAY "
                        f"(BLOSUM={blosum_score:.4f}, eşik={threshold})"
                    )
                else:
                    hit["verified"] = False
                    hit["filter_note"] = (
                        f"BLOSUM62 REDDEDİLDİ: {blosum_score:.4f} < {threshold}"
                    )
                    logger.info(
                        f"  {gene}/{hit['protein_id']}: RED "
                        f"(BLOSUM={blosum_score:.4f} < {threshold})"
                    )

            verified_results[gene] = verified

        return verified_results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    dl = DoubleLayerFilter()
    print("DoubleLayerFilter başarıyla yüklendi.")
    print(f"Referanslar: { {k: len(v) for k, v in dl.references.items()} }")
