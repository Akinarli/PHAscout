"""
PHAscout HMM Tarayıcı Modülü (Katman 1: Geniş Ağ)
====================================================
pyhmmer kullanarak protein dizilerini PFAM HMM profillerine karşı tarar.
Gathering Threshold (GA) ve E-value filtrelemesi uygular.

Bu modül pipeline'ın ilk katmanıdır: Geniş ağ atarak tüm olası
PHA gen adaylarını yakalar. False positive'ler Katman 2'de (double_layer.py)
elenir.

Kullanım:
    from phascout.detection.hmm_scanner import HMMScanner
    scanner = HMMScanner()
    hits = scanner.scan(protein_records)
"""

import os
import logging
import urllib.request
import pyhmmer

from phascout.config import (
    PFAM_PROFILES,
    PFAM_HMM_DIR,
    HMM_EVALUE_THRESHOLD,
)

logger = logging.getLogger(__name__)


def download_pfam_hmm(pfam_id: str, output_dir: str) -> str:
    """
    InterPro API'den tek bir PFAM HMM profilini indirir.

    Args:
        pfam_id: PFAM accession (örn. "PF07167")
        output_dir: Kayıt klasörü.

    Returns:
        str: İndirilen HMM dosyasının yolu.
    """
    os.makedirs(output_dir, exist_ok=True)
    hmm_path = os.path.join(output_dir, f"{pfam_id}.hmm")

    if os.path.exists(hmm_path):
        logger.debug(f"PFAM profili zaten mevcut: {pfam_id}")
        return hmm_path

    url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{pfam_id}?annotation=hmm"
    logger.info(f"PFAM profili indiriliyor: {pfam_id}")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PHAscout"})
        with urllib.request.urlopen(req, timeout=30) as response:
            data = response.read()

            # InterPro API gzip döndürebilir
            if data[:2] == b'\x1f\x8b':
                import gzip
                data = gzip.decompress(data)

            with open(hmm_path, "wb") as f:
                f.write(data)

        logger.info(f"PFAM profili indirildi: {pfam_id} -> {hmm_path}")
        return hmm_path
    except Exception as e:
        logger.error(f"PFAM profili indirilemedi ({pfam_id}): {e}")
        return None


class HMMScanner:
    """
    pyhmmer tabanlı HMM tarayıcı.

    PFAM profillerini yükler ve protein dizilerini bu profillere karşı
    tarayarak aday genleri tespit eder.
    """

    def __init__(self):
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.hmm_profiles = {}  # gen_adı -> [HMM objeleri]
        self._load_profiles()

    def _load_profiles(self):
        """PFAM HMM profillerini diskten yükle veya indir."""
        os.makedirs(PFAM_HMM_DIR, exist_ok=True)

        for gene_name, pfam_ids in PFAM_PROFILES.items():
            self.hmm_profiles[gene_name] = []

            for pfam_id in pfam_ids:
                hmm_path = os.path.join(PFAM_HMM_DIR, f"{pfam_id}.hmm")

                # Dosya yoksa indir
                if not os.path.exists(hmm_path):
                    hmm_path = download_pfam_hmm(pfam_id, PFAM_HMM_DIR)
                    if hmm_path is None:
                        logger.warning(
                            f"{gene_name} için {pfam_id} profili yüklenemedi, atlanıyor."
                        )
                        continue

                # pyhmmer ile HMM yükle
                try:
                    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                        hmm = hmm_file.read()
                        self.hmm_profiles[gene_name].append(hmm)
                        logger.debug(f"HMM yüklendi: {gene_name} <- {pfam_id}")
                except Exception as e:
                    logger.warning(f"HMM okunamadı ({hmm_path}): {e}")

        loaded = sum(len(v) for v in self.hmm_profiles.values())
        logger.info(f"Toplam {loaded} HMM profili yüklendi ({len(self.hmm_profiles)} gen).")

    def scan(self, records: list) -> dict:
        """
        Protein dizilerini tüm PFAM profillerine karşı tara.

        Args:
            records: Biopython SeqRecord listesi.

        Returns:
            dict: Gen adı -> hit listesi.
            Her hit: {
                'protein_id': str,
                'gene': str,
                'pfam_id': str,
                'evalue': float,
                'score': float,  # bit score
                'description': str,
                'sequence': str,
                'seq_length': int,
            }
        """
        if not records:
            return {}

        # SeqRecord'ları pyhmmer DigitalSequence'a çevir
        digital_seqs = []
        seq_map = {}  # pyhmmer name -> orijinal SeqRecord

        for rec in records:
            seq_str = str(rec.seq).replace("*", "").replace("X", "A")
            try:
                ds = pyhmmer.easel.TextSequence(
                    name=rec.id.encode(),
                    sequence=seq_str,
                ).digitize(self.alphabet)
                digital_seqs.append(ds)
                seq_map[rec.id] = rec
            except Exception as e:
                logger.debug(f"Dizi dönüştürülemedi ({rec.id}): {e}")

        if not digital_seqs:
            logger.warning("Hiçbir dizi HMM taramasına uygun değil.")
            return {}

        logger.info(f"{len(digital_seqs)} dizi, HMM taramasına hazır.")

        # Her gen için tarama yap
        results = {}

        for gene_name, hmms in self.hmm_profiles.items():
            results[gene_name] = []

            for hmm in hmms:
                # pyhmmer.hmmsearch ile tarama
                try:
                    all_hits = list(pyhmmer.hmmsearch(
                        [hmm],
                        digital_seqs,
                        E=HMM_EVALUE_THRESHOLD,
                    ))

                    for top_hits in all_hits:
                        for hit in top_hits:
                            if hit.evalue <= HMM_EVALUE_THRESHOLD:
                                # pyhmmer version uyumu: str veya bytes olabilir
                                hit_name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
                                original_rec = seq_map.get(hit_name)
                                seq_str = str(original_rec.seq) if original_rec else ""

                                hmm_name = hmm.name.decode() if isinstance(hmm.name, bytes) else str(hmm.name) if hmm.name else "unknown"
                                hit_desc = hit.description.decode() if isinstance(hit.description, bytes) else str(hit.description) if hit.description else ""

                                hit_info = {
                                    "protein_id": hit_name,
                                    "gene": gene_name,
                                    "pfam_id": hmm_name,
                                    "evalue": hit.evalue,
                                    "score": hit.score,
                                    "description": hit_desc,
                                    "sequence": seq_str,
                                    "seq_length": len(seq_str),
                                }
                                results[gene_name].append(hit_info)

                except Exception as e:
                    logger.error(f"HMM tarama hatasi ({gene_name}): {e}")

            # En iyi hit'i en başa koy (en düşük E-value)
            results[gene_name].sort(key=lambda x: x["evalue"])

            if results[gene_name]:
                best = results[gene_name][0]
                logger.info(
                    f"  {gene_name}: {len(results[gene_name])} aday bulundu. "
                    f"En iyi: {best['protein_id']} (E={best['evalue']:.2e})"
                )
            else:
                logger.info(f"  {gene_name}: Aday bulunamadı.")

        return results

    def get_best_hits(self, scan_results: dict) -> dict:
        """
        Her gen için en iyi (en düşük E-value) hit'i döndür.

        Args:
            scan_results: scan() çıktısı.

        Returns:
            dict: gen_adı -> en iyi hit dict veya None.
        """
        best = {}
        for gene_name, hits in scan_results.items():
            if hits:
                best[gene_name] = hits[0]
            else:
                best[gene_name] = None
        return best


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    from phascout.input.ncbi_datasets import fetch_proteome

    acc = sys.argv[1] if len(sys.argv) > 1 else "GCF_000009285.1"
    print(f"Test: {acc}")

    records = fetch_proteome(acc)
    print(f"Proteom: {len(records)} protein")

    scanner = HMMScanner()
    results = scanner.scan(records)
    best = scanner.get_best_hits(results)

    print("\n=== En İyi Adaylar ===")
    for gene, hit in best.items():
        if hit:
            print(f"  {gene}: {hit['protein_id']} (E={hit['evalue']:.2e}, Score={hit['score']:.1f})")
        else:
            print(f"  {gene}: Bulunamadı")
