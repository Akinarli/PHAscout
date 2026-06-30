"""
PHAscout NCBI Datasets Modülü
==============================
NCBI Assembly accession numarasından (GCF_/GCA_) protein FASTA dosyasını
indiren modül. NCBI Datasets v2 REST API kullanır.

Kullanım:
    from phascout.input.ncbi_datasets import fetch_proteome
    records = fetch_proteome("GCF_000009285.1")  # Cupriavidus necator H16
"""

import os
import io
import json
import logging
import zipfile
import requests
from Bio import SeqIO

from phascout.config import NCBI_DATASETS_BASE_URL, NCBI_API_KEY, NCBI_TIMEOUT

logger = logging.getLogger(__name__)


def fetch_proteome(accession: str, output_dir: str = None) -> list:
    """
    NCBI Datasets API v2 ile verilen assembly accession'dan
    protein FASTA dosyasını indirir ve SeqRecord listesi döndürür.

    Args:
        accession: NCBI Assembly accession (GCF_XXXXXXXXX.X veya GCA_XXXXXXXXX.X)
        output_dir: İndirilen FASTA'nın kaydedileceği klasör (opsiyonel).
                    None ise geçici klasör kullanılır.

    Returns:
        list[SeqRecord]: Biopython SeqRecord nesnelerinin listesi.

    Raises:
        ValueError: Geçersiz accession formatı.
        ConnectionError: NCBI'a bağlantı başarısız.
        FileNotFoundError: İndirilen pakette protein dosyası bulunamadı.
    """
    # Accession validasyonu
    accession = accession.strip()
    if not (accession.startswith("GCF_") or accession.startswith("GCA_")):
        raise ValueError(
            f"Geçersiz accession: '{accession}'. "
            f"GCF_ veya GCA_ ile başlamalıdır. Örnek: GCF_000009285.1"
        )

    logger.info(f"NCBI'dan proteom indiriliyor: {accession}")

    # NCBI Datasets v2 API - genome download endpoint
    url = (
        f"{NCBI_DATASETS_BASE_URL}/genome/accession/{accession}/download"
        f"?include_annotation_type=PROT_FASTA"
    )

    headers = {"Accept": "application/zip"}
    if NCBI_API_KEY:
        headers["api-key"] = NCBI_API_KEY

    try:
        response = requests.get(url, headers=headers, timeout=NCBI_TIMEOUT, stream=True)
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        if response.status_code == 404:
            raise FileNotFoundError(
                f"Accession '{accession}' NCBI'da bulunamadı. "
                f"Lütfen doğruluğunu kontrol edin."
            )
        raise ConnectionError(f"NCBI API hatası: {e}")
    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"NCBI'a bağlantı başarısız: {e}")

    # ZIP dosyasını aç ve protein FASTA'yı bul
    zip_data = io.BytesIO(response.content)
    records = []

    try:
        with zipfile.ZipFile(zip_data) as zf:
            # ZIP içindeki protein FASTA dosyasını ara
            protein_file = None
            for name in zf.namelist():
                if name.endswith("protein.faa"):
                    protein_file = name
                    break

            if protein_file is None:
                raise FileNotFoundError(
                    f"'{accession}' için protein FASTA dosyası bulunamadı. "
                    f"Bu assembly'de protein annotasyonu olmayabilir."
                )

            logger.info(f"Protein dosyası bulundu: {protein_file}")

            # FASTA dosyasını oku
            with zf.open(protein_file) as fasta_handle:
                fasta_text = fasta_handle.read().decode("utf-8")
                fasta_io = io.StringIO(fasta_text)
                records = list(SeqIO.parse(fasta_io, "fasta"))

            logger.info(f"{len(records)} protein dizisi okundu.")

            # İsteğe bağlı olarak diske kaydet
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                out_path = os.path.join(output_dir, f"{accession}_proteins.faa")
                with open(out_path, "w", encoding="utf-8") as f:
                    f.write(fasta_text)
                logger.info(f"Proteom kaydedildi: {out_path}")

    except zipfile.BadZipFile:
        raise ConnectionError(
            f"NCBI'dan gelen veri bozuk. Lütfen tekrar deneyin."
        )

    if not records:
        raise FileNotFoundError(
            f"'{accession}' için protein dizisi bulunamadı."
        )

    return records


def fetch_proteome_and_gff(accession: str) -> tuple:
    accession = accession.strip()
    if not (accession.startswith("GCF_") or accession.startswith("GCA_")):
        raise ValueError(f"Geçersiz accession: '{accession}'. GCF_ veya GCA_ ile başlamalıdır.")

    logger.info(f"NCBI'dan proteom ve GFF3 indiriliyor: {accession}")

    url = (
        f"{NCBI_DATASETS_BASE_URL}/genome/accession/{accession}/download"
        f"?include_annotation_type=PROT_FASTA,GENOME_GFF"
    )

    headers = {"Accept": "application/zip"}
    if NCBI_API_KEY:
        headers["api-key"] = NCBI_API_KEY

    try:
        response = requests.get(url, headers=headers, timeout=NCBI_TIMEOUT, stream=True)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"NCBI API hatası: {e}")

    zip_data = io.BytesIO(response.content)
    records = []
    gff_data = {}

    try:
        with zipfile.ZipFile(zip_data) as zf:
            protein_file = None
            gff_file = None
            for name in zf.namelist():
                if name.endswith("protein.faa"):
                    protein_file = name
                elif name.endswith("genomic.gff"):
                    gff_file = name

            if protein_file:
                with zf.open(protein_file) as fasta_handle:
                    fasta_text = fasta_handle.read().decode("utf-8")
                    records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
                    
            if gff_file:
                with zf.open(gff_file) as gff_handle:
                    for line in gff_handle:
                        line = line.decode("utf-8")
                        if line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 9 and parts[2] == "CDS":
                            contig = parts[0]
                            start = int(parts[3])
                            end = int(parts[4])
                            strand = parts[6]
                            attributes = parts[8]
                            
                            prot_id = None
                            for attr in attributes.split(";"):
                                if attr.startswith("protein_id="):
                                    prot_id = attr.split("=")[1]
                                    break
                                elif attr.startswith("Name=WP_") or attr.startswith("Name=NP_"):
                                    prot_id = attr.split("=")[1]
                                    
                            if prot_id:
                                gff_data[prot_id] = {
                                    "contig": contig,
                                    "start": start,
                                    "end": end,
                                    "strand": strand
                                }
    except zipfile.BadZipFile:
        raise ConnectionError("NCBI'dan gelen veri bozuk (Bad Zip).")

    return records, gff_data


def get_organism_info(accession: str) -> dict:
    """
    NCBI Datasets API v2 ile verilen assembly accession'dan
    organizma bilgilerini (isim, taksonomi ID vb.) çeker.

    Args:
        accession: NCBI Assembly accession.

    Returns:
        dict: Organizma bilgileri.
    """
    url = f"{NCBI_DATASETS_BASE_URL}/genome/accession/{accession}/dataset_report"

    headers = {"Accept": "application/json"}
    if NCBI_API_KEY:
        headers["api-key"] = NCBI_API_KEY

    try:
        response = requests.get(url, headers=headers, timeout=NCBI_TIMEOUT)
        response.raise_for_status()
        data = response.json()
    except Exception as e:
        logger.warning(f"Organizma bilgisi alınamadı: {e}")
        return {"organism_name": "Bilinmeyen", "taxon_id": None}

    # API yanıtından organizma bilgisini çıkar
    try:
        reports = data.get("reports", [])
        if reports:
            org = reports[0].get("organism", {})
            return {
                "organism_name": org.get("organism_name", "Bilinmeyen"),
                "taxon_id": org.get("tax_id", None),
                "assembly_level": reports[0].get("assembly_info", {}).get("assembly_level", "Bilinmiyor"),
                "accession": accession,
            }
    except (KeyError, IndexError):
        pass

    return {"organism_name": "Bilinmeyen", "taxon_id": None, "accession": accession}


if __name__ == "__main__":
    # Hızlı test
    import sys
    acc = sys.argv[1] if len(sys.argv) > 1 else "GCF_000009285.1"
    print(f"Test: {acc}")

    info = get_organism_info(acc)
    print(f"Organizma: {info}")

    records = fetch_proteome(acc)
    print(f"Toplam protein: {len(records)}")
    if records:
        print(f"İlk protein: {records[0].id} ({len(records[0].seq)} aa)")
