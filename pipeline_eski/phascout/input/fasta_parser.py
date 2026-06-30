"""
PHAscout FASTA Ayrıştırıcı Modülü
===================================
Kullanıcının yüklediği veya NCBI'dan indirilen FASTA dosyalarını
ayrıştırıp, pipeline'a hazır SeqRecord listesine dönüştüren modül.

İki modlu çalışır:
  1. Proteom modu: Tam proteom (yüzlerce/binlerce protein)
  2. Tekli protein modu: Tek bir protein dizisi
"""

import os
import io
import logging
from Bio import SeqIO

logger = logging.getLogger(__name__)


class FastaInput:
    """
    FASTA girdisini yöneten sınıf.

    Attributes:
        records (list[SeqRecord]): Ayrıştırılmış protein dizileri.
        mode (str): 'proteome' veya 'single_protein'.
        source (str): Girdinin kaynağı (dosya yolu veya 'text_input').
    """

    def __init__(self):
        self.records = []
        self.mode = None
        self.source = None

    @classmethod
    def from_file(cls, filepath: str) -> "FastaInput":
        """
        Disk üzerindeki bir FASTA dosyasından oku.

        Args:
            filepath: FASTA dosyasının mutlak yolu.

        Returns:
            FastaInput: Ayrıştırılmış girdi nesnesi.

        Raises:
            FileNotFoundError: Dosya bulunamadı.
            ValueError: Dosya geçerli bir FASTA formatında değil.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Dosya bulunamadı: {filepath}")

        instance = cls()
        instance.source = filepath

        try:
            raw_records = list(SeqIO.parse(filepath, "fasta"))
            records = []
            for rec in raw_records:
                if not rec.seq: continue
                seq_str = str(rec.seq).replace("*", "")
                x_count = seq_str.count("X")
                if len(seq_str) > 0 and (x_count / len(seq_str)) > 0.05:
                    logger.warning(f"{rec.id} cok fazla belirsiz rezidu iceriyor (%{x_count/len(seq_str)*100:.1f} X), atlandi.")
                    continue
                seq_str = seq_str.replace("X", "")
                if not seq_str: continue
                rec.seq = rec.seq.__class__(seq_str)
                records.append(rec)
            instance.records = records
        except Exception as e:
            raise ValueError(f"FASTA ayrıştırma hatası: {e}")

        if not instance.records:
            raise ValueError(
                f"Dosyada geçerli FASTA kaydı bulunamadı: {filepath}"
            )

        instance._determine_mode()
        logger.info(
            f"{len(instance.records)} dizi okundu ({instance.mode} modu), "
            f"kaynak: {filepath}"
        )
        return instance

    @classmethod
    def from_text(cls, fasta_text: str) -> "FastaInput":
        """
        Metin olarak yapıştırılan FASTA verisinden oku.
        (Streamlit arayüzünde text_area'dan gelen girdi için.)

        Args:
            fasta_text: FASTA formatında metin.

        Returns:
            FastaInput: Ayrıştırılmış girdi nesnesi.
        """
        if not fasta_text or not fasta_text.strip():
            raise ValueError("Boş FASTA girdisi.")

        instance = cls()
        instance.source = "text_input"

        fasta_io = io.StringIO(fasta_text.strip())
        records = []
        for rec in SeqIO.parse(fasta_io, "fasta"):
            if not rec.seq:
                continue

            seq_str = str(rec.seq).replace("*", "")
            
            # X karakteri kontrolu
            x_count = seq_str.count("X")
            if len(seq_str) > 0 and (x_count / len(seq_str)) > 0.05:
                logger.warning(f"{rec.id} cok fazla belirsiz rezidu iceriyor (%{x_count/len(seq_str)*100:.1f} X), atlandi.")
                continue
                
            seq_str = seq_str.replace("X", "")
            if not seq_str:
                continue
            
            rec.seq = rec.seq.__class__(seq_str)
            records.append(rec)
        instance.records = records

        if not instance.records:
            raise ValueError(
                "Girdi geçerli bir FASTA formatında değil. "
                "Her dizi '>' ile başlamalıdır."
            )

        instance._determine_mode()
        logger.info(
            f"{len(instance.records)} dizi okundu ({instance.mode} modu), "
            f"kaynak: text_input"
        )
        return instance

    @classmethod
    def from_records(cls, records: list, source: str = "ncbi_api") -> "FastaInput":
        """
        Biopython SeqRecord listesinden oluştur.
        (ncbi_datasets.py'den gelen kayıtlar için.)

        Args:
            records: SeqRecord listesi.
            source: Kaynak tanımı.

        Returns:
            FastaInput: Ayrıştırılmış girdi nesnesi.
        """
        if not records:
            raise ValueError("Boş SeqRecord listesi.")

        instance = cls()
        instance.records = records
        instance.source = source
        instance._determine_mode()
        logger.info(
            f"{len(instance.records)} dizi yüklendi ({instance.mode} modu), "
            f"kaynak: {source}"
        )
        return instance

    def _determine_mode(self):
        """
        Girdi boyutuna göre çalışma modunu belirle.
        - 1 dizi → single_protein (sadece tespit + sınıflandırma, skor yok)
        - 2+ dizi → proteome (tam pipeline)
        """
        if len(self.records) == 1:
            self.mode = "single_protein"
        else:
            self.mode = "proteome"

    def get_summary(self) -> dict:
        """
        Girdi hakkında özet bilgi döndür.

        Returns:
            dict: Özet istatistikler.
        """
        if not self.records:
            return {"count": 0, "mode": None}

        lengths = [len(rec.seq) for rec in self.records]
        return {
            "count": len(self.records),
            "mode": self.mode,
            "source": self.source,
            "avg_length": sum(lengths) / len(lengths),
            "min_length": min(lengths),
            "max_length": max(lengths),
        }

    def __len__(self):
        return len(self.records)

    def __repr__(self):
        return (
            f"FastaInput(records={len(self.records)}, "
            f"mode='{self.mode}', source='{self.source}')"
        )


if __name__ == "__main__":
    # Hızlı test
    import sys

    if len(sys.argv) > 1:
        fi = FastaInput.from_file(sys.argv[1])
    else:
        test_fasta = (
            ">test_protein_1\n"
            "MATGKGAAASTQEGKSQPFKVTPGPFDPATWLEWSRQWQGTEGNGHAAASGIPGLDALAG\n"
            ">test_protein_2\n"
            "MFPIDIRPDKLTQEMLDYSRKLGQGMENLLNAEAIDTGVSPKQAVYSEDKLVLYRYDRPE\n"
        )
        fi = FastaInput.from_text(test_fasta)

    print(fi)
    print(fi.get_summary())
