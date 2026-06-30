"""
PHAscout Gen Çağırıcı (Prodigal)
=================================
Anote edilmemiş (ham) genom/kontig girdisi için gen çağırma katmanı.

PHAscout'un çekirdeği protein dizileri üzerinde çalışır. Kullanıcı anote
edilmemiş bir NÜKLEOTİT FASTA'sı (genom/kontig) verdiğinde, bu modül Prodigal
ile gen çağırıp protein dizilerini üretir ve aynı zamanda Prodigal'in protein
başlıklarına gömdüğü koordinatlardan (`# start # end # strand`) bir
gff_data (protein_id -> {contig, start, end, strand}) çıkarır. Böylece ham
genom girdisinde de operon/sinteni analizi çalışabilir.

Tasarım:
  - Nükleotit mi protein mi olduğu, dizilerin alfabesinden otomatik anlaşılır
    (ACGTUN oranı yüksekse nükleotit). Protein dizilerinde A/C/G/T/N de birer
    amino asittir, ama oranları ~%25'tir; eşik 0.9 ikisini net ayırır.
  - Kısa/parçalı girdilerde Prodigal eğitilemez; toplam uzunluk <100 kb ise
    `-p meta` (önceden eğitilmiş), aksi halde `-p single` kullanılır ve
    `single` başarısız olursa `meta`'ya düşülür.
"""

import os
import shutil
import logging
import tempfile
import subprocess
from Bio import SeqIO

logger = logging.getLogger(__name__)

# Nükleotit alfabesi (RNA dahil). Protein dizilerinde de bu harfler bulunur
# (A,C,G,T,N,U birer amino asit ya da belirsiz) ama oranları düşüktür.
_NUC_CHARS = set("ACGTUNacgtun")
_NUC_FRACTION_THRESHOLD = 0.9
_META_MODE_BP_CEILING = 100_000  # bu uzunluğun altında Prodigal -p single eğitilemez


def _records_from(path=None, text=None):
    if path:
        return list(SeqIO.parse(path, "fasta"))
    import io
    return list(SeqIO.parse(io.StringIO(text), "fasta"))


def is_nucleotide_input(path=None, text=None, sample=5):
    """İlk birkaç kaydın alfabesine bakarak girdinin nükleotit olup olmadığını döndür.

    Boş/okunamaz girdide False döner (protein varsayımı; mevcut akış bozulmaz).
    """
    try:
        records = _records_from(path=path, text=text)
    except Exception:
        return False
    checked = nucleotide = 0
    for rec in records[:sample]:
        s = str(rec.seq)
        if not s:
            continue
        frac = sum(c in _NUC_CHARS for c in s) / len(s)
        checked += 1
        if frac >= _NUC_FRACTION_THRESHOLD:
            nucleotide += 1
    return checked > 0 and nucleotide == checked


def _parse_prodigal_faa(faa_path):
    """Prodigal protein FASTA'sını oku; (records, gff_data) döndür.

    Prodigal başlığı:  >{seqid}_{geneidx} # {start} # {end} # {strand} # {attrs}
    strand: 1 -> '+', -1 -> '-'. contig = id'nin son '_geneidx' öncesi kısmı.
    """
    records = []
    gff_data = {}
    for rec in SeqIO.parse(faa_path, "fasta"):
        pid = rec.id
        parts = [p.strip() for p in rec.description.split("#")]
        if len(parts) >= 4:
            try:
                start = int(parts[1])
                end = int(parts[2])
                strand = "+" if parts[3] == "1" else "-"
                contig = pid.rsplit("_", 1)[0] if "_" in pid else pid
                gff_data[pid] = {
                    "contig": contig, "start": start, "end": end, "strand": strand,
                }
            except (ValueError, IndexError):
                pass
        rec.description = pid  # başlığı sadeleştir (downstream temiz id kullansın)
        records.append(rec)
    return records, gff_data


def call_genes(nucleotide_fasta_path=None, fasta_text=None):
    """Nükleotit girdiden Prodigal ile protein dizileri + koordinat çıkar.

    Returns:
        (records: list[SeqRecord], gff_data: dict[str, dict])

    Raises:
        RuntimeError: Prodigal PATH'te yoksa veya gen çağrılamazsa.
    """
    prodigal = shutil.which("prodigal")
    if not prodigal:
        raise RuntimeError(
            "Prodigal PATH'te bulunamadı. Ham (anote edilmemiş) genom girdisi "
            "için Prodigal gereklidir: `conda install -c bioconda prodigal`."
        )

    tmpdir = tempfile.mkdtemp(prefix="phascout_prodigal_")
    try:
        in_path = os.path.join(tmpdir, "input.fna")
        if nucleotide_fasta_path:
            shutil.copyfile(nucleotide_fasta_path, in_path)
        else:
            with open(in_path, "w", encoding="utf-8") as fh:
                fh.write(fasta_text or "")

        total_bp = sum(len(r.seq) for r in SeqIO.parse(in_path, "fasta"))
        faa_path = os.path.join(tmpdir, "proteins.faa")
        mode = "single" if total_bp >= _META_MODE_BP_CEILING else "meta"

        def _run(p_mode):
            return subprocess.run(
                [prodigal, "-i", in_path, "-a", faa_path, "-p", p_mode, "-q"],
                capture_output=True, text=True,
            )

        proc = _run(mode)
        if (proc.returncode != 0 or not os.path.exists(faa_path)) and mode == "single":
            logger.warning("Prodigal -p single başarısız; -p meta'ya düşülüyor.")
            proc = _run("meta")

        if proc.returncode != 0 or not os.path.exists(faa_path):
            raise RuntimeError(f"Prodigal gen çağrısı başarısız: {proc.stderr.strip()[:400]}")

        records, gff_data = _parse_prodigal_faa(faa_path)
        if not records:
            raise RuntimeError("Prodigal hiçbir gen çağıramadı (girdi çok kısa/bozuk olabilir).")

        logger.info(
            f"Prodigal ({mode}): {total_bp} bp girdiden {len(records)} protein çağrıldı."
        )
        return records, gff_data
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if len(sys.argv) > 1:
        path = sys.argv[1]
        print(f"Nükleotit mi? {is_nucleotide_input(path=path)}")
        recs, gff = call_genes(nucleotide_fasta_path=path)
        print(f"{len(recs)} protein, {len(gff)} koordinat kaydı")
        for r in recs[:3]:
            print(f"  {r.id}: {len(r.seq)} aa  coords={gff.get(r.id)}")
