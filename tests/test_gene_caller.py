"""
gene_caller testleri: nükleotit tespiti, Prodigal başlık koordinat ayrıştırma
ve (Prodigal kuruluysa) uçtan uca gen çağırma.

Prodigal PATH'te yoksa uçtan uca test atlanır (birim testler yine koşar).

Çalıştırma:  python -m pytest tests/test_gene_caller.py
"""

import os
import sys
import shutil

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phascout.input.gene_caller import (
    is_nucleotide_input,
    _parse_prodigal_faa,
    call_genes,
)

PROTEIN_FASTA = (
    ">p1\nMATGKGAAASTQEGKSQPFKVTPGPFDPATWLEWSRQWQGTEGNGHAAASGIPGLDALAG\n"
    ">p2\nMFPIDIRPDKLTQEMLDYSRKLGQGMENLLNAEAIDTGVSPKQAVYSEDKLVLYRYDRPE\n"
)
NUCLEOTIDE_FASTA = (
    ">contig1\n" + "ATGCATGCATGCGGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGC\n" * 4
)


def test_protein_not_detected_as_nucleotide():
    assert is_nucleotide_input(text=PROTEIN_FASTA) is False


def test_nucleotide_detected():
    assert is_nucleotide_input(text=NUCLEOTIDE_FASTA) is True


def test_empty_input_is_not_nucleotide():
    assert is_nucleotide_input(text="") is False


def test_parse_prodigal_header_coords(tmp_path):
    # Prodigal başlık biçimi: >{seqid}_{idx} # start # end # strand # attrs
    faa = tmp_path / "p.faa"
    faa.write_text(
        ">NODE_1_2 # 150 # 900 # 1 # ID=1_2;partial=00\nMAAAGGG\n"
        ">NODE_1_3 # 1200 # 1800 # -1 # ID=1_3;partial=00\nMCCCDDD\n"
    )
    records, gff = _parse_prodigal_faa(str(faa))
    assert len(records) == 2
    assert gff["NODE_1_2"] == {"contig": "NODE_1", "start": 150, "end": 900, "strand": "+"}
    assert gff["NODE_1_3"]["strand"] == "-"
    # başlık sadeleşmeli (downstream temiz id kullanır)
    assert records[0].description == "NODE_1_2"


@pytest.mark.skipif(shutil.which("prodigal") is None, reason="Prodigal kurulu değil")
def test_call_genes_end_to_end(tmp_path):
    # Tek kodonlu geri-çeviri ile sentetik bir ORF; Prodigal meta modunda çağırmalı.
    codon = {"M": "ATG", "A": "GCG", "G": "GGC", "K": "AAA", "L": "CTG",
             "E": "GAA", "V": "GTG", "T": "ACC", "S": "AGC", "R": "CGC",
             "D": "GAC", "P": "CCG", "F": "TTC", "I": "ATC", "Q": "CAG"}
    aa = "MAKLEGVTSRDAKLEGVTSRDPFIQ" * 6  # ~150 aa, hiç stop kodonu yok
    orf = "".join(codon[c] for c in aa) + "TAA"
    spacer = "GCGCGCGCGCATCGATCGATCG" * 4
    contig = spacer + orf + spacer
    fna = tmp_path / "g.fna"
    fna.write_text(">c1\n" + contig + "\n")

    records, gff = call_genes(nucleotide_fasta_path=str(fna))
    assert records, "Prodigal en az bir gen çağırmalı"
    assert all(r.id in gff for r in records)
    r0 = records[0]
    assert gff[r0.id]["contig"] == "c1"
    assert gff[r0.id]["start"] >= 1 and gff[r0.id]["end"] > gff[r0.id]["start"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
