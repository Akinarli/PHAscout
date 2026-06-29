import unittest
import urllib.request
import json
from Bio import SeqIO
from io import StringIO
from phascout.detection.phac_validator import PhaCValidator

class TestPhaCValidator(unittest.TestCase):
    def setUp(self):
        self.validator = PhaCValidator()

    def get_uniprot_seq(self, acc):
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout'})
        try:
            resp = urllib.request.urlopen(req)
            data = resp.read().decode('utf-8')
            rec = next(SeqIO.parse(StringIO(data), "fasta"))
            return str(rec.seq)
        except Exception:
            return None

    def test_class_I_independent(self):
        # A0A1B1YJ44 (Halomonas bluephagenesis PhaC - independent)
        seq = self.get_uniprot_seq("A0A1B1YJ44")
        if not seq:
            self.skipTest("Could not fetch sequence from UniProt")
            
        result = self.validator.validate_triad_hmm(seq, "Class_I")
        self.assertTrue(result["triad_found"], f"Triad should be found. Notes: {result['notes']}")
        self.assertTrue(result["box_found"], "PhaC box should be found.")
        self.assertEqual(result["triad_residues"]["Cys"]["residue"], "C")
        self.assertEqual(result["triad_residues"]["Asp"]["residue"], "D")
        self.assertEqual(result["triad_residues"]["His"]["residue"], "H")

    def test_false_positive_rejection(self):
        # P0A953 (E. coli FabG - totally unrelated)
        seq = self.get_uniprot_seq("P0A953")
        if not seq:
            self.skipTest("Could not fetch sequence from UniProt")
            
        result = self.validator.validate_triad_hmm(seq, "Class_I")
        self.assertFalse(result["is_functional"], "FabG should not pass the PhaC functional test.")

if __name__ == '__main__':
    unittest.main()
