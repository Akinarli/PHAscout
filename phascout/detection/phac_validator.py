import re
import pyhmmer
from phascout.config import CATALYTIC_TRIAD, PHAC_BOX_REGEX, PHAC_CLASS_PROFILES

class PhaCValidator:
    """
    PhaC enziminin islevselligini HMM Domain Alignment koordinatlari uzerinden dogrular.
    Katalitik triad (Cys-Asp-His) ve PhaC Box motifi aranir.
    """
    def __init__(self, background=None):
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.background = background if background else pyhmmer.plan7.Background(self.alphabet)
        
        # HMM modellerini onceden yukle
        self.class_hmms = {}
        for cls, path in PHAC_CLASS_PROFILES.items():
            try:
                with pyhmmer.plan7.HMMFile(path) as f:
                    self.class_hmms[cls] = next(f)
            except Exception as e:
                pass # HMM henuz uretilmediyse atla
                
        # Referans HMM pozisyonlari (Merkez kolonlar)
        self.ref_positions = {
            "Class_I":   {"Cys": 319, "Asp": 447, "His": 477},
            "Class_II":  {"Cys": 296, "Asp": 424, "His": 453},
            "Class_III": {"Cys": 149, "Asp": 302, "His": 331},
            "Class_IV":  {"Cys": 149, "Asp": 302, "His": 331},
        }

    def validate_triad_hmm(self, sequence: str, phac_class: str) -> dict:
        result = {
            "triad_found": False,
            "triad_residues": {},
            "box_found": False,
            "box_match": None,
            "is_functional": False,
            "notes": []
        }
        
        hmm = self.class_hmms.get(phac_class)
        if not hmm:
            result["notes"].append("HATA: HMM yuklenemedi.")
            return result
            
        ref_cols = self.ref_positions.get(phac_class, {})
        seq_clean = sequence.replace("*", "").replace("X", "")
        
        try:
            digital_seq = pyhmmer.easel.TextSequence(
                name=b"query",
                sequence=seq_clean
            ).digitize(self.alphabet)
            
            hits = list(pyhmmer.hmmsearch([hmm], [digital_seq], background=self.background, E=1e10, T=0))
            if not hits or not hits[0] or not hits[0][0].domains:
                result["notes"].append("HMM domain alignment alinamadi.")
                return result
                
            domain = hits[0][0].domains[0]
            aln = domain.alignment
            
            found_residues = {}
            target_col = aln.target_from  # 1-indexed
            hmm_col = aln.hmm_from  # 1-indexed
            
            # Tolerans penceresi (HMM guncellemelerinde kolonlar biraz kayabilir)
            TOLERANCE = 25
            
            for hmm_char, target_char in zip(aln.hmm_sequence, aln.target_sequence):
                is_hmm_match = (hmm_char != '.')
                is_target_match = (target_char != '-')
                
                if is_target_match:
                    if is_hmm_match:
                        # Katalitik pencerede miyiz?
                        for res_name, ref_col in ref_cols.items():
                            if abs(hmm_col - ref_col) <= TOLERANCE:
                                if target_char.upper() == CATALYTIC_TRIAD[res_name]:
                                    # En yakin olani veya ilk buldugunu kaydet
                                    if res_name not in found_residues:
                                        found_residues[res_name] = {
                                            "hmm_col": hmm_col,
                                            "target_pos": target_col,
                                            "residue": target_char.upper()
                                        }
                    target_col += 1
                    
                if is_hmm_match:
                    hmm_col += 1
                    
            triad_ok = True
            for res_name, exp_aa in CATALYTIC_TRIAD.items():
                found = found_residues.get(res_name, {})
                actual = found.get("residue", "?")
                if actual != exp_aa:
                    triad_ok = False
                    result["notes"].append(f"{res_name}: bulunamadi (HMM tolerans penceresinde {exp_aa} yok).")
                else:
                    result["notes"].append(f"{res_name}: {exp_aa} dogrulandi (pozisyon {found.get('target_pos')}).")
                    
            result["triad_found"] = triad_ok
            result["triad_residues"] = found_residues
            
        except Exception as e:
            result["notes"].append(f"Hizalama hatasi: {e}")
            
        # PhaC Box kontrolu
        box_match = re.search(PHAC_BOX_REGEX, sequence)
        if box_match:
            result["box_found"] = True
            result["box_match"] = box_match.group()
        else:
            result["notes"].append("PhaC Box motifi bulunamadi.")
            
        result["is_functional"] = result["triad_found"] and result["box_found"]
        return result

    def full_analysis(self, sequence: str) -> dict:
        result = {
            "phac_confirmed": False,
            "best_class": "Unknown",
            "best_score": 0.0,
            "confidence": 0.0,
            "triad_found": False,
            "box_found": False,
            "is_functional": False,
            "notes": [],
            "all_scores": {}
        }
        
        seq_clean = sequence.replace("*", "").replace("X", "A")
        
        try:
            digital_seq = pyhmmer.easel.TextSequence(
                name=b"query", sequence=seq_clean
            ).digitize(self.alphabet)
        except Exception:
            return result
            
        # 1. HMM tabanli siniflandirma
        best_score = 0
        best_class = "Unknown"
        
        for cls, hmm in self.class_hmms.items():
            try:
                hits = list(pyhmmer.hmmsearch([hmm], [digital_seq], background=self.background))
                if hits and len(hits[0]) > 0:
                    score = hits[0][0].score
                    result["all_scores"][cls] = score
                    if score > best_score and score > 20: # Minimum cutoff
                        best_score = score
                        best_class = cls
            except Exception:
                pass
                
        result["best_class"] = best_class
        result["best_score"] = best_score
        
        if best_class == "Unknown":
            result["notes"].append("Hicbir PhaC sinifina uygun degil.")
            return result
            
        result["phac_confirmed"] = True
        
        # 2. Triad dogrulama
        triad_res = self.validate_triad_hmm(sequence, best_class)
        result["triad_found"] = triad_res.get("triad_found", False)
        result["triad_residues"] = triad_res.get("triad_residues", {})
        result["box_found"] = triad_res.get("box_found", False)
        result["box_match"] = triad_res.get("box_match")
        result["is_functional"] = triad_res.get("is_functional", False)
        result["notes"].extend(triad_res.get("notes", []))
        
        return result
