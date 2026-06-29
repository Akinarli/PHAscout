import re
import pyhmmer
from phascout.config import (
    CATALYTIC_TRIAD,
    CATALYTIC_HMM_COLUMNS,
    TRIAD_TOLERANCE,
    PHAC_BOX_REGEX,
    PHAC_CLASS_PROFILES,
)


class PhaCValidator:
    """
    PhaC enziminin islevselligini HMM Domain Alignment koordinatlari uzerinden dogrular.

    Katalitik triad (Cys-Asp-His) kontrolu, her sinif icin deneysel olarak
    dogrulanmis ve HMM match-state kolonlarina esleştirilmiş referans
    pozisyonlara (config.CATALYTIC_HMM_COLUMNS) dayanir. Her katalitik kalinti
    icin, beklenen kolona EN YAKIN hizalanmis kalinti secilir ve dar bir
    tolerans (config.TRIAD_TOLERANCE) penceresinde dogru amino asit aranir.
    PhaC Box motifi, katalitik Cys'in kendisine sabitlenir (box icindeki Cys =
    katalitik nukleofil).
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
            except Exception:
                pass  # HMM henuz uretilmediyse atla

        # Sinif-spesifik katalitik HMM kolonlari (config'den)
        self.ref_positions = CATALYTIC_HMM_COLUMNS

    def validate_triad_hmm(self, sequence: str, phac_class: str) -> dict:
        result = {
            "triad_found": False,
            "triad_residues": {},
            "box_found": False,
            "box_match": None,
            "is_functional": False,
            "notes": [],
        }

        hmm = self.class_hmms.get(phac_class)
        if not hmm:
            result["notes"].append("HATA: HMM yuklenemedi.")
            return result

        ref_cols = self.ref_positions.get(phac_class, {})
        if not ref_cols:
            result["notes"].append(f"HATA: {phac_class} icin katalitik kolon tanimi yok.")
            return result

        seq_clean = sequence.replace("*", "").replace("X", "")

        try:
            digital_seq = pyhmmer.easel.TextSequence(
                name=b"query", sequence=seq_clean
            ).digitize(self.alphabet)

            hits = list(pyhmmer.hmmsearch([hmm], [digital_seq], background=self.background, E=1e10, T=0))
            if not hits or not hits[0] or not hits[0][0].domains:
                result["notes"].append("HMM domain alignment alinamadi.")
                return result

            domain = hits[0][0].domains[0]
            aln = domain.alignment

            # Hizalama boyunca yuru: her hedef kalinti icin HMM kolonunu izle.
            # Her katalitik kalinti icin, referans kolona en yakin (tolerans
            # icindeki) eslesmeyi sec.
            best_match = {}  # res_name -> (col_distance, hmm_col, target_pos, residue)
            target_col = aln.target_from  # 1-indexed
            hmm_col = aln.hmm_from         # 1-indexed

            for hmm_char, target_char in zip(aln.hmm_sequence, aln.target_sequence):
                is_hmm_match = (hmm_char != ".")
                is_target_res = (target_char != "-")

                if is_target_res and is_hmm_match:
                    for res_name, ref_col in ref_cols.items():
                        dist = abs(hmm_col - ref_col)
                        if dist <= TRIAD_TOLERANCE:
                            prev = best_match.get(res_name)
                            if prev is None or dist < prev[0]:
                                best_match[res_name] = (
                                    dist, hmm_col, target_col, target_char.upper()
                                )

                if is_target_res:
                    target_col += 1
                if is_hmm_match:
                    hmm_col += 1

            # Triad degerlendirmesi: her kalinti icin, secilen pozisyonda dogru
            # amino asit var mi?
            found_residues = {}
            triad_ok = True
            cys_target_pos = None
            for res_name, exp_aa in CATALYTIC_TRIAD.items():
                match = best_match.get(res_name)
                if match and match[3] == exp_aa:
                    found_residues[res_name] = {
                        "hmm_col": match[1],
                        "target_pos": match[2],
                        "residue": match[3],
                    }
                    result["notes"].append(
                        f"{res_name}: {exp_aa} dogrulandi (kolon {match[1]}, pozisyon {match[2]})."
                    )
                    if res_name == "Cys":
                        cys_target_pos = match[2]
                else:
                    triad_ok = False
                    if match:
                        result["notes"].append(
                            f"{res_name}: kolon {match[1]}'de {match[3]} bulundu, "
                            f"{exp_aa} bekleniyordu."
                        )
                    else:
                        result["notes"].append(
                            f"{res_name}: katalitik kolon (~{ref_cols[res_name]}) "
                            f"+/-{TRIAD_TOLERANCE} penceresinde hizalanmis kalinti yok."
                        )

            result["triad_found"] = triad_ok
            result["triad_residues"] = found_residues

            # PhaC Box: katalitik Cys'e SABITLENMIS lipase box motifi.
            # Box icindeki Cys, triad'da bulunan katalitik Cys ile ayni
            # kalinti olmalidir.
            if cys_target_pos is not None:
                # Box motifi G-x-C-x-G-[GA]: katalitik Cys merkez (indeks +2).
                # Cys'in 2 once / 3 sonrasini iceren pencereyi test et.
                start = cys_target_pos - 3  # 1-indexed Cys -> 0-indexed window start
                window = seq_clean[max(0, start):cys_target_pos + 3]
                m = re.search(PHAC_BOX_REGEX, window)
                if m and "C" in m.group():
                    # box icindeki Cys gercekten katalitik Cys mi?
                    cys_idx_in_window = max(0, start)
                    box_cys_abs = cys_idx_in_window + m.start() + m.group().index("C")
                    if box_cys_abs + 1 == cys_target_pos:
                        result["box_found"] = True
                        result["box_match"] = m.group()
                    else:
                        result["notes"].append(
                            "PhaC Box bulundu ancak katalitik Cys'e sabitlenemedi."
                        )
                else:
                    result["notes"].append(
                        "Katalitik Cys cevresinde PhaC Box (G-x-C-x-G-G) motifi yok."
                    )
            else:
                result["notes"].append("Katalitik Cys bulunamadigi icin PhaC Box kontrol edilemedi.")

        except Exception as e:
            result["notes"].append(f"Hizalama hatasi: {e}")

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
            "all_scores": {},
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
                    if score > best_score and score > 20:  # Minimum cutoff
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
