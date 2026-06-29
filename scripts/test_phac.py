import sys
sys.path.insert(0, ".")
import logging
logging.basicConfig(level=logging.INFO)

from phascout.detection.phac_validator import PhaCValidator

v = PhaCValidator()

# Cupriavidus necator H16 PhaC (P23608 - Class I)
seq = (
    "MATGKGAAASTQEGKSQPFKVTPGPFDPATWLEWSRQWQGTEGNGHAAASGIPGLDALAG"
    "VKIAPAQLGDIQQRYMKDFSALWQAMAEGKAEATGPLHDRRFAGDAWRTNLPYRFAAAFY"
    "LLNARALTELADAVEADAKTRQRIRFAISQWVDAMSPANFLATNPEAQRLLIESGGESLR"
    "AGVRNMMEDLTRGKISQTDESAFEVGRNVAVTEGAVVFENEYFQLLQYKPLTDKVHARPL"
    "LMVPPCINKYYILDLQPESSLVRHVVEQGHTVFLVSWRNPDASMAGSTWDDYIEHAAIRA"
    "IEVARDISGQDKINVLGFCVGGTIVSTALAVLAARGEHPAASVTLLTTLLDFADTGILDV"
    "FVDEGHVQLREATLGGGAGAPCALLRGLELANTFSFLRPNDLVWNYVVDNYLKGNTPVPF"
    "DLLFWNGDATNLPGPWYCWYLRHTYLQNELKVPGKLTVCGVPVDLASIDVPTYIYGSRED"
    "HIVPWTAAYASTALLANKLRFVLGASGHIAGVINPPAKNKRSHWTNDALPESPQQWLAGA"
    "IEHHGSWWPDWTAWLAGQAGAKRAAPANYGNARYRAIEPAPGRYVKAKA"
)

r = v.full_analysis(seq)

print(f"\nSinif: {r['best_class']} (skor: {r['best_score']:.1f})")
print(f"Guven: {r['confidence']:.1f}")
print(f"Triad: {r['triad_found']}")
print(f"Box: {r['box_found']} ({r['box_match']})")
print(f"Fonksiyonel: {r['is_functional']}")
print(f"Onaylandi: {r['phac_confirmed']}")
for note in r["notes"]:
    print(f"  -> {note}")
