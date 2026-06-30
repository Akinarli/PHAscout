# Kör Koşu #01 — Baseline Sonuçları (N=15)

**Tarih:** 2026-06-30
**Kod durumu:** phaJ ko-polimer düzeltmesinden ÖNCE (baseline).
**Önemli:** Bu 15 genom artık YAKILDI (held-out sayılmaz). Sonuçlar phaJ
düzeltmesinin ÖNCESİNE aittir; düzeltme sonrası doğrulama TAZE bir sette yapılmalıdır.

## Tespit (üretici mi?)

| | Değer |
|---|---|
| TP / FN / FP / TN | 11 / 0 / 0 / 2 |
| Çekimser (belirsiz) | 2 (1 pozitif, 1 negatif) |
| Duyarlılık (recall) | 1.000 [%95 GA 0.741–1.000], n=11 |
| Özgüllük | 1.000 [%95 GA 0.342–1.000], n=2 — **istatistiksel olarak anlamsız** |

Tespit "kusursuz" görünüyor ama güvenilmez: N küçük, özgüllük n=2, ve en zor 2 vaka
(Archaea üreticisi + TAG negatifi) "belirsiz"e kaçtı.

## Tip doğruluğu (yalnızca gerçek üreticiler): 8/11 = 0.727 [%95 GA 0.434–0.903]

Karışıklık (satır=gerçek → sütun=tahmin):
- MCL → MCL: 2 ✓
- SCL → SCL: 6 ✓
- SCL → MCL: 1 ✗ (P. extremaustralis)
- SCL → SCL-co-MCL: 2 ✗ (P. sacchari, B. thuringiensis)

## Teşhis edilen başarısızlık modları

1. **phaJ ko-polimer over-claim (ANA BULGU, DÜZELTİLDİ):** phaJ tek başına SCL-class
   organizmayı SCL-co-MCL'e yükseltiyordu. 3HHx ancak MCL-yetenekli (Class II) sentazla
   polimerleşir; Class I/III/IV SCL-spesifiktir. B. thuringiensis için Class IV + phaJ ile
   SCL-co-MCL demek biyokimyasal olarak imkansızdı. → `pha_type.py` düzeltildi, `tests/test_pha_type.py`.
2. **Çift sentaz (P. extremaustralis, AÇIK):** Class I (PHB) + Class II makinesi birlikte;
   araç route-complete olan Class II'yi seçip MCL dedi, baskın PHB fenotipiyle çelişti.
   İfade verisi olmadan çözümü zor.
3. **Arkeal kapsam boşluğu (Haloferax, AÇIK):** Bakteriyel HMM'ler arkeal phaA/route
   genlerini kaçırdı → gerçek üreticide "belirsiz". phaE bulup sınıfı Class_I demesi de tutarsız.
4. **Gizli özgüllük riski (R. opacus TAG, AÇIK):** TAG üreticisinde fonksiyonel Class_II phaC
   bulundu; yalnızca rota-yokluğu yanlış-pozitifi engelledi.
