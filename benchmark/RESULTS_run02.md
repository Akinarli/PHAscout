# Kör Koşu #02 — phaJ/PHBV düzeltmesi + dekontaminasyon sonrası (N=14)

**Tarih:** 2026-06-30
**Kod durumu:** phaJ ko-polimer düzeltmesi + PHBV/3HV `delta` düzeltmesi SONRASI.
**Veri durumu:** BLOSUM phaA/phaB referansları dekontamine edildi (Synechocystis 6803
ve S. meliloti 1021 kendi referanslarından çıkarıldı). *A. vinosum* DSM 180,
Class III'ü tanımlayan dizi olduğu (ve katalitik kolonların ondan türetildiği) için
referansta tutuldu; bunun yerine **benchmark'tan çıkarıldı** (yakılmış genom, demir
kural #5). Bu yüzden N: 15 → 14.

> **UYARI:** Bu 14 genom artık YAKILDI — held-out sayılmaz. Düzeltmelerin gerçek
> etkisi yalnızca TAZE, hiç dokunulmamış bir sette doğrulanabilir. Aşağıdaki sayılar
> bir regresyon kontrolüdür, bir genelleme iddiası DEĞİLDİR.

## Tespit (üretici potansiyeli var mı?)

| | Değer |
|---|---|
| TP / FN / FP / TN | 9 / 0 / 0 / 2 |
| Çekimser (belirsiz) | 3 (2 pozitif, 1 negatif) |
| Duyarlılık (recall) | 1.000 [%95 GA 0.701–1.000], n=9 |
| Kesinlik (precision) | 1.000, n=9 |
| Özgüllük | 1.000 [%95 GA 0.342–1.000], n=2 — **istatistiksel olarak anlamsız** |

Tespit yine "kusursuz" görünüyor ama güvenilmez: N küçük, negatif seti n=2, ve 3 zor
vaka "belirsiz"e kaçtı (aşağıya bak). Bu sayılar headline doğruluk olarak ASLA
sunulmamalı.

## Tip doğruluğu (yalnızca gerçek üreticiler): 8/9 = 0.889 [%95 GA 0.565–0.980]

Karışıklık (satır=gerçek → sütun=tahmin):
- MCL → MCL: 2 ✓
- SCL → SCL: 6 ✓
- SCL → MCL: 1 ✗ (*P. extremaustralis* — çift sentaz)

## run#01'e göre değişen (düzeltmelerin etkisi)

1. **PHBV/ko-polimer over-claim GİDERİLDİ.** run#01'de iki vaka SCL → SCL-co-MCL
   yanlışı veriyordu (*P. sacchari*, *B. thuringiensis*); ikisi de artık doğru **SCL**.
   Tip hatası 3 → 1'e düştü. (`pha_type.py` + `pathway_engine.py` `delta` düzeltmesi.)
2. **Veri sızıntısı ortaya çıktı — *Synechocystis* 6803 SCL → belirsiz.** run#01'de
   SCL deniyordu, ÇÜNKÜ kendi phaA/phaB dizileri referans setteydi (ezberleme).
   Referanslardan çıkarılınca araç phaA/phaB'yi doğrulayamadı ve dürüstçe çekimser
   kaldı. *Synechocystis* bilinen bir PHB (SCL) üreticisidir; bu artık dürüst bir
   "bilmiyorum"dur (gizli bir doğru-pozitif değil). **Bu, sızıntının somut kanıtıdır.**
3. ***S. meliloti* 1021 SCL kaldı** — kendi referansı çıkarılmasına rağmen phaA/phaB
   diğer Class I referanslarıyla eşleşti; bu tespit sağlamdı, sızıntıya bağlı değildi.

## Açık başarısızlık modları (değişmedi)

- **Çift sentaz (*P. extremaustralis*):** Class I (PHB) + Class II makinesi; araç
  route-complete Class II'yi seçip MCL dedi, baskın PHB fenotipiyle çelişti. İfade
  verisi olmadan çözümü zor.
- **Arkeal/dekontaminasyon kapsam boşluğu (*Haloferax*, *Synechocystis*):** monomer
  rotası genleri (phaA/phaB) doğrulanamayınca gerçek üreticide "belirsiz". Doğru
  davranış (uydurma yerine çekimserlik) ama kapsam eksikliğini gösterir.
- **Gizli özgüllük riski (*R. opacus* TAG):** TAG üreticisinde fonksiyonel Class_II
  phaC bulundu; yalnızca rota-yokluğu (çekimserlik) yanlış-pozitifi engelledi.

## Sonraki adım

Düzeltmeleri TAZE bir held-out sette (bu 14 + yakılmış eski genomlar hariç) doğrula.
Negatif seti genişlet (şu an n=2 — özgüllük ölçülemiyor). Arkeal/monomer-rota
kapsamını artır.
