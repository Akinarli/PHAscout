# PHAscout — Proje Planı

## 1. Dürüst kapsam (omurga)

**Ne yapar:** Bir bakteri genomu (NCBI GCF/GCA accession veya proteom FASTA) verildiğinde,
PHA biyosentez genlerini bulur, PhaC'yi sınıflandırır (I–IV) ve katalitik triad'ı doğrular;
ardından bilinen biyokimyaya göre organizmanın hangi **PHA tipini üretme potansiyeli**
taşıdığını raporlar.

**Ne yapmaz / iddia etmez:** "Bu suş PHA üretir" demez. Genomda gen bulunması üretimi
garanti etmez (ifade düzeyi, karbon kaynağı, büyüme koşulları, pseudogenler devrede).
Çıktı bir **potansiyel/gen-donanımı** raporudur, üretim hükmü değil.

**Niş (literatürle):** Mevcut phaC sınıflandırma yöntemleri "sadece iyi çalışılmış cinslerde
işliyor" ve "dört sınıfı birden hedefleyemiyor" (Microbial Cell Factories 2025, dejenere primer
çalışması). PHAscout bu açığı **in silico**, filumlar arası, dört sınıfı da kapsayan,
triad-doğrulamalı bir araçla kapatır. Yöntem yeniliği değil; sağlam + erişilebilir entegrasyon.

## 2. Boru hattı

1. **Girdi:** GCF → proteom (+ GFF, operon analizi için).
2. **Katman 1 (geniş HMM ağı):** aday genler — phaC, phaA, phaB, phaG, phaJ, phaE, phaR, phaP.
3. **PhaC çekirdeği (SAĞLAM):** sınıf I–IV ataması + katalitik triad (kesin HMM kolonlarında
   Cys-Asp-His) + lipase box doğrulaması. Tüm phaC adayları taranır, en güçlü fonksiyonel olan seçilir.
4. **Yardımcı gen onayı:** HMM + double-layer + **operon konteksti** (phaB↔FabG ayrımı için sinteni).
5. **Yorum (`pha_type`):** sınıf + gen donanımı → PHA tipi potansiyeli (none/SCL/MCL/SCL-co-MCL) + rota.
6. **Rapor:** kanıt tablosu + güven + açık belirsizlik etiketleri + görsel kart.

## 3. PhaC sınıfı → PHA tipi eşlemesi (literatür: Rehm 2003; Steinbüchel)

| Sınıf | Alt birim | Temel PHA tipi |
|---|---|---|
| Class I | PhaC (homodimer) | SCL-PHA (C3–C5): P(3HB), yan genlerle ko-polimer |
| Class II | PhaC1/PhaC2 | MCL-PHA (C6–C14) |
| Class III | PhaC + PhaE | SCL-PHA |
| Class IV | PhaC + PhaR | SCL-PHA |

Yardımcı genlerle rota netleşir:
- **phaA + phaB** → asetil-CoA'dan 3HB (şekerden, SCL).
- **phaG** (Class II) → de novo yağ asidi sentezinden MCL (şekerden).
- **phaJ** → β-oksidasyondan MCL/ko-polimer (yağ asidinden); Class I + phaJ → P(3HB-co-3HHx) (*Aeromonas*).

## 4. Çıktının dürüst çerçevesi (en kritik değişiklik)

- `produces_pha: EVET/HAYIR` damgası **kaldırıldı**. Yerine `pha_potential`:
  `{ tip: none|SCL|MCL|SCL-co-MCL, ürünler, rotalar, gen-kanıtı, uyarı }`.
- Her gen için **dürüst etiket**: "phaB: operon-destekli" / "phaB: aday — FabG'den ayrılamadı".
- Dil her yerde "potansiyel/kapasite", asla "üretir".

## 5. Görselleştirme

1. **Genom rapor kartı** (ana çıktı): organizma + PHA tipi potansiyeli + PhaC sınıf/güven + triad +
   gen donanımı + aktif rota + belirsizlik uyarısı.
2. Yolak şeması (substrat → genler → PHA tipi; var/eksik adımlar renkli).
3. Gen donanım matrisi (✓/✗ + güven).

## 6. Varyasyon dayanıklılığı (organizmadan organizmaya mutasyon)

- **PhaC fonksiyon çağrısı sağlam:** katalitik triad saflaştırıcı seçilim altında korunur;
  kesin kolonlara çapalıyoruz (arkeden Bacillus'a çalıştığı doğrulandı). Katalitik kalıntı
  mutasyona uğrarsa triad kontrolü bunu "fonksiyonel değil" diye yakalar.
- **Yardımcı genler kırılgan:** tek-BLOSUM-eşiği + Proteobakteri-yanlı referanslar, ıraksak gerçek
  pozitifleri kaçırır (ör. B. megaterium phaA: 0.662 vs eşik 0.716). Çözüm yol haritasında (madde 8).

## 7. Değerlendirme (dürüst)

- 20-genomluk set (10 poz + 10 neg) artık **tuning hedefi değil**, yalnızca hızlı smoke-test.
- Gerçek performans için: dokunulmamış, daha büyük, çoklu-filum, held-out genom seti — bir kez ölç.
- Sınıf→tip eşlemesi literatüre dayanır (yöntem değil, biyokimya).

## 8. Yol haritası (öncelik sırası)

- [x] PhaC çekirdeği: kesin katalitik kolonlar, ±2 tolerans, box katalitik Cys'e çapalı.
- [x] Referans dekontaminasyonu + kalibrasyon (PhaB↔FabG F1 0.50 → 0.92).
- [x] phaJ PFAM düzeltmesi (PF00767 → PF01575).
- [x] Tüm phaC adaylarını tarama (best-E-value tuzağı düzeltildi).
- [ ] **Çıktıyı "PHA potansiyeli" diline çevir** (bu PR).
- [ ] Operon'u phaB↔FabG için birinci-sınıf kanıt yap.
- [ ] Yardımcı gen referanslarını filogenetik olarak genişlet (Firmicutes/Archaea/Actinobacteria).
- [ ] Tek-BLOSUM-eşiğini kalibre profil HMM'lere geçir.
- [ ] Görsel rapor kartını gerçek pipeline çıktısından üret.
- [ ] Büyük, bağımsız held-out benchmark ile bir kez dürüst ölçüm.
