# PHAscout — Dürüst Benchmark (Adım A)

Bu klasör, PHAscout'un GERÇEK doğruluğunu ölçmek içindir. 20-genom smoke
testinin aksine, buradaki tek amaç **acımasız dürüstlük**: aracın nerede
çalıştığını ve — daha önemlisi — NEREDE PATLADIĞINI görmek.

## Demir kurallar (ihlal = benchmark çöp olur)

1. **Etiketler ıslak-laboratuvardan gelir, anotasyondan DEĞİL.**
   "NCBI class I PHA synthase diyor" bir kanıt değildir — bu, kendi
   gen-bulucumuzu başka bir gen-bulucuyla kıyaslamaktır (döngüsellik).
   Geçerli kanıt: GC-MS / NMR / FTIR ile karakterize edilmiş polimer,
   Nile-red / Sudan-black boyama, gravimetrik birikim %, TEM granülleri.
   `evidence_method` sütunu bu yüzden zorunludur ve `annotation`/boş olamaz.

2. **Suş-eşleşmesi.** `accession`, etiketin geldiği çalışmadaki SUŞA ait
   olmalı. "B. megaterium PHB üretir" yetmez — hangi suş, o suşun genomu mu?

3. **Atıf zorunlu.** Her satırda DOI veya PMID. Atıfsız satır koşturulmaz.

4. **TEST SETİNDE AYAR YASAK.** Genomları bir kez koştur, sonuca bak, hata
   analizini yap. Eşik/parametre oynatıp AYNI sette tekrar koşmak overfitting'tir
   — kaçtığımız tiyatronun ta kendisi. İyileştirmeden sonra FARKLI, taze bir
   held-out sette doğrula.

5. **Yakılmış genomlar geçersiz.** Geliştirmede kullanılmış her genom
   (Halomonas alkalicola, P. putida KT2440, V. cholerae, S. pyogenes ve eski
   20-genom seti) held-out SAYILMAZ. Sadece daha önce hiç dokunulmamış genomlar.

## Dosyalar

- `dataset.csv` — etiketli held-out genomlar (sen doldurursun; şema aşağıda).
- `run.py` — her accession için TAM pipeline'ı bir kez koşturur, ham raporu
  `cache/` altına yazar (tekrar koşunca NCBI'a gitmez, "tek koşu" disiplinini
  cache zorlar), `predictions.csv` üretir.
- `metrics.py` — `dataset.csv` + `predictions.csv`'i birleştirip dürüst metrik
  basar: tespit karışıklık matrisi, tip doğruluğu, filum/sınıf kırılımı, Wilson
  %95 GA, ve negatif-setin zayıflığına dair açık uyarılar.

## `dataset.csv` şeması

| sütun | açıklama |
|-------|----------|
| `accession` | GCF_/GCA_ (suş-eşleşmiş) |
| `organism` | tür adı |
| `strain` | kesin suş |
| `phylum` | filum (yanlılık analizi için) |
| `synthase_class` | beklenen sınıf (Class_I..IV / unknown) — biliniyorsa |
| `label_produces` | `yes` / `no` (ıslak-lab fenotip) |
| `label_type` | `SCL` / `MCL` / `SCL-co-MCL` / `none` |
| `carbon_source` | çalışmada kullanılan substrat |
| `evidence_method` | GC-MS / NMR / FTIR / Nile_red / gravimetric / TEM (annotation YASAK) |
| `citation` | DOI veya PMID (zorunlu) |
| `notes` | serbest not |

## Çalıştırma

```
python -m benchmark.run --dataset benchmark/dataset.csv
python -m benchmark.metrics
```

## Beklenti

İlk sayı muhtemelen düşük olacak — özellikle özgüllük (zayıf negatifler) ve
model-organizma dışı taksonlarda. Bu bir başarısızlık değil, ÖLÇÜMÜN KENDİSİDİR.
Çıkan hataları teşhis edip iteratif çözeceğiz.
