"""
PHAscout Streamlit Web Arayuzu
================================
Kullanici dostu web arayuzu.
Calistirmak icin: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import json
import logging
from phascout.pipeline import PHAscoutPipeline

# Streamlit sayfa ayarlari
st.set_page_config(
    page_title="PHAscout | PHA Uretici Bakteri Taramasi",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Arka plan loglarini gizle (sadece UI uzerinden gosterecegiz)
logging.getLogger("phascout").setLevel(logging.CRITICAL)


@st.cache_resource
def get_pipeline():
    return PHAscoutPipeline()


def main():
    st.title("🧬 PHAscout")
    st.subheader("Bakteriyel Polihidroksialkanoat (PHA) Üreticisi Tespit Aracı")
    
    st.markdown("""
    PHAscout, genom dizilerinden (proteomlardan) yola çıkarak bir bakterinin PHA plastiği 
    üretip üretemeyeceğini, hangi karbon kaynaklarını kullanabileceğini ve hangi yolaklara 
    sahip olduğunu **Yapay Zeka Destekli HMM Modelleri** ve **Çift Katmanlı BLOSUM62 Filtresi** 
    ile yüksek doğrulukla tahmin eder.
    """)
    
    st.sidebar.header("Taramayı Başlat")
    input_method = st.sidebar.radio(
        "Girdi Yöntemi",
        ["NCBI Accession", "FASTA Metni Yükle"]
    )
    
    pipeline = get_pipeline()
    
    if input_method == "NCBI Accession":
        accession = st.sidebar.text_input("NCBI Assembly Accession (Örn: GCF_000009285.1)", "")
        if st.sidebar.button("Analizi Başlat", type="primary"):
            if not accession:
                st.sidebar.error("Lütfen bir accession numarası girin.")
            else:
                run_analysis(pipeline, accession=accession.strip())
                
    elif input_method == "FASTA Metni Yükle":
        fasta_text = st.sidebar.text_area("Protein FASTA Metni", height=200)
        if st.sidebar.button("Analizi Başlat", type="primary"):
            if not fasta_text:
                st.sidebar.error("Lütfen FASTA formatında metin girin.")
            else:
                run_analysis(pipeline, fasta_text=fasta_text)


def run_analysis(pipeline, accession=None, fasta_text=None):
    with st.spinner("PHAscout analiz yapıyor... Lütfen bekleyin (NCBI indirmesi 1-2 dakika sürebilir)."):
        try:
            report = pipeline.run(accession=accession, fasta_text=fasta_text)
            display_report(report)
        except Exception as e:
            st.error(f"Bir hata oluştu: {str(e)}")


def display_report(report):
    st.success("Analiz başarıyla tamamlandı!")
    
    # --- ORGANİZMA VE SONUÇ ÖZETİ ---
    org = report.get("organism", {})
    summ = report.get("summary", {})
    
    st.header(f"🦠 Organizma: {org.get('organism_name', 'Bilinmiyor')}")
    st.caption(f"Accession: {org.get('accession', 'N/A')} | Seviye: {org.get('assembly_level', 'N/A')}")
    
    # 3'lü Metrik Kartları
    c1, c2, c3 = st.columns(3)
    
    with c1:
        if summ.get("produces_pha"):
            st.success("✅ PHA Üreticisi (Pozitif)")
        else:
            st.error("❌ PHA Üreticisi Değil (Negatif)")
            
    with c2:
        phac_class = summ.get("phac_class")
        if phac_class and phac_class != "None":
            st.info(f"🧬 Polimeraz: {phac_class.replace('_', ' ')}")
        else:
            st.warning("🧬 Polimeraz: Bulunamadı")
            
    with c3:
        score = summ.get("heuristic_score", 0)
        max_score = summ.get("heuristic_max", 92)
        tier = summ.get("heuristic_tier", "")
        
        if "Yüksek" in tier:
            st.success(f"🏆 Potansiyel: {score}/{max_score} ({tier})")
        elif "Orta" in tier:
            st.info(f"⭐ Potansiyel: {score}/{max_score} ({tier})")
        else:
            st.warning(f"⚪ Potansiyel: {score}/{max_score} ({tier})")

    st.divider()

    # --- DETAYLI ANALİZ SEKME YAPISI ---
    tab1, tab2, tab3, tab4 = st.tabs([
        "🧬 Genetik Profil", 
        "⚙️ PhaC (Polimeraz) Analizi", 
        "🔄 Metabolik Yolaklar", 
        "📊 Ham JSON Raporu"
    ])
    
    # Sekme 1: Genetik Profil
    with tab1:
        st.subheader("Bakteriyel PHA Sentez Genleri")
        genes = report.get("genes", {})
        details = genes.get("details", {})
        
        # Tablo verisi hazırla
        table_data = []
        for gene in ["phaC", "phaA", "phaB", "phaJ", "phaG", "phaP", "phaR", "phaE"]:
            g_det = details.get(gene, {})
            is_found = g_det.get("detected", False)
            table_data.append({
                "Gen Adı": gene,
                "Durum": "✅ Bulundu" if is_found else "❌ Yok",
                "Protein ID": g_det.get("protein_id", "-"),
                "E-value": f"{g_det.get('evalue'):.1e}" if g_det.get("evalue") is not None else "-",
                "BLOSUM Filtresi": g_det.get("filter_note", "Gerekmiyor / Hesaplanmadı")
            })
            
        df = pd.DataFrame(table_data)
        st.dataframe(df, use_container_width=True, hide_index=True)
        
    # Sekme 2: PhaC Analizi
    with tab2:
        st.subheader("Polihidroksialkanoat Sentaz (PhaC) Detayları")
        phac = report.get("phac_analysis", {})
        
        if phac.get("class"):
            c1, c2 = st.columns(2)
            with c1:
                st.metric("Sınıflandırma", phac.get('class'))
                st.metric("Güven Skoru", f"{phac.get('confidence', 0):.1f}")
            with c2:
                func = "✅ Aktif (Fonksiyonel)" if phac.get("functional") else "❌ Pasif"
                st.metric("Enzim Durumu", func)
                st.metric("PhaC Box Motifi", phac.get("box_match", "Bulunamadı"))
            
            st.markdown("#### Biyoinformatik Doğrulama (Triad & Motif)")
            for note in phac.get("notes", []):
                if "FONKSİYONEL" in note:
                    st.success(note)
                else:
                    st.info(note)
        else:
            st.error("Bu organizmada geçerli bir PhaC geni bulunamamıştır.")
            
    # Sekme 3: Metabolik Yolaklar
    with tab3:
        st.subheader("Aktif PHA Sentez Yolakları")
        pathways = report.get("pathways", [])
        
        active_pw = [pw for pw in pathways if pw["active"]]
        if active_pw:
            for pw in active_pw:
                with st.expander(f"✅ {pw['name']} (AKTİF)", expanded=True):
                    st.write(f"**Ürün Eğilimi:** {pw['product_tendency']}")
                    st.write(f"**Uyumlu Karbon Kaynakları:** {', '.join(pw['carbon_sources'])}")
                    st.write(f"**Açıklama:** {pw['note']}")
        else:
            st.warning("Bu organizmada çalışan bir PHA sentez yolağı tespit edilemedi.")
            
        st.subheader("İnaktif Yolaklar")
        inactive_pw = [pw for pw in pathways if not pw["active"]]
        for pw in inactive_pw:
            with st.expander(f"❌ {pw['name']} (Eksik)"):
                st.write(f"**Eksik Genler:** {', '.join(pw['missing_genes'])}")
                
    # Sekme 4: Ham JSON
    with tab4:
        st.download_button(
            label="Raporu İndir (JSON)",
            data=json.dumps(report, indent=2, ensure_ascii=False),
            file_name=f"{org.get('accession', 'report')}_phascout.json",
            mime="application/json"
        )
        st.json(report, expanded=False)

if __name__ == "__main__":
    main()
