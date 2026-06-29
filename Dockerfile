FROM python:3.12-slim

LABEL maintainer="PHAscout Project"
LABEL description="PHAscout - PHA Üretici Bakteri Tarama Aracı"

# Sistem bağımlılıkları
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Çalışma dizini
WORKDIR /app

# Python bağımlılıkları
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Proje dosyaları
COPY . .

# PFAM veritabanı ve HMM profilleri konteyner içinde kalır
VOLUME ["/app/data"]

# Varsayılan komut: Streamlit arayüzü
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]


