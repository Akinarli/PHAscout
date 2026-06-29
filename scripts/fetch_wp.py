import urllib.request

def fetch_protein(accession):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        fasta = response.read().decode('utf-8')
        return fasta

print(fetch_protein("WP_045991792.1"))
