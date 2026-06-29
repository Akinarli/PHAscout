import urllib.request
import json
import urllib.parse
import csv
import time
import random

TAXA_TARGETS = {
    'SCL': [('Cupriavidus', 10), ('Halomonas', 10), ('Bacillus', 5), ('Azotobacter', 5), ('Chromobacterium', 5)],
    'MCL': [('Pseudomonas', 25)],
    'SCL-co-MCL': [('Aeromonas', 10)],
    'None': [('Escherichia', 10), ('Staphylococcus', 10), ('Yersinia', 5), ('Streptococcus', 5)]
}

def get_genomes_for_taxon(taxon, count):
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{urllib.parse.quote(taxon)}/dataset_report?filters.assembly_level=complete_genome'
    results = []
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            reports = data.get('reports', [])
            # Shuffle to get diverse strains
            random.shuffle(reports)
            for rep in reports:
                acc = rep.get('accession')
                org = rep.get('organism', {}).get('organism_name', taxon)
                if acc:
                    results.append({'Accession': acc, 'Species': org})
                if len(results) >= count:
                    break
    except Exception as e:
        print(f"Error fetching {taxon}: {e}")
    return results

def build_dataset():
    dataset = []
    for expected_type, taxon_list in TAXA_TARGETS.items():
        for taxon, count in taxon_list:
            print(f"Fetching {count} genomes for {taxon} ({expected_type})...")
            genomes = get_genomes_for_taxon(taxon, count)
            for g in genomes:
                g['Expected_Type'] = expected_type
                g['Expected_Class'] = 'Unknown'
                g['Notes'] = f"{taxon} {expected_type} producer"
                dataset.append(g)
            time.sleep(1)
            
    with open('benchmark_dataset_100.csv', 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Accession', 'Species', 'Expected_Type', 'Expected_Class', 'Notes'])
        writer.writeheader()
        writer.writerows(dataset)
    print(f"Dataset built with {len(dataset)} organisms.")

if __name__ == '__main__':
    build_dataset()
