import urllib.request
import json
import urllib.parse
import csv
import os
import time

def get_best_gcf(species_name):
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{urllib.parse.quote(species_name)}/dataset_report?filters.reference_only=true'
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            if 'reports' in data and len(data['reports']) > 0:
                # Prefer RefSeq
                for rep in data['reports']:
                    if rep['accession'].startswith('GCF_'):
                        return rep['accession']
                return data['reports'][0]['accession']
    except:
        pass
        
    url2 = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{urllib.parse.quote(species_name)}/dataset_report'
    try:
        req2 = urllib.request.Request(url2, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req2) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            if 'reports' in data and len(data['reports']) > 0:
                for rep in data['reports']:
                    if rep['assembly_info']['assembly_level'] == 'Complete Genome' and rep['accession'].startswith('GCF_'):
                        return rep['accession']
                for rep in data['reports']:
                    if rep['accession'].startswith('GCF_'):
                        return rep['accession']
                return data['reports'][0]['accession']
    except:
        pass
    return None

def process_table():
    species_list = []
    with open('table3_output.txt', 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Halomonas ') or line.startswith('Candidatus Halomonas '):
                species_list.append(line)
                
    # Remove duplicates
    species_list = list(set(species_list))
    print(f"Bulunan tur sayisi: {len(species_list)}")
    
    dataset = []
    for sp in species_list:
        print(f"Ara: {sp}...", end='')
        acc = get_best_gcf(sp)
        if acc:
            print(f" Bulundu ({acc})")
            dataset.append({
                'Accession': acc,
                'Species': sp,
                'Expected_Type': 'Unknown',
                'Expected_Class': 'Unknown',
                'Notes': 'Table3'
            })
        else:
            print(" YOK")
        time.sleep(0.5)
        
    with open('table3_dataset.csv', 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Accession', 'Species', 'Expected_Type', 'Expected_Class', 'Notes'])
        writer.writeheader()
        writer.writerows(dataset)
        
    print(f"\n{len(dataset)} genom bulundu. Test basliyor...")

if __name__ == '__main__':
    process_table()
