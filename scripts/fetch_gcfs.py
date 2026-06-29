import urllib.request
import json
import urllib.parse

def get_gcf(query):
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{urllib.parse.quote(query)}/dataset_report?filters.reference_only=true&filters.assembly_level=complete_genome'
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            if 'reports' in data and len(data['reports']) > 0:
                print(f'{query}: {data["reports"][0]["accession"]}')
                return
            else:
                url2 = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{urllib.parse.quote(query)}/dataset_report'
                req2 = urllib.request.Request(url2, headers={'Accept': 'application/json'})
                with urllib.request.urlopen(req2) as resp2:
                    data2 = json.loads(resp2.read().decode('utf-8'))
                    if 'reports' in data2 and len(data2['reports']) > 0:
                        print(f'{query}: {data2["reports"][0]["accession"]}')
                        return
    except Exception as e:
        print(f'{query} error: {e}')
    print(f'{query}: NOT FOUND')

get_gcf('Cupriavidus necator H16')
get_gcf('Pseudomonas putida KT2440')
get_gcf('Priestia megaterium DSM 319')
get_gcf('Aeromonas caviae')
get_gcf('Escherichia coli K-12')
get_gcf('Staphylococcus aureus')
