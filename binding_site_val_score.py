from ccdc.protein import Protein
import os
import requests


## get all the entries for a specific ligand from pdb

base_url = 'https://www.ebi.ac.uk/pdbe/api/'
ligand_url = 'pdb/compound/in_pdb/'
validation_url = 'validation/summary_quality_scores/entry/'

def url_response(url):
    """
    Getting JSON response from URL
    :param url: String
    :return: JSON
    """
    r = requests.get(url=url)
    # Status code 200 means 'OK'
    if r.status_code == 200:
        json_result = r.json()
        return json_result
    else:
        print(r.status_code, r.reason)
        return None


def run_search(ligand_id):
    """

    Check pdbe search api documentation for more detials
    :param pdbe_search_term: String
    :return: JSON
    """
    # This constructs the complete query URL
    full_query = base_url + ligand_url+ ligand_id
    print(full_query)
    response = url_response(full_query)
    # if 'response' in response.keys() and 'docs' in response['response']:
    #    return response['response']['docs']
    # else:
    #    return None
    return response

sti_entries = run_search('STI')

list_of_pdb_entries = list(sti_entries.values())[0]
print(list_of_pdb_entries)


## download pdb files fromn list of entries

def fetch_pdb(entry_id):
    req = requests.get('http://files.rcsb.org/download/{}.pdb'.format(entry_id), '{}.pdb'.format(entry_id))
    filename = '{}.pdb'.format(entry_id)
    with open(filename, 'wb') as f:
        f.write(req.content)


for entry in list_of_pdb_entries:
    print(entry)
    fetch_pdb(entry)

## create protein objects  and dictionary of binding objects from list of entries

directory = r'C:\Users\amukhopadhyay\Documents\test_scripts'

list_of_protein_objs = []
binding_site_dict = {}

for filename in os.listdir(directory):
    if filename.endswith(".pdb"):
        print(filename)
        pdb_id = filename.split(".")[0]
        print(pdb_id)
        protein_from_entry = Protein.from_file(filename)
        print(len(protein_from_entry.chains))
        list_of_protein_objs.append(protein_from_entry)
        list_of_ligands = protein_from_entry.ligands
        list_of_ligands_identifiers = [ligands.identifier for ligands in list_of_ligands]
        print(list_of_ligands_identifiers)
        sti_list = (list(filter(lambda x: 'STI' in x, list_of_ligands_identifiers)))

        sti_indices = [i for i, s in enumerate(list_of_ligands_identifiers) if 'STI' in s]

        print(sti_indices)
        first_index_of_sti = sti_indices[0]
        binding_site = Protein.BindingSiteFromMolecule(protein_from_entry, list_of_ligands[first_index_of_sti], 6.)
        #print(binding_site.residues)
        binding_site_dict[pdb_id] = binding_site




# find global score of pdb entry list


def run_val_search(pdb_entry):
    """
    Check pdbe search api documentation for more detials
    :param pdbe_search_term: String
    :return: JSON
    """
    # This constructs the complete query URL
    full_query = base_url + validation_url+ pdb_entry
    print(full_query)
    val_score = url_response(full_query)

    return val_score

val_score_dict = {}
for entry in list_of_pdb_entries:
    val_dict = run_val_search(entry)
    for k,v in val_dict.items():
        val_score_dict[entry] = v['overall_quality']
print(val_score_dict)

#find the entry with best val score

best_entry = max(val_score_dict, key = val_score_dict.get)
print(best_entry)

