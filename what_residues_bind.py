import requests
import os
import pandas as pd
from ccdc.protein import Protein
import argparse
import logging
logging.basicConfig(level=logging.INFO)

class Pdbesearch:
    def __init__(self, ligand_id):
        self.search_url = 'https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/'
        self.search_options = '&wt=json&rows=100000'
        self.ligand_id = ligand_id


    def url_response(self, url):
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

    def run_search(self):

        full_query = self.search_url + self.ligand_id
        print(full_query)
        response = self.url_response(full_query)
        return response



def fetch_pdb(entry_id):
    req = requests.get('http://files.rcsb.org/download/{}.pdb'.format(entry_id), '{}.pdb'.format(entry_id))
    
    filename = '{}.pdb'.format(entry_id)
    with open(filename, 'wb') as f:
                       f.write(req.content)

def create_plot(list_of_aa):
    residue_df = pd.DataFrame(list_of_aa, columns=['aa'])

    residue_df.count()

    pt = pd.value_counts(residue_df['aa']).plot.bar()
    fig = pt.get_figure()
    fig.savefig("output.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--het_code', '--het', help='het_code pdb')
    args = parser.parse_args()

    P = Pdbesearch(args.het_code)
    sti_entries = P.run_search()
    list_of_entries = list(sti_entries.values())[0]

    for entry in list_of_entries:
        print(entry)
        fetch_pdb(entry)
    logging.info(list_of_entries)

    directory = r'C:\Users\amukhopadhyay\Documents\test_scripts'

    list_of_sti_binding_residues = []
    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            protein_from_entry = Protein.from_file(filename)
            list_of_ligands = protein_from_entry.ligands
            list_of_ligands_identifiers = [ligands.identifier for ligands in list_of_ligands]
            sti_list = (list(filter(lambda x: 'STI' in x, list_of_ligands_identifiers)))
            sti_indices = [i for i, s in enumerate(list_of_ligands_identifiers) if 'STI' in s]

            for ligand_sti in sti_indices:
                list_of_sti_binding_residues.append(((Protein.BindingSiteFromMolecule(protein_from_entry, list_of_ligands[ligand_sti], 6.).residues)))
        list_of_residues = [item for t in list_of_sti_binding_residues for item in t]
    list_of_aa = []
    for i in list_of_residues:
        list_of_aa.append((i.identifier[2:5]))

    create_plot(list_of_aa)


if __name__ == "__main__":
    main()