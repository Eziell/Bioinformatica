##Functions
from Bio import SeqIO
import json
import requests
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def getProteinArray(uniprot_file):   
        parkinson_sequences = SeqIO.parse(open(uniprot_file),'fasta')
        #retrieving protein names and IDs
        
        
        protein_info = []
        for fasta in parkinson_sequences:
                name, sequence, description= fasta.id, str(fasta.seq), str(fasta.description)
                id = name.split("|")
                entry =  [id[1], id[2], description, sequence]
                protein_info.append(entry)
        return protein_info

def get_Full_Description(uniprot_id):

    #Get the Uniprot data
    response = requests.get('http://www.uniprot.org/uniprot/' + uniprot_id + '.txt')
    data = response.text
    
    return data

def xml_get_NucleotideVar(uniprot_id):
    
    #Get the Uniprot data in xml
    response = requests.get('http://www.uniprot.org/uniprot/' + uniprot_id + '.xml')
    data = response.text
    
def get_NucleotideVar(uniprot_id):

    #Get the Uniprot data
    response = requests.get('http://www.uniprot.org/uniprot/' + uniprot_id + '.txt')
    data = response.text
    data_lines = data.splitlines()
    
    # Look for the RP identifiers followed by NUCLEOTIDE SEQUENCE
    rp_ids = set()
    for line in data_lines:
        if line.startswith('RP') and "VARIANT" in line:
            rp_ids.add(line)
    
    return "\n".join(rp_ids)

def get_full_annotation(uniprot_id):
    distinct_terms = set()
    
    # First get the Uniprot data
    response = requests.get('http://www.uniprot.org/uniprot/' + uniprot_id + '.txt')
    data = response.text
    data_lines = data.splitlines()
    
    # Look for the GO identifiers in the response and save them in a set
    go_ids = set()
    for line in data_lines:
        if line.startswith('DR   GO'):
            line_prefix, go_id, go_term_full, evidence_code = line.split(';')
            go_id = go_id.strip()
            go_ids.add(go_id)
    
    # Construct the URL for the QuickGO Request with all those GO identifiers
    the_url = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/'
    the_url += ','.join(go_ids)
    the_url += '/ancestors?relations=is_a'
    
    # Make the request and parse the results
    response = requests.get(the_url)
    anc_data = json.loads(response.text)
    results = anc_data['results']
    for record in results:
        distinct_terms.update(record['ancestors'])
                                        
    return distinct_terms

def get_term_info(go_id):
    the_url = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/' + go_id
    response = requests.get(the_url)
    data = json.loads(response.text)
    
    record = data['results'][0] # Get the first result, because we are only requesting one!
    
    # Return a dictionary that contains the information of this GO term
    go_info = {
        'go_id': record['id'],
        'aspect': record['aspect'],
        'name': record['name'],
        'definition': record['definition']['text'],
        'obsolete': record['isObsolete'],
    }
    
    return go_info

def heat_plot_analysis(protein_label, protein_terms, show_graph):
    protein_annotation_count = []
    for protein1 in protein_terms.items():
            pac_line = []
            for protein2 in protein_terms.items():
                    l1 = len(protein1)
                    l2 = len(protein2)
                    
                    if(l1 > l2):
                            count_ref = l1
                    else:
                            count_ref = l2
                    p = set(protein1)&set(protein2)
                    pac_line.append(round((len(p)*1.0/count_ref)*100))
            
            protein_annotation_count.append(pac_line)

    fig, ax = plt.subplots()
    im = ax.imshow(protein_annotation_count)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(protein_label)))
    ax.set_yticks(np.arange(len(protein_label)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(protein_label)
    ax.set_yticklabels(protein_label)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.

    ax.set_title("Relative number of common annotations:")
    fig.tight_layout()
    
    if(show_graph == True):
        plt.show()
    else:
        output = []
        header = "\t\t"+"\t\t".join(protein_label)
        output.append(header)
        for index, protein in enumerate(protein_annotation_count):
            output.append(protein_label[index] + "\t" + "\t\t".join(str(e) for e in protein))

    
    return output
            
            