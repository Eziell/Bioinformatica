
import requests
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import uniRetrieval

#opening proteins.tsv for writting procedures.
f = open('proteins.tsv', "w")

OMIM="blast.fasta"


protein_info = uniRetrieval.getProteinArray(OMIM)
print(protein_info)

#writting focus proteins to file
f.write(str("Retrieving protein names and IDs from "+ OMIM + "\n"))
f.write(str("ProteinName\tProteinID\n"))
for protein in protein_info:
        f.write(protein[1] + "\t" + protein[0] + "\n")

#GO term counting 
protein_terms = dict()
protein_label = []
all_terms = []
for protein in protein_info:
        #get the annotations
        protein_id = protein[0]
        annotations = uniRetrieval.get_full_annotation(protein_id)
        protein_terms[protein_id] = annotations
        protein_label.append(protein_id)
        all_terms.extend(annotations)

#Heatmap
f.write("By looking at the correlation plot, there is a clear divide between 2 groups of proteins\n")
corr_plot = uniRetrieval.heat_plot_analysis(protein_label, protein_terms, False)
f.writelines(["%s\n" % item for item in corr_plot])

#################################################################
term_counts = {term: all_terms.count(term) for term in all_terms}
distinct_terms = list(term_counts)


# Get common GO terms
f.write(str("\nCommon GO terms ("+str(len(protein_info))+" proteins):\n"))
f.write(str("GO ID:\tAspect:\tName:\n"))
for go_id in distinct_terms:
    if term_counts[go_id] == 11:
        go_info = uniRetrieval.get_term_info(go_id)
        f.write(str(go_id+ '\t'+ go_info['aspect']+ '\t'+ go_info['name'] + '\n'))
###########################

# Get unique GO terms
f.write(str("\nUnique GO terms ("+str(len(protein_info))+" proteins):\n"))
f.write(str("GO ID:\tAspect:\tName:\n"))
for go_id in distinct_terms:
    if term_counts[go_id] == 1:
        
        go_info = uniRetrieval.get_term_info(go_id)
        f.write(str(go_id+ '\t'+ go_info['aspect']+ '\t'+ go_info['name'] + '\n'))
###########################