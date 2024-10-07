import re
import os
import time
import requests
import pandas as pd
import datetime
import gzip
import json
 
from openpyxl import load_workbook
from pathlib import PureWindowsPath, Path
from download_proteomes import download_proteome, log_proteomes

referece_proteomes_table_path = Path(PureWindowsPath('proteomes_AND_proteome_type_1_2023_11_28.tsv'))
proteome_table = pd.read_csv(referece_proteomes_table_path, sep = "\t")
proteome_table = proteome_table.set_index('Organism Id')
log_file = Path(PureWindowsPath('database_log.xlsx'))

## Extract ID, AC (accession), DT (date, last version), DE (description) full and short name (recname first), GN (gene), OH (host organism if applicable)
## Extrac CC (comments) - Caution - Cofactor - Developmental Stg - Disease - Disruption phenotype - Domain - Interactions
## Extract CC - Misc - Pathway - Pharmaceutical - Similarity - Subcellular location - Subunit - Tissue specificity - toxic dose 
## Extract DR - GO - KEGG - REACTOME
## Extract KW (keywords)

def grab_ID_info(line):
    out_dict = dict()
    info = re.split('[ ]+', line)
    out_dict['entry_ID'] = info[1]
    out_dict['review_status'] = info[2].split(';')[0]
    out_dict['seq_len'] = info[3]
    return out_dict

def grab_accession(line, secondary_accessions = False, out_dict = dict()):
    info = re.split('[ ;]+', line)
    accession = info[1:-1]
    if secondary_accessions:
        if 'secondary_accession' not in out_dict.keys():
            out_dict['secondary_accession'] = accession
        else:
            out_dict['secondary_accession'] = out_dict['secondary_accession'] + accession
    else:
        out_dict['primary_accession'] = accession[0]
    if accession[1:]:
        out_dict['secondary_accession'] = accession[1:]
    return out_dict

def grab_date(line):
    out_dict = dict()
    info = re.split('[ ]{3}', line)
    info = re.split(', ', info[1])
    out_dict['version_date'] = info[0]
    out_dict['version_number'] = info[1].split('\n')[0]
    return out_dict

def grab_description(file, line):
    main_flags = ['RecName', 'SubName']
    out_dict = dict()
    before_contains = True
    while(line.startswith('DE   ')):
        if line.startswith('DE   Contains:'):
            before_contains = False

        if before_contains:
            for main_flag in main_flags:
                if line.startswith(f'DE   {main_flag}: Full='):
                    fullname = line.split(f'DE   {main_flag}: Full=')[1]
                    fullname = fullname.split(' {')[0].split(';')[0]
                    out_dict[f'{main_flag}_full'] = fullname      
        bline = file.readline()
        line = bline.decode("utf-8")
    return(out_dict, file, line)
        
def get_all_features(proteome_path):
    all_features = []
    with open(proteome_path) as proteome_fasta:
        for line in proteome_fasta:
            if line.startswith('>'):  
                all_features = all_features + [line.split('|')[1]]
    return all_features

def grab_comments(file, line):
    ## Extract CC - Caution - Cofactor - Developmental Stg - Disease - Disruption phenotype - Domain - Interactions
    ## Extract CC - Misc - Pathway - Pharmaceutical - Similarity - Subcellular location - Subunit - Tissue specificity - toxic dose 
    main_flags = ['CAUTION', 'COFACTOR', 'DEVELOPMENTAL STAGE', 'DISEASE', 'DISRUPTION PHENOTYPE', 'DOMAIN']
    main_flags = main_flags + ['FUNCTION', 'INTERACTION', 'MISCELLANEOUS', 'PATHWAY', 'PHARMACEUTICAL', 'SIMILARITY']
    main_flags = main_flags + ['SUBCELLULAR LOCATION', 'SUBUNIT', 'TISSUE SPECIFICITY', 'TOXIC DOSE']
    out_dict = dict()
    comment = ''
    flag_key = ''

    while(line.startswith('CC  ')):
       if line.startswith(f'CC   -!- '):
            if flag_key != '':
                out_dict[flag_key] = comment
            for main_flag in main_flags:
                if line.startswith(f'CC   -!- {main_flag}: '):
                    comment = line.split(f'CC   -!- {main_flag}: ')[1].split('\n')[0]
                    flag_key = main_flag
       elif line.startswith('CC       '):
            new_text = line.split('CC       ')[1].split('\n')[0]
            comment = f'{comment} {new_text}'
       bline = file.readline()
       line = bline.decode("utf-8")
    out_dict[flag_key] = comment
    return(out_dict, file, line)

def parse_request(request):
    with gzip.open(request.raw) as file:
        new_dict = dict()
        all_annotations = dict()
        new_annotation = dict()
        bline = file.readline()
        line = bline.decode("utf-8")
        gn_done = False
        while bline:
            get_new_line = True
            if line.startswith('//'):
                all_annotations[accession] = new_annotation
                new_annotation = dict()
                gn_done = False
            elif line.startswith('ID   '):
                new_dict = grab_ID_info(line)
                new_annotation['Identification (ID)'] = new_dict  
                new_dict = dict()
            elif line.startswith('AC   '):
                if new_dict:
                    new_dict = grab_accession(line, secondary_accessions = True, out_dict = new_dict)
                else:
                    new_dict = grab_accession(line, out_dict = dict())
                accession = new_dict['primary_accession']
                new_annotation['Accession (AC)'] = new_dict  
            elif line.startswith('DT   ') and ', entry version' in line:
                new_dict = grab_date(line)
                new_annotation['Date (DT)'] = new_dict  
            elif line.startswith('DE   '):
                new_dict, file, line = grab_description(file, line)
                new_annotation['Description (DE)'] = new_dict 
                ## grab_description gives the next line
                get_new_line = False  
            elif line.startswith('GN   Name=') and not gn_done:
                gn_done = True
                new_dict = dict()
                new_dict['official_uniprot_gene'] = line.split('GN   Name=')[1].split(';')[0].split(' {')[0]
                new_annotation['Gene name (GN)'] = new_dict
            elif line.startswith('CC   -!- '):
                new_dict, file, line = grab_comments(file, line)
                new_annotation['Comments (CC)'] = new_dict
                get_new_line = False
            elif line.startswith('DR   GO; '):
                go_annotation = line.split('DR   GO; ')[1].split('; ')
                new_dict = dict()
                new_dict['GO_code'] = go_annotation[0]
                new_dict['GO_evidence'] = go_annotation[2].split('\n')[0]
                if 'Gene Ontology annotations (DR)' not in new_annotation.keys():
                    new_annotation['Gene Ontology annotations (DR)'] = dict()
                new_annotation['Gene Ontology annotations (DR)'][go_annotation[1]] = new_dict
            elif line.startswith('DR   KEGG; '):
                new_annotation['KEGG identifier (DR)'] = line.split('DR   KEGG; ')[1].split(';')[0]
            elif line.startswith('DR   Reactome; '):
                reactome_annotation = line.split('DR   Reactome; ')[1].split('; ')
                new_dict = dict()
                new_dict['Reactome_code'] = reactome_annotation[0]
                if 'Reactome annotations (DR)' not in new_annotation.keys():
                    new_annotation['Reactome annotations (DR)'] = dict()
                new_annotation['Reactome annotations (DR)'][reactome_annotation[1]] = new_dict
            elif line.startswith('KW   '):
                new_keywords = line.split('KW   ')[1].split('\n')[0].split(';')
                new_keywords = [x for x in new_keywords if x != '']
                if 'Keywords (KW)' in new_annotation.keys():
                    new_annotation['Keywords (KW)'] = new_annotation['Keywords (KW)'] + new_keywords 
                else:
                    new_annotation['Keywords (KW)'] = new_keywords 
            if get_new_line:
                bline = file.readline()
                line = bline.decode("utf-8")
    return all_annotations

def get_annotations(taxid):
    print(f'Working on taxa {taxid}')
    url_pattern = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes'
    lineage = proteome_table['Taxonomic lineage'][taxid].split(', ')
    proteome_id = proteome_table['Proteome Id'][taxid]
    ## Viruses show first in lineage, all others show second.
    if lineage[0] == 'Viruses':
        lineage = 'Viruses'
    else:
        lineage = lineage[1]
    annotations_path = Path(PureWindowsPath(f'{proteome_id}_taxaid_{taxid}_annotations.json'))
    all_annotations = dict()
    if not (annotations_path.exists()):
        print(f'Fetching annotations for taxid {taxid}')
        url = f'{url_pattern}/{lineage}/{proteome_id}/{proteome_id}_{taxid}.dat.gz'
        try:
            request = requests.get(url, stream = True)
            start_time = time.time()
            all_annotations = parse_request(request)
            annotations_exist = True
        except:
            annotations_exist = False
            print(f'Annotations not found, likely {taxid} removed from reference proteomes')
        print(time.time()-start_time)

        additional_annotations = dict()
        try:
            url = f'{url_pattern}/{lineage}/{proteome_id}/{proteome_id}_{taxid}_additional.dat.gz'
            request = requests.get(url, stream = True)
            print('Fetching additional annotations')
            start_time = time.time()
            additional_annotations = parse_request(request)
            print(time.time()-start_time)
        except:
            print('No additional annotation to fetch')
        if additional_annotations:
            assert(len(set(additional_annotations.keys()).intersection(set(all_annotations.keys()))) == 0)
            for new_accession in additional_annotations.keys():
                all_annotations[new_accession] = additional_annotations[new_accession]

        if annotations_exist:
            with open(annotations_path, "w") as outfile: 
                json.dump(all_annotations, outfile) 
        
    else: 
        print(f'Annotation path {annotations_path} already exists')
        annotations_exist = True

    ## Verify consistency with proteome
    if annotations_exist:
        with open(annotations_path) as json_file:
            all_annotations = json.load(json_file)
        annotated_features = all_annotations.keys()

        print('Verifying consistency of annotations file and FASTA')
        proteome_path = Path(PureWindowsPath(f'{proteome_id}_taxaid_{taxid}.FASTA'))
        all_features = get_all_features(proteome_path)
        consistency_check = [x in annotated_features for x in all_features]
        if not all(consistency_check):
            print(f'Downloading updated proteome for taxid {taxid}')
            download_proteome(taxid)
            log_proteomes([taxid], [datetime.datetime.now()], [proteome_path], log_file)
            all_features = get_all_features(proteome_path)
            consistency_check = [x in annotated_features for x in all_features]
        assert(all(consistency_check))
        print('Passes consistency check')
    return f'{proteome_id}_{taxid}', annotations_exist

all_changed_proteomes = []
for taxid in list(proteome_table.index):
    proteome_taxa_id, annotations_exist = get_annotations(taxid)
    if not annotations_exist:
        all_changed_proteomes = all_changed_proteomes + [proteome_taxa_id]

with open('all_changed_proteomes.txt', 'w') as file:
    for changed_proteome in all_changed_proteomes:
        file.write(f'{changed_proteome}\n')

#     request = requests.get(url, stream = True)
#     start_time = time.time() 


# with gzip.open(request.raw) as file:
#     start_time = time.time() 
#     last_pos = 0
#     new_dict = dict()
#     all_annotations = dict()
#     new_annotation = dict()
#     bline = file.readline()
#     line = bline.decode("utf-8")
#     gn_done = False
#     while bline:
#         get_new_line = True
#         if line.startswith('//'):
#             all_annotations[accession] = new_annotation
#             new_annotation = dict()
#             gn_done = False
#         elif line.startswith('ID   '):
#             new_dict = grab_ID_info(line)
#             new_annotation['Identification (ID)'] = new_dict  
#             new_dict = dict()
#         elif line.startswith('AC   '):
#             if new_dict:
#                 new_dict = grab_accession(line, secondary_accessions = True, out_dict = new_dict)
#             else:
#                 new_dict = grab_accession(line, out_dict = dict())
#             accession = new_dict['primary_accession']
#             new_annotation['Accession (AC)'] = new_dict  
#         elif line.startswith('DT   ') and ', entry version' in line:
#             new_dict = grab_date(line)
#             new_annotation['Date (DT)'] = new_dict  
#         elif line.startswith('DE   '):
#             new_dict, file, line = grab_description(file, line)
#             new_annotation['Description (DE)'] = new_dict 
#             ## grab_description gives the next line
#             get_new_line = False  
#         elif line.startswith('GN   Name=') and not gn_done:
#             gn_done = True
#             new_dict = dict()
#             new_dict['official_uniprot_gene'] = line.split('GN   Name=')[1].split(';')[0].split(' {')[0]
#             new_annotation['Gene name (GN)'] = new_dict
#         elif line.startswith('CC   -!- '):
#             new_dict, file, line = grab_comments(file, line)
#             new_annotation['Comments (CC)'] = new_dict
#             get_new_line = False
#         elif line.startswith('DR   GO; '):
#             go_annotation = line.split('DR   GO; ')[1].split('; ')
#             new_dict = dict()
#             new_dict['GO_code'] = go_annotation[0]
#             new_dict['GO_evidence'] = go_annotation[2].split('\n')[0]
#             if 'Gene Ontology annotations (DR)' not in new_annotation.keys():
#                 new_annotation['Gene Ontology annotations (DR)'] = dict()
#             new_annotation['Gene Ontology annotations (DR)'][go_annotation[1]] = new_dict
#         elif line.startswith('DR   KEGG; '):
#             new_annotation['KEGG identifier (DR)'] = line.split('DR   KEGG; ')[1].split(';')[0]
#         elif line.startswith('DR   Reactome; '):
#             reactome_annotation = line.split('DR   Reactome; ')[1].split('; ')
#             new_dict = dict()
#             new_dict['Reactome_code'] = reactome_annotation[0]
#             if 'Reactome annotations (DR)' not in new_annotation.keys():
#                 new_annotation['Reactome annotations (DR)'] = dict()
#             new_annotation['Reactome annotations (DR)'][reactome_annotation[1]] = new_dict
#         elif line.startswith('KW   '):
#             new_keywords = line.split('KW   ')[1].split('\n')[0].split(';')
#             new_keywords = [x for x in new_keywords if x != '']
#             if 'Keywords (KW)' in new_annotation.keys():
#                 new_annotation['Keywords (KW)'] = new_annotation['Keywords (KW)'] + new_keywords 
#             else:
#                 new_annotation['Keywords (KW)'] = new_keywords 
#         if get_new_line:
#             bline = file.readline()
#             line = bline.decode("utf-8")
    
#     with open("sample.json", "w") as outfile: 
#         json.dump(all_annotations, outfile)
#     print(time.time()-start_time)

        
        
