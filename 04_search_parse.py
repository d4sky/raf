import pandas as pd

from _utils import *

def Parse_search_txt(txt_path, searchID, pdbTable):
    pdbchTable_ID = 0
    blastTable_ID = 0
    print(txt_path)
    fr = open(txt_path,'r')
    for line in fr:
      items = line.split()
      if len(items) == 7:
        pdb_id = items[0].upper()
        pdb_ch = items[1].upper()
        seq    = items[2]
        score  = items[3]
        beg = items[4]
        end = items[5]
        match = items[6]
        
        pdbchTable_ID += 1
        pdbTable['ID'].append(pdbchTable_ID)
        pdbTable['PDB'].append(pdb_id)
        pdbTable['chains'].append(pdb_ch)
        pdbTable['search_id'].append(f"search:{searchID}")

        pdbTable['score'].append(score)
        pdbTable['query'].append(seq)
        pdbTable['match'].append(match)
        pdbTable['sbjct'].append(match)
        pdbTable['beg'].append(beg)
        pdbTable['end'].append(end)
        
        length = int(end) - int(beg)
        
    fr.close()    

output_dir = "shorts_search/"
dom_prot = "raf"

pdbColumns = ['ID', 'PDB', 'chains', 'search_id', 'score', 'query', 'match', 'sbjct', 'beg', 'end']
pdbchTable = {column:[] for column in pdbColumns}

zoznam = Get_files_from_dir_by_extension(output_dir, "txt")
for i, txt_file in enumerate(zoznam):
  print(f"FILE {txt_file}")
  Parse_search_txt(txt_file, dom_prot, pdbchTable)
print("Finished Parsing")

df_pdbch = pd.DataFrame(pdbchTable)
df_pdbch.to_excel('./' + f'{dom_prot}_shorts.xlsx')

