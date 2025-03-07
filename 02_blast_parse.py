import sys
import pandas as pd

from Bio.Blast import NCBIXML

sys.path.append('../')
from _utils import *

class PdbChain():
  def __init__(self, inp_string):
    self.pdb    = {}
    self.chains = {}
    print(inp_string)
    items = inp_string.split('|')
    print(items)
    for ii, item in enumerate(items):
      if (ii%2)==1 and ii < 3: 
        pdbID = item.strip()
        if pdbID not in self.pdb: self.pdb[pdbID] = []
        
        if (ii+1) > len(items):
          print("There seems to be no Chain in for PDB:" + pdbID)
          self.pdb[pdbID].append("N@")
        else:
          if (len(items[ii+1].strip()))==0:
            print("There seems to be no Chain in for PDB:" + pdbID)
            self.pdb[pdbID].append("N@")
          else:
            self.pdb[pdbID].append(items[ii+1].strip()[0])
            
    for pdbID in self.pdb:
      if   len(self.pdb[pdbID])>0:
        self.chains[pdbID] = ','.join(self.pdb[pdbID])
      else:
        print("Something is wrong")
        
  def Report(self):
    out_string = ""
    for pdbID in self.pdb:
      if len(out_string)==0: out_string = out_string + pdbID
      else:                  out_string = out_string + ";" + pdbID
      
      out_string = out_string + ":" + self.chains[pdbID]
    return out_string

def Parse_blast_xml(fname, searchID, headers, blTable, pdbTable):
  blastTable_ID = 0
  pdbchTable_ID = 0
  
  result_handle = open(fname)
  blast_records = NCBIXML.parse(result_handle)
  for blast_record in blast_records:
    for ai, alignment in enumerate(blast_record.alignments):
      pdb_chain   = PdbChain(alignment.title)
      pdbC_string = pdb_chain.Report()
        
      for hsp in alignment.hsps:  
        blastTable_ID = blastTable_ID + 1
        blTable['ID_blast_search'].append(blastTable_ID)
        blTable['SH3'].append(searchID)
        blTable['title'].append(alignment.title)
        blTable['PDB_Chains'].append(pdbC_string)
        pdbItems = pdbC_string.split(':')
        
        blTable['length'].append(alignment.length)
        for header in headers: 
          blTable[header].append(getattr(hsp, header))
          
        # INSERTING PDBs
        for pdb_id in pdb_chain.pdb:
          pdbchTable_ID += 1
          pdbTable['ID'].append(pdbchTable_ID)
          pdbTable['PDB'].append(pdb_id)
          pdbTable['chains'].append(pdb_chain.chains[pdb_id])
          pdbTable['search_id'].append(blastTable_ID)

          pdbTable['score'].append(getattr(hsp, 'score'))
          pdbTable['query'].append(getattr(hsp, 'query'))
          pdbTable['match'].append(getattr(hsp, 'match'))
          pdbTable['sbjct'].append(getattr(hsp, 'sbjct'))
          pdbTable['beg'].append(getattr(hsp, 'sbjct_start'))
          pdbTable['end'].append(getattr(hsp, 'sbjct_end'))

blastColumns = ['ID_blast_search', 'SH3', 'title', 'PDB_Chains', 'length']
blastHeaders = ['score', 'bits', 'expect', 'identities', 'positives', 'gaps', 'align_length', 'query', 'query_start', 'query_end', 'match', 'sbjct', 'sbjct_start', 'sbjct_end']
blastTable = {column:[] for column in blastColumns + blastHeaders}

pdbColumns = ['ID', 'PDB', 'chains', 'search_id', 'score', 'query', 'match', 'sbjct', 'beg', 'end']
pdbchTable = {column:[] for column in pdbColumns}

blast_xml_dir = "pdb_blast/"

zoznam = Get_files_from_dir_by_extension(blast_xml_dir, "xml")
for i, xml_file in enumerate(zoznam):
  print(f"FILE {xml_file}")
  prefix = xml_file.split('/')[-1].replace(".xml", "")
  Parse_blast_xml(xml_file, prefix, blastHeaders, blastTable, pdbchTable)
  print()

print("Finished Parsing")

df_blast = pd.DataFrame(blastTable)
df_blast.to_excel('./' + 'searches.xlsx')

df_pdbch = pd.DataFrame(pdbchTable)
df_pdbch.to_excel('./' + 'matches.xlsx')


