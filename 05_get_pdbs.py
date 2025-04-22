import sys 

from _utils import *

#sys.path.append('C:/Python/')
#from nutils import PandasDF

def Add_matches(rowsDF):
  print(rowsDF)
  for index, row in rowsDF.iterrows():
    pdbID = row["PDB"]
    chain = row["chains"]
    score = row["score"]
    pMatches.Add_item(pdbID, chain, {"score": score})

pMatches = Matches("PdbID", "Chain")

Table1 = PandasDF("raf_matches.xlsx")
matches1 = Table1.Select_rows(['PDB', 'chains', 'score', 'query', 'match', 'sbjct', 'beg', 'end'], {'score': ('>', 100.0)})
Add_matches(matches1)

Table2 = PandasDF("raf_shorts.xlsx")
matches2 = Table2.Select_rows()
Add_matches(matches2)

got_pdbs = pMatches.Get_values(['score'])

unitedTable = PandasDF("raf_united.xlsx")
unitedTable.Add_columns(["PdbID", "Chain", "score"])
unitedTable.Append_rows(got_pdbs[1:])
unitedTable.Save()

pdb_out_dir = "pdbs/"
Check_and_create_dir(pdb_out_dir)

pdb_url = 'https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb'
cif_url = 'https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF'

for uniquePDB in pMatches:
  print(uniquePDB)
  if Download_pdb_entry(uniquePDB, pdb_out_dir, pdb_url, cif_url) != 'downloaded':
    break