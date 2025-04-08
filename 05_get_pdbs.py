from _utils import *

def Process_outputs(table, threshold=0.0):
  for index, row in table.iterrows():
    score = row['score']
    if score > threshold:
      #print(row)
      query = row['query']
      match = row['match']
      sbjct = row['sbjct']
      beg   = row['beg']
      end   = row['end']
      mr = MatchRecord(score, query, match, sbjct, beg, end)
      #print(mr)
      pMatches.Add_item(row['PDB'], row['chains'], mr)

pMatches = Matches()

Table1 = PandasDF("raf_matches.xlsx")
Process_outputs(Table1, 20.0)

Table2 = PandasDF("raf_shorts.xlsx")
Process_outputs(Table2, 0.0)

#print(pMatches.Line_formated(col_spaces={"key1": 10, "key2": 5, "item": 20}))
'''
print(pMatches.Formated(
    record_fields=["score", "query", "sbjct", "beg", "end"],
    col_spaces={"PdbID": 6, "Chain": 4, "score": 6, "beg": 5, "end": 5}
))
'''
pdb_out_dir = "pdbs/"
Check_and_create_dir(pdb_out_dir)

pdb_url = 'https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb'
cif_url = 'https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF'

for uniquePDB in pMatches:
  print(uniquePDB)
  Download_pdb_entry(uniquePDB, pdb_out_dir, pdb_url, cif_url)
  break