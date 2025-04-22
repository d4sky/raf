import sys, os

from _utils import *

unitedTable = PandasDF("raf_united.xlsx")
unitedTable.Add_columns(['seq'])

localPdb_path = "F:/pdbs/infos/"

pdbTables = {}
for index, row in unitedTable.iterrows():
  pdbID = row['PdbID']
  chain = row["Chain"]
  sada = int(pdbID[0])
  print(pdbID, chain, sada)
  if not sada in pdbTables:
    print(f"PDB Table for set '{sada}' is not yet IN and is being downloaded now...")
    filename = f"s{sada}_chains.xlsx"  
    tablePath = os.path.join(localPdb_path, filename)  
    pdbTable = PandasDF(tablePath)
    pdbTable["pdb_id"] = pdbTable["pdb_id"].str.upper()
    pdbTable["chain"]  = pdbTable["chain"].str.upper()
    pdbTables[sada] = pdbTable
    print("... downloading finished")
      
  actPdbTable = pdbTables[sada]
  foundRows = actPdbTable[(actPdbTable["pdb_id"] == pdbID.upper()) & (actPdbTable["chain"] == chain.upper())]
  if len(foundRows) > 0:
    if len(foundRows) > 1:
      print(f"Found more than one sequence for PDB '{pdbID}' with Chain '{chain}'. Taken first one.")  
    foundRow = actPdbTable[(actPdbTable["pdb_id"] == pdbID.upper()) & (actPdbTable["chain"] == chain.upper())].iloc[0]
    foundSeq = foundRow['seq']
    print("Found sequence")
    print(foundSeq)
    if foundSeq != '':
      if len(str(row['seq'])) == 0 or str(row['seq']) == 'nan':
        unitedTable.Update_row(index, {'seq': foundSeq})
    else:
      print(f"No sequence found for PDB '{pdbID}' with Chain '{chain}'")  
    print()
  else:
    print(f"No sequence found'")  
 
unitedTable.Save() 
sys.exit()