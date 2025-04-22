import sys, os

from _utils import *

#pMatches = Matches("PdbID", "Chain", ["score", "query", "match", "sbjct", "beg", "end"])
pMatches = Matches("PdbID", "Chain", ["query", "beg", "end"])

ref_seq_path = "seqs/RAF1_HUMAN_P04049.fas"
unitedTable = PandasDF("raf_united.xlsx")

_, ref_seq = Read_fasta(ref_seq_path)

thr = 400
open_gap_penalty = -5.9
extend_gap_penalty = -1.0

for index, row in unitedTable.iterrows():
  if row['score'] > thr:
    pdbID   = row['PdbID']
    chain   = row["Chain"]
    act_seq = row["seq"]
    if act_seq != '' and pd.notna(act_seq):
      tmp_alignment = Muscle_without_files(act_seq, ref_seq, open_gap_penalty, extend_gap_penalty)
      as1 = ''.join([aa for aa in tmp_alignment[0]])
      as2 = ''.join([aa for aa in tmp_alignment[1]])
      os1, os2 = Process_pairwise_SA(as1, as2)
      #pMatches.Add_item(pdbID, chain, {"query": os1, "beg": beg1, "end": end1})
      pMatches.Add_item(pdbID, chain, {"query": os1})

#print(pMatches.As_lines())
#print(pMatches.As_tree())
#print(pMatches.As_json())

pMatches.Write('json', "matches.json")
    
