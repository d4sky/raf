import sys, os, re
import subprocess

from difflib import SequenceMatcher

from _utils import *

def Process_seqs(seqs):
  # Creating the dictionary of comparing sequences
  out_seqs = {}
  for sID, entry in seqs.items():
    header, seq = entry
    items = sID.split("mol:")
    if len(items) > 0:
      if len(seq) >= lenTol:
        out_seqs[items[0]] = seq
  return out_seqs

def Run_shorts_search(inpf, outf):
  ref_prot = Read_fasta_huge(inpf)
  if len(ref_prot):
    pID = next(iter(ref_prot))  # Get the first Protein
    header, ref_seq = ref_prot[pID]
  else:
    print("Only one Ref protein is feasible")
    sys.exit()

  out_lines = []
  for sID, short_seq in comp_seqs.items():
    if short_seq in ref_seq:
      begends = [(m.start(), m.end()) for m in re.finditer(short_seq, ref_seq)]
      if len(begends) > 0:
        #out_lines.append(f"{sID} {short_seq} {begends[0][0]} {begends[0][1]}")
        #out_lines.append(f"{sID.replace('_',' ')} {short_seq} 1000.0 {begends[0][0]} {begends[0][1]}")
        matching_sequence = ref_seq[begends[0][0]:begends[0][1]]
        out_lines.append(f"{sID.replace('_',' ')} {short_seq} 1000.0 {begends[0][0]} {begends[0][1]} {matching_sequence}")
        if (len(begends)) > 1:
          print(f"More then one occurence of '{short_seq}' found!")
    #else:
    matcher = SequenceMatcher(None, ref_seq, short_seq)
    if matcher.ratio() > 0.03:
      match = matcher.find_longest_match(0, len(ref_seq), 0, len(short_seq))
      matching_sequence = ref_seq[ match.a: match.a + match.size]

      if match.size > 5:
        print(sID, matcher.ratio(), match.a, match.size)
        out_lines.append(f"{sID.replace('_',' ')} {short_seq} {matcher.ratio()} {match.a} {match.a + match.size} {matching_sequence}")
      
      '''
      if "SVQIVYK" in short_seq:
        print(matcher.ratio(), match.a, match.size)
        sys.exit()
      '''
           
  with open(outf, "w") as file:
    file.writelines(line + "\n" for line in out_lines)
  
print("STARTING")
pdb_short_seqs_path = "pdb_seqres_u20.fas"
inp_seq_dir = "seqs/"

lenTol = 5 # Smalles lenght of sequence to be considered

output_dir = "shorts_search/"

pdb_seqs  = Read_fasta_huge(pdb_short_seqs_path)
comp_seqs = Process_seqs(pdb_seqs)

Check_and_create_dir(output_dir)

file_list = Get_files_from_dir_by_extension(inp_seq_dir, "fas")
for i, inp_seq_file in enumerate(file_list):
  print(inp_seq_file)
  seqID = inp_seq_file.split('/')[-1].replace(".fas", "")
  out_search_file = output_dir + seqID + ".txt"
  print(i+1, inp_seq_file, out_search_file)
  Run_shorts_search(inp_seq_file, out_search_file)
