import sys, os
from Bio.Blast.Applications import NcbiblastpCommandline

sys.path.append('../')
from _utils import *

def Run_blast(db_dir, inpf, outf):
  myblast = NcbiblastpCommandline(query=inpf, db=db_dir, evalue=0.001, outfmt=5, out=outf)
  myblast()

#pdbaa_dir = "pdbaa_20230822"
#pdbaa_dir = "pdbaa_20250307"
##################################
##### !!!!!!!!!!!!!!!
##### Needs to be run from within PDB AA database-directory. e.i:
##### 1. change directory to pdbaa_20230822
##### 2. from there run the script as "../01_blast_run.py"
##### !!!!!!!!!!!!!!!
##################################

pdbaa_dir = "pdbaa"
inp_seq_dir = "../seqs/"
blast_xml_dir = "../pdb_blast/"

file_list = Get_files_from_dir_by_extension(inp_seq_dir, "fas")
for i, inp_seq_file in enumerate(file_list):
  print(inp_seq_file)
  seqID = inp_seq_file.split('/')[-1].replace(".fas", "")
  out_xml_file = blast_xml_dir + seqID + ".xml"
  print(i+1, inp_seq_file, out_xml_file)
  Run_blast(pdbaa_dir, inp_seq_file, out_xml_file)
