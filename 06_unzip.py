import os, sys
import urllib.request

import gzip

from _utils import *

def unzip_gz_file(directory, inp_file):
    out_file = inp_file[:-3]  # Remove the .gz extension
    
    if out_file.endswith(".ent"):
      out_file = out_file.replace(".ent", ".pdb")

    if out_file.startswith("pdb"):
      out_file = out_file[3:]

    inp_path = os.path.join(directory, inp_file)
    out_path = os.path.join(directory, out_file)

    #print(inp_path, out_path)
    '''
    with gzip.open(input_file, 'rb') as f_in:
      with open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"Unzipped: {input_file} to {output_file}")
    return True
    '''
    
    try:
      # Unzip the file 
      with gzip.open(inp_path, 'rt') as f_in:  # Open the input file in text mode
        with open(out_path, 'wt') as f_out:  # Open the output file in text mode
          for line in f_in:
            f_out.write(line)    
                
      print(f"Unzipped: '{inp_path}' to '{out_path}'")
      return True
      
    except:
      print(f"Unzipping of '{inp_path}' to '{out_path}' failed!")
      return False

pdb_dir = "pdbs/"

file_names = os.listdir(pdb_dir)
for file_name in file_names:
  unzip_gz_file(pdb_dir, file_name)

sys.exit()