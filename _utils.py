import os
import pathlib
from collections import defaultdict

import urllib.request 
import pandas as pd

def Get_files_from_dir_by_extension(wdir, ext, filterIN = ""):
  out_files = []
  for file_name in os.listdir(wdir):
    file_path = wdir + file_name
    if os.path.isfile(file_path):
      if file_name.endswith(ext): 
        if filterIN == "" or filterIN in file_path:
          out_files.append(file_path)
          
  return out_files

def Check_and_create_dir(inp_path):
    pathlib.Path(inp_path).mkdir(parents=True, exist_ok=True) 

def Read_fasta_huge(file_name, with_spaces = True, what_items = [1], sep = '|'):
  print("Reading sequences from '" + file_name + "'")
  out={}
  try:
    first_test = 0
    fr = open(file_name,'r')
    si = 0   #counter for sequences
    sid="QUATSCH"
    seq=""
    head=""
            
    for line in fr:
      line = line.strip()
      if len(line)>0:
        if (line[0]=='>'):
          first_test = 1
          '''  should be initiated with another number if the sequence does not start with 1 '''
          if (si>0): 
            if len(seq)>0:	    
              if sid in out:
                print("Sequence ID '" + sid + "' is already in the pool. Trying another ID")
              else:
                out[sid]=(head,seq)
                
              seq=""              
          si = si + 1
    	
          items=line[1:].split('|')
          if with_spaces:  head = line
          else:            head = line.replace(' ','_')
          
          #if len(items)>1:
          if len(items) > max(what_items):
            '''
            if with_spaces:  sid = items[what_item]
            else:            sid = items[what_item].replace(' ','_')
            '''
            sid = sep.join([items[idx] for idx in what_items]) 
            if not with_spaces: sid = sid.replace(' ','_')
          else:
            if with_spaces:  sid = line[1:]
            else:            sid = line[1:].replace(' ','_')
  	  
        else:
          if first_test == 0:
            print("FASTA file '" + file_name.strip() + "' does not begin with '>'")
            return False
          else:  
            #!!! CHECK FOR GAPs '-' !!!
            seq = seq + line  
     	      
    if len(seq)>0:	    
      out[sid]=(head,seq)
	  
    fr.close()
    print("...Finished")
    print() 
    return out

  except IOError:
    print("File '" + file_name.strip() + "' not found")
    return False
  
def Download_pdb_entry(rID, out_dir, from_pdb_url, from_cif_url):
  def Actual_download(inp_url, out_path):
    try:
      urllib.request.urlretrieve(inp_url, out_path)
      print(f"Downloaded '{out_path}' successfully.")
      return True
    except Exception as err:
      print(f"Failed to retrieve PDB IDs: {rID} {err}")
      print("URL:", inp_url)
      return False
  
  pdbEntry = f"pdb{rID.lower()}.ent.gz"
  cifEntry = f"{rID.lower()}.cif.gz"
  
  url_for_pdb = f"{from_pdb_url}/{rID[1:3].lower()}/{pdbEntry}"
  url_for_cif = f"{from_cif_url}/{rID[1:3].lower()}/{cifEntry}"
  out_pdb_path = os.path.join(out_dir, pdbEntry)
  out_cif_path = os.path.join(out_dir, cifEntry)
  
  if os.path.exists(out_pdb_path):
    print(f"File to '{out_pdb_path}' already downloaded.")
    return 'downloaded'
  else:
    if os.path.exists(out_cif_path):
      print(f"File to '{out_cif_path}' already downloaded.")
      return 'downloaded'
    else:
      if Actual_download(url_for_pdb, out_pdb_path):
        return 'pdb'
      else:
        if Actual_download(url_for_cif, out_cif_path):
          return 'cif'
        else:
          return 'err'
  
class PandasDF(pd.DataFrame):
  def __init__(self, file_path=""):
    if file_path and os.path.exists(file_path):
      print(f"Reading Table from {file_path}")
      df = pd.read_excel(file_path, index_col=None, usecols=lambda x: x != 0).dropna(how='all')
      print("... finished!")
    else:
      print(f"File name {file_path} not found!")
      df = pd.DataFrame()  # Initialize an empty DataFrame
      print("Created empty DataFrame!")
    
    super().__init__(df)  # Initialize the parent class (pd.DataFrame)
    self.__finalize__(df)  # Ensure new objects keep class attributes
    
    self.path = file_path
    
  def Save(self, new_path = ''):  
    if new_path:
      out_path = new_path
    else:
      out_path = self.path

    if out_path:          
      try:
        print(f"GOING TO SAVE {out_path}")  
        #print(self)
        self.to_excel(out_path, index=False)
        print(f"DataFrame saved successfully to {out_path}")
      except PermissionError:
        new_path = modify_filename(out_path)
        self.to_excel(new_path, index=False)
        print(f"Error: The file {out_path} is open. Saved to a new file: {new_path}")
      except OSError as e:
        print(f"Error: Unable to save the file to {out_path}. Details: {e}")  
        
  def Add_columns(self, columns, with_save = True):
    if isinstance(columns, dict):
      # Add new columns with default values from the dictionary
      for col_name, default_value in columns.items():
        if col_name not in self.columns:
          self[col_name] = default_value
        else:
          print(f"Column with name '{col_name}' already exists!")
    elif isinstance(columns, list):
      # Add new columns with None values
      for col_name in columns:
        if col_name not in self.columns:
          self[col_name] = None
        else:
          print(f"Column with name '{col_name}' already exists!")
    else:
      raise TypeError("The 'columns' parameter must be either a dictionary or a list.")

    if with_save:
      self.Save()
  
  def Append_rows(self, inp_rows, with_save = True):
    if isinstance(inp_rows, list):
      # Inserting values according to order in the list
      Ncol = len(self.columns)
      if len(inp_rows) > 0:
        empty_cols = Ncol - len(inp_rows[0])
        if empty_cols > 0:
          rows = [row + [None] * empty_cols for row in inp_rows]
          #rows = [row + [''] * empty_cols for row in inp_rows]
        else:
          rows = [row[:Ncol] for row in inp_rows]
        
        excel_List = self.values.tolist()
        #print("BEFORE ", len(excel_List))
        for row in rows:
          excel_List.append(row)
            
        #print("AFTER ", len(excel_List))
        #self = pd.DataFrame(excel_List, columns=self.columns) # DOES NOT WORK !!!!!
        #self[:] = pd.DataFrame(excel_List, columns=self.columns) # DOES NOT WORK !!!!!
        #self.__init__(pd.DataFrame(excel_List, columns=self.columns)) # DOES NOT WORK !!!!!
        self._update_inplace(pd.DataFrame(excel_List, columns=self.columns))
        #print("APPENDED DATA FRAME")
        #print(self)
      else:
        print("Input Rows are empty")
    elif isinstance(inp_rows, dict):
      # Inserting values from dictionary, with keys as column names
      print("DICT to be implemented yet")
    else:
      print("Unsupported DataType/Format for Rows")
    
    if with_save:
      self.Save()
      
class MatchRecord():
  def __init__(self, score=0.0, query="", match="", sbjct="", beg=0, end=0):
    self.score = score
    self.query = query
    self.match = match
    self.sbjct = sbjct
    self.beg   = beg
    self.end   = end
      
class Matches(defaultdict):
  def __init__(self):
    super().__init__(lambda: defaultdict(list))
    
  def Add_item(self, key1, key2, row_record):
    self[key1][key2].append(row_record)  

  def Formated(self, record_fields=None, col_spaces=None, def_spaces=15):
    """
    record_fields: list of attribute names from MatchRecord to include (e.g. ["score", "query"])
    col_spaces: dict specifying column widths (keys can be "key1", "key2", or any record field)
    def_spaces: default spacing if not specified in col_spaces
    """
    if record_fields is None:
      record_fields = ["score", "query", "match", "sbjct", "beg", "end"]
    if col_spaces is None:
      col_spaces = {}

    lines = []
    # optional header
    headers = ["PdbID", "Chain"] + record_fields
    header_line = " ".join([f"{h:<{col_spaces.get(h, def_spaces)}}" for h in headers])
    lines.append(header_line)
    lines.append("-" * len(header_line))
    for key1 in self:
        for key2 in self[key1]:
            for record in self[key1][key2]:
                values = [key1, key2] + [getattr(record, field, "") for field in record_fields]
                formatted = []
                for h, val in zip(headers, values):
                    width = col_spaces.get(h, def_spaces)
                    formatted.append(f"{str(val)[:width]:<{width}}")
                lines.append(" ".join(formatted))

    return "\n".join(lines)
    
  def Line_formated(self, col_spaces={}, def_spaces=20):
    lines = []
    for key1 in self:
      for key2 in self[key1]:
        for item in self[key1][key2]:
          row = {
            "key1": key1,
            "key2": key2,
            "item": item
          }
          formatted = []
          for col in ["key1", "key2", "item"]:
            width = col_spaces.get(col, def_spaces)
            value = str(row[col])[:width]
            formatted.append(f"{value:<{width}}")
          lines.append(" ".join(formatted))

    return "\n".join(lines)
    
  def Tree_formated(self):
    output = []
    for key1 in self:
      output.append(f"[{key1}]")
      for key2 in self[key1]:
        items = self[key1][key2]
        output.append(f"  {key2}: {len(items)} item(s)")
        for item in items:
          output.append(f"    - {item}")
    return "\n".join(output)