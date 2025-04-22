import os
import pathlib
from collections import defaultdict

import random
import string
import subprocess

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import urllib.request 
import pandas as pd

import json

from dataclasses import dataclass

def to_regular_dict(d):
  if isinstance(d, defaultdict):
    d = {k: to_regular_dict(v) for k, v in d.items()}
  elif isinstance(d, list):
    d = [to_regular_dict(i) for i in d]
  return d

def modify_filename(filename, insLength=3):
    """
    Modify the filename according to the given rules.
    """
    if len(filename) < 5:
        S1 = ''
        S2 = filename
    else:
        S1 = filename[:-5]
        S2 = filename[-5:]

    inS = ''.join(random.choice(string.ascii_letters) for _ in range(insLength))
    if S1: inS = '_' + inS

    newfilename = S1 + inS + S2
    return newfilename
    
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

def Get_terminal_gaps(string):
  first = -1
  last = -1
  i = 0
  while i < len(string):
    if string[i] != '-': 
      first = i
      break
    i += 1
     
  i = len(string) - 1
  while i > 0:
    if string[i] != '-': 
      last = i
      break
    i -= 1
    
  return first + 1, last + 1 

def Muscle_without_files(seq1, seq2, gapOpen = -4.90, gapExtend = -0.00, muscle_exe= r"C:\python\execs\muscle3.8.31_i86win32.exe"):
    ##############################
    ## WITHOU INP/OUT FILES
    #seqs = [SeqRecord(Seq(m.bobject[chain].seq)),SeqRecord(Seq(sequences[protein][1]))]
    seqs = [SeqRecord(Seq(seq1)), SeqRecord(Seq(seq2))]
    #args = [muscle_exe,'-gapopen', '-4.90','-gapextend', '-0.00']
    args = [muscle_exe,'-gapopen', f'{gapOpen}','-gapextend', f'{gapExtend}']
    muscle  = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
        
    SeqIO.write(seqs, muscle.stdin, "fasta")
    muscle.stdin.close()
    aligned = AlignIO.read(muscle.stdout, 'fasta')
    muscle.stdout.close()

    return aligned
    
def Read_fasta(seq_path):
  h = ''
  s = ''
  with open(seq_path) as f:
    lines = [line.strip() for line in f]
    #lines = f.readlines()

  for line in lines:
    if line.startswith('>'):
      h = line[1:]
    else:
      s = s + line  
      
  return h, s

def Process_pairwise_SA(seq1, seq2):
  os1, os2 = '', ''
  if len(seq1) == len(seq2):
    for si, aa1 in enumerate(seq1):
      aa2 = seq2[si]      
      if    aa2 == '-':
        pass
      else:
        os1 = os1 + aa1  
        os2 = os2 + aa2  
  else:
    print("Sequences are not equally long")
  
  return os1, os2

class PandasDF(pd.DataFrame):
  def __init__(self, file_path=""):
    if file_path and os.path.exists(file_path):
      print(f"Reading Table from {file_path}")
      df = pd.read_excel(file_path, index_col=None, usecols=lambda x: x != 0).dropna(how='all')
      print("... finished!")
    else:
      df = pd.DataFrame()  # Initialize an empty DataFrame
      if file_path != "":
        print(f"File name {file_path} not found!")
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

  def Update_row(self, idx, update_row):
    for key, value in update_row.items():
      if key in self.columns:
        self.at[idx, key] = value
      else:
        print(f"There is no Column named '{key}'")        
        
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
      
  def Select_rows(self, out_columns = [], where_columns = {}):
    if out_columns:
      for wcol in out_columns:
        if not wcol in self.columns:
          print(f"There is no column '{wcol}' in Table specified as SELECT fields")
          return []
      
    if where_columns:
      for wcol in where_columns.keys():
        if not wcol in self.columns:
          print(f"There is no column '{wcol}' in Table usen in WHERE fields")
          return []
          
      taken_indexes = []          
      for index, row in self.iterrows():
        takeIt = 1
        for wcol, clause in where_columns.items():
          sign, value = clause
          if   sign == '=':
            if row[wcol] == value:
              takeIt = takeIt*1
            else:
              takeIt = takeIt*0
          elif sign == '<':
            if row[wcol] < value:
              takeIt = takeIt*1
            else:
              takeIt = takeIt*0
          elif sign == '>':
            if row[wcol] > value:
              takeIt = takeIt*1
            else:
              takeIt = takeIt*0
        if takeIt:
          taken_indexes.append(index)
    else:
       taken_indexes = [index for index, row in self.iterrows()]
    
    if out_columns == []:
      out_columns = list(self.columns)
    
    returning_rows = []
    for idx in taken_indexes:
      ret_row = [self.at[idx, wcol] for wcol in out_columns]
      returning_rows.append(ret_row)
    
    #return returning_rows    
    returningDF = PandasDF()
    returningDF.Add_columns(out_columns)
    returningDF.Append_rows(returning_rows)
    
    # Return the new DataFrame
    return returningDF

def default_factory():
  return defaultdict(list, )
      
class Matches(defaultdict):
  def __init__(self, key1, key2, iFields = ["score", "query", "match", "sbjct", "beg", "end"]):
    #super().__init__(lambda: defaultdict(list))

    super().__init__(default_factory)
    
    self.key1 = key1 
    self.key2 = key2 
    self.fields = iFields

    self.min_field = ''
    self.max_field = ''
    
    self.min_value = 0
    self.max_value = 0

  def iterator(self):
    # Iterate through each key1, key2, and record
    for key1 in self:
      for key2 in self[key1]:
        for record in self[key1][key2]:
          yield record  # Yielding each record

  def Populate_beg(self):
    for record in self.iterator():
      '''
      for ai, aa in enumerate(record['query']):
        if aa != '-': break
      fi = ai
      '''
      fi = next(ai for ai, aa in enumerate(record['query']) if aa != '-')
      record['beg'] = fi

  def Not_in_fields_msq(self, field):
    print(f"Key '{field}' is not one of these:")
    print(self.fields)

  def Set_min(self, field):
    if field in self.fields:
      self.min_field = field
    else:
      self.Not_in_fields_msq(rkey)

  def Set_max(self, field):
    if field in self.fields:
      self.max_field = field
    else:
      self.Not_in_fields_msq(rkey)
      
  def Set_min_max(self, min_field, max_field):
    self.Set_min(min_field)    
    self.Set_max(max_field)

  def calc_MinMax(self):
    if self.min_field and self.max_field:  
      '''
      begs = [record.beg for record in self.iterator()]
      ends = [record.end for record in self.iterator()]
      self.minX = min(begs)
      self.maxX = max(ends)
      '''
      mines = [record[self.min_field] for record in self.iterator()]
      maxes = [record[self.max_field] for record in self.iterator()]
      self.min_value = min(mines)
      self.max_value = max(maxes)

    else:
      print("Min and/or Max do(es) not seem to be set|")
    
  def Add_item(self, key1, key2, record_dict):
    row_record = {}
    for rkey, rval in record_dict.items():
      if rkey in self.fields:
        row_record[rkey] = rval
      else:
        self.Not_in_fields_msq(rkey)
        
      for ckey in self.fields:
        if ckey not in record_dict:
          row_record[ckey] = '' # For all items in record an empty string is as Default Value. Could be extended for various default values

    self[key1][key2].append(row_record)  

  def Get_values(self, record_fields=None, with_keys = True):
    if   record_fields is None:
      #record_fields = list(vars(MatchRecord()).keys())
      record_fields = self.fields

    out_lines = []
    if with_keys:
      #out_lines.append(["PdbID", "Chain"] + record_fields)
      # Taken from initialized attributes 
      out_lines.append([self.key1, self.key2] + record_fields)
    else:
      if record_fields:  
        out_lines.append(record_fields)
    for key1 in self:
      for key2 in self[key1]:
        for record in self[key1][key2]:
          if with_keys:
            #values = [key1, key2] + [getattr(record, field, "") for field in record_fields]
            values = [key1, key2] 
          else:
            #values = [getattr(record, field, "") for field in record_fields]
            values = []
          values = values + [record[field] for field in record_fields]
          if values: 
            out_lines.append(values)

    return out_lines

  def As_lines(self, which_fields=None, col_spaces=None, def_spaces=15, also_keys = True):
    if which_fields:
      for field in which_fields:
        if not field in self.fields:
          self.Not_in_fields_msq(field)
          return ''
        
    got_lines = self.Get_values(which_fields, also_keys)   
      
    if col_spaces is None:
      col_spaces = {}

    lines = []
    # optional header
    headers = got_lines[0]
    header_line = " ".join([f"{h:<{col_spaces.get(h, def_spaces)}}" for h in headers])
    lines.append(header_line)
    lines.append("-" * len(header_line))
    for values in got_lines:
      formatted = []
      for h, val in zip(headers, values):
        width = col_spaces.get(h, def_spaces)
        formatted.append(f"{str(val)[:width]:<{width}}")
      lines.append(" ".join(formatted))

    return "\n".join(lines)
    
  def As_tree(self):
    output = []
    for key1 in self:
      output.append(f"{key1}")
      for key2 in self[key1]:
        records = self[key1][key2]
        output.append(f" {key2}")
        for record in records:
          for ikey, ival in record.items():
            output.append(f"  {ival}")
    return "\n".join(output)
    
  '''
  def As_aligned(self, head_fields={}):
    # head_fields - dictionary with keys as field name to be showen at the beginning of each lines with values for the number of space for each field at output
    #self.calc_MinMax()
    lines = []
    for pdbID in self:
      for Chain in self[pdbID]:
        for record in self[pdbID][Chain]:
          #line = f"{pdbID} {Chain:{5}} {record.beg:{5}} {record.end:{5}} "
          line = f"{pdbID} {Chain:{5}} "
          for field, fSpaces in head_fields.items():
            line = line + f"{record[field]:{fSpaces}} "
          #line = line + '-'*(record.beg-1) + record.query + '-'*(self.maxX - record.end)
          #!!!!
          #line = line + '-'*(record['beg'] - 1) + record['query'] + '-'*(self.max_value - record['end'])
          line = line + record['query'] 
          lines.append(line)
          
    #return lines
    return "\n".join(lines)
  '''
  
  def As_aligned(self, head_fields={}):
    # head_fields - dictionary with keys as field name to be showen at the beginning of each lines with values for the number of space for each field at output
    #self.calc_MinMax()
    lines = []
    '''
    for pdbID in self:
      for Chain in self[pdbID]:
        for record in self[pdbID][Chain]:
          #line = f"{pdbID} {Chain:{5}} {record.beg:{5}} {record.end:{5}} "
          line = f"{pdbID} {Chain:{5}} "
          for field, fSpaces in head_fields.items():
            line = line + f"{record[field]:{fSpaces}} "
          #line = line + '-'*(record.beg-1) + record.query + '-'*(self.maxX - record.end)
          #!!!!
          #line = line + '-'*(record['beg'] - 1) + record['query'] + '-'*(self.max_value - record['end'])
          line = line + record['query'] 
          lines.append(line)
    '''
    
    got_records = self.Get_values(['beg', 'query'])
    #for record in got_records:
    #  print(record[0], record[2])
    
    sorted_recs = sorted(got_records[1:], key=lambda x: x[2])      
    for record in sorted_recs:
      line = f"{record[0]} {record[1]:{5}} {record[3]}"
      lines.append(line)
          
    #return lines
    return "\n".join(lines)

          
  def As_json(self):
    data = {
      "key1": self.key1,
      "key2": self.key2,
      "fields": self.fields,
      "min_field": self.min_field,
      "max_field": self.max_field,
      "min_value": self.min_value,
      "max_value": self.max_value,
      "data": to_regular_dict(dict(self))  # convert nested defaultdicts
    }
    out_str = json.dumps(data, indent=2)
        
    return out_str
    
  def Write(self, as_what, file_path = ''):
    if file_path:
      output = ''
      if   as_what == "lines":
        output = self.As_lines()
      elif as_what == "tree":
        output = self.As_tree()
      elif as_what == "aligned":
        output = self.As_aligned()
      elif as_what == "json":
        output = self.As_json()
      else:
        print(f"No valid option '{as_what}'")  
      
      if output:
        with open(file_path, 'w') as f:
          f.write(output)    
          
    else:
      print("No path to file is specified")
      
  @classmethod
  def from_json(cls, file_path):
    with open(file_path, 'r') as f:
      data = json.load(f)
    
    obj = cls(data["key1"], data["key2"], data["fields"])
    obj.min_field = data.get("min_field", "")
    obj.max_field = data.get("max_field", "")
    obj.min_value = data.get("min_value", 0)
    obj.max_value = data.get("max_value", 0)

    # Convert back nested dict structure
    for k1, subdict in data["data"].items():
      for k2, records in subdict.items():
        obj[k1][k2] = records  # list of dicts
    
    return obj
    