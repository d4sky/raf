import os

def Get_files_from_dir_by_extension(wdir, ext, filterIN = ""):
  out_files = []
  for file_name in os.listdir(wdir):
    file_path = wdir + file_name
    if os.path.isfile(file_path):
      if file_name.endswith(ext): 
        if filterIN == "" or filterIN in file_path:
          out_files.append(file_path)
          
  return out_files