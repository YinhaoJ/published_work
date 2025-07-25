import pandas as pd
import numpy as np
from tqdm import tqdm, trange

def load_and_process_SS_data(file_path):
    # Load txt file
    with open(file_path, 'r', encoding='utf8') as file:
        lines = file.readlines()
        
    # Skip the first 9 rows
    lines = lines[9:]
    
    # Convert into DataFrame
    df = pd.DataFrame([line.split() for line in lines])
    
    # Import data into DataFrame and drop unnecessary data
    df.columns = ["residue#", "chain_identifier", "empty", "frame", "secondary_structure"]
    SS = df.drop(['chain_identifier', 'empty'], axis=1)
    SS['residue#'] = pd.to_numeric(SS['residue#'])
    
    # Calculate start_residue and end_residue
    start_residue = SS['residue#'].min()
    end_residue = SS['residue#'].max()
    
    # Sort for secondary structure data
    Lr = []
    for i in trange(start_residue, end_residue + 1):
        Lr_temp = SS.loc[SS["residue#"] == i]['secondary_structure'].reset_index(drop=True)
        Lr.append(Lr_temp)

    # Convert into DataFrame
    ss = pd.DataFrame(Lr)
    ss = ss.set_index(pd.Index(range(start_residue, end_residue + 1)))
    
    return ss

ss_dict={}
def SimplifySecStructure(ssl):
  """
  This function takes in a list ssl of standard dssp secondary structure letters (G,H,I,E,B,T,C,S)
  and returns a simplified list where each element is assigned to the broader category of helix, beta or coil (H,E,C).
  So it basically maps H to H; E to E; H,I,B,T,C and S to C.
  """
  ss_dict={}
  for el in ['H']:ss_dict[el]='H'
  for el in ['E']:ss_dict[el]='E'
  for el in ['T','C','S','I','G','B']:ss_dict[el]='C'

  return [ss_dict[ss] for ss in ssl]

def SimplifySecStructureList(ss_list):
  """
  This function takes a list of lists of secondary structures as obtained from AnalyzeSecondaryStructure
  and reassigns each element to the broader categories of Helix, Extended and Coil (see SimplifySecStructure)
  Inputs:
  ss_list: A list of lists of secondary structures as obtained from AnalyzeSecondaryStructure, see SecStrucColorMap.
  Outpu:
  ss_list: Simplified list of lists of secondary structures
  """
  return [SimplifySecStructure(ssl) for ssl in ss_list]

def process_and_map_data(ss_LYZ_500K, ss_label_mapping, start=0, end=10001, step=10, row_limit=129):
    """
    Processes and maps secondary structure data from ss_LYZ_500K.
    
    Args:
        ss_LYZ_500K (list): List of lists containing secondary structure data.
        ss_label_mapping (dict): Mapping of simplified secondary structures to final labels.
        start (int): Starting index.
        end (int): Ending index.
        step (int): Step size for iteration.
        row_limit (int): Maximum number of rows to process.
        
    Returns:
        np.ndarray: A 2D numpy array representation of the processed data.
    """
    
    # List to store processed records
    ss_dict1 = []

    for i in range(start, end, step):
        # Simplify and flatten the structure
        ss_dict_temp2 = []
        ss_dict_temp = SimplifySecStructureList(ss_LYZ_500K[i])
        
        # Flattening the list of lists
        for sublist in ss_dict_temp:
            ss_dict_temp2.extend(sublist)
        
        # Mapping labels
        ss_dict_temp2 = [ss_label_mapping.get(item, item) for item in ss_dict_temp2]
        
        ss_dict1.append(ss_dict_temp2)
        
    # Convert into DataFrame
    ss_full = pd.DataFrame(ss_dict1)

    # Ensure the index aligns with the desired range
    ss_full.index = range(start, end, step)
    
    # Transpose and convert to a numpy array
    ss_map2 = ss_full.transpose().to_numpy()

    return ss_map2
