'''
Path manage 10-03-2014
'''

from os import makedirs, listdir, path
from re import search



globals()["dir_initial"] = "/home/borrel/Yue_project/"



def result ( dir_in="" ):
    """
    Create result directory
    args: Directory in dataSet
    return: path
    """
    
    dir = dir_initial + "result/" 
    try : makedirs( dir, mode=0777 )
    except : pass
    if dir_in != "" : 
        dir_in_dataSet = dir + dir_in + "/"
        try : makedirs( dir_in_dataSet, mode=0777 )
        except : pass
        return dir_in_dataSet
    
    return dir


def alignmentOutput ( dir_in="" ):
    """
    Create result directory
    args: Directory in dataSet
    return: path
    """
    
    dir = dir_initial + "alignment/" 
    try : makedirs( dir, mode=0777 )
    except : pass
    if dir_in != "" : 
        dir_in_dataSet = dir + dir_in + "/"
        try : makedirs( dir_in_dataSet, mode=0777 )
        except : pass
        return dir_in_dataSet
    
    return dir



def dataset ( dir_in="" ):
    """
    Create result directory
    args: Directory in dataSet
    return: path
    """
    
    dir = dir_initial + "dataset/" 
    try : makedirs( dir, mode=0777 )
    except : pass
    if dir_in != "" : 
        dir_in_dataSet = dir + dir_in + "/"
        try : makedirs( dir_in_dataSet, mode=0777 )
        except : pass
        return dir_in_dataSet
    
    return dir



def generatePath (path_directory):
    
    try : makedirs( path_directory, mode=0777 )
    except : pass
    
    return path_directory



def findPDBRef(p_dataset_folder) : 
    
    l_filesin = listdir(p_dataset_folder)
    
    name_ref = path.dirname(p_dataset_folder).split ("/")[-1]
    
    for filein in l_filesin : 
        if search("^" + name_ref, filein) and search(".pdb" , filein): 
            
            return p_dataset_folder + filein
    
    return 0

def findPDBQueryDataset (p_dataset_folder):
    
    l_filesin = listdir(p_dataset_folder)
    
    name_ref = path.dirname(p_dataset_folder).split ("/")[-1]
    
    l_out = []
    for filein in l_filesin : 
        if search("^" + name_ref, filein) : 
            continue
        elif search(".pdb" , filein) and len(filein.split ("_")[0]) == 4: 
            l_out.append (p_dataset_folder + filein)
            
    
    return l_out    
 
 
    
def findPDBQueryTransloc (p_result):
    
    try : l_filesin = listdir(p_result)
    except : return []    

    l_out = []
    for filein in l_filesin : 
        if search("^CX_", filein) : 
            l_out.append (p_result + filein)
            
    
    return l_out     
    
    
def findligandRef(p_dataset_folder, name_lig) :     
    
    l_filesin = listdir(p_dataset_folder)
    name_ref = path.dirname(p_dataset_folder).split ("/")[-1]
    for filein in l_filesin : 
        if search( name_ref, filein) and  search( name_lig, filein) and not search("subref", filein): 
            return p_dataset_folder + filein
    return 0
    
 
def findSubstructRef (p_dataset_folder, substruct): 
    
    l_out = []
    l_filesin = listdir(p_dataset_folder)
    name_ref = path.dirname(p_dataset_folder).split ("/")[-1]
    for filein in l_filesin : 
        if search(name_ref, filein) and search("^subref_", filein) and search(".pdb", filein): 
            l_out.append (p_dataset_folder + filein)
        
    return l_out
    
    
    
    
def findMatrix(p_lig_ref, p_lig, substruct) : 
    
    begin_name = p_lig_ref.split ("/")[-1][4:-4]
#     print begin_name
    end_name = p_lig.split ("/")[-1][4:-4]
#     print end_name
    
    return alignmentOutput( substruct + "/" + begin_name + "__" + end_name) + "matrix.out"
    
    
def findListSmileFile(substruct) : 
    
    l_out = []
    pr_result = result(substruct)
    
    l_file = listdir(pr_result)
    
    for name_file in l_file : 
        if search("smile", name_file) and search(".txt$", name_file): 
            l_out.append (pr_result + name_file)
    
            
    return l_out
    
    
def findFileBS(pr_result, PDB_query) : 
    
    l_out = []
    
    l_file_res = listdir(pr_result)
    
    for file_result in l_file_res : 
        if search ("BS", file_result)  and search (PDB_query, file_result): 
            l_out.append (pr_result + file_result)
    
    return l_out 
            
    
def findligandQuery(pr_dataset, name_ligand, PDB_query) : 
    
    l_filin = listdir(pr_dataset)
    
    for filin in l_filin : 
        if search (name_ligand, filin) and search (PDB_query, filin) and search (".pdb", filin) : 
            return pr_dataset + filin
    
    

def findSubstructFind(pr_result, name_ligand, PDB_query, substruct) : 
    
    l_filin = listdir(pr_result )
    
    for filin in l_filin : 
        if search (name_ligand, filin) and search ("^substituent", filin) and search (".pdb", filin) and search (PDB_query, filin) and search (substruct, filin) : 
            return pr_result + filin    
    
    
    
    
def findFamilyFile (name_lig): 
    
    pr_dataset = dataset (name_lig)
    return pr_dataset + "family_PDB.txt"
    
     
                
            
