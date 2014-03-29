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
        if search("^" + name_ref, filein) : 
            
            return p_dataset_folder + filein
    
    return 0
    
    
def findligandRef(p_dataset_folder, substruct) :     
    
    l_filesin = listdir(p_dataset_folder)
    name_ref = path.dirname(p_dataset_folder).split ("/")[-1]
    for filein in l_filesin : 
        if search( name_ref, filein) and  search( substruct, filein): 
            return p_dataset_folder + filein
    return 0
    
 
def findSubstructRef (p_dataset_folder, substruct): 
    
    l_filesin = listdir(p_dataset_folder)
    name_ref = path.dirname(p_dataset_folder).split ("/")[-1]
    for filein in l_filesin : 
        if search(name_ref, filein) and search("^subref_", filein): 
            return p_dataset_folder + filein
        
    return 0
    
    
    
    
def findMatrix(p_lig_ref, p_lig) : 
    
    begin_name = p_lig_ref.split ("/")[-1][4:-4]
#     print begin_name
    end_name = p_lig.split ("/")[-1][4:-4]
#     print end_name
    
    return alignmentOutput( begin_name + "__" + end_name) + "matrix.out"
    
    
def findListSmileFile(substruct) : 
    
    l_out = []
    pr_result = result(substruct)
    
    l_file = listdir(pr_result)
    
    for name_file in l_file : 
        if search("smile.txt", name_file) : 
            l_out.append (pr_result + name_file)
    
            
    return l_out
    
    
def findFileBS(pr_result, PDB_query) : 
    
    l_out = []
    
    l_file_res = listdir(pr_result)
    
    for file_result in l_file_res : 
        if search ("BS", file_result)  and search (PDB_query, file_result): 
            l_out.append (pr_result + file_result)
    
    return l_out 
            
            
            
            
            
    
    
    
    
    
    print    file_ref[0:-4]   
        
    
    
    


