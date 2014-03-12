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
    print name_ref, "****"
    
    for filein in l_filesin : 
        if search("^" + name_ref, filein) : 
            print "IN"
            print filein
            
            return p_dataset_folder + filein
    
    return 0
    
    
    
    
    
    
    
    
    


