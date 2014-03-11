'''
Path manage 10-03-2014
'''

from os import makedirs



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




def generatePath (path_directory):
    
    try : makedirs( path_directory, mode=0777 )
    except : pass
    
    return path_directory


