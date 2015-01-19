from os import system, listdir, makedirs
from re import search, compile, findall
from os.path import isfile, exists, getsize
from shutil import copy
import parsePDB




def formatFilePDB(dir_PDB):
    """Manage global PDB database
    manage PDB database, decompress the file and move and rename file with .pdb extension"""
    listRepertory = listdir(dir_PDB)
    
    for subRepertory in listRepertory :
        filesPDB = listdir(dir_PDB + subRepertory)   
        for file in filesPDB :
            pathFile = dir_PDB + subRepertory + '/' + file
            cmdDecompress = "gunzip " + pathFile
            system(cmdDecompress)  ####run decompress
            pathFile = pathFile[0:-3]
            namePDBFile = pathFile[-8:-4] + '.pdb'
            pathFileOut = pathFile[0:-14] + namePDBFile
            cmdMove = "mv " + pathFile + ' ' + pathFileOut
            system(cmdMove)  ####run move file
    
    for repertoryPDB in listRepertory:
        repertoryPDB = dir_PDB + repertoryPDB
        cmdRemove = "rm -r " + repertoryPDB
        system(cmdRemove)  ####run remove repertory



def searchLigands(pr_init, pr_PDB):
    '''search ligands in PDB database
    out : list of ligands with PDB files associated'''
    
    d_stock = {}
    
    # control file exist
    if exists(pr_init + "resultLigandInPDB") and getsize(pr_init + "resultLigandInPDB") != 0: 
        return pr_init + "resultLigandInPDB"
    
    # import list PBD from file .dat
    l_p_PDB = retriveListPDB(pr_PDB)
    
    for p_PDB in l_p_PDB:
        namePDB = p_PDB[-8:-4]
        l_ligand = []
        try : 
            filinPDB = open (p_PDB, "r")
            l_linesPDB = filinPDB.readlines()
            filinPDB.close ()
        except : continue

        for linePDB in l_linesPDB:
            if(search ("^HETATM", linePDB)):
                atom = parsePDB.lineCoords(linePDB)
                if not atom["resName"] in l_ligand:
                    l_ligand.append(atom["resName"])

        d_stock[namePDB] = l_ligand

    filout = open (pr_init + "resultLigandInPDB", "w")
    
    for pdb in d_stock.keys ():
        filout.write(pdb + "\t")
        if d_stock[pdb] != []:
            filout.write(" ".join(d_stock[pdb]))
        filout.write("\n")
        
    filout.close()
    
    return pr_init + "resultLigandInPDB"



def retriveListPDB (pr_PDB):
    
    list_files = listdir(pr_PDB)
        
    l_p_pdb = []
    for file_PDB in list_files : 
        if file_PDB[-4:] == ".pdb" : 
            l_p_pdb.append (pr_PDB + file_PDB)
            
    return l_p_pdb


