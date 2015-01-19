from os import system, listdir, makedirs
from re import search, compile, findall
from os.path import isfile
from urllib import urlretrieve
from shutil import copy




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
