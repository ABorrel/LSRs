"""
BORREL Alexandre
09-2012
"""

import runOtherProg
import pathDirectory
from os import remove, system
from re import search, sub
import parsePDB
import writePDBfile

def superposeApoHolo (dictionary_dataset, name_directory, pocket_retireve):
    """
    Superpose Apo and holo and run pocket detection on apo form
    return: list with apo superposed and holo superposed files
    """
    
    
#     print dictionary_dataset
    path_file_RMSD = pathDirectory.descriptor(name_directory) + "RMSD_apo_holo"
    fileRMSD = open (path_file_RMSD, "w")
    
    list_PDB = dictionary_dataset.keys ()
#    print list_PDB
    for PDB_ID in list_PDB :
        if dictionary_dataset[PDB_ID]['Type structure'] == "apo structure" : 
            path_file_superpose = runOtherProg.runTMalign( pathDirectory.dataSet(name_directory + "/" + dictionary_dataset[PDB_ID]['PDB holo'][0] )  + "protein.pdb",pathDirectory.dataSet(name_directory + "/" + PDB_ID) + "protein.pdb", pathDirectory.descriptor(name_directory + "/" + pocket_retireve +"/" + PDB_ID))
#            path_file_superpose = runOtherProg.runMMalign(pathDirectory.dataSet(name_directory + "/" + PDB_ID) + "protein.pdb", pathDirectory.dataSet(name_directory + "/" + dictionary_dataset[PDB_ID]['PDB holo'][0] )  + "protein.pdb", pathDirectory.descriptor(name_directory + "/" + PDB_ID))

#            print path_file_superpose
            RMSD = retrieveRMSDFileTMalign (path_file_superpose[-1])
            fileRMSD.write (PDB_ID + "\t" + dictionary_dataset[PDB_ID]['PDB holo'][0] + "\t" + RMSD + "\n")
            
#             remove (path_file_superpose[0])
#             remove (path_file_superpose[1])
#             remove (path_file_superpose[2])
#             remove (path_file_superpose[-1])
    
#             path_files_ApoHolo_superposed = diviseSuperposeFile (path_file_superpose[-2])# list of file holo and apo
    
 
 
 
 
def retrieveRMSDFileTMalign (path_file_RMSD) : 
    """
    Retrieve RMSD in TMalign out file
    args: -> path file RMSD
    return: value of RMSD
    """
    
    filin = open (path_file_RMSD, "r")
    filin_read = filin.read()
    filin.close ()
    
    part_file = filin_read.split ("RMSD=")[1]
    RMSD = part_file.split (",")[0]
    try : RMSD = RMSD.replace (" ", "")
    except : pass
    
    return RMSD
    
    
    

def diviseSuperposeFile (path_file_superpose) : 
    """
    Divise files with apo and holo structure
    arg: - path file superpose all atoms
    return: - path file apo
            - path file holo
    """
    path_file_apo = path_file_superpose + "_apo"
    path_file_holo = path_file_superpose + "_holo"
    
    filin = open (path_file_superpose, "r")
    list_lines = filin.readlines ()
    
    filout = open (path_file_apo, "w")
    
    for lines in list_lines [:-1]: 
        if search("^ATOM", lines) : 
            filout.write (lines)
        elif search("^TER", lines):
            filout.close ()
            filout = open (path_file_holo, "w")
    
    return [path_file_apo, path_file_holo]
            
    
    
def manageTMalign (path_protein ) : 
    
    list_atoms = parsePDB.loadCoordSectionPDB(path_protein)
    dico_residues = parsePDB.arrangeResidues(list_atoms)
    list_res = dico_residues.keys()
    list_res.sort ()
    
    filout = open (path_protein, "w")
    for resID in list_res : 
        for atom in dico_residues[resID] : 
            writePDBfile.coordinateStructure(atom, "ATOM", filout)
            
    
    filout.close ()
    return path_protein
    
    
    
def applyTranslocMatrix (path_protein, path_matrix) :
    """
    return list of atoms
    """
    # format matrix
    matrix_transloc = formatMatrix(path_matrix)
    
    # apply matrix
    return applyMatrixProt (path_protein, matrix_transloc)
    
    
 
 
 
    
def formatMatrix(path_file_matrix):

    dico_matrix = {}
    filin = open (path_file_matrix, "r")
    list_lines = filin.readlines ()
    filin.close ()

    m = 1
    for line_file in list_lines[2:5] :

        line_format = sub("[ ]{2,}", " ", line_file.strip())
        line_format = line_format.split (" ")
        dico_matrix["t" + str(m)] = float(line_format[1])
        dico_matrix["u" + str(m) + "1"] = float(line_format[2])
        dico_matrix["u" + str(m) + "2"] = float(line_format[3])
        dico_matrix["u" + str(m) + "3"] = float(line_format[4])
        m = m + 1

    return dico_matrix


def applyMatrixProt (path_file_PDB, matrix_transloc) :

    list_atom = parsePDB.loadCoordSectionPDB(path_file_PDB, section = "ATOM")
    for atom in list_atom : 
        atomx = matrix_transloc["t1"] + matrix_transloc["u11"] * atom["x"] + matrix_transloc["u12"] * atom["y"] + matrix_transloc["u13"] * atom["z"]
        atomy = matrix_transloc["t2"] + matrix_transloc["u21"] * atom["x"] + matrix_transloc["u22"] * atom["y"] + matrix_transloc["u23"] * atom["z"]
        atomz = matrix_transloc["t3"] + matrix_transloc["u31"] * atom["x"] + matrix_transloc["u32"] * atom["y"] + matrix_transloc["u33"] * atom["z"]
        atom["x"] = atomx
        atom["y"] = atomy
        atom["z"] = atomz

    return list_atom
        
def applyMatrixLigand (l_atoms, matrix_transloc) :

    for atom in l_atoms : 
        atomx = matrix_transloc["t1"] + matrix_transloc["u11"] * float(atom["x"]) + matrix_transloc["u12"] * float(atom["y"]) + matrix_transloc["u13"] *float( atom["z"])
        atomy = matrix_transloc["t2"] + matrix_transloc["u21"] * atom["x"] + matrix_transloc["u22"] * atom["y"] + matrix_transloc["u23"] * atom["z"]
        atomz = matrix_transloc["t3"] + matrix_transloc["u31"] * atom["x"] + matrix_transloc["u32"] * atom["y"] + matrix_transloc["u33"] * atom["z"]
        atom["x"] = atomx
        atom["y"] = atomy
        atom["z"] = atomz

        