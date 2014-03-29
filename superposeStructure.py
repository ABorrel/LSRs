"""
BORREL Alexandre
09-2012
"""

from re import sub
import parsePDB
import writePDBfile


 
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


def applyMatrixProt (l_atom, p_matrix) :
    
    matrix_transloc = formatMatrix(p_matrix)

    for atom in l_atom : 
        atomx = matrix_transloc["t1"] + matrix_transloc["u11"] * atom["x"] + matrix_transloc["u12"] * atom["y"] + matrix_transloc["u13"] * atom["z"]
        atomy = matrix_transloc["t2"] + matrix_transloc["u21"] * atom["x"] + matrix_transloc["u22"] * atom["y"] + matrix_transloc["u23"] * atom["z"]
        atomz = matrix_transloc["t3"] + matrix_transloc["u31"] * atom["x"] + matrix_transloc["u32"] * atom["y"] + matrix_transloc["u33"] * atom["z"]
        atom["x"] = atomx
        atom["y"] = atomy
        atom["z"] = atomz

        
def applyMatrixLigand (l_atoms, p_matrix ) :
    
    matrix_transloc = formatMatrix(p_matrix)
    for atom in l_atoms : 
        atomx = matrix_transloc["t1"] + matrix_transloc["u11"] * float(atom["x"]) + matrix_transloc["u12"] * float(atom["y"]) + matrix_transloc["u13"] *float( atom["z"])
        atomy = matrix_transloc["t2"] + matrix_transloc["u21"] * atom["x"] + matrix_transloc["u22"] * atom["y"] + matrix_transloc["u23"] * atom["z"]
        atomz = matrix_transloc["t3"] + matrix_transloc["u31"] * atom["x"] + matrix_transloc["u32"] * atom["y"] + matrix_transloc["u33"] * atom["z"]
        atom["x"] = atomx
        atom["y"] = atomy
        atom["z"] = atomz

        