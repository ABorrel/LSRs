from os import path
import parsePDB
import writePDBfile








def removeChain (path_protein_PDB) : 
    
    path_directory = path.dirname(path_protein_PDB)
    path_filout = path_directory + "/pralign.pdb"
    
    list_atom = parsePDB.loadCoordSectionPDB(path_protein_PDB)
    
    for atom in list_atom : 
        atom["chainID"] = ""
    
    writePDBfile.coordinateSection(path_filout, list_atom, "ATOM")
    
    return path_filout
