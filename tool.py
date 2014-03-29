from os import path
import parsePDB
import writePDBfile





def removeChain (path_protein_PDB, p_dir_out) : 
    
    name_file = path.split(path_protein_PDB)[-1]
    path_filout = p_dir_out + name_file
    
    list_atom = parsePDB.loadCoordSectionPDB(path_protein_PDB, section="ATOM")
    
    for atom in list_atom : 
        atom["chainID"] = ""
    
    writePDBfile.coordinateSection(path_filout, list_atom, "ATOM")
    
    return path_filout



def transformAA (aa):
    
    aa = aa.upper()
    dico_code = {"S":"SER", "T":"THR", "N":"ASN", "Q":"GLN", "E":"GLU", "D":"ASP", "K":"LYS", "R":"ARG", "H":"HIS", "M":"MET", "C":"CYS", "W":"TRP", "F":"PHE", "Y":"TYR", "A":"ALA", "V":"VAL", "L":"LEU", "I":"ILE", "P":"PRO", "G":"GLY"}
    
    if len (aa) == 1 : 
        return dico_code[aa]
    else :
        for aa_one in dico_code.keys ():
            if dico_code[aa_one] == aa : 
                return  aa_one
            

            

