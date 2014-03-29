from copy import deepcopy, copy
from re import search

import writePDBfile
import runOtherSoft
import superposeStructure
import parsePDB


def retrieveAtomCentral (lig_parsed, substruct):
    
    l_out = []
    for atom_lig in lig_parsed : 
        if atom_lig["element"] == "P" : 
            if atom_lig["resName"] == substruct : 
                l_out.append (atom_lig)
    return l_out    




def searchNeighborAtom (lig_ref, lig_superimposed, p_complex, substruct, p_lig, p_matrix, p_dir_result, thresold_substruct = 6.0, thresold_binding = 6.0):
    
    l_atom_interest = retrieveAtomCentral (lig_ref, substruct)
    
    # Step1 -> search substructure
    l_atom_substruct = []
    for atom_interest in l_atom_interest : 
        for atom_lig_superimposed in lig_superimposed : 
#             print parsePDB.distanceTwoatoms(atom_interest, atom_lig_superimposed), atom_interest["element"], atom_lig_superimposed["element"]
            if parsePDB.distanceTwoatoms(atom_interest, atom_lig_superimposed) <= thresold_substruct : 
                out = copy(atom_lig_superimposed)
                l_atom_substruct.append (out)
    
    if l_atom_substruct == [] : 
        print "Not substruct find"
        return 0
#     else : 
#         print len (l_atom_substruct), "*******"
    
    p_file_substruct = p_dir_result +  "substruct_" + p_lig.split ("/")[-1]
#     print p_file_substruct, "check l27 neighbor"
    writePDBfile.coordinateSection(p_file_substruct, l_atom_substruct, "HETATM", p_lig.split ("/")[-1], connect_matrix = 1)
    
    # Step2 -> convert to smile -> review to smart if a find a good software
    smile_code = runOtherSoft.babelConvertPDBtoSMILE(p_file_substruct)
    #print smile_code
    
    #if search ("P",smile_code) : 
    #    print smile_code
    # Step3 -> find the binding site 
    
    # 1. search in PDB file
    l_atom_complex = parsePDB.loadCoordSectionPDB(p_complex, "ATOM")
    # 2.apply rotated matrix on protein
    superposeStructure.applyMatrixProt(l_atom_complex, p_matrix)
    p_file_cx = p_dir_result +  "CX_" + p_lig.split ("/")[-1]
    writePDBfile.coordinateSection(p_file_cx, l_atom_complex, "ATOM", p_lig.split ("/")[-1], connect_matrix = 0)
    
    
    l_atom_binding_site = []
    for atom_substruct in l_atom_substruct : 
        for atom_complex in l_atom_complex : 
#             print parsePDB.distanceTwoatoms (atom_substruct, atom_complex), atom_substruct["element"], atom_complex["element"]
            if parsePDB.distanceTwoatoms (atom_substruct, atom_complex) <= thresold_binding :
                l_atom_binding_site.append (deepcopy(atom_complex)) 
    # 3. retrieve complet residue
    l_atom_res = parsePDB.getResidues(l_atom_binding_site, l_atom_complex)
    #print len (l_atom_res)
    
    # 4. write binding site
    p_binding = p_dir_result +  "BS_" + p_lig.split ("/")[-1]
    writePDBfile.coordinateSection(p_binding, l_atom_res, "ATOM", p_binding, connect_matrix = 0)
    
    
    return smile_code


