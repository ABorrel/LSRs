from re import search 
from os import listdir, path
from math import sqrt
from copy import deepcopy

import pathManage
import parseShaep
import parsePDB
import writePDBfile


def selectSmileCode (p_file_smile, minimal_length_smile = 3):
    """
    selection: -> length 
               -> separated by .
    """
    p_filout = p_file_smile[0:-4] + ".filter"
    filout = open (p_filout, "w")
    
    d_smile = {}
    
    filin = open (p_file_smile, "r")
    l_line_smile = filin.readlines ()
    filin.close ()
    
    for line_smile in l_line_smile : 
        
        l_element = line_smile.strip ().split ("\t")
        l_pdb = l_element [2].split (" ")
        l_ligand = l_element[3].split (" ")
        smile = l_element[0]
        
        # length 
        if len (smile) < minimal_length_smile : 
            continue
        elif search("\.", smile) : 
            l_smile = smile.split (".")
            for smile_separated in l_smile : 
                if len (smile_separated) < minimal_length_smile : 
                    continue
                elif not smile_separated in d_smile.keys () : 
                    d_smile[smile_separated] = {}
                    d_smile[smile_separated]["ligand"] = l_ligand
                    d_smile[smile_separated]["PDB"] = l_pdb
                else :
                    for lig in l_ligand : 
                        if not lig in d_smile[smile_separated]["ligand"] : 
                            d_smile[smile_separated]["ligand"].append (lig)
                    for pdb in l_pdb : 
                        if not pdb in d_smile[smile_separated]["PDB"] : 
                            d_smile[smile_separated]["PDB"].append (pdb)
        else : 
            if not smile in d_smile.keys () : 
                d_smile[smile] = {}
                d_smile[smile]["ligand"] = l_ligand
                d_smile[smile]["PDB"] = l_pdb
            else :
                for lig in l_ligand : 
                    if not lig in d_smile[smile]["ligand"] : 
                        d_smile[smile]["ligand"].append (lig)
                for pdb in l_pdb : 
                    if not pdb in d_smile[smile]["PDB"] : 
                        d_smile[smile]["PDB"].append (pdb)
        
    for smile_code in d_smile.keys () : 
        filout.write (str (smile_code) + "\t" + str (len (d_smile[smile_code]["PDB"])) + "\t" + " ".join (d_smile[smile_code]["PDB"]) + "\t" + " ".join(d_smile[smile_code]["ligand"]) + "\n")
            
    filout.close ()
        
    return p_filout
    

def globalShaepStat (substruct):
    
    pr_result = pathManage.result(substruct)
    
    p_filout = pr_result + "shaep_global.txt"
    filout = open (p_filout, "w")
    filout.write ("best_similarity\tshape_similarity\tESP_similarity\n")
    
    l_folder = listdir(pr_result)
    
    
    for ref_folder in l_folder  :
        if not path.isdir(pr_result + ref_folder + "/") : continue
        l_file_result = listdir(pr_result + ref_folder + "/")
        for file_result in l_file_result : 
            if search(".hit", file_result) :
                d_shaep_parsed = parseShaep.parseOutputShaep(pr_result + ref_folder + "/" + file_result) 
                if d_shaep_parsed != {} : 
                    filout.write (ref_folder + "_" + file_result[10:-4] + "\t" + str(d_shaep_parsed["best_similarity"]) + "\t" + str(d_shaep_parsed["shape_similarity"]) + "\t" + str(d_shaep_parsed["ESP_similarity"]) + "\n")
    filout.close ()
                
        
        
def computeRMSDBS (p_ref, p_BS, pr_result) :
    
    p_filout_pdb = pr_result + p_BS.split ("/")[-1][0:-4] + "_" + p_ref.split ("/")[-1]
    filout_pdb = open (p_filout_pdb, "w")
    
    l_atomBS_parsed = parsePDB.loadCoordSectionPDB(p_BS, "ATOM")
    l_pdb_parsed = parsePDB.loadCoordSectionPDB(p_ref)
    
    writePDBfile.coordinateSection(filout_pdb, l_pdb_parsed, recorder = "ATOM")
    writePDBfile.coordinateSection(filout_pdb, l_atomBS_parsed, recorder = "ATOM")
    
    l_BS_ref = []
    for atomBS_parsed in l_atomBS_parsed :
#         print  atomBS_parsed 
        d_max = 100.0 
        for atom_ref in l_pdb_parsed :
            if atomBS_parsed["resName"] ==  atom_ref["resName"] and atomBS_parsed["name"] ==  atom_ref["name"] : 
#                 print "@@@@@@@@@@@"
                d = parsePDB.distanceTwoatoms(atomBS_parsed, atom_ref)
#                 print d
                if d < d_max : 
                    d_max = d
                    res_temp = atom_ref
        
        if d_max != 100.0 : 
            l_BS_ref.append (deepcopy(res_temp))
        else : 
            return []
    
    filout_pdb.close ()
    l_RMSD = RMSDTwoList (l_atomBS_parsed, l_BS_ref, )
    return l_RMSD
    
    
    
def RMSDTwoList (l_atom1, l_atom2) : 
    
    nb_ca = 0.0
    diff_position_all = 0.0
    diff_position_ca = 0.0
    
    if len (l_atom1) != len (l_atom2) or len (l_atom2) == 0 : 
        print "ERROR - RMSD: list length different or null"
        return []
    else : 
        i = 0
        while i < len (l_atom1): 
            if l_atom1[i]["name"] != l_atom2[i]["name"] : 
                print l_atom1[i]["name"] , l_atom2[i]["name"]
                print "ERROR"
                return []
            else : 
                diff_position_all = diff_position_all + (l_atom1[i]["x"] - l_atom2[i]["x"]) * (l_atom1[i]["x"] - l_atom2[i]["x"]) + (l_atom1[i]["y"] - l_atom2[i]["y"]) * (l_atom1[i]["y"] - l_atom2[i]["y"]) + (l_atom1[i]["z"] - l_atom2[i]["z"]) * (l_atom1[i]["z"] - l_atom2[i]["z"])
                
                if l_atom1[i]["name"] == "CA" : 
                    diff_position_ca = diff_position_ca + (l_atom1[i]["x"] - l_atom2[i]["x"]) * (l_atom1[i]["x"] - l_atom2[i]["x"]) + (l_atom1[i]["y"] - l_atom2[i]["y"]) * (l_atom1[i]["y"] - l_atom2[i]["y"]) + (l_atom1[i]["z"] - l_atom2[i]["z"]) * (l_atom1[i]["z"] - l_atom2[i]["z"])
                    nb_ca = nb_ca + 1
            i = i + 1
    
    return [sqrt(diff_position_all / len (l_atom1)), sqrt (diff_position_ca / nb_ca), len (l_atom1)]
                
            
            
            
            
            
                
    
    
    
    
    
    
    
    





