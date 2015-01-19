from copy import deepcopy, copy
from re import search

import writePDBfile
import runOtherSoft
import superposeStructure
import parsePDB


def retrievePi (substruct_parsed):
    
    l_out = []
    for atom_substruct in substruct_parsed : 
        if atom_substruct["element"] == "P" or atom_substruct["element"] == "O": 
            out = copy(atom_substruct)
            l_out.append (out)
    if l_out == [] : 
        print "ERROR l15 -- neighborSearch"
    return l_out




def searchNeighborAtom(substruct_parsed, lig_query_parsed, struct_type, log_file, thresold_superimposed_ribose = 2.5, thresold_superimposed_pi = 3 ) : 
    
    l_atom_substituate = []
    if struct_type == "ribose" : 
        for atom_substruct in substruct_parsed :
            for atom_query in lig_query_parsed : 
                if parsePDB.distanceTwoatoms(atom_substruct, atom_query) <= thresold_superimposed_ribose : 
                    out = copy(atom_query)
                    if not out in l_atom_substituate : 
                        l_atom_substituate.append (out)
        
    else : 
        l_atom_interest = retrievePi (substruct_parsed)
        for atom_interest in l_atom_interest :
            for atom_query in lig_query_parsed : 
                if parsePDB.distanceTwoatoms(atom_interest, atom_query) <= thresold_superimposed_pi : 
                    out = copy(atom_query)
                    if not out in l_atom_substituate : 
                        l_atom_substituate.append (out)
    
    # control out empty
    if l_atom_substituate == [] : 
        log_file.write ("[Not substituate] -> " + substruct_parsed[0] ["resName"]+ struct_type + "\n")
        return []
    else : 
        return l_atom_substituate
        
        
    
