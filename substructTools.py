from copy import deepcopy
# from sys import exit

import parsePDB


def retrieveSubstruct (l_lig_parsed, name_ligand):
    
    d_out = {}
    l_ribose = []
    l_pi1 = []
    l_pi2 = []
    l_pi3 = []
    
    if name_ligand == "POP" : 
        for lig_parsed in l_lig_parsed : 
            if lig_parsed["name"] == "P1" or \
            lig_parsed["name"] == "O1" or \
            lig_parsed["name"] == "O2" or \
            lig_parsed["name"] == "O3" or \
            lig_parsed["name"] == "O": 
                l_pi1.append (deepcopy(lig_parsed))
            if lig_parsed["name"] == "P2" or \
            lig_parsed["name"] == "O4" or \
            lig_parsed["name"] == "O5" or \
            lig_parsed["name"] == "O6" or \
            lig_parsed["name"] == "O": 
                l_pi2.append (deepcopy(lig_parsed))
    
    else : 
        for lig_parsed in l_lig_parsed : 
            
            if lig_parsed["name"] == "O5'" or \
            lig_parsed["name"] == "C5'" or \
            lig_parsed["name"] == "C4'" or \
            lig_parsed["name"] == "O4'"or \
            lig_parsed["name"] == "C3'" or \
            lig_parsed["name"] == "O3'" or \
            lig_parsed["name"] == "C2'"or \
            lig_parsed["name"] == "O2'"or \
            lig_parsed["name"] == "C1'": 
                l_ribose.append (deepcopy(lig_parsed))
            
            elif  lig_parsed["name"] == "PA" or \
            lig_parsed["name"] == "O1A" or \
            lig_parsed["name"] == "O2A" or \
            lig_parsed["name"] == "O3A" : 
                print lig_parsed["name"]
                l_pi1.append (deepcopy(lig_parsed))
            
            elif  lig_parsed["name"] == "PB" or \
            lig_parsed["name"] == "O1B" or \
            lig_parsed["name"] == "O2B" or \
            lig_parsed["name"] == "O3B" : 
                l_pi2.append (deepcopy(lig_parsed))
            
            elif  lig_parsed["name"] == "PG" or \
            lig_parsed["name"] == "O1G" or \
            lig_parsed["name"] == "O2G" or \
            lig_parsed["name"] == "O3G" : 
                l_pi3.append (deepcopy(lig_parsed))
    
    if l_ribose != [] : 
        d_out["ribose"] = l_ribose
    if l_pi1 != [] : 
        d_out["pi1"] = l_pi1
    if l_pi2 != [] : 
        d_out["pi2"] = l_pi2
    if l_pi3 != [] : 
        d_out["pi3"] = l_pi3   
                
    return d_out
    
    
