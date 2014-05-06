from copy import deepcopy
# from sys import exit



def retrieveSubstruct (l_lig_parsed, name_ligand):
    
    l_out = []
    if name_ligand == "AMP" : 
        for lig_parsed in l_lig_parsed : 
            if lig_parsed["name"] == "O3P" or lig_parsed["name"] == "O2P" or lig_parsed["name"] == "O1P" or lig_parsed["name"] == "O5'" or lig_parsed["name"] == "P": 
                l_out.append (deepcopy(lig_parsed))
    elif name_ligand == "ADP":
        for lig_parsed in l_lig_parsed : 
            if lig_parsed["name"] == "O1B" or lig_parsed["name"] == "O2B" or lig_parsed["name"] == "O3B" or lig_parsed["name"] == "PB" or lig_parsed["name"] == "O3A" or lig_parsed["name"] == "O1A" or lig_parsed["name"] == "O2A" or lig_parsed["name"] == "PA" or lig_parsed["name"] == "O5'" :
                l_out.append (deepcopy(lig_parsed))
    elif name_ligand == "POP" : 
        l_out = l_lig_parsed
    elif name_ligand == "ATP" : 
        for lig_parsed in l_lig_parsed : 
            if lig_parsed["name"] == "O1A" or lig_parsed["name"] == "O2A" or lig_parsed["name"] == "O3A" or lig_parsed["name"] == "PA" or lig_parsed["name"] == "O1B" or lig_parsed["name"] == "O2B" or lig_parsed["name"] == "O3B" or lig_parsed["name"] == "PB" or lig_parsed["name"] == "O1G" or lig_parsed["name"] == "O2G" or lig_parsed["name"] == "O3G" or lig_parsed["name"] == "PG" or lig_parsed["name"] == "O5'" :
                l_out.append (deepcopy(lig_parsed))
                
    return l_out
    
    
    
