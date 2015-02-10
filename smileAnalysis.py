from re import search
import pathManage
import parsePDB
import writePDBfile


def smallLSR (smile, thresold_small = 3) : 

    smile = smile.upper ()
    smile = smile.replace ("+", "")
    smile = smile.replace ("-", "")
    smile = smile.replace ("[", "")
    smile = smile.replace ("]", "")
    smile = smile.replace ("(", "")
    smile = smile.replace (")", "")
    smile = smile.replace ("=", "")

    l_part = smile.split (".")

    for part in l_part : 
        if len (part) > thresold_small : 
            return 0

    return 1



def searchP (smile):
    
    if search ("P", smile) : 
        return 1
    else : 
        return 0



def searchB (smile):
    
    
    if not search ("B", smile) : 
        return 0
    else : 
        nb_char = len(smile)
        i = 0 
        while (i < nb_char) : 
            if smile[i] == "B" and (i + 1) != nb_char: 
                # check if not BE, Br, Bk, Bh
                if not smile[i + 1] in ["e", "r", "k", "h"] : 
                    print smile
                    return 1
            i = i + 1
    return 0


def searchBe (smile) : 

    if search ("Be", smile) : 
        return 1
    else : 
        return 0

    
def searchMetal (smile):
    """Pb, MAJ on metals and the ionic form is sometimes associated with other groups"""
    l_metal = ["B", "F", "I", "K", "V", "W", "Y","AG", "AL", "AR" ,"AU", "BA", "BE", "BR", "CA","CD","CE","CF","CL","CO","CR","CS","CU","EU","FE","GA","GD","HE","HF","HG","IN","IR","KR","LA","LI","LU" ,"MG","MN" ,"MO" ,"NA","ND","NE","NI","OS","PB","PD","PR","PT","RB","RE","RU","SB","SE","SI","SM","SR","TA","TB","TE","TL","XE","YB","ZN","ZR"]
    
    smile = smile.upper ()
    smile = smile.replace ("+","")
    smile = smile.replace ("-","")
    smile = smile.replace ("[", "")
    smile = smile.replace ("]", "")
    
    for metal in l_metal : 
        if search("\." + metal + "\.", smile) : 
            return metal
        elif search ("\." + metal + "$", smile) : 
            return metal
        elif search ("^" + metal + "\." , smile) : 
            return metal
    return 0

def searchRing (smile):
    """Search only number now"""
    
    if search ("1", smile) : 
        return 1 
    else : 
        return 0
    
    

# def countlenRing (smile):
#    
#    c = 0
#    size = len (smile)
#    i = 0
#    ring_open = 0
#    
#    while i < size : 
#        if ring_open == 0 : 
#            if smile[i] == "1" and smile[i-1].upper () == "C":
#                c = c + 1
#                ring_open = 1
#        elif ring_open == 1 : 
#            if smile[i] == "(" : 
#                ring_open = 2
#            elif smile[i].upper () == "C" : 
#                c = c + 1
#            elif smile[i] == "1": 
#                return c
#            elif smile[i] != "=" and smile[i] != "[" and smile[i] != "]" and smile[i] != "#" and smile[i] != "@" and smile[i] != "H": 
#                return 99
#        elif ring_open == 2 : 
#            if smile[i] == ")" : 
#                ring_open = 1
#        
#        i = i + 1
#    
#    return 1       


def searchSulfonyl (smile): 
    
    smile = smile.replace ("[", "")
    smile = smile.replace ("]", "")
    if search("O=S\(=O\)", smile) or search ("S\(=O\)\(=O\)", smile) or  search("S\(=O\)O", smile ) or search("O=S=O", smile) or search ("\(=O\)S\(=O\)", smile) or search ("S\(O\)O", smile): 
        return 1
    else : 
        return 0


def searchS (smile):
    if search ("S", smile) : 
        return 1
    return 0


def searchCl (smile):
    if search ("Cl", smile) : 
        return 1
    return 0



def searchF (smile):
    
    
    if not search ("F", smile) : 
        return 0
    else : 
        nb_char = len(smile)
        i = 0 
        while (i < nb_char) : 
            if smile[i] == "F" and (i + 1) != nb_char: 
                # check if not Fe,..
                if not smile[i + 1] in ["e", "r", "k", "h", "m", "n"] : 
                    print smile
                    return 1
            i = i + 1
    return 0


def searchBr (smile) : 
    if search ("Br", smile) : 
        return 1
    return 0

def searchNO2 (smile): 
    
    smile = smile.replace ("[", "")
    smile = smile.replace ("]", "")

    if search("O=N\(=O\)", smile) or search ("N\(=O\)\(=O\)", smile) or  search("N\(=O\)O", smile ) or search("O=N=O", smile) or search ("\(=O\)N\(=O\)", smile) or search ("N\(O\)O", smile): 
        print smile
        return 1
    else : 
        return 0

    
             
def searchCON (smile):
    
    if search ("NC\(=O\)", smile) or search ("C\(=O\)N", smile) : 
        return 1
    else : 
        return 0

def searchCarboxy (smile) :
    if search ("O=C\(O\)", smile) or search ("C\(=O\)O", smile) : 
        return 1
    else : 
        return 0

def searchCandOandN (smile) : 
    
    for at in smile : 
        if search("[a-z,A-Z]", at) : 
            if at.upper() != "C" and at.upper() != "O" and at.upper() != "N" : 
                return 0
    return 1


def searchCandO (smile) : 
    
    for at in smile : 
        if search("[a-z,A-Z]", at) : 
            if at.upper() != "C" and at.upper() != "O" : 
                return 0
    return 1

def searchCandN (smile) : 
    
    for at in smile : 
        if search("[a-z,A-Z]", at) : 
            if at.upper() != "C" and at.upper() != "N" : 
                return 0
    return 1

 
 
def searchConly (smile) : 
    
    for at in smile : 
        if search("[a-z,A-Z]", at) : 
            if at.upper() != "C" : 
                return 0
    return 1
 

def searchReplacement (smile, PDB_query, PDB_ref, name_ligand, in_cycle = 0) : 
    
    metal_find  = searchMetal (smile)
    if metal_find != 0 : 
        p_dir_dataset = pathManage.dataset(name_ligand + "/" + PDB_ref)
        l_PDB_query = pathManage.findPDBQueryDataset(p_dir_dataset)
        for p_query in l_PDB_query : 
            if search (PDB_query, p_query): 
                p_PDB_query = p_query
                break
        if "p_PDB_query" in locals() : 
            l_atom_parsed = parsePDB.loadCoordSectionPDB(p_query)
            l_ions_PDB = parsePDB.retrieveListIon(l_atom_parsed)
            
            if metal_find in l_ions_PDB : 
                l_atom_ion = parsePDB.retrieveLigand(l_atom_parsed, metal_find)
                filout = open (p_dir_dataset + str(metal_find) + "_" + p_query.split("/")[-1], "w")
                for atom_ion in l_atom_ion : 
                    writePDBfile.coordinateSection(filout, atom_ion, recorder = "HETATM", header = str(metal_find), connect_matrix = 0)
                filout.close ()
                return "metal", metal_find
    


    if in_cycle == 0:
        if searchRing(smile) == 1 : 
            return "cycle",""      
    if searchP(smile) == 1 : 
        return "P", ""
    elif searchB(smile) == 1 : 
        return "B",""
    elif searchF (smile) == 1 :
        return "F", ""
    elif searchCl (smile) == 1 :
        return "Cl", ""
    elif searchBr (smile) == 1 :
        return "Br", ""
    elif searchBe (smile) == 1 : 
        return "Be", ""
    elif searchNO2 (smile) == 1 : 
        return "NO2", ""
    elif searchSulfonyl(smile) == 1: 
        return "SO2",""
    elif searchS (smile) == 1 :
        return "S", ""
    elif searchCON (smile) == 1 : 
        return "CON",""
    elif searchCarboxy (smile) == 1 : 
        return "COO",""
    elif searchConly(smile) == 1 :
        return "onlyC", ""
    elif searchCandO (smile) == 1 : 
        return "C+O", ""
    elif searchCandN (smile) == 1 : 
        return "C+N", "" 
    elif searchCandOandN (smile) == 1 :
        return "C+O+N", ""

    return "other"  ,""

 
