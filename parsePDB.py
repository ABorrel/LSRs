"""
BORREL Alexandre
04-2012
Analysis PDB file
"""
# global module
from re import search, sub
from math import sqrt, acos, asin, cos, sin, degrees
from copy import deepcopy

# personal
import tool


def retrieveListLigand ( path_file_PDB, option_not_ligand = 0 , debug=0 ):
    """
    Retrieve list ligand in PDB file, only ligands, remove metals
    args: - path file
    return: list ligands
    """
    
    list_not_ligand = ["SO4", "FE2", "GOL", "FES", "PO4", "MSE", "DMS", "URE", "FMT", "TRS", "NCO"]
    
    filin = open ( path_file_PDB, "r" )
    list_line_pdb = filin.readlines ()
    filin.close()
    list_out = []
    for line_pdb in list_line_pdb : 
        if search ( "^HETATM", line_pdb ) : 
            ligand = line_pdb[17:21].replace ( " ", "" )
            if debug : print ligand
            if not ligand in list_out and ligand != "HOH"  :
                if len (ligand) > 2 : # remove metals just 2 letters
                    if not ligand in list_not_ligand : 
                        list_out.append ( ligand )
    return list_out 
    


def lineCoords (line):
    """Parsing line of coordinate PDB File
    in: line
    out: dictionnary atom"""

    atom = {}
    atom["type"] = line[0:6].replace (" ", "")
    try :atom["serial"] = int(line[6:11].replace (" ", ""))
    except :line[6:11].replace (" ", "")
    atom["name"] = line[12:16].replace (" ", "")
    atom["char"] = line[16]
    atom["resName"] = line[17:20].replace (" ", "")
    try : atom["chainID"] = str(line[21])
    except : print line
    try : atom["resSeq"] = int (line[22:26].replace (" ", ""))
    except : atom["resSeq"] = 0
    atom["iCode"] = str(line[26])
    atom["x"] = float (line[30:38].replace (" ", ""))
    atom["y"] = float (line[38:46].replace (" ", ""))
    atom["z"] = float (line[46:54].replace (" ", ""))
    atom["element"] = line[76:78].replace (" ", "")
    # pqr without element
    if atom["element"] == "" :
        if type (atom["name"][0]) is int :
            atom["element"] = atom["name"][1]
        else :
            atom["element"] = atom["name"][0] 
    
    atom["charge"] = line[78:80].replace (" ", "")
    atom["occupancy"] = line[54:60].replace (" ", "")
    atom["tempFactor"] = line[60:66].replace (" ", "")
    
    atom["connect"] = []
    return atom


def loadCoordSectionPDB (path_PDB_file, section = "", remove_H = 0, debug = 0):
    """
    Retrieve every atom in cordinate section. If it is NMR complex
    retrieve only first model
    
    """
    
    list_atom = []
    filin = open (path_PDB_file, "r")
    list_line_PDB = filin.readlines()
    filin.close ()
    
    
    for line_PDB in list_line_PDB :
        #End model
        if search ("^ENDMDL", line_PDB) : 
            break
        
        if section == "" : 
            if search ("^ATOM", line_PDB) or search ("^HETATM", line_PDB) : 
                atom = lineCoords(line_PDB)
                if remove_H != 0 : 
                    if atom["element"] == "H" : 
                        continue
                list_atom.append (atom)
        else : 
            if search ("^" + section, line_PDB)  : 
                atom = lineCoords(line_PDB)
                if remove_H != 0 : 
                    if atom["element"] == "H" : 
                        continue
                list_atom.append (lineCoords(line_PDB))
                
#    if debug : 
#        print "TEST"            
#        print len (list_atom)
#    
    return list_atom
    




def retrieveLigand  (list_atom_parsed, name_ligand, extend = 2, debug = 0):
    """
    Retrieve list of ligand in PDB
    args: -> PDB parsed
          -> name ligand
    return: list of ligand with atoms
    """
    
    if debug : 
        print "-control ligand select-"
        print name_ligand, "NAME ligand"
        print "l144 - parse PDB"
        print "****************"
        
    name_ligand = name_ligand.upper ()
    list_atom_ligand = []
    for element in list_atom_parsed : 
        if element ["resName"] == name_ligand : 
            list_atom_ligand.append (deepcopy(element))
            
    for atomLigand in list_atom_ligand:
        atomLigand["connect"].append(atomLigand["serial"])
        for atom_ligand_connect in list_atom_ligand:
            distance = distanceTwoatoms(atomLigand, atom_ligand_connect)
            if distance < extend and distance != 0:
                if not atom_ligand_connect["serial"] in atomLigand["connect"]:
                    atomLigand["connect"].append(atom_ligand_connect["serial"])
    
    # Check exotic bound (SE-C...) with special length of bound
    if checkConnectMatrix(list_atom_ligand) == 0 :
        if debug : 
            print "- exotic ligand -, recursive distance"
            print "l165 - parse PDB"
            print "Threshold: ", extend
            print "--------------------"
        if extend <= 4.0 : 
            return retrieveLigand (list_atom_parsed, name_ligand, extend + 0.1, debug = 0)
    
    return separateByLigand (list_atom_ligand)

def checkConnectMatrix (list_atom_ligand):
    
    list_nb_atom = []
    list_atom_serial_by_ligand = retrieveListAtomID (list_atom_ligand)
    for atom_serial in list_atom_serial_by_ligand : 
        list_nb_atom.append(len(atom_serial)), "len of ligand"
        
    list_nb_atom = list(set(list_nb_atom))
    if len(list_nb_atom) != 1 : 
        return 0
    return 1
    
    

def checkLigandHooked (PDB_parsed, list_atom_ligand_parsed):
    
    for atom_ligan in list_atom_ligand_parsed : 
        for atom_pdb in PDB_parsed : 
            if not atom_pdb["resName"] == atom_ligan["resName"] : 
                if (distanceTwoatoms(atom_pdb, atom_ligan)) < 1.5 : 
                    return 1
    return 0




def separateByLigand (l_atom_ligand, debug = 0) :
    """
    Separate list atoms ligand with same name or ID by atomic position
    args: -> list atoms ligand
    return: -> list of list with several atom by ligand
    """
    if debug : 
        print "--- check ligand ---"
        print  "NB atoms: ", len (l_atom_ligand)
        print "First atoms: ", l_atom_ligand[0]
        print "l209 parsePDB"
        print "--------------------"
    
    ####################################
    # ligand separed by ID and num res #
    ####################################
    if l_atom_ligand == [] :
        return []
    
    
    l_code_ligand = []
    for atom_ligand in l_atom_ligand : 
        code_ligand = atom_ligand["chainID"] + "_" + str(atom_ligand["resSeq"])
        if not code_ligand in l_code_ligand : 
            l_code_ligand.append (code_ligand)
    
    if len (l_code_ligand) == 1 : 
        return [l_atom_ligand]
    
    else : 
        # control number of atom but eg DEX case where same ligand without number of atoms
#         if len (l_atom_ligand) % len (l_code_ligand) != 0 : 
#             print "- ligand several ligand without same number of atoms"
#             print "l229 parsePDB"
#             print "----------------"
        
        l_out = []
        # separe by ligand, list of list
        for code_ligand in l_code_ligand : 
            l_atom_ligand_out = []
            for atom_ligand in l_atom_ligand : 
                code_atom = atom_ligand["chainID"] + "_" + str(atom_ligand["resSeq"])
                if code_atom == code_ligand : 
                    l_atom_ligand_out.append (atom_ligand)
            l_out.append (l_atom_ligand_out)
        
        return l_out
                
    
#     # First -> try with ligand name -> append condition case of all ligand have the same ID
#     d_resSeq = {}
#     for atom in l_atom_ligand : 
#         res_seq = atom["resSeq"]
#         # check resSeq empty
#         if res_seq == "" : 
#             continue
#         
#         # build dictionnary with resSeq -> keys
#         if not res_seq in d_resSeq.keys () : 
#             d_resSeq[res_seq] = []
#         d_resSeq[res_seq].append (atom)
#     
#     l_lig_atom = d_resSeq.values ()
#     l_nb_atom = []
#     for lig_atom in l_lig_atom : 
#         l_nb_atom.append(len(lig_atom))
#     l_nb_atom = list(set(l_nb_atom))
#     if len(l_nb_atom) == 1 :
#         return l_lig_atom
#     
#     # case with connect matrix
#     l_works_serial2atom = retrieveListAtomID (l_atom_ligand)
#     list_nb_atom = []
#     for atom_serial in l_works_serial2atom : 
#         list_nb_atom.append(len(atom_serial))
#     list_nb_atom = list(set(list_nb_atom))
#     
#     if len(list_nb_atom) != 1 : 
#         print "ERROR, retrieve ligand" 
#         print list_nb_atom
#         
#     if debug : 
#         print l_works_serial2atom, "List of list serial ligand"
#         for atom_serial in l_works_serial2atom : 
#             print len(atom_serial), "len of ligand"
#     
#     while len (l_atom_ligand) != 0 :
# #        if debug : print l_atom_ligand, "LIST ATOM LIGAND"
#         for list_atom_serial in l_works_serial2atom :
# #            if debug : print list_atom_serial, "LIST SERIAL"
#             if l_atom_ligand[0]["serial"] in list_atom_serial :
#                 list_atom_serial.remove (l_atom_ligand[0]["serial"]) 
#                 list_atom_serial.append (l_atom_ligand[0])
#                 del l_atom_ligand[0]
#                 
#                 
#     return l_works_serial2atom


    
def retrieveListAtomID (list_atom):
    """
    Retrieve list ID atoms
    args: -> list atoms
    return: list serial
    """
    
    list_out = []
    list_atom_temp = deepcopy(list_atom)
    #except : return []
    
    while len (list_atom_temp) != 0 : 
        list_out.append (retrieveSerialLigand(list_atom_temp))
    
    return list_out
    
    
def retrieveSerialLigand (list_atom):
    """
    Retrieve serial atom by ligand
    args: -> list atoms
    return: -> list serial atoms
    """
    
    list_out = list_atom[0]["connect"]
    del list_atom[0]
    nb_atom = len (list_atom)
    
    validate = 100
    while validate != 0 :
        nb_atom_temp = nb_atom 
        i = 0
        while i < nb_atom :
            if list_atom[i]["serial"] in list_out :
                list_out = list_out + list_atom[i]["connect"]
                del list_atom[i]
                nb_atom = nb_atom - 1
            else :
                i = i + 1
        validate = nb_atom - nb_atom_temp 
    return list(set(list_out))


def arrangeResidues(list_atoms) : 
    dico_res = {}
    for atom in list_atoms : 
        resID = atom["resSeq"]
        if not resID in dico_res.keys () : 
            dico_res[int(resID)] = []
        dico_res[resID].append (atom)
    
    return dico_res
        
        
        
def recountRes (l_at_parsed):       
     
    c_at = 0
    c_temp = -1
    for atom_parsed in l_at_parsed:
        if atom_parsed["resSeq"] != c_temp : 
            c_temp = atom_parsed["resSeq"]
            c_at = c_at + 1
            atom_parsed["resSeq"] = c_at
        else : 
            atom_parsed["resSeq"] = c_at
        
    
 
def retrieveSeq (p_PDB) : 
    
    print p_PDB, "###################"
    l_atom =  loadCoordSectionPDB(p_PDB, "ATOM")
    print l_atom
    
    l_temp = []
    s_out = ""
    for atom in l_atom : 
        res_ID = str (atom["resSeq"]) + "_" + str (atom["resName"]) + "_" + str (atom["chainID"])
        print res_ID
        if not res_ID in l_temp : 
            l_temp.append (res_ID)
            s_out = s_out + str(tool.transformAA (atom["resName"]))
            
            
    return s_out
 
def buildMatrixConnect (l_atom_parsed):
    
    for atomLigand in l_atom_parsed:
        atomLigand["connect"] = []

    for atomLigand in l_atom_parsed:
        atomLigand["connect"].append(atomLigand["serial"])
        for atomPDB in l_atom_parsed:
            distance = distanceTwoatoms(atomLigand, atomPDB)
            if distance < 1.8 and distance != 0.1 and distance != "ERROR":
                if not atomPDB["serial"] in atomLigand["connect"]:
                    atomLigand["connect"].append(atomPDB["serial"])

def getResidues(l_atom_binding, l_atom_complex_parsed) : 
    """
    """
    dico_code = {"S":"SER", "T":"THR", "N":"ASN", "Q":"GLN", "E":"GLU", "D":"ASP", "K":"LYS", "R":"ARG", "H":"HIS", "M":"MET", "C":"CYS", "W":"TRP", "F":"PHE", "Y":"TYR", "A":"ALA", "V":"VAL", "L":"LEU", "I":"ILE", "P":"PRO", "G":"GLY"}
    # -> retrieve residue
    memory_res = []
    l_res = []
    for atom_bs in l_atom_binding : 
        # case where consider all atom, het and atom
        if not atom_bs["resName"] in dico_code.values () : 
            l_res.append (atom_bs)
            continue
        for atom_complex in l_atom_complex_parsed : 
            if atom_bs["x"] == atom_complex["x"] and atom_bs["y"] == atom_complex["y"] and atom_bs ["z"] == atom_complex["z"] : 
                chain_id = atom_complex["chainID"]
                res_num = atom_complex["resSeq"]
                res_code = str(chain_id) + "_" + str(res_num)
                if not res_code in memory_res : 
                    memory_res.append (res_code)
                    for atom_complex2 in l_atom_complex_parsed : 
                        if atom_complex2["chainID"] == chain_id and atom_complex2["resSeq"] == res_num : 
                            l_res.append (atom_complex2)
    
    return l_res  
 
 
 
############## 
## calcul   ##
##############
 

def distanceTwoatoms(atom1, atom2):  ##############to review
    '''calculate distance of 2 atoms
    in : - atom1 structure
         - atom2 structure
    out : distance -> float
          100 if impossible calcul'''

    try:
        x1 = float(atom1['x'])
        x2 = float(atom2['x'])
        xd = x2 - x1

        y1 = float(atom1['y'])
        y2 = float(atom2['y'])
        yd = y2 - y1

        z1 = float(atom1['z'])
        z2 = float(atom2['z'])
        zd = z2 - z1

        return sqrt(xd * xd + yd * yd + zd * zd)
    except:
        return "ERROR"

 
    
def scalar(vector1, vector2):
    
    x = vector1[0] * vector2[0]
    y = vector1[1] * vector2[1]
    z = vector1[2] * vector2[2]

    return x + y + z


def angleVector(pointD, pointCentral, pointG):
    
    vectorNC1 = vectoriel(pointCentral, pointD)
    vectorNC2 = vectoriel(pointCentral, pointG)
    normeNC1 = normeVector(pointCentral, pointD)
    normeNC2 = normeVector(pointCentral, pointG)

    scalarNC1NC2 = scalar(vectorNC1, vectorNC2) 

    try :
        alpha = degrees(acos(scalarNC1NC2 / (normeNC1 * normeNC2)))
        return alpha
    except :
        return 0

def normeVectoriel (listVectoriel1, listVectoriel2):
    
    terme1 = pow(listVectoriel1[1] * listVectoriel2[2] - listVectoriel1[2] * listVectoriel2[1], 2)
    terme2 = pow(listVectoriel1[2] * listVectoriel2[0] - listVectoriel1[0] * listVectoriel2[2], 2)
    terme3 = pow(listVectoriel1[0] * listVectoriel2[1] - listVectoriel1[1] * listVectoriel2[0], 2)
    
    return sqrt (terme1 + terme2 + terme3)

def vectoriel (N, C1):
    
    
    NC1x = C1["x"] - N["x"]
    NC1y = C1["y"] - N["y"]
    NC1z = C1["z"] - N["z"]
    
    return [NC1x, NC1y, NC1z]


def normeVector (point1, point2):
    
    terme1 = pow(point1["x"] - point2["x"], 2)
    terme2 = pow(point1["y"] - point2["y"], 2)
    terme3 = pow(point1["z"] - point2["z"], 2)
    
    return sqrt(terme1 + terme2 + terme3)


def resolution(p_PDB):
    """Retrieve by PDB file the resolution if X-ray methods
    in : name of pdb file
    out : resolution -> format string"""

    filin = open (p_PDB, "r")
    fileLines =filin.readlines()
    filin.close () 

    for line in fileLines:
        if search("^REMARK   2 RESOLUTION", line):
            line = sub('[ ]{2,}', ' ', line)
            try:
                resolution = line.split(" ")[3].replace (" ", "")
            except:
                resolution = "NA"
            return resolution
    return "NA"



def nameProtein (p_filin) : 
    
    filin = open (p_filin, "r")
    fileLines =filin.readlines()
    filin.close () 

    for line in fileLines:
        if search('^COMPND   2', line):
            
            try: protein_name = line.strip().split(': ')[1].strip(';')
            except : return "NA"
            return protein_name 
            
    return "NA"    
    
    
def UniProtID (p_filin):   
    
    filin = open (p_filin, "r")
    fileLines =filin.readlines()
    filin.close () 

    for line in fileLines:
        if search('^DBREF', line):
            try : uniprot_id = line.split()[7]
            except : uniprot_id = "NA"
            return uniprot_id 
    return "NA"
    
def keywords (p_filin): 
    filin = open (p_filin, "r")
    fileLines =filin.readlines()
    filin.close () 

    for line in fileLines:
        if search('^KEYWDS', line):
            try : kwords = line.strip()
            except : kwords = "NA"
            
            return kwords 
        
    return "NA"


def retrieveListIon(l_atom_parsed) : 
    
    l_metal = ["B", "F", "I", "K", "V", "W", "Y","AG", "AL", "AR" ,"AU", "BA", "BE", "BR", "CA","CD","CE","CF","CL","CO","CR","CS","CU","EU","FE","GA","GD","HE","HF","HG","IN","IR","KR","LA","LI","LU" ,"MG","MN" ,"MO" ,"NA","ND","NE","NI","OS","PB","PD","PR","PT","RB","RE","RU","SB","SE","SI","SM","SR","TA","TB","TE","TL","XE","YB","ZN","ZR"]
    l_out = []
    
    for atom_parsed in l_atom_parsed : 
        if atom_parsed["resName"] in l_metal and not atom_parsed["resName"] in l_out: 
            l_out.append (atom_parsed["resName"])
    
    return l_out
        
        
        
def computeBS (p_protein, p_ligand, thresold = 4.50, option_onlyATOM = 0):
    
    
    l_atom_BS = []
    l_atom_pr = loadCoordSectionPDB(p_protein)
    l_atom_lig = loadCoordSectionPDB(p_ligand)
    
    for atom_pr in l_atom_pr : 
        for atom_lig in l_atom_lig : 
            if distanceTwoatoms(atom_pr, atom_lig) <= thresold : 
                l_atom_BS.append (atom_pr)
    
    
    l_res_BS = getResidues(l_atom_BS, l_atom_pr)
    
    return l_res_BS
    
    
          
    
    
    
    

