"""
BORREL Alexandre
04-2012
Analysis PDB file
"""
# global module
from re import search
from math import sqrt
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
    try :atom["serial"] = int(line[6:11].replace (" ", ""))
    except :line[6:11].replace (" ", "")
    atom["name"] = line[12:16].replace (" ", "")
    atom["char"] = line[16]
    atom["resName"] = line[17:20].replace (" ", "")
    atom["chainID"] = str(line[21])
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


def loadCoordSectionPDB (path_PDB_file, section = "", debug = 1):
    """
    Retrieve every atom in cordiante section. If it is NMR complex
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
                list_atom.append (lineCoords(line_PDB))
        else : 
            if search ("^" + section, line_PDB)  : 
                list_atom.append (lineCoords(line_PDB))
                
#    if debug : 
#        print "TEST"            
#        print len (list_atom)
#    
    return list_atom
    



def distanceTwoatoms(atom1, atom2):##############to review
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
        return 100



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
 
 
 
 
    
#pdb_parsed = loadCoordSectionPDB ("1B8O.pdb") 
#print pdb_parsed
#list_ligand = retrieveLigand  (pdb_parsed, "IMH")
#print "----------"
#print list_ligand
#writePDBfile.coordinateSection ("test", list_ligand,"HETATM")
