from re import search
import parsePDB

def coordinateSection (f_write, listAtom, recorder, header = "", connect_matrix = 0):
    """
    Write list atom in PDB file
    in: list atoms, name of file out
    out: write file
    
    arg -> header = 0 not header in file
    """
    
    
    if type (f_write) == str : 
        filout = open(f_write, "w")
    elif  type (f_write) == file  : 
        filout = f_write
    else : 
        print "======ERROR========"
        print "flux writing ERROR"
        print "l17 - write PDB file"
        print "==================="
        return 0
        
    if header != 0 : 
        filout.write ("HEADER " + str (header) + "\n")
    for atom in listAtom : 
        coordinateStructure(atom, recorder, filout)
    
    if connect_matrix :
        if  listAtom [0]["connect"] == [] : 
            parsePDB.buildMatrixConnect (listAtom)
        for atom in listAtom : 
            connect(atom, filout)
    filout.write("END\n")
    
    if type (f_write) == str : 
        filout.close ()
    else : 
        pass
    
    return f_write








def coordinateStructure(atom, recorder, fileWrite):
    """Write coordinate lines and format structure
    in: atom with coordinate, type of recorder, file write
    out: write line in file"""
    
    # check \n PDB not formating
    for key_PDB in atom.keys () : 
        try : atom[key_PDB] = atom[key_PDB].replace ("\n", "")
        except : pass
    
    
    try : atomSerial = atom["serial"]
    except :atomSerial = ""
    try : nameAtom = atom["name"]
    except :nameAtom = ""
    try : resName = atom["resName"]
    except :resName = ""
    try : ChainID = atom["chainID"]
    except :ChainID = ""
    try : iCode = atom["iCode"]
    except :iCode = ""
    try : resNumber = atom["resSeq"]
    except : resNumber = ""
    try : char = atom["char"]
    except : char = ""
    try :  x = atom["x"]
    except : x = ""
    try : y = atom["y"]
    except : y = ""
    try : z = atom["z"]
    except :z = ""
    try : element = atom["element"]
    except : element = ""
    try : charge = atom["charge"]
    except : charge = ""
    try : occupancy = atom["occupancy"]
    except : occupancy = ""
    try : tempFactor = atom["tempFactor"]
    except : tempFactor = ""
    
    coordinate(recorder, atomSerial, nameAtom, char, resName, ChainID, resNumber, iCode, x, y, z, occupancy, tempFactor, element, charge, fileWrite)


def coordinate(recorder, serialAtom, nameAtom, char, resName, ChainID, resNumber, iCode, x, y, z, occupancy, tempFactor, element, charge, fileWrite):
    """Write coordinate section line
    in: type recorder, serial atom, name atom (PDB database format), ligand or amino acid where atom is. chain ID, serial atom in structure, ID code, x, y, z coordinate, value of occupancy, temperature factor, charge, flux write
    out: write in flux """
    
    line = formatRecorder(recorder, 6) + formate(serialAtom, 5) + " " + formate(nameAtom, 4) + formate(char, 1) + formate(resName, 3) + " " + formate(ChainID, 1) + formate(resNumber, 4) + formate(iCode, 1) + "   " + formate(formatCoord(x), 8) + formate(formatCoord(y), 8) + formate(formatCoord(z), 8) + formate(occupancy, 6) + formate(tempFactor, 6) + "          " + formate(element, 2) + formate(charge, 2)
    
    fileWrite.write(line)
    fileWrite.write("\n")


def connect(atom, fileWrite) : 
    """write in file the connect section for one atom
    in: atom, flux write
    out: write in flux write"""
    
    
    line = "CONECT"
    for serial in atom["connect"] : 
        line = line + formate(serial, 5)

    fileWrite.write(line + "\n")

def formatRecorder(recorder, lengthMax): 
    
    while len(recorder) < lengthMax :
        recorder = str(recorder) + " "
    
    return recorder 
        
    
def formate(inVariable, nbCarac):
    """formate string for write
    in: string or value integers or float
    out: string formated for write"""
    
    inVariable = str(inVariable)
    
    lenStr = len(inVariable)
    
    caracAppend = nbCarac - lenStr
    
    if caracAppend == 0 : 
        return inVariable
    elif caracAppend > 0 :
        for i in range(0, caracAppend) : 
            inVariable = " " + inVariable 
        
        return inVariable
    else:
        return (inVariable[0:nbCarac - 1])


def formatCoord(float):
    """Format float in coordinate section
    in: float
    out: string formated"""
    
    return str("%.3f" % float)    
