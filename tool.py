from os import path
from re import search
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
            

            
def parseLigandPDBList (p_filin):
    """Load file result, PDB associated at ligand
    in: Name of file
    out: list ligand with PDB associated"""


    fileOpen = open(p_filin, "r")

    lineFile = fileOpen.readlines()
    outFile = {}

    for line in lineFile:
        line = line.split("\t")
        PDB = line[0]

        listLigand = line[1]
        listLigand = listLigand.split("\n")[0]
        listLigand = listLigand.split(" ")
        if listLigand == []:
            continue

        for ligand in listLigand:
            if ligand != "":
                try:
                    outFile[ligand].append(PDB)
                except:
                    outFile[ligand] = []
                    outFile[ligand].append(PDB)

    return outFile
    
    
def fusionchainfasta(p_fasta) : 
    
    filin = open (p_fasta, "r")
    l_line_fasta = filin.readlines ()
    filin.close ()
    seq = ""
    comment = l_line_fasta [0][0:5]
    for line_fasta in l_line_fasta : 
        if not search("^>", line_fasta) :
            seq = seq + line_fasta.strip ()
    
    # format 80
    width= 80
    seq_write = "\n".join( [seq[i:i+width] for i in xrange(0,len(seq),width)] )
    
    filout = open (p_fasta, "w")
    filout.write (comment + "\n" + seq_write)
    filout.close ()
            

def closeDicoFile (d_filout_lig) : 
    
    for k in d_filout_lig.keys () : 
        try : d_filout_lig[k].close ()
        except : pass
        

def NameFamily (name_family) : 

    if name_family == "phosphatase":
        return "PT"
    elif name_family == "phosphorylase":
        return "PL"
    elif name_family == "kinase":
        return "KS" 
    elif name_family == "helicase":
        return "HL" 
    elif name_family == "myosine":
        return "MY" 
    elif name_family == "transferase":
        return "TF" 
    elif name_family == "synthetase":
        return "SY" 
    elif name_family == "hydrolase":
        return "HD" 
    elif name_family == "ligase":
        return "LG"
    elif name_family == "ATPase":
        return "AT" 
    elif name_family == "carboxylase":
        return "CB" 
    elif name_family == "actin":
        return "AC" 
    elif name_family == "HSP":
        return "HP" 

    return "OT"    
    
    
    
    
    
    
    

