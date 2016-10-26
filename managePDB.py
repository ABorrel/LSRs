from os import system, listdir, path
import parsePDB


def formatPDBDatabase(pr_PDB):
    """Manage global PDB database
    manage PDB database, decompress the filePDB and move and rename filePDB with .pdb extension"""
    
    l_folders = listdir(pr_PDB)
    for folder in l_folders :
        l_filePDB = listdir(pr_PDB + folder)   
        for filePDB in l_filePDB :
            p_filePDB = pr_PDB + folder + '/' + filePDB
            cmd_gunzip = "gunzip " + p_filePDB
            system(cmd_gunzip)  ####run decompress
            p_filePDB = p_filePDB[0:-3]
            namePDB = p_filePDB[-8:-4] + '.pdb'
            p_filout = p_filePDB[0:-14] + namePDB
            cmd_mv = "mv " + p_filePDB + ' ' + p_filout
            system(cmd_mv)  ####run move filePDB
    
    
    for folder in l_folders:
        folder = pr_PDB + folder
        cmd_rm = "rm -r " + folder
        system(cmd_rm)  ####run remove repertory



def retriveListPDB (pr_PDB):
    
    l_files = listdir(pr_PDB)
        
    l_p_pdb = []
    for file_PDB in l_files : 
        if file_PDB[-4:] == ".pdb" : 
            l_p_pdb.append (pr_PDB + file_PDB)
            
    return l_p_pdb


def searchLigands(prPDB, prresult):
    '''search ligands in PDB database
    out : list of ligands with PDB files associated'''

    print "Start Search Ligand In PDB file"
    pfilout = prresult + "resultLigandInPDB"
    # control file exist
    if path.exists(pfilout) and path.getsize(pfilout) != 0:
        return pfilout

    l_PDB = retriveListPDB(prPDB)

    filout = open(pfilout, "w")

    for PDBid in l_PDB:
        llig = parsePDB.retrieveListLigand(prPDB + PDBid.lower() + ".pdb")
        if llig != []:
            filout.write(PDBid + "\t" + " ".join(llig) + "\n")
        else:
            continue

    filout.close()
    return pfilout