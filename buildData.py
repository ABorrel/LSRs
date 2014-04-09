import pathManage
import tool
import downloadFile
import runOtherSoft
import parsePDB
import writePDBfile
import parseWater
from os import path




# rebuild the dataset
def builtDatasetGlobal (p_list_ligand, substruct, thresold_RX = 2.5, thresold_IDseq = 90, thresold_blast = 1e-4, verbose = 1 ):
    
    # directory with dataset
    p_dir_dataset = pathManage.dataset(substruct)
    # directory with result
    p_dir_result = pathManage.result(substruct + "/datasetBuilding")
    
    # first extract reference
    d_dataset = extractReference (p_list_ligand, p_dir_dataset, p_dir_result, substruct)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # select reference
    # remove RX and same chain
    p_dir_align = pathManage.result(substruct + "/datasetBuilding/aligmentRef")
    filterReferenceByOne (d_dataset, p_dir_align, thresold_RX = 2.5)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # conserve only unique protein
    filterGlobalDataset (d_dataset, p_dir_align)
    
    if verbose : toolViewStructDataset (d_dataset)

    
    

def extractReference (p_list_ligand, p_dir_dataset, p_dir_result, substruct):    
    
    # struct reference
    d_dataset = {}
    
    # retrieve list of ligand in PDB
    d_ligand = tool.parseLigandPDBList (p_list_ligand)
    
    # download PDB and fasta associated
    l_p_PDB = []
    l_p_fasta = []
    for PDB_ID in d_ligand[substruct][0:50] :
        PDB_ID = PDB_ID.upper() 
        p_pdb = downloadFile.importPDB(PDB_ID, p_dir_dataset, dir_by_PDB = 1, debug = 1)
        p_fasta = downloadFile.importFasta(PDB_ID, p_dir_dataset, dir_by_PDB = 1, debug = 1)
        
        if p_pdb != 0 and p_fasta != 0 : 
            l_p_pdb_chain = separeByChain (p_pdb)
            l_p_fasta_chain = separeChainFasta(p_fasta)
            d_dataset[PDB_ID] = {}
            d_dataset[PDB_ID] ["p_pdb"] = p_pdb
            d_dataset[PDB_ID] ["p_fasta"] = p_fasta
            d_dataset[PDB_ID] ["p_pdb_chain"] = l_p_pdb_chain
            d_dataset[PDB_ID] ["p_fasta_chain"] = l_p_fasta_chain
            d_dataset[PDB_ID] ["conserve"] = 1
        
    # plot resolution
    p_file_RX = p_dir_result + "resolution_ref.txt"
    file_RX = open (p_file_RX, "w")
    for PDB_ID in d_dataset.keys () : 
        RX = parsePDB.resolution(d_dataset[PDB_ID]["p_pdb"])
        d_dataset[PDB_ID] ["RX"] = float(RX)
        file_RX.write (PDB_ID + "\t" + str (RX) + "\n") 
    file_RX.close ()
    
    runOtherSoft.Rhistogram (p_file_RX, "RX_ref_no_filter")
    return d_dataset
    

def filterReferenceByOne (d_dataset, pr_aligmenent_water, thresold_RX = 2.5):
    
    for PDB_ref in d_dataset.keys () : 
        if float(d_dataset[PDB_ref]["RX"]) >= thresold_RX : 
            d_dataset[PDB_ref] ["conserve"] = 0
        else : 
            # filter redundency by PDB // run PDB
            result_align = calculIdenticWaterCross(d_dataset[PDB_ref]["p_fasta_chain"], pr_aligmenent_water)
            # select best by PDB ID -> remove chain if 100% sequence identity
            selectBestFilePDBFasta (d_dataset[PDB_ref], result_align)
            
            

def separeByChain (path_PDB_file):
    """
    separe PDB file by chain
    """
    file_PDB_parsed = parsePDB.loadCoordSectionPDB(path_PDB_file, section="")
    
    l_file = []
    file_open_write = {}
    for atom_PDB in file_PDB_parsed : 
        chain = atom_PDB["chainID"]
        if not chain in file_open_write.keys () : 
            file_open_write [chain] = open(path_PDB_file[0:-4] + "_" + chain + ".pdb", "w")
            l_file.append (path_PDB_file[0:-4] + "_" + chain + ".pdb")
            writePDBfile.coordinateStructure(atom_PDB, file_open_write [chain] )
        else : 
            writePDBfile.coordinateStructure(atom_PDB, file_open_write [chain] )
    
    # close files
    for chain in file_open_write.keys () : 
        file_open_write[chain].close ()
        
    return l_file
        

def separeChainFasta (p_fasta, debug = 0) : 
    """
    Separe fasta by chain to run water for retrieve identity
    """
    l_filout = []
    
    if not path.exists(p_fasta) : 
        print "ERROR -> File fasta not exist"
        return []
    else : 
        filin = open (p_fasta, "r")
        element_file = filin.read ()
        if debug : print element_file
        filin.close ()
        element_file = element_file.split (">")
        if debug :print element_file
        
        for seq_chain in element_file[1:] : 
            chain = seq_chain[5]
            if debug : print chain, "Chain fasta retrieve"
            path_filout = p_fasta[0:-6] + "_" + str (chain) + ".fasta"
            l_filout.append (path_filout)
            filout = open (path_filout,"w")
            filout.write (">" + seq_chain)
            filout.close ()
            
    return l_filout




def calculIdenticWaterCross(l_p_fasta, pr_alignement) : 
    """
    Calcul identity structure
    args: -> list of fasta files
    return: -> dictionary with identity
    """
    if len (l_p_fasta) == 1 : 
        return []

 
    dico_out = {}
    
    if len (l_p_fasta) == 2 : 
        first_key = l_p_fasta[0].split("/")[-1][0:-6] 
        dico_out[first_key] = {}
        if not path.exists( pr_alignement + first_key + "_" + l_p_fasta[1].split("/")[-1][0:-6] + ".water") : 
            path_file_water = runOtherSoft.water(l_p_fasta[0],  l_p_fasta[1], pr_alignement + first_key + "_" + l_p_fasta[1].split("/")[-1][0:-6] + ".water")
        else : 
            path_file_water = pr_alignement + first_key + "_" + l_p_fasta[1].split("/")[-1][0:-6] + ".water"
        dico_out[first_key] [l_p_fasta[1].split("/")[-1][0:-6]] = parseWater.waterFile(path_file_water)[3]
    else :
        nb_fasta = len (l_p_fasta)
        i = 0
        while i < nb_fasta-1 : 
            first_key = l_p_fasta[i].split("/")[-1][0:-6]
            dico_out[first_key] = {}
            j = i + 1
            while j < nb_fasta : 
                if not path.exists( pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".water") :
                    path_file_water = runOtherSoft.water(l_p_fasta[i], l_p_fasta[j], pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".water")
                dico_out[first_key] [l_p_fasta[j].split("/")[-1][0:-6]] = parseWater.waterFile(pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".water")[3]
                j = j + 1
            i = i + 1
    
    return dico_out    



def selectBestFilePDBFasta (d_PDB, result_align) : 
    
    # case with only one PDB
    d_PDB["best"] = {}
    if result_align == [] : 
        d_PDB["best"]["PDB"] = d_PDB["p_pdb"]
        d_PDB["best"]["fasta"] = d_PDB["p_fasta"]
    
    else :
        for primary_key in result_align.keys () : 
            for secondary_key in result_align[primary_key] .keys (): 
                if result_align[primary_key][secondary_key] != "100.0%" : 
                    d_PDB["best"]["PDB"] = d_PDB["p_pdb"]
                    d_PDB["best"]["fasta"] = d_PDB["p_fasta"]
                    return 
        
        d_PDB["best"]["PDB"] = d_PDB["p_pdb_chain"][0]
        d_PDB["best"]["fasta"] = d_PDB["p_fasta_chain"][0]
        
    
    
def filterGlobalDataset (d_dataset, p_dir_align) : 
    
    # retrieve list of sequences
    l_all_fasta = []
    
    for PDB in d_dataset.keys () : 
        if d_dataset[PDB]["conserve"] == 1 : 
            l_all_fasta.append (d_dataset[PDB]["best"]["fasta"])
             
    # alignmenet
    result_align = calculIdenticWaterCross(l_all_fasta, p_dir_align)
    
    # remove same protein
    for primary_key in result_align.keys () : 
            for secondary_key in result_align[primary_key] .keys (): 
                if result_align[primary_key][secondary_key] == "100.0%" : 
                    if d_dataset[primary_key[0:4]]["conserve"] == 0 or d_dataset[secondary_key[0:4]]["conserve"] == 0 : 
                        continue
                    else : 
                        if d_dataset[primary_key[0:4]]["RX"] < d_dataset[secondary_key[0:4]]["RX"] : 
                            d_dataset[secondary_key[0:4]]["conserve"] = 0
                        else : 
                            d_dataset[primary_key[0:4]]["conserve"] = 0
                        
    
    
def toolViewStructDataset (d_dataset):    
    
    print "+++++++++++++++++++++++++++++++++++++"
    for k in d_dataset.keys () : 
        print k, d_dataset[k]["conserve"], d_dataset[k]["RX"]
    print "--------------------------------------"    
    
    
   