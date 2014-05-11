import pathManage
import tool
import downloadFile
import runOtherSoft
import parsePDB
import writePDBfile
import parseWater
import RunBlast
import downloadFile
from os import path, removedirs, remove, listdir
from shutil import rmtree
import analysis




# rebuild the dataset
def builtDatasetGlobal (p_list_ligand, substruct, thresold_RX = 2.5, thresold_IDseq = 90, thresold_blast = 1e-4, verbose = 1 ):
    
    # directory with dataset
    p_dir_dataset = pathManage.dataset(substruct)
    # directory with result
    p_dir_result = pathManage.result(substruct + "/datasetBuilding")
    
    # first extract reference
    d_dataset = extractReference (p_list_ligand, p_dir_dataset, p_dir_result, substruct)
    
    # file with name and family
    analysis.familyPDBRef (d_dataset, p_dir_dataset + "family_PDB.txt")
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # select reference
    # remove RX and same chain
    p_dir_align = pathManage.result(substruct + "/datasetBuilding/aligmentRef")
    filterReferenceByOne (d_dataset, p_dir_align, substruct, thresold_RX = 2.5)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # conserve only unique protein
    filterGlobalDataset (d_dataset, p_dir_align)
    
    if verbose : toolViewStructDataset (d_dataset)

    # run blast by sequence conserved 
    p_dir_blast = pathManage.result(substruct + "/datasetBuilding/blast")
    RunBlast.globalRun (d_dataset, p_dir_blast)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # filter by e-value and RX
    filterBlastResult (d_dataset, p_dir_dataset,substruct, thresold_RX = 2.5, thresold_blast = 1e-4)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # clean folder dataset
    cleanFolderDataset (d_dataset, p_dir_dataset)
    
    

def extractReference (p_list_ligand, p_dir_dataset, p_dir_result, substruct):    
    
    # struct reference
    d_dataset = {}
    
    # retrieve list of ligand in PDB
    d_ligand = tool.parseLigandPDBList (p_list_ligand)
    
    # download PDB and fasta associated
    l_p_PDB = []
    l_p_fasta = []
    for PDB_ID in d_ligand[substruct] :
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
        try : d_dataset[PDB_ID] ["RX"] = float(RX)
        except : d_dataset[PDB_ID] ["RX"] = 100.0
        file_RX.write (PDB_ID + "\t" + str (RX) + "\n") 
    file_RX.close ()
    
    runOtherSoft.Rhistogram (p_file_RX, "RX_ref_no_filter")
    return d_dataset
    

def filterReferenceByOne (d_dataset, pr_aligmenent_water, substruct, thresold_RX = 2.5):
    
    for PDB_ref in d_dataset.keys () : 
        if float(d_dataset[PDB_ref]["RX"]) >= thresold_RX : 
            d_dataset[PDB_ref] ["conserve"] = 0
        else : 
            # filter redundency by PDB // run PDB
            result_align = calculIdenticWaterCross(d_dataset[PDB_ref]["p_fasta_chain"], pr_aligmenent_water)
            # select best by PDB ID -> remove chain if 100% sequence identity
            selectBestFilePDBFasta (d_dataset[PDB_ref], result_align, substruct)
            
            

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
                if path.exists( pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".water") : # case DNA sequence in fasta
                    dico_out[first_key] [l_p_fasta[j].split("/")[-1][0:-6]] = parseWater.waterFile(pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".water")[3]
                j = j + 1
            i = i + 1
    
    return dico_out    



def selectBestFilePDBFasta (d_PDB, result_align, substruct) : 
    
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
                    # fusion fasta for blast remove header
                    tool.fusionchainfasta(d_PDB["p_fasta"])
                    return 
        
        # control substructure presented // case identity 100%
        nb_chain = len (d_PDB["p_pdb_chain"])
        i = 0
        while i < nb_chain : 
            l_ligand = parsePDB.retrieveListLigand(d_PDB["p_pdb_chain"][i])
            if substruct in l_ligand : 
                d_PDB["best"]["PDB"] = d_PDB["p_pdb_chain"][i]
                d_PDB["best"]["fasta"] = d_PDB["p_fasta_chain"][i]
                return
            else : 
                i = i + 1
        
    
    
def filterGlobalDataset (d_dataset, p_dir_align) : 
    
    # retrieve list of sequences
    l_all_fasta = []
    
    for PDB in d_dataset.keys () : 
        if d_dataset[PDB]["conserve"] == 1 : 
            if not "fasta" in d_dataset[PDB]["best"] : 
                d_dataset[PDB]["conserve"] = 0
                continue
            else : 
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
                        


def filterBlastResult (d_dataset, p_dir_dataset, sustruct, thresold_RX = 2.5, thresold_blast = 1e-4, debug = 0) : 
    
    """
    Filter resolution PDB
    Filter evalue
    """
    
    for pdb_ref in d_dataset.keys () : 
        if d_dataset[pdb_ref]["conserve"] == 0 : continue
        for pdb_blast_chain in d_dataset[pdb_ref]["align"].keys () : 
            # filter e.value
            if debug == 1 : print d_dataset[pdb_ref]["align"][pdb_blast_chain], pdb_ref, thresold_blast
            # remove thresold and reference cleanner
            if d_dataset[pdb_ref]["align"][pdb_blast_chain] <= thresold_blast and not pdb_blast_chain[0:4] in d_dataset.keys (): 
                # dowload PDB files
                pdb_blast = pdb_blast_chain[0:4]
                p_pdb_blast = downloadFile.importPDB(pdb_blast, p_dir_dataset + pdb_ref + "/", dir_by_PDB=0)
                separeByChain(p_pdb_blast)
                try : RX = parsePDB.resolution(p_pdb_blast)
                except : RX = 100.0
                l_ligand = parsePDB.retrieveListLigand(p_pdb_blast)
                if debug == 1 : print p_pdb_blast, l_ligand, RX
                # remove apo forms and remove not substiuant
                if l_ligand == [] or sustruct in l_ligand: 
                    continue
                
                print 
                
                # case RMN structure
                try : RX = float(RX)
                except: continue
                
                if debug == 1 :print "----", RX, thresold_RX, "----"
                if float(RX) <= thresold_RX : 
                    if not "blast" in d_dataset[pdb_ref].keys () : 
                        d_dataset[pdb_ref]["blast"] = [pdb_blast_chain]
                    else : 
                        d_dataset[pdb_ref]["blast"].append (pdb_blast_chain)
                    if debug == 1 :print d_dataset[pdb_ref]["blast"], pdb_ref
            
def cleanFolderDataset (d_dataset, p_dir_dataset) : 
    print d_dataset.keys ()
    
    for PDB_ref in  d_dataset.keys () : 
        if d_dataset[PDB_ref]["conserve"] == 0 or not "blast" in d_dataset[PDB_ref].keys (): 
            rmtree(p_dir_dataset + PDB_ref + "/")
            print "RM1", d_dataset[PDB_ref]["conserve"]
            continue
        else : 
            pr_ref = p_dir_dataset + PDB_ref + "/"
            l_file_dir_ref = listdir(p_dir_dataset + PDB_ref)
            for filein_rep in l_file_dir_ref : 
                
                print "*******"
                print pr_ref +  filein_rep 
                print d_dataset[PDB_ref]["best"]["PDB"], d_dataset[PDB_ref]["best"]["fasta"], d_dataset[PDB_ref]["blast"]
                print "*********"


                if pr_ref +  filein_rep in [d_dataset[PDB_ref]["best"]["PDB"], d_dataset[PDB_ref]["best"]["fasta"]] : 
                    continue
                elif filein_rep[0:-4] in d_dataset[PDB_ref]["blast"] : 
                    continue
                else : 
                    remove(pr_ref + filein_rep)
                
    
    
def toolViewStructDataset (d_dataset):    
    
    print "+++++++++++++++++++++++++++++++++++++"
    for k in d_dataset.keys () : 
#         print k, d_dataset[k]["conserve"], d_dataset[k]["RX"]
        
        if "best" in d_dataset[k] : 
            print "@@@@@@@@@@"
            print d_dataset[k]["best"]
            if "blast" in d_dataset[k]:
                print "****"
                print d_dataset[k]["blast"]
                print "****"
    print "--------------------------------------"    
    
    
   
