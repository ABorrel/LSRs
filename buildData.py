import pathManage
import tool
import downloadFile
import runOtherSoft
import parsePDB
import writePDBfile
import parseEMBOSS
import RunBlast
from os import path, remove, listdir
from re import search
from shutil import rmtree
import analysis




# rebuild the dataset
def builtDatasetGlobal (p_list_ligand, ligand_ID, thresold_RX = 2.5, thresold_blast = 1e-4, l_ligand_out = [], verbose = 1 ):
    
    # directory with dataset
    p_dir_dataset = pathManage.dataset(ligand_ID)
    # directory with result
    p_dir_result = pathManage.result(ligand_ID + "/datasetBuilding")
    
    # first extract reference
    d_dataset = extractReference (p_list_ligand, p_dir_dataset, p_dir_result, ligand_ID)
    
    # file with name and family
    analysis.familyPDBRef (d_dataset, p_dir_dataset + "family_PDB.txt")
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # select reference
    # remove RX and same chain
    p_dir_align = pathManage.result(ligand_ID + "/datasetBuilding/aligmentRef")
    filterReferenceByOne (d_dataset, p_dir_align, ligand_ID, thresold_RX = thresold_RX)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # conserve only unique protein
    filterGlobalDataset (d_dataset, p_dir_align)
    
    if verbose : toolViewStructDataset (d_dataset)

    # run blast by sequence conserved 
    p_dir_blast = pathManage.result(ligand_ID + "/datasetBuilding/blast")
    RunBlast.globalRun (d_dataset, p_dir_blast)
    
    if verbose : toolViewStructDataset (d_dataset)
    
    # filter by e-value and RX
    filterBlastResult (d_dataset, p_dir_dataset,ligand_ID, thresold_RX = thresold_RX, thresold_blast = thresold_blast, l_ligand_out = l_ligand_out)
    
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
        p_pdb = downloadFile.importPDB(PDB_ID, p_dir_dataset, dir_by_PDB = 1, debug = 1, dbPDB = "/home/borrel/PDB/" )
        p_fasta = downloadFile.importFasta(PDB_ID, p_dir_dataset, dir_by_PDB = 1, debug = 1, fastaGlobal = "/home/borrel/Yue_project/pdb_seqres.txt")
        
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
    

def filterReferenceByOne (d_dataset, pr_aligmenent_needle, substruct, thresold_RX = 2.5):
    
    for PDB_ref in d_dataset.keys () : 
        if float(d_dataset[PDB_ref]["RX"]) >= thresold_RX : 
            d_dataset[PDB_ref] ["conserve"] = 0
        else : 
            # filter redundency by PDB // run PDB
            result_align = calculIdenticNeedleCross(d_dataset[PDB_ref]["p_fasta_chain"], pr_aligmenent_needle)
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




def calculIdenticNeedleCross(l_p_fasta, pr_alignement) : 
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
        if not path.exists( pr_alignement + first_key + "_" + l_p_fasta[1].split("/")[-1][0:-6] + ".needle") : 
            path_file_needle = runOtherSoft.needle(l_p_fasta[0],  l_p_fasta[1], pr_alignement + first_key + "_" + l_p_fasta[1].split("/")[-1][0:-6] + ".needle")
        else : 
            path_file_needle = pr_alignement + first_key + "_" + l_p_fasta[1].split("/")[-1][0:-6] + ".needle"
        dico_out[first_key] [l_p_fasta[1].split("/")[-1][0:-6]] = parseEMBOSS.embossFile(path_file_needle)[3]
    else :
        nb_fasta = len (l_p_fasta)
        i = 0
        while i < nb_fasta-1 : 
            first_key = l_p_fasta[i].split("/")[-1][0:-6]
            dico_out[first_key] = {}
            j = i + 1
            while j < nb_fasta : 
                if not path.exists( pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".needle") :
                    path_file_needle = runOtherSoft.needle(l_p_fasta[i], l_p_fasta[j], pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".needle")
                if path.exists( pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".needle") : # case DNA sequence in fasta
                    dico_out[first_key] [l_p_fasta[j].split("/")[-1][0:-6]] = parseEMBOSS.embossFile(pr_alignement + first_key + "_" + l_p_fasta[j].split("/")[-1][0:-6] + ".needle")[3]
                j = j + 1
            i = i + 1
    
    
#     for i in dico_out.keys () : 
#         print dico_out[i]
    
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
		# print "999999",d_PDB
                d_PDB["best"]["PDB"] = d_PDB["p_pdb_chain"][i]
                d_PDB["best"]["fasta"] = d_PDB["p_fasta_chain"][i]
                return
            else : 
                i = i + 1
        
    
    
def filterGlobalDataset (d_dataset, p_dir_align) : 
    
    # retrieve list of sequences
    l_all_fasta = []
    
    for PDB in d_dataset.keys () : 
        d_dataset[PDB]["keep"] = PDB
        if d_dataset[PDB]["conserve"] == 1 : 
            if not "fasta" in d_dataset[PDB]["best"] : 
                d_dataset[PDB]["conserve"] = 0
                d_dataset[PDB]["keep"] = PDB
                continue
            else : 
                l_all_fasta.append (d_dataset[PDB]["best"]["fasta"])
             
    # alignmenet
    result_align = calculIdenticNeedleCross(l_all_fasta, p_dir_align)
    
    # remove same protein
    for primary_key in result_align.keys () : 
            for secondary_key in result_align[primary_key] .keys (): 
                if result_align[primary_key][secondary_key] == "100.0%" : 
                    if d_dataset[primary_key[0:4]]["conserve"] == 0 or d_dataset[secondary_key[0:4]]["conserve"] == 0 : 
                        continue
                    else : 
                        if d_dataset[primary_key[0:4]]["RX"] < d_dataset[secondary_key[0:4]]["RX"] : 
                            d_dataset[secondary_key[0:4]]["conserve"] = 0
                            d_dataset[secondary_key[0:4]]["keep"] = primary_key[0:4]
                            d_dataset[primary_key[0:4]]["keep"] = primary_key[0:4]
                        else : 
                            d_dataset[primary_key[0:4]]["conserve"] = 0
                            d_dataset[secondary_key[0:4]]["keep"] = secondary_key[0:4]
                            d_dataset[primary_key[0:4]]["keep"] = secondary_key[0:4]
    
    filout_bilan = open (p_dir_align + "global_align_selected.txt", "w")
    
    for pdb_ref in d_dataset.keys () : 
        print pdb_ref
        print d_dataset[pdb_ref].keys ()
        
        try : filout_bilan.write (pdb_ref + "\t" + str(d_dataset[pdb_ref]["best"]["PDB"] ) + "\t" + str(d_dataset[pdb_ref]["keep"] ) + "\t" + str(d_dataset[pdb_ref]["conserve"]) + "\n")
        except : filout_bilan.write( pdb_ref + "\t" + str("------") + "\t" + str(d_dataset[pdb_ref]["keep"]) + "\t" + str(d_dataset[pdb_ref]["conserve"]) + "\n")
            
    filout_bilan.close ()
    
    
                        


def filterBlastResult (d_dataset, p_dir_dataset, substruct, thresold_RX = 2.5, thresold_blast = 1e-4, l_ligand_out = [], debug = 1) : 
    
    """
    Filter resolution PDB
    Filter evalue
    """
    
    for pdb_ref in d_dataset.keys () : 
        if d_dataset[pdb_ref]["conserve"] == 0 : continue
        for pdb_blast_chain in d_dataset[pdb_ref]["align"].keys () : 
            # filter e.value
            if debug == 1 : print d_dataset[pdb_ref]["align"][pdb_blast_chain], pdb_ref, thresold_blast
            # remove thresold and reference cleanner, remove if alignement on same protein
            if d_dataset[pdb_ref]["align"][pdb_blast_chain] <= thresold_blast and pdb_blast_chain [0:4].upper() != pdb_ref: 
                if debug == 1 :print "CONTROL", pdb_blast_chain, d_dataset[pdb_ref]["align"][pdb_blast_chain]
                # dowload PDB files
                pdb_blast = pdb_blast_chain[0:4]
                p_pdb_blast = downloadFile.importPDB(pdb_blast, p_dir_dataset + pdb_ref + "/", dir_by_PDB=0, dbPDB = "/home/borrel/PDB/")
                if p_pdb_blast == 0 : 
                    continue
                
#                 print p_pdb_blast, pdb_blast, "************"
                l_queries_by_chain = separeByChain(p_pdb_blast) # divise chain in reference folder
                try : RX = parsePDB.resolution(p_pdb_blast)
                except : RX = 100.0
                l_ligand = parsePDB.retrieveListLigand(p_pdb_blast)
                if debug == 1 : print p_pdb_blast, l_ligand, RX, "control blast\n" + " ".join(l_queries_by_chain) + "\n" 
                
                # remove apo forms and remove not substiuant + list ligand do not considered
                if l_ligand == [] or substruct in l_ligand: 
                    remove (p_pdb_blast)
                    for queries_by_chain in l_queries_by_chain : 
                        remove (queries_by_chain)
                    continue
                
                for ligand_out in l_ligand_out : 
                    if ligand_out in l_ligand: 
                        try : remove (p_pdb_blast)
                        except : 
                            if debug: print "CONTROL File", p_pdb_blast
                        for queries_by_chain in l_queries_by_chain : 
                            try : remove (queries_by_chain)
                            except : 
                                if debug: print "CONTROL File", " ".join(l_queries_by_chain)
                    continue
                
                # case RMN structure
                try : RX = float(RX)
                except: continue
                
                if debug == 1 :print "----", RX, thresold_RX, "----"
                if float(RX) <= thresold_RX : 
                    if not "blast" in d_dataset[pdb_ref].keys () : 
                        d_dataset[pdb_ref]["blast"] = [pdb_blast_chain.upper ()]
                    else : 
                        # control if other chain in list
                        flag_in = 0
                        for PDB_in in d_dataset[pdb_ref]["blast"] : 
                            if PDB_in.split ("_")[0] == pdb_blast_chain.upper ().split ("_")[0] :
                                flag_in = 1
                        if flag_in == 0 : 
                            d_dataset[pdb_ref]["blast"].append (pdb_blast_chain.upper ())
                            
                    if debug == 1 :print d_dataset[pdb_ref]["blast"], pdb_ref
        
                    
                    
            
def cleanFolderDataset (d_dataset, p_dir_dataset, debug = 1) : 
    
    for PDB_ref in  d_dataset.keys () : 
        if debug : 
            print "*******"
            print  p_dir_dataset + PDB_ref 
            try : print d_dataset[PDB_ref]["best"]["PDB"], d_dataset[PDB_ref]["best"]["fasta"]
            except : pass
            try: print d_dataset[PDB_ref]["blast"]
            except: pass
            print "*********"
        
        if d_dataset[PDB_ref]["conserve"] == 0 or not "blast" in d_dataset[PDB_ref].keys (): 
            rmtree(p_dir_dataset + PDB_ref + "/")
            print "RM1", d_dataset[PDB_ref]["conserve"]
            continue
        else : 
            pr_ref = p_dir_dataset + PDB_ref + "/"
            l_file_dir_ref = listdir(pr_ref)
            for filein_rep in l_file_dir_ref : 
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
            print "BEST", d_dataset[k]["best"]
            if "blast" in d_dataset[k]:
                print "****"
                print "BLAST", d_dataset[k]["blast"]
                print "****"
    print "--------------------------------------"    
    
    
   
