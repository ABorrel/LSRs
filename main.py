import pathManage
import parsePDB
import runOtherSoft
import writePDBfile
import superposeStructure
import neighborSearch
import ionSearch
import substructTools
import analysis
import parseTMalign
import buildData
import parseShaep
import tool
import arrangeResult
import managePDB
import refClassification
import ligandSimilarity
import classifResults

from os import listdir, path, remove, rename, system
from re import search
from copy import deepcopy, copy




# dataset constuction


# PDB extraction
def downloadPDB (pr_in):

    pathManage.generatePath(pr_in)
    managePDB.formatPDBDatabase(pr_in)



# step 1
# first step -> preparation
# - extract the ligand


def datasetPreparation (ligand_ID, clean = 1):

    p_dir_dataset = pathManage.dataset(ligand_ID)
    l_folder = listdir(p_dir_dataset)
    indent = 0

    for ref_folder in l_folder  :
        # file include in dataset folder
        if len (ref_folder) != 4:
            continue
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        indent = indent + 1
        print ref_folder, indent

        # clean repertory -> only PDB ref and PDB
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        if clean == 1 : 
            for pdbfile in l_pdbfile : 
                p_file_pdb = p_dir_dataset + ref_folder + "/" + pdbfile
                if not search (".pdb", pdbfile ) or search ("subref", pdbfile) or len (pdbfile.split("_")[0]) == 3: 
                    remove (p_file_pdb)

        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        for pdbfile in l_pdbfile : 
            p_file_pdb = p_dir_dataset + ref_folder + "/" + pdbfile
            # extract ligand in PDB
            l_ligand = parsePDB.retrieveListLigand(p_file_pdb)
#             print l_ligand
            if l_ligand == []:
                continue
            else:
                l_atom_pdb_parsed = parsePDB.loadCoordSectionPDB(p_file_pdb)
                for name_ligand in l_ligand : 
                    l_lig_parsed = parsePDB.retrieveLigand(l_atom_pdb_parsed, name_ligand)
                    if l_lig_parsed == [] : 
                        continue
                    p_filout_ligand = p_dir_dataset + ref_folder + "/" + name_ligand + "_" + path.split(p_file_pdb)[1]
                    writePDBfile.coordinateSection(p_filout_ligand , l_lig_parsed[0], "HETATM", header=0 , connect_matrix = 1)

        # ligand_ID write for shaep
#         print p_dir_dataset + ref_folder + "/"
        p_lig_ref = pathManage.findligandRef(p_dir_dataset + ref_folder + "/", ligand_ID)
        if p_lig_ref == 0:

            continue
#         print p_lig_ref
        lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref)
        d_l_atom_substruct = substructTools.retrieveSubstruct(lig_ref_parsed, ligand_ID)
        # case with AMP without phosphate
        if d_l_atom_substruct == {}:
            continue
        # write ligand_ID
        for subs in d_l_atom_substruct.keys ():
            p_filout_substruct = p_dir_dataset + ref_folder + "/subref_" +  subs + "_" + ref_folder + ".pdb"
            writePDBfile.coordinateSection(p_filout_substruct , d_l_atom_substruct [subs], "HETATM", header=0 , connect_matrix = 1)

    return 1


# step 2
# second step -> surperimposed all protein on reference
# - tmaling
# - generated new folder alignement with the TMalign ouput


def applyTMAlign (substruct):

    p_dir_dataset = pathManage.dataset(substruct)
    l_folder = listdir(p_dir_dataset)

    for ref_folder in l_folder:
        if len (ref_folder) != 4:
            continue
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        p_pdb_ref = pathManage.findPDBRef(p_dir_dataset + ref_folder + "/")


        for pdbfile in l_pdbfile:
            # try if PDB not ligand
            if len(pdbfile.split ("_")[0]) != 4 or not search (".pdb", pdbfile):
                continue
            # same alignment
            elif p_dir_dataset + ref_folder + "/" + pdbfile == p_pdb_ref:
                continue
            else:
                p_file_pdb = p_dir_dataset + ref_folder + "/" + pdbfile
                p_dir_align = pathManage.alignmentOutput(substruct + "/" + p_pdb_ref.split ("/")[-1][:-4] + "__" + p_file_pdb.split ("/")[-1][:-4])

                # superimpose
                runOtherSoft.runTMalign(p_file_pdb, p_pdb_ref, p_dir_align)
    return 1


# step 3
# Third step -> generate the SMART code
# - apply rotated matrix on the ligand - done
# - retrieve substruct at 4a
# - convert in smart
# - generated pdb files for each reference the ligand superimposition - done
# - generated a list of SMART by global phosphate


def retrieveSubstructSuperimposed (name_lig, thresold_BS = 4.5, thresold_superimposed_ribose = 2.5, thresold_superimposed_pi = 3, thresold_shaep = 0.4):

    # ouput
    p_dir_dataset = pathManage.dataset(name_lig)
    p_dir_result = pathManage.result(name_lig )
    l_folder_ref = listdir(p_dir_dataset)

    # log control
    p_log = open(p_dir_result + "log_superimposed.txt", "w")

    # control extraction
    d_control = {}
    d_control["pr ref"] = 0
    d_control["lig query"] = 0
    d_control["subref"] = {}
    d_control["subref empty"] = {}
    d_control["out sheap"] = {}
    filout_control = open (p_dir_result + "quality_extraction.txt", "w")

    # stock smile code
    d_smile = {}

    # sheap control
    d_filout_sheap = {}
    d_filout_sheap ["list"] = [p_dir_result + "shaep_global.txt"]
    d_filout_sheap["global"] = open (p_dir_result + "shaep_global.txt", "w") 
    d_filout_sheap["global"].write ("name\tbest_similarity\tshape_similarity\tESP_similarity\n")

    for ref_folder in l_folder_ref :
        # control folder reference name
        if len (ref_folder) != 4 : 
            p_log.write ("[ERROR folder] -> " + ref_folder + "\n")
            continue

        # reference
        p_lig_ref = pathManage.findligandRef(p_dir_dataset + ref_folder + "/", name_lig)
        try:
            lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref, "HETATM")
#             print len (lig_ref_parsed)
        except:
            p_log.write ("[ERROR ligand ref] -> " + p_lig_ref + "\n")
            continue

        #control
        d_control["pr ref"] = d_control["pr ref"] + 1

        # output by reference
        p_dir_result_ref = pathManage.result(name_lig + "/" + ref_folder)
        d_filout_superimposed = {}
        d_filout_superimposed["global"] = open (p_dir_result_ref + "all_ligand_aligned.pdb", "w")
        d_filout_superimposed["sheap"] = open (p_dir_result_ref + "all_ligand_aligned_" + str (thresold_shaep)  + ".pdb", "w")
        
        
        
        # write lig ref -> connect matrix corrrect in all reference and all sheap
        writePDBfile.coordinateSection(d_filout_superimposed["global"], lig_ref_parsed, "HETATM", connect_matrix = 1)
        writePDBfile.coordinateSection(d_filout_superimposed["sheap"], lig_ref_parsed, "HETATM", connect_matrix = 1)
        
        # inspect folder dataset
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        for pdbfile in l_pdbfile : 
            # no ligand file
            if len (pdbfile.split ("_")) == 1 : 
                continue
            pdbfile = pdbfile[:-4] # remove extention
            
            if len(pdbfile.split ("_")[0]) == 3  and len(pdbfile.split ("_")[1]) == 4 and pdbfile.split ("_")[1] != ref_folder:
                p_lig = p_dir_dataset + ref_folder + "/" + pdbfile  + ".pdb"
                if p_lig_ref != p_lig : 
                    # pass case where ligand replace same ligand -> does not need run
                    if pdbfile.split ("_")[0] == name_lig : 
                        p_log.write ("[REMOVE] -> same ligand substituate")
                        continue
                    
                    # parsed ligand query
                    lig_parsed = parsePDB.loadCoordSectionPDB(p_lig, "HETATM")

                    # find matrix of rotation
                    p_matrix = pathManage.findMatrix(p_lig_ref, p_lig, name_lig)
                    # control file matrix exist
                    if not path.exists(p_matrix) : 
                        p_log.write ("[ERROR] -> Matrix transloc " + p_lig_ref + " " + p_lig + " " + name_lig + "\n")
                        continue
                    
                    # control
                    d_control["lig query"] = d_control["lig query"] + 1
                    
                    # find the path of complex used
                    p_complex = p_dir_dataset + ref_folder + "/" + p_lig.split ("/")[-1][4:]
                    
                    # ligand rotated -> change the referentiel
                    superposeStructure.applyMatrixLigand(lig_parsed, p_matrix)
                    
                    
                    # use substruct
                    l_p_substruct_ref = pathManage.findSubstructRef (pathManage.dataset(name_lig) + ref_folder + "/" , name_lig)
                    for p_substruct_ref in l_p_substruct_ref : 
                        # ribose or phosphate
                        struct_type = p_substruct_ref.split ("_")[-2]
                        substruct_parsed = parsePDB.loadCoordSectionPDB(p_substruct_ref, "HETATM")
                        
                        l_atom_substituate = neighborSearch.searchNeighborAtom(substruct_parsed, lig_parsed, struct_type, p_log, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi)    
                        # control find 
                        if len (l_atom_substituate) == 0 :  
                            if not struct_type in d_control["subref empty"].keys () : 
                                d_control["subref empty"][struct_type] = 1
                            else : 
                                d_control["subref empty"][struct_type] = d_control["subref empty"][struct_type] + 1
                            continue
                        
                        else : 
                            if not struct_type in d_control["subref"].keys () : 
                                d_control["subref"][struct_type] = 1
                            else : 
                                d_control["subref"][struct_type] = d_control["subref"][struct_type] + 1
                            
                            # write PDB file, convert smile
                            p_substituate_pdb = p_dir_result_ref + "substituent_" + pdbfile.split ("_")[0] + "_" + pdbfile.split ("_")[1] + "_" + struct_type + ".pdb"
                            writePDBfile.coordinateSection(p_substituate_pdb, l_atom_substituate, recorder="HETATM", header=0, connect_matrix = 1)
    
                            # sheap reference on part of ligand
                            p_sheap = runOtherSoft.runShaep (p_substruct_ref, p_substituate_pdb, p_substituate_pdb[0:-4] + ".hit", clean = 0)
                            val_sheap = parseShaep.parseOutputShaep (p_sheap)
                            if val_sheap == {} : 
                                p_log.write ("[ERROR] -> ShaEP " + p_substituate_pdb + " " + p_substruct_ref + "\n")
                                
                                if not struct_type in d_control["out sheap"].keys () :
                                    d_control["out sheap"][struct_type] = 1
                                else : 
                                    d_control["out sheap"][struct_type] = d_control["out sheap"][struct_type] + 1
                                continue
                            
                            # control thresold sheap
                            if not struct_type in d_filout_sheap.keys () : 
                                d_filout_sheap[struct_type] = {}
                                d_filout_sheap[struct_type] = open (p_dir_result + "shaep_global_" + struct_type + ".txt", "w")
                                d_filout_sheap[struct_type].write ("name\tbest_similarity\tshape_similarity\tESP_similarity\n")
                                d_filout_sheap["list"].append (p_dir_result + "shaep_global_" + struct_type + ".txt") # to improve with python function
                            
                            # write value in ShaEP control
                            d_filout_sheap[struct_type].write (ref_folder + "_" +  str(pdbfile.split ("_")[1]) + "_" + struct_type + "_" + str (pdbfile.split ("_")[0]) + "\t" + str(val_sheap["best_similarity"]) + "\t" + str(val_sheap["shape_similarity"]) + "\t" + str(val_sheap["ESP_similarity"]) + "\n")
                            d_filout_sheap["global"].write (ref_folder + "_" +  str(pdbfile.split ("_")[1]) + "_" + struct_type + "_" + str (pdbfile.split ("_")[0]) + "\t" + str(val_sheap["best_similarity"]) + "\t" + str(val_sheap["shape_similarity"]) + "\t" + str(val_sheap["ESP_similarity"]) + "\n")
                            
                            # rename file substituent with shaEP value
                            rename(p_substituate_pdb, p_substituate_pdb[:-4] + "_" + str (val_sheap["best_similarity"]) + ".pdb")
                            # rename and change the file name
                            p_substituate_pdb = p_substituate_pdb[:-4] + "_" + str (val_sheap["best_similarity"]) + ".pdb"
                            
                            # write all substruct in global file
                            writePDBfile.coordinateSection(d_filout_superimposed["global"], lig_parsed, recorder= "HETATM", header = str(p_lig.split ("/")[-1]) + "_" + str (val_sheap["best_similarity"]) ,  connect_matrix = 1)
                            
                            # control sheap thresold    
                            if float(val_sheap["best_similarity"]) >= thresold_shaep  : 
                                
                                # write subligand superimposed selected in global files
                                writePDBfile.coordinateSection(d_filout_superimposed["sheap"], lig_parsed, recorder= "HETATM", header = str(p_lig.split ("/")[-1]) + "_" + str (val_sheap["best_similarity"]) ,  connect_matrix = 1)
                                
                                ############
                                # write BS #
                                ############
                                # not only protein superimposed -> also ion and water
                                l_atom_complex = parsePDB.loadCoordSectionPDB(p_complex)
                                superposeStructure.applyMatrixProt(l_atom_complex, p_matrix)
                                p_file_cx = p_dir_result_ref +  "CX_" + p_lig.split ("/")[-1]
                                # write CX
                                writePDBfile.coordinateSection(p_file_cx, l_atom_complex, recorder="ATOM", header= p_lig.split ("/")[-1], connect_matrix = 0)
    
                                # search atom in BS
                                l_atom_binding_site = []
                                for atom_complex in l_atom_complex : 
                                    for atom_substruct in lig_parsed : 
                                        if parsePDB.distanceTwoatoms (atom_substruct, atom_complex) <= thresold_BS :
                                            if not atom_complex in l_atom_binding_site : 
                                                l_atom_binding_site.append (deepcopy(atom_complex))
                                
                                # 3. retrieve complet residue
                                l_atom_BS_res = parsePDB.getResidues(l_atom_binding_site, l_atom_complex)
                                                
                                # 4. write binding site
                                p_binding = p_dir_result_ref +  "BS_" + p_lig.split ("/")[-1]
                                writePDBfile.coordinateSection(p_binding, l_atom_BS_res, "ATOM", p_binding, connect_matrix = 0)
                                
                                # smile code substituate analysis                    
                                # Step smile -> not conversion if shaep not validate 
                                smile_find = runOtherSoft.babelConvertPDBtoSMILE(p_substituate_pdb)
                                if not struct_type in d_smile.keys ()  :
                                    d_smile[struct_type] = {}
                                    d_smile[struct_type][smile_find] = {}
                                    d_smile[struct_type][smile_find]["count"] = 1
                                    d_smile[struct_type][smile_find]["PDB"] = [pdbfile.split ("_")[1]]
                                    d_smile[struct_type][smile_find]["ligand"] = [pdbfile.split ("_")[0]]
                                    d_smile[struct_type][smile_find]["ref"] = [ref_folder]
                                else : 
                                    if not smile_find in d_smile[struct_type].keys () : 
                                        d_smile[struct_type][smile_find] = {}
                                        d_smile[struct_type][smile_find]["count"] = 1
                                        d_smile[struct_type][smile_find]["PDB"] = [pdbfile.split ("_")[1]]
                                        d_smile[struct_type][smile_find]["ligand"] = [pdbfile.split ("_")[0]] 
                                        d_smile[struct_type][smile_find]["ref"] = [ref_folder]
                                    else : 
                                        d_smile[struct_type][smile_find]["count"] = d_smile[struct_type][smile_find]["count"] + 1
                                        d_smile[struct_type][smile_find]["PDB"].append (pdbfile.split ("_")[1])
                                        d_smile[struct_type][smile_find]["ligand"].append (pdbfile.split ("_")[0])
                                        d_smile[struct_type][smile_find]["ref"].append (ref_folder)

                            else : 
                                if not struct_type in d_control["out sheap"].keys () : 
                                    d_control["out sheap"][struct_type] = 1
                                else : 
                                    d_control["out sheap"][struct_type] = d_control["out sheap"][struct_type] + 1

        tool.closeDicoFile (d_filout_superimposed)

    # sheap control    
    tool.closeDicoFile (d_filout_sheap)
    for p_file_sheap in d_filout_sheap["list"] : 
        runOtherSoft.RhistogramMultiple (p_file_sheap)    
        
            
    # write list of smile
    for substruct in d_smile.keys () : 
        p_list_smile = pathManage.result(name_lig) + "list_" + substruct + "_" + str (thresold_shaep) + "_smile.txt"
        filout_smile = open (p_list_smile, "w")
        for smile_code in d_smile[substruct].keys () : 
            l_lig = d_smile[substruct][smile_code]["ligand"]
            l_PDB = d_smile[substruct][smile_code]["PDB"]
            l_ref = d_smile[substruct][smile_code]["ref"]
            filout_smile.write (str (smile_code) + "\t" + str (d_smile[substruct][smile_code]["count"]) + "\t" + " ".join (l_PDB) + "\t" + " ".join (l_ref) + "\t" + " ".join(l_lig) + "\n")
        filout_smile.close ()
    p_log.close ()
    
    # control
    filout_control.write ("NB ref: " + str(d_control["pr ref"]) + "\n")
    filout_control.write ("Ligand query: " + str(d_control["lig query"]) + "\n")
    for k in d_control["subref"].keys () :
        filout_control.write ("LSR " + str (k) + ": " + str(d_control["subref"][k]) + "\n")
    for k in d_control["subref empty"].keys () :
        filout_control.write ("NB LSR empty " + str (k) + ": " + str(d_control["subref empty"][k]) + "\n")
    for k in d_control["out sheap"].keys () :
        filout_control.write ("LSR out by sheap " + str (k) + ": " + str(d_control["out sheap"][k]) + "\n")
    
    filout_control.write ("**********************\n\n")
    for k in d_control["subref"].keys () :
        filout_control.write ("LSR keep" + str (k) + ": " + str(d_control["subref"][k] - d_control["out sheap"][k]) + "\n")
    
    filout_control.close ()
    
    return 1

# step 4 
# search in the close environment if metal is here
# compute distance and angles

def ionIdentification (name_ligand):
    """
    step 4 
    search in the close environment if metal is here
    compute distance and angles
    """
    
    
    # in folder
    p_dir_dataset = pathManage.dataset(name_ligand)
    p_filout = pathManage.result(name_ligand) + "ionsAnalysis.txt"
    ionSearch.analyseIons (p_dir_dataset, name_ligand, p_filout)




# step 6
# Analysis  
# - smile, filtering
# - shaep on substructure and ligand

def analysisSmile (substruct):

    l_p_smile = pathManage.findListSmileFile(substruct) 
    for p_smile in l_p_smile : 
        analysis.selectSmileCode(p_smile, minimal_length_smile = 4)
    return 1



def analysisShaep (substruct):
    analysis.globalShaepStat(substruct)
    return 1

def analysisBS (name_lig, ID_seq = '0.0', debug = 1):
    
    pr_result = pathManage.result(name_lig)
    pr_out = pathManage.result(name_lig + "/sameBS")
    
    # log files
    p_log_file = pr_out + "log.txt"
    filout_log = open (p_log_file, "w")

    # dictionnar with files
    d_file_BS = {}
    d_file_BS["global"] = open (pr_out + name_lig + "_", "w")
    d_file_BS["global"].write ("name_bs\tRMSD_prot\tRMSD_BS_ca\tRMSD_BS_all\tD_max\tl_at_BS\tidentic\n")
    d_file_BS["summary"] = open (pr_out + "summary.txt", "w")
    pr_dataset = pathManage.dataset(name_lig)
     
     
    l_folder_ref = listdir(pr_result)
    nb_BS = 0
    nb_BS_filtered = 0
    nb_same_BS = 0  
    for PDB_ref in l_folder_ref  :
        if debug : print PDB_ref
        if len (PDB_ref) != 4 : 
            continue
         
        p_pdb_ref = pathManage.findPDBRef(pr_dataset + PDB_ref + "/")
        l_p_query = pathManage.findPDBQueryTransloc (pathManage.result(name_lig) + PDB_ref + "/")
        
        if debug : print l_p_query
        for p_query in l_p_query : 
            
            # read TM Align
            if debug : print p_query.split ("/")[-1][7:-4]
            
            p_TMalign =  pathManage.alignmentOutput(name_lig) + p_pdb_ref.split ("/")[-1][0:-4] + "__" + p_query.split ("/")[-1][7:-4] + "/RMSD"
            try : score_align = parseTMalign.parseOutputTMalign(p_TMalign)
            except : 
                filout_log.write ("ERROR TM align " + p_TMalign + "\n")
                continue
            nb_BS = nb_BS + 1
            
            if score_align["IDseq"] >= ID_seq : 
                nb_BS_filtered = nb_BS_filtered + 1
                
                l_p_substruct_ref = pathManage.findSubstructRef (pr_dataset + PDB_ref + "/", name_lig)
                
                # sub BS
                for p_substruct_ref in l_p_substruct_ref : 
                    struct_substitued = p_substruct_ref.split ("_")[-2]
                    
                    # write header
                    if not struct_substitued in d_file_BS.keys () : 
                        d_file_BS[struct_substitued] = open (pr_out + name_lig + "_" + struct_substitued + "_", "w")
                        d_file_BS[struct_substitued].write ("name_bs\tRMSD_prot\tRMSD_BS_ca\tRMSD_BS_all\tD_max\tl_at_BS\tidentic\n")
                        
                    RMSD_bs = analysis.computeRMSDBS (p_pdb_ref, p_query, p_substruct_ref, pr_out)
                    if RMSD_bs != [] : 
                        d_file_BS[struct_substitued].write (p_substruct_ref.split("/")[-1][0:-4] +  "_*_" + p_query.split ("/")[-1][0:-4] + "\t" + str(score_align["RMSD"]) + "\t" + str(RMSD_bs[1]) + "\t" + str(RMSD_bs[0]) + "\t" + str(RMSD_bs[2]) + "\t" + str(RMSD_bs[-2]) + "\t" + str(RMSD_bs[-1]) + "\n")
                      
  

                p_ligand_ref = pathManage.findligandRef(pr_dataset + PDB_ref + "/", name_lig)
                RMSD_bs_lig = analysis.computeRMSDBS (p_pdb_ref, p_query, p_ligand_ref, pr_out)
                if RMSD_bs_lig != [] : 
                    d_file_BS["global"].write (p_ligand_ref.split("/")[-1][0:-4] +  "_*_" + p_query.split ("/")[-1][0:-4] + "\t" + str(score_align["RMSD"]) + "\t" + str(RMSD_bs_lig[1]) + "\t" + str(RMSD_bs_lig[0]) + "\t" + str(RMSD_bs_lig[2]) + "\t" + str(RMSD_bs_lig[-2]) + "\t" + str(RMSD_bs_lig[-1]) + "\n")
                    if RMSD_bs_lig [-1] == 1 : 
                        nb_same_BS = nb_same_BS + 1


    # write summary
    d_file_BS["summary"].write ("BS global: " + str (nb_BS) + "\n")
    d_file_BS["summary"].write ("BS - IDseq " + str (ID_seq) + "%: " +  str (nb_BS_filtered) + "\n")
    d_file_BS["summary"].write ("BS - same atom number: " + str (nb_same_BS) + "\n")
    
    filout_log.close ()
                    
    
    # close files and run histograms                
    for k_dico in d_file_BS.keys () : 
        p_file = d_file_BS[k_dico].name
        d_file_BS[k_dico].close ()
        if name_lig == "ATP" : 
            runOtherSoft.RhistogramRMSD(p_file, max_RMSD = 5.0)
        elif name_lig == "ADP" : 
            runOtherSoft.RhistogramRMSD(p_file, max_RMSD = 4.0)
        elif name_lig == "AMP" : 
            runOtherSoft.RhistogramRMSD(p_file, max_RMSD = 4.0)
        else : 
            runOtherSoft.RhistogramRMSD(p_file, max_RMSD = 3.5)
             

        
    return 1
                


# step 5
# manage results  
# - table result
# - folder tree

def manageResult (l_ligand, name_final, l_out = []):
    
    pr_result = pathManage.result("final_" + name_final)
    # remove the folder 
#     pr_pi = pathManage.result("final/phosphates")
#     pr_ribose = pathManage.result("final/ribose")
    
    
    for name_lig in l_ligand : 
        l_p_smile = pathManage.findListSmileFile(name_lig)
        p_file_famile = pathManage.findFamilyFile (name_lig)
        for p_smile in l_p_smile : 
            if search("ribose", p_smile) and  search(".txt", p_smile) and search("smile", p_smile): 
                arrangeResult.globalArrangement(pr_result, p_smile, p_file_famile, name_lig, l_out) 
            elif search("smile", p_smile) and search(".txt", p_smile) : 
                arrangeResult.globalArrangement(pr_result, p_smile, p_file_famile, name_lig, l_out) 
        
    return 1


 


#################
# RUN MAIN !!!! #
#################

#######
# PDB #
#######

# Folder including the PDB
pPDB = "/home/borrel/PDB/" # need to change
# Folder including results
presult = "/home/borrel/Yue_project/" # need to change

# to download the all PDB database
# downloadPDB (pPDB)
managePDB.searchLigands(pPDB, presult)

# constante
thresold_RX = 2.7
thresold_BS = 4.5
thresold_blast = 1e-100
thresold_superimposed_ribose = 2.5
thresold_superimposed_pi = 2.5
thresold_IDseq = 100
thresold_shaep = 0.2
l_ligand_out = ["AMP", "ADP", "ATP", "TTP", "DCP", "DGT", "DTP", "DUP", "ACP", "AD9", "NAD", "AGS", "UDP", "POP", "APC", "CTP", "AOV", "ANP", "GDP", "GTP", "ANP"]


###########
#   AMP   #
###########

# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "AMP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("AMP")
# applyTMAlign ("AMP")
# ionIdentification ("AMP")
# retrieveSubstructSuperimposed ("AMP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisSmile ("AMP")
# analysisBS ("AMP")

###########
#   ADP   #
###########

# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "ADP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("ADP")
# applyTMAlign ("ADP")
# ionIdentification ("ADP")
# retrieveSubstructSuperimposed ("ADP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisBS ("ADP")
# analysisSmile ("ADP")


##########
#  POP   #
##########
#
# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "POP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("POP")
# applyTMAlign ("POP")
# ionIdentification ("POP")
# retrieveSubstructSuperimposed ("POP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisBS ("POP")
# analysisSmile ("POP")


###########
#   ATP   #
###########

# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "ATP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("ATP")
# applyTMAlign ("ATP")
# ionIdentification ("ATP")
# retrieveSubstructSuperimposed ("ATP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisBS ("ATP")
# analysisSmile ("ATP")


#######################################
# classification of reference protein #
#######################################
prdataset = "/home/borrel/Yue_project/dataset/"
# arrangeResult.controlResult (["AMP", "ADP", "POP"])
# refClassification.classifRefProtein (prdataset, ["AMP", "ADP", "ATP", "POP"])
#
#
#######################
# CLASSIFICATION LSRs #
#######################

name_folder_final = "withoutLig"
# manageResult (["AMP", "ADP", "POP", "ATP"], name_folder_final, l_ligand_out)
# arrangeResult.qualityExtraction (["AMP", "ADP", "POP", "ATP"], name_folder_final, p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB", thresold_sheap = thresold_shaep)
#arrangeResult.countingSubstituent(name_folder_final)


###################################################
# AFFINITY AND INTERACTIONS BY PROTEIN REFERENCE  #
###################################################

# folder final
pr_classif = pathManage.result("final_" + name_folder_final) + "Pi_LSR"
ligandSimilarity.analyseLGDProximity(pr_classif)


#########################################
# ANALYSE CLASSIFICATION BASED ON SHEAP #
#########################################

#classifResults.SheapScoreToClass(pr_classif)



######################
# ANALYSE REFERENCE  #
######################

# analyse enantiomer
# arrangeResult.enantiomer(["AMP", "ADP", "ATP"], name_folder_final)
# analyse the superimposition of ligand references
# arrangeResult.superpositionAllRef(["AMP", "ADP", "POP", "ATP"], name_folder_final)
