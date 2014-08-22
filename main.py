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

from os import listdir, path, remove, rename
from re import search
from copy import deepcopy, copy




# dataset constuction



# step 1
# first step -> preparation
# - extract the ligand


def datasetPreparation (ligand_ID, clean = 1):
    
    
    p_dir_dataset = pathManage.dataset(ligand_ID)
    
    l_folder = listdir(p_dir_dataset)
    indent = 0
    
    for ref_folder in l_folder  :
        # file include in dataset folder
        if len (ref_folder) != 4 : 
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
            if l_ligand == []  : 
                continue
            else : 
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
        if p_lig_ref == 0 : 
            
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
    
    
    for ref_folder in l_folder  :
        if len (ref_folder) != 4 : 
            continue
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        p_pdb_ref = pathManage.findPDBRef(p_dir_dataset + ref_folder + "/")
        
        
        for pdbfile in l_pdbfile : 
            # try if PDB not ligand
            if len(pdbfile.split ("_")[0]) != 4 or not search (".pdb", pdbfile): 
                continue
            # same alignment
            elif p_dir_dataset + ref_folder + "/" + pdbfile == p_pdb_ref : 
                continue
            else : 
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
    
    # append thresold in name file ???
    
    # ouput
    p_dir_dataset = pathManage.dataset(name_lig)
    p_dir_result = pathManage.result(name_lig )
    l_folder_ref = listdir(p_dir_dataset)
    
    # log control
    p_log = open(p_dir_result + "log_superimposed.txt", "w")
    
    
    d_smile = {}
    
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
        print p_lig_ref
        try : 
            lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref, "HETATM")
            print len (lig_ref_parsed)
        except : 
            p_log.write ("[ERROR ligand ref] -> " + p_lig_ref + "\n")
            continue
        
        # outup by reference
        p_dir_result_ref = pathManage.result(name_lig + "/" + ref_folder)
        d_filout_superimposed = {}
        d_filout_superimposed["global"] = open (p_dir_result_ref + "all_ligand_aligned.pdb", "w")
        d_filout_superimposed["sheap"] = open (p_dir_result_ref + "all_ligand_aligned_" + str (thresold_shaep)  + ".pdb", "w")
        
        
        
        # write lig ref -> connect matrix corrrect
        writePDBfile.coordinateSection(d_filout_superimposed["global"], lig_ref_parsed, "HETATM", connect_matrix = 1)
        
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
                    
                    lig_parsed = parsePDB.loadCoordSectionPDB(p_lig, "HETATM")

                    # find matrix of rotation
                    p_matrix = pathManage.findMatrix(p_lig_ref, p_lig, name_lig)
                    # control file matrix exist
                    if not path.exists(p_matrix) : 
                        p_log.write ("[ERROR] -> Matrix transloc " + p_lig_ref + " " + p_lig + " " + name_lig + "\n")
                        continue
                    
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
                            continue
                        else : 
                            # write PDB file, convert smile
                            p_substituate_pdb = p_dir_result_ref + "substituent_" + pdbfile.split ("_")[0] + "_" + pdbfile.split ("_")[1] + "_" + struct_type + ".pdb"
                            writePDBfile.coordinateSection(p_substituate_pdb, l_atom_substituate, recorder="HETATM", header=0, connect_matrix = 1)
    
                            # Step2 -> convert to smile -> review to smart if a find a good software
                            smile_find = runOtherSoft.babelConvertPDBtoSMILE(p_substituate_pdb)
                    
                            # step3 apply sheap---> try without mol2 --- convertion ---
#                             p_substruct_ref_mol2 = runOtherSoft.babelPDBtoMOL2 (p_substruct_ref)
#                             p_subs_query_mol2 = runOtherSoft.babelPDBtoMOL2 (p_substituate_pdb)
                    
                            p_sheap = runOtherSoft.runShaep (p_substruct_ref, p_substituate_pdb, p_substituate_pdb[0:-4] + ".hit", clean = 1)
                            val_sheap = parseShaep.parseOutputShaep (p_sheap)
                            if val_sheap == {} : 
                                p_log.write ("[ERROR] -> ShaEP " + p_substituate_pdb + " " + p_substruct_ref + "\n")
                                continue
                            # control thresold sheap
                            if not struct_type in d_filout_sheap.keys () : 
                                d_filout_sheap[struct_type] = {}
                                d_filout_sheap[struct_type] = open (p_dir_result + "shaep_global_" + struct_type + ".txt", "w")
                                d_filout_sheap[struct_type].write ("name\tbest_similarity\tshape_similarity\tESP_similarity\n")
                                d_filout_sheap["list"].append (p_dir_result + "shaep_global_" + struct_type + ".txt")
                            
                            # write value in ShaEP control
                            d_filout_sheap[struct_type].write (ref_folder + "_" +  str(pdbfile.split ("_")[1]) + "_" + struct_type + "_" + str (pdbfile.split ("_")[0]) + "\t" + str(val_sheap["best_similarity"]) + "\t" + str(val_sheap["shape_similarity"]) + "\t" + str(val_sheap["ESP_similarity"]) + "\n")
                            d_filout_sheap["global"].write (ref_folder + "_" +  str(pdbfile.split ("_")[1]) + "_" + struct_type + "_" + str (pdbfile.split ("_")[0]) + "\t" + str(val_sheap["best_similarity"]) + "\t" + str(val_sheap["shape_similarity"]) + "\t" + str(val_sheap["ESP_similarity"]) + "\n")
                            
                            # rename file substituent with shaEP value
                            rename(p_substituate_pdb, p_substituate_pdb[:-4] + "_" + str (val_sheap["best_similarity"]) + ".pdb")
                            
                            if val_sheap["best_similarity"] >= thresold_shaep  : 
                                # write subligand superimposed selected in global files
                                writePDBfile.coordinateSection(d_filout_superimposed["sheap"], lig_parsed, recorder= "HETATM", header = str(p_lig.split ("/")[-1]) + "_" + str (val_sheap["best_similarity"]) ,  connect_matrix = 1)
                                
                                ############
                                # write BS #
                                ############
                                # not only protein superimposed -> also ion and water
                                l_atom_complex = parsePDB.loadCoordSectionPDB(p_complex)
                                superposeStructure.applyMatrixProt(l_atom_complex, p_matrix)
                                p_file_cx = p_dir_result_ref +  "CX_" + p_lig.split ("/")[-1]
                                writePDBfile.coordinateSection(p_file_cx, l_atom_complex, recorder="ATOM", header= p_lig.split ("/")[-1], connect_matrix = 0)
    
                                # search atom in BS
                                l_atom_binding_site = []
                                for atom_substruct in lig_parsed : 
                                    for atom_complex in l_atom_complex : 
                                        if parsePDB.distanceTwoatoms (atom_substruct, atom_complex) <= thresold_BS :
                                            l_atom_binding_site.append (deepcopy(atom_complex)) 
                                            # 3. retrieve complet residue
                                            l_atom_res = parsePDB.getResidues(l_atom_binding_site, l_atom_complex)
                                            # 4. write binding site
                                            p_binding = p_dir_result_ref +  "BS_" + p_lig.split ("/")[-1]
                                            writePDBfile.coordinateSection(p_binding, l_atom_res, "ATOM", p_binding, connect_matrix = 0)
                                
                            # smile code substituate analysis                    
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
    return 1

# step 4 
# search in the close environment if metal is here
# compute distance and angles

def ionIdentification (substruct):
    # in folder
    p_dir_dataset = pathManage.dataset(substruct)
    l_folder_ref = listdir(p_dir_dataset)
    
    # output
    p_filout = pathManage.result(substruct) + "ionsAnalysis.txt"
    filout = open (p_filout, "w")
    filout.write ("PDB\tIon\tD1\tD2\tAngle\tAt1\tAt2\n")
    
    for ref_folder in l_folder_ref  :
        if len (ref_folder) != 4 : 
            continue
        p_complex = pathManage.findPDBRef(p_dir_dataset + ref_folder + "/")

        ionSearch.analyseIons (p_dir_dataset + ref_folder + "/", substruct, filout)
    
    filout.close ()




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

def analysisSameBS (substruct, ID_seq = '1.000'):
    
    pr_result = pathManage.result(substruct + "/sameBS")
    
    d_file_sameBS = {}
    d_file_sameBS["global"] = open (pr_result + "RMSD_BS.txt", "w")
    d_file_sameBS["global"].write ("name_bs\tRMSD_prot\tRMSD_BS_ca\tRMSD_BS_all\tD_max\tl_at_BS\tidentic\n")
    pr_dataset = pathManage.dataset(substruct)
    
    
    l_folder_ref = listdir(pr_dataset)
    
    for ref_folder in l_folder_ref  :
        if len (ref_folder) != 4 : 
            continue
        l_reffile = listdir(pr_dataset + ref_folder + "/")
        
        p_pdb_ref = pathManage.findPDBRef(pr_dataset + ref_folder + "/")
        
        for file_ref in l_reffile : 
#             print file_ref, p_pdb_ref.split ("/")[-1]
            if len(file_ref.split("_")[0]) != 4 or file_ref == p_pdb_ref.split ("/")[-1] or search(".fasta", file_ref): 
#                 print file_ref, p_pdb_ref.split ("/")[-1], "*************"
                continue
            else : 
                p_TMalign =  pathManage.alignmentOutput(substruct) + p_pdb_ref.split ("/")[-1][0:-4] + "__" + file_ref[0:-4] + "/RMSD"
                try : score_align = parseTMalign.parseOutputTMalign(p_TMalign)
                except : continue
#                 print score_align
#                 print p_TMalign
                
                if score_align["IDseq"] >= ID_seq : 
                    
                    l_p_substruct_ref = pathManage.findSubstructRef (pr_dataset + ref_folder + "/", substruct)
                    l_p_query = pathManage.findPDBQueryTransloc (pathManage.result(substruct) + ref_folder + "/")
                    
                    for p_query in l_p_query : 
                        for p_substruct_ref in l_p_substruct_ref : 
                            struct_substitued = p_substruct_ref.split ("_")[-2]
                            
                            if not struct_substitued in d_file_sameBS.keys () : 
                                d_file_sameBS[struct_substitued] = open (pr_result + struct_substitued + "_RMSD_BS.txt", "w")
                                d_file_sameBS[struct_substitued].write ("name_bs\tRMSD_prot\tRMSD_BS_ca\tRMSD_BS_all\tD_max\tl_at_BS\tidentic\n")
                        
                            RMSD_bs = analysis.computeRMSDBS (p_pdb_ref, p_query, p_substruct_ref, pr_result)
                            if RMSD_bs != [] : 
                                d_file_sameBS["global"].write (p_substruct_ref.split("/")[-1][0:-4] +  "_*_" + p_query.split ("/")[-1][0:-4]  + "\t" + str(score_align["RMSD"]) + "\t" + str(RMSD_bs[1]) + "\t" + str(RMSD_bs[0]) + "\t" + str(RMSD_bs[2]) + "\t" + str(RMSD_bs[-2]) + "\t" + str(RMSD_bs[-1]) + "\n")
                                d_file_sameBS[struct_substitued].write (p_substruct_ref.split("/")[-1][0:-4] +  "_*_" + p_query.split ("/")[-1][0:-4] + "\t" + str(score_align["RMSD"]) + "\t" + str(RMSD_bs[1]) + "\t" + str(RMSD_bs[0]) + "\t" + str(RMSD_bs[2]) + "\t" + str(RMSD_bs[-2]) + "\t" + str(RMSD_bs[-1]) + "\n")
                                
    
    tool.closeDicoFile(d_file_sameBS)
    return 1
                
        
        
        




# step 5
# manage results  
# - table result
# - folder tree

def manageResult (l_ligand):
    
    pr_result = pathManage.result("final")
    pr_pi = pathManage.result("final/phosphates")
    pr_ribose = pathManage.result("final/ribose")
    
    
    for name_lig in l_ligand : 
        l_p_smile = pathManage.findListSmileFile(name_lig)
        p_file_famile = pathManage.findFamilyFile (name_lig)
        for p_smile in l_p_smile : 
            if search("ribose", p_smile) and  search("txt", p_smile) and search("smile", p_smile): 
                arrangeResult.globalArrangement(pr_ribose, p_smile, p_file_famile, name_lig) 
            elif search("smile", p_smile) and search(".txt", p_smile) : 
                arrangeResult.globalArrangement(pr_pi, p_smile, p_file_famile, name_lig) 
        
    return 1


 




#################
# RUN MAIN !!!! #
#################

# constante
thresold_RX = 2.7
thresold_BS = 4.5
thresold_blast = 1e-100
thresold_superimposed_ribose = 2.5
thresold_superimposed_pi = 3
thresold_IDseq = 100
thresold_shaep = 0.2

### AMP ###
###########


# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "AMP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("AMP")
# applyTMAlign ("AMP")
# ionIdentification ("AMP")
# retrieveSubstructSuperimposed ("AMP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisSmile ("AMP")
# analysisSameBS ("AMP")


### ADP ###
###########

# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "ADP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("ADP")
# applyTMAlign ("ADP")
# ionIdentification ("ADP")
# retrieveSubstructSuperimposed ("ADP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisSameBS ("ADP")
# analysisSmile ("ADP")
# 
# 
# ### POP ###
# ###########
# # # 
# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "POP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("POP")
# applyTMAlign ("POP")
# ionIdentification ("POP")
# retrieveSubstructSuperimposed ("POP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisSameBS ("POP")
# analysisSmile ("POP")
# 
# 
# ### ATP ###
# ###########
# 
# buildData.builtDatasetGlobal(p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB" , ligand_ID = "ATP", thresold_RX = thresold_RX, thresold_blast = thresold_blast, verbose = 1)
# datasetPreparation ("ATP")
# applyTMAlign ("ATP")
# retrieveSubstructSuperimposed ("ATP", thresold_BS = thresold_BS, thresold_superimposed_ribose = thresold_superimposed_ribose, thresold_superimposed_pi = thresold_superimposed_pi, thresold_shaep = thresold_shaep)
# analysisSameBS ("ATP")
# analysisSmile ("ATP")
# 
# 
# 
# manageResult (["AMP", "POP", "ADP", "ATP"])
# arrangeResult.controlResult (["AMP", "ADP", "ATP", "POP"])

arrangeResult.qualityExtraction (["AMP", "ADP", "ATP", "POP"], p_list_ligand = "/home/borrel/Yue_project/resultLigandInPDB",  thresold_sheap=thresold_shaep)
