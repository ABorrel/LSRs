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

from os import listdir, path
from re import search




# dataset constuction



# step 1
# first step -> preparation
# - extract the ligand


def datasetPreparation (substruct):
    
    p_dir_dataset = pathManage.dataset(substruct)
    
    l_folder = listdir(p_dir_dataset)
    
    
    for ref_folder in l_folder  :
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        
        
        for pdbfile in l_pdbfile : 
            p_file_pdb = p_dir_dataset + ref_folder + "/" + pdbfile
            
            # extract ligand in PDB
            if not search (".pdb", pdbfile ) or len (pdbfile.split ("_")[0])==3: 
                continue
            l_ligand = parsePDB.retrieveListLigand(p_file_pdb)
            if l_ligand == []  : 
                continue
            else : 
                l_atom_pdb_parsed = parsePDB.loadCoordSectionPDB(p_file_pdb)
                for name_ligand in l_ligand : 
                    l_lig_parsed = parsePDB.retrieveLigand(l_atom_pdb_parsed, name_ligand)
                    p_filout_ligand = p_dir_dataset + ref_folder + "/" + name_ligand + "_" + path.split(p_file_pdb)[1]
                    writePDBfile.coordinateSection(p_filout_ligand , l_lig_parsed[0], "HETATM", name_ligand + "_" + p_file_pdb , connect_matrix = 1)
                    
        
        # substruct write for shaep
        print p_dir_dataset + ref_folder + "/"
        p_lig_ref = pathManage.findligandRef(p_dir_dataset + ref_folder + "/", substruct)
        print p_lig_ref
        lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref)
        d_l_atom_substruct = substructTools.retrieveSubstruct(lig_ref_parsed, substruct)
        # case with AMP without phosphate
        if d_l_atom_substruct == {}: 
            continue
        # write substruct
        for subs in d_l_atom_substruct.keys (): 
            p_filout_substruct = p_dir_dataset + ref_folder + "/subref_" +  subs + "_" + ref_folder + ".pdb"
            writePDBfile.coordinateSection(p_filout_substruct , d_l_atom_substruct [subs], "HETATM", name_ligand + "_" + p_file_pdb , connect_matrix = 1)
        
    return 1


# step 2
# second step -> surperimposed all protein on reference
# - tmaling
# - generated new folder alignement with the TMalign ouput


def applyTMAlign (substruct):
    
    p_dir_dataset = pathManage.dataset(substruct)
    l_folder = listdir(p_dir_dataset)
    
    
    for ref_folder in l_folder  :
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        p_pdb_ref = pathManage.findPDBRef(p_dir_dataset + ref_folder + "/")
        
        
        for pdbfile in l_pdbfile : 
            # try if PDB not ligand
            if len(pdbfile.split ("_")[0]) != 4 or not search (".pdb", pdbfile): 
                continue
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


def retrieveSubstructSuperimposed (name_lig, thresold_shaep  = 0.0, thresold_substruct = 3.5, thresold_binding = 3.5):
    
    # append thresold in name file ???
    
    # ouput
    p_dir_dataset = pathManage.dataset(name_lig)
    l_folder_ref = listdir(p_dir_dataset)
    d_smile = {}
    
    
    
    
    for ref_folder in l_folder_ref  :
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
#         print ref_folder
        
        # reference
        p_lig_ref = pathManage.findligandRef(p_dir_dataset + ref_folder + "/", name_lig)
        try : lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref, "HETATM")
        except : continue
        
        # outup by reference
        p_dir_result = pathManage.result(name_lig + "/" + ref_folder)
        d_filout_lig = {}
        d_filout_lig["global"] = open (p_dir_result + "all_ligand_aligned_" + str (thresold_shaep) + ".pdb", "w")
        
        
        
        # write lig ref
        writePDBfile.coordinateSection(d_filout_lig["global"], lig_ref_parsed, "HETATM", connect_matrix = 1)
        
#         print p_lig_ref
        
        for pdbfile in l_pdbfile : 
#             print pdbfile[4:10], "****"
            if len(pdbfile.split ("_")[0]) == 3  and len(pdbfile.split ("_")[1]) == 4 and pdbfile.split ("_")[1] != ref_folder:
                p_lig = p_dir_dataset + ref_folder + "/" + pdbfile
                if p_lig_ref != p_lig : 
                    lig_parsed = parsePDB.loadCoordSectionPDB(p_lig, "HETATM")

                    # find matrix of rotation
                    p_matrix = pathManage.findMatrix(p_lig_ref, p_lig, name_lig)
                    # control filce exist
                    if not path.exists(p_matrix) : 
                        continue
                    
                    # find the path of complex used
                    p_complex = p_dir_dataset + ref_folder + "/" + p_lig.split ("/")[-1][4:]
                    
                    # ligand rotated -> change the referentiel
                    superposeStructure.applyMatrixLigand(lig_parsed, p_matrix)
                    
                    
                    
                    
                    # retrieve lig substructure
                    smile_find, p_file_substruct = neighborSearch.searchNeighborAtom(lig_ref_parsed, lig_parsed, p_complex, name_lig, p_lig, p_matrix, p_dir_result, thresold_substruct = 3.5, thresold_binding = 3.5)    
                    if smile_find == 0 and p_file_substruct == 0 : 
                        continue
                    
                    # control sheap // run sheap
                    # run sheap
                    # search ref in dataset
                    l_p_substruct_ref = pathManage.findSubstructRef (pathManage.dataset(name_lig) + ref_folder + "/" , name_lig)
# # # # # # # # # # #                     print l_p_substruct_ref, "?????"
                    
                    
                    for p_substruct_ref in l_p_substruct_ref : 
                        struct_substitued = p_substruct_ref.split ("_")[-2]
                        
                        # need convert in mol2
                        p_substruct_ref = runOtherSoft.babelPDBtoMOL2 (p_substruct_ref)
                        p_subs_file = runOtherSoft.babelPDBtoMOL2 (p_file_substruct)
                        
                        # run shaep
                        p_sheap = runOtherSoft.runShaep (p_substruct_ref, p_subs_file, p_subs_file[0:-5] + ".hit", clean = 1)
                        val_sheap = parseShaep.valueShapeSimilarity (p_sheap)
                        if val_sheap < thresold_shaep  : 
                            continue
                        
                        # Count the smile found + remove not sheap
                        writePDBfile.coordinateSection(d_filout_lig["global"], lig_parsed, "HETATM", p_lig, connect_matrix = 1)
                        
                        if not struct_substitued in d_smile.keys (): 
                            
                            # control superimpostion -> write global alignement
                            d_filout_lig[struct_substitued] = open (p_dir_result + "all_ligand_aligned_" + str (thresold_shaep) + "_" + struct_substitued + ".pdb","w")
                            writePDBfile.coordinateSection(d_filout_lig[struct_substitued], lig_parsed, "HETATM", connect_matrix = 1)
                            
                            # append smile structure
                            
                            d_smile[struct_substitued] = {}
                            d_smile[struct_substitued][smile_find] = {}
                            d_smile[struct_substitued][smile_find]["count"] = 1
                            d_smile[struct_substitued][smile_find]["PDB"] = [pdbfile.split ("_")[1]]
                            d_smile[struct_substitued][smile_find]["ligand"] = [pdbfile.split ("_")[0]]
                            d_smile[struct_substitued][smile_find]["ref"] = [ref_folder]
                        else : 
                            if not struct_substitued in d_filout_lig.keys () : 
                                d_filout_lig[struct_substitued] = open (p_dir_result + "all_ligand_aligned_" + str (thresold_shaep) + "_" + struct_substitued + ".pdb","w")
                            writePDBfile.coordinateSection(d_filout_lig[struct_substitued], lig_parsed, "HETATM", connect_matrix = 1)
                            if not smile_find in d_smile[struct_substitued].keys () : 
                                d_smile[struct_substitued][smile_find] = {}
                                d_smile[struct_substitued][smile_find]["count"] = 1
                                d_smile[struct_substitued][smile_find]["PDB"] = [pdbfile.split ("_")[1]]
                                d_smile[struct_substitued][smile_find]["ligand"] = [pdbfile.split ("_")[0]] 
                                d_smile[struct_substitued][smile_find]["ref"] = [ref_folder]
                            else : 
                                d_smile[struct_substitued][smile_find]["count"] = d_smile[struct_substitued][smile_find]["count"] + 1
                                d_smile[struct_substitued][smile_find]["PDB"].append (pdbfile.split ("_")[1])
                                d_smile[struct_substitued][smile_find]["ligand"].append (pdbfile.split ("_")[0])
                                d_smile[struct_substitued][smile_find]["ref"].append (ref_folder)
                            
        tool.closeDicoFile (d_filout_lig)
            
    # write list of smile
    for substruct in d_smile.keys () : 
        
        p_list_smile = pathManage.result(name_lig) + "list_" + substruct + "_" + str (thresold_shaep) + "_smile.txt"
        p_list_smile_P = pathManage.result(name_lig) + "list_" + substruct + "_" + str (thresold_shaep) + "_smileP.txt"
        
        filout_smile = open (p_list_smile, "w")
        filout_smile_P = open (p_list_smile_P, "w")
        for smile_code in d_smile[substruct].keys () : 
            l_lig = d_smile[substruct][smile_code]["ligand"]
            l_PDB = d_smile[substruct][smile_code]["PDB"]
            l_ref = d_smile[substruct][smile_code]["ref"]
            print smile_code
            if search ("P", str(smile_code)) : 
                filout_smile_P.write (str (smile_code) + "\t" + str (d_smile[substruct][smile_code]["count"]) + "\t" + " ".join (l_PDB) + "\t" + " ".join (l_ref) + "\t" + " ".join(l_lig) + "\n")
            else : 
                filout_smile.write (str (smile_code) + "\t" + str (d_smile[substruct][smile_code]["count"]) + "\t" + " ".join (l_PDB) + "\t" + " ".join (l_ref) + "\t" + " ".join(l_lig) + "\n")
        
        filout_smile.close ()
        filout_smile_P.close ()
    
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

### AMP ###
###########


buildData.builtDatasetGlobal("/home/borrel/Yue_project/resultLigandInPDB", "AMP")
datasetPreparation ("AMP")
applyTMAlign ("AMP")
ionIdentification ("AMP")
retrieveSubstructSuperimposed ("AMP", thresold_shaep = 0.4)
analysisShaep ("AMP")
analysisSmile ("AMP")
analysisSameBS ("AMP")


### ADP ###
###########

buildData.builtDatasetGlobal("/home/borrel/Yue_project/resultLigandInPDB", "ADP")
datasetPreparation ("ADP")
applyTMAlign ("ADP")
ionIdentification ("ADP")
retrieveSubstructSuperimposed ("ADP", thresold_shaep = 0.4)
analysisShaep ("ADP")
analysisSameBS ("ADP")
analysisSmile ("ADP")


### POP ###
###########
# # 
buildData.builtDatasetGlobal("/home/borrel/Yue_project/resultLigandInPDB", "POP")
datasetPreparation ("POP")
applyTMAlign ("POP")
ionIdentification ("POP")
retrieveSubstructSuperimposed ("POP", thresold_shaep = 0.4)
analysisShaep ("POP")
analysisSameBS ("POP")
analysisSmile ("POP")


### ATP ###
###########

buildData.builtDatasetGlobal("/home/borrel/Yue_project/resultLigandInPDB", "ATP")
datasetPreparation ("ATP")
applyTMAlign ("ATP")
# ionIdentification ("ATP") ----> to do
retrieveSubstructSuperimposed ("ATP", thresold_shaep = 0.4)
analysisShaep ("ATP")
analysisSameBS ("ATP")
analysisSmile ("ATP")



manageResult (["POP", "ATP", "AMP", "ADP"])


