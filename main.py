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
        l_atom_substruct = substructTools.retrieveSubstruct(lig_ref_parsed, substruct)
        # case with AMP without phosphate
        if l_atom_substruct == [] : 
            continue
        # write substruct
        p_filout_substruct = p_dir_dataset + ref_folder + "/subref_" + ref_folder + ".pdb"
        writePDBfile.coordinateSection(p_filout_substruct , l_atom_substruct, "HETATM", name_ligand + "_" + p_file_pdb , connect_matrix = 1)
        
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


def retrieveSubstructSuperimposed (substruct):
    
    # ouput
    p_dir_dataset = pathManage.dataset(substruct)
    l_folder_ref = listdir(p_dir_dataset)
    d_smile = {}
    p_list_smile = pathManage.result(substruct ) + "list_smile.txt"
    p_list_smile_P = pathManage.result(substruct) + "list_smile_P.txt"
    
    
    
    for ref_folder in l_folder_ref  :
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
#         print ref_folder
        
        # reference
        p_lig_ref = pathManage.findligandRef(p_dir_dataset + ref_folder + "/", substruct)
        lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref, "HETATM")
        
        
        # outup by reference
        p_dir_result = pathManage.result(substruct + "/" + ref_folder)
        p_lig_align = p_dir_result + "all_ligand_aligned.pdb"
        filout_lig = open (p_lig_align, "w")
        
        # write lig ref
        writePDBfile.coordinateSection(filout_lig, lig_ref_parsed, "HETATM", p_lig_align, connect_matrix = 1)
        
#         print p_lig_ref
        
        for pdbfile in l_pdbfile : 
#             print pdbfile[4:10], "****"
            if len(pdbfile.split ("_")[0]) == 3  and len(pdbfile.split ("_")[1]) == 4 and pdbfile.split ("_")[1] != ref_folder:
                p_lig = p_dir_dataset + ref_folder + "/" + pdbfile
                if p_lig_ref != p_lig : 
                    lig_parsed = parsePDB.loadCoordSectionPDB(p_lig, "HETATM")

                    # find matrix of rotation
                    p_matrix = pathManage.findMatrix(p_lig_ref, p_lig, substruct)
                    
                    # find the path of complex used
                    p_complex = p_dir_dataset + ref_folder + "/" + p_lig.split ("/")[-1][4:]
                    
                    # ligand rotated -> change the referentiel
                    superposeStructure.applyMatrixLigand(lig_parsed, p_matrix)
                    
                    # control superimpostion
                    writePDBfile.coordinateSection(filout_lig, lig_parsed, "HETATM", p_lig, connect_matrix = 1)
                    
                    # retrieve lig substructure
                    smile_find = neighborSearch.searchNeighborAtom(lig_ref_parsed, lig_parsed, p_complex, substruct, p_lig, p_matrix, p_dir_result, thresold_substruct = 3.5, thresold_binding = 3.5)    
                    
                    # Count the smile found
                    if not smile_find in d_smile.keys () : 
                        d_smile[smile_find] = {}
                        d_smile[smile_find]["count"] = 1
                        d_smile[smile_find]["PDB"] = [pdbfile.split ("_")[1]]
                        d_smile[smile_find]["ligand"] = [pdbfile.split ("_")[0]]
                    else : 
                        d_smile[smile_find]["count"] = d_smile[smile_find]["count"] + 1
                        d_smile[smile_find]["PDB"].append (pdbfile.split ("_")[1])
                        d_smile[smile_find]["ligand"].append (pdbfile.split ("_")[0])
        filout_lig.close ()
            
    # write list of smile
    filout_smile = open (p_list_smile, "w")
    filout_smile_P = open (p_list_smile_P, "w")
    for smile_code in d_smile.keys () : 
        l_lig = list(set(d_smile[smile_code]["ligand"]))
        l_PDB = list(set(d_smile[smile_code]["PDB"]))
        print smile_code
        if search ("P", str(smile_code)) : 
            filout_smile_P.write (str (smile_code) + "\t" + str (d_smile[smile_code]["count"]) + "\t" + " ".join (l_PDB) + "\t" + " ".join(l_lig) + "\n")
        else : 
            filout_smile.write (str (smile_code) + "\t" + str (d_smile[smile_code]["count"]) + "\t" + " ".join (l_PDB) + "\t" + " ".join(l_lig) + "\n")
    
    filout_smile.close ()
    
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




# step 5
# shaep run for every reference
# - input 1 Phosphate with oxygen
# - input 2 substructure 

def applyShaep (substruct):
    
    p_dir_result = pathManage.result(substruct)
    l_file_folder = listdir(p_dir_result)
    
    
    
    for folder_ref in l_file_folder  :
        PDB_ID_ref = folder_ref
        
        try : l_dir_ref = listdir(p_dir_result + PDB_ID_ref + "/")
        except : continue
        
        for file_ref in l_dir_ref : 
            
            if search("^substruct", file_ref) and search("pdb", file_ref) : 
                p_subs_file = p_dir_result + PDB_ID_ref + "/" + file_ref
#                 print p_subs_file
#                 print PDB_ID_ref, "******"
                
                # reference substructure
                try : p_substruct_ref = pathManage.findSubstructRef (pathManage.dataset(substruct) + PDB_ID_ref + "/" , substruct)
                except : continue
                
#                 print p_substruct_ref
                
                
                # run shaep
                runOtherSoft.runShaep(p_substruct_ref, p_subs_file,p_subs_file[0:-4] + ".hit")
                
    return 1



# step 6
# Analysis  
# - smile, filtering
# - shaep on substructure and ligand

def analysisSmile (substruct):

    l_p_smile = pathManage.findListSmileFile(substruct) 
    for p_smile in l_p_smile : 
        analysis.selectSmileCode(p_smile, minimal_length_smile = 3)
    

    return 1



def analysisShaep (substruct):
    analysis.globalShaepStat(substruct)
    return 1

def analysisSameBS (substruct):
    
    pr_result = pathManage.result(substruct + "/sameBS")
    p_file_result = pr_result + "RMSD_BS.txt"
    filout_res = open (p_file_result, "w")
    filout_res.write ("name_bs\tRMSD_prot\tRMSD_BS_ca\tRMSD_BS_all\tl_BS\n")
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
                score_align = parseTMalign.parseOutputTMalign(p_TMalign)
#                 print score_align
#                 print p_TMalign
                
                if score_align["IDseq"] == '1.000' : 
                    l_p_BS = pathManage.findFileBS(pathManage.result(substruct + "/" + ref_folder), file_ref[0:-4])
                    for p_BS in l_p_BS : 
                        print p_BS
                        RMSD_bs = analysis.computeRMSDBS (p_pdb_ref, p_BS, pr_result)
                        if RMSD_bs != [] : 
                            filout_res.write (p_BS.split ("/")[-1][0:-4] + "_" + p_pdb_ref.split ("/")[-1][0:-4] + "\t" + str(score_align["RMSD"]) + "\t" + str(RMSD_bs[1]) + "\t" + str(RMSD_bs[0]) + "\t" + str(RMSD_bs[-1]) + "\n")
    
    
    filout_res.close ()
    return 1
                
        
        
        




# step 5
# manage results  
# - table result
# - folder tree

def manageResult ():
    
    return 0


 




#################
# RUN MAIN !!!! #
#################

### AMP ###
###########


# buildData.builtDatasetGlobal("/home/borrel/Yue_project/resultLigandInPDB", "AMP")
# datasetPreparation ("AMP")
# applyTMAlign ("AMP")
# retrieveSubstructSuperimposed ("AMP")
# ionIdentification ("AMP")
# applyShaep ("AMP")
# analysisShaep ("AMP")
analysisSameBS ("AMP")
analysisSmile ("AMP")



### ADP ###
###########
 
# datasetPreparation ("ADP")
# applyTMAlign ("ADP")
# retrieveSubstructSuperimposed ("ADP")
# ionIdentification ("ADP")
# applyShaep ("ADP")
# analysisShaep ("ADP")
# analysisSameBS ("ADP")
# analysisSmile ("ADP")



### POP ###
###########

 
# datasetPreparation ("POP")
# applyTMAlign ("POP")
# retrieveSubstructSuperimposed ("POP")
# ionIdentification ("POP")
# applyShaep ("POP")
# analysisShaep ("POP")
# analysisSameBS ("POP")
# analysisSmile ("POP")

