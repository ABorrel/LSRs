import smileAnalysis
import pathManage
import parsePDB
import analysis
import superposeStructure
import writePDBfile
import tool
import runOtherSoft
import substructTools
import superimpose

from os import makedirs, listdir, path
from shutil import copy2
from re import search
import numpy as np
from copy import deepcopy

def globalArrangement (pr_orgin, p_smile, p_family, name_ligand, l_ligand_out):
    
#     print "--------"
#     print pr_orgin
#     print p_smile
#     print p_family
#     print name_ligand
#     print "--------"
    
    
    subst = p_smile.split ("_")[-3]
    
    filin = open (p_smile, "r")
    l_line_smile = filin.readlines ()
    filin.close()
    
    for line_smile in l_line_smile : 
        
        # search substructure
#         print line_smile
        l_PDB_query = line_smile.split ("\t")[-3].split (" ")
#         print l_PDB_query
        l_PDB_ref = line_smile.split ("\t")[-2].split (" ")
        l_ligand = line_smile.strip().split ("\t")[-1].split (" ")
        
        # search replacement
        smile = line_smile.split ("\t")[0]
        
        # search if LSR is small -> thresold < 3
        small_LSR = smileAnalysis.smallLSR (smile) 
        if subst == "ribose" :  
            if small_LSR == 1 : 
                first_folder == "ribose_small"
            else : 
                first_folder = "ribose"
        else : 
            if small_LSR == 1 : 
                first_folder = "Pi_small"
            else : 
                first_folder = "Pi"
        
        
        print smile, l_PDB_query, l_PDB_ref, l_ligand, subst, small_LSR
        replacement, metal = smileAnalysis.searchReplacement (smile, l_PDB_query[0], l_PDB_ref[0], name_ligand)
        
        # case with cycle -> search replacement 2
        if replacement == "cycle" : 
            replacement2, metal = smileAnalysis.searchReplacement (smile, l_PDB_query[0], l_PDB_ref[0], name_ligand, in_cycle = 1)
            replacement = replacement + "/" + replacement2 # new folder

        # case metal
        if replacement == "metal" : 
            print metal, l_PDB_query, l_PDB_ref, name_ligand
        
        len_find = len (l_PDB_ref)
        i = 0
        while i < len_find : 
            
            # exclusion of ligand out
            if l_ligand[i] in l_ligand_out : 
                i = i + 1
                continue
            
            
            family = analysis.findFamily(l_PDB_ref[i], p_family)
            
            # folder reference
            pr_dataset = pathManage.dataset(name_ligand + "/" + l_PDB_ref[i])
            
            PDB_ref = pathManage.findPDBRef(pr_dataset)
            p_ligand_ref = pathManage.findligandRef(pr_dataset , name_ligand)
            l_frag_ref = pathManage.findSubstructRef(pr_dataset, name_ligand)
            for f_ref in l_frag_ref :
                if search (subst, f_ref) : 
                    p_frag_ref = f_ref
                    break
            
            # folder_query
            pr_result = pathManage.result(name_ligand + "/" + l_PDB_ref[i])
            l_protein_tranloc = pathManage.findPDBQueryTransloc(pr_result)
            for p_t in l_protein_tranloc : 
                if search (l_ligand[i], p_t) and search (l_PDB_query[i], p_t) : 
                    p_protein_query = p_t
                    break
                
            if replacement != "metal" : 
                p_lig_query = pathManage.findligandQuery(pr_dataset , l_ligand[i], l_PDB_query[i])
            else : 
                p_lig_query = pathManage.findligandQuery(pr_dataset ,metal, l_PDB_query[i])
            # need apply transloc matrix
            matrix_transloc = pathManage.findMatrix(p_ligand_ref, p_lig_query, name_ligand)
            lig_query_parsed = parsePDB.loadCoordSectionPDB(p_lig_query)
            try : superposeStructure.applyMatrixLigand(lig_query_parsed, matrix_transloc)
            except : 
                i = i + 1
                continue
            
            
            p_lig_substituate = pathManage.findSubstructFind(pr_result, l_ligand[i], l_PDB_query[i], subst)
            l_p_BS = pathManage.findFileBS(pr_result, l_PDB_query[i])
            for BS in l_p_BS : 
                if search (l_ligand[i], BS) : 
                    p_BS = BS
                    break
            
            
#             print pr_final
#             print "***************"
#             print PDB_ref
#             print p_ligand_ref
#             print p_frag_ref
#             print "----"
#             print p_protein_query
#             print p_lig_query
#             print p_lig_substituate
#             print p_BS
#             print "**************"
            
            pr_final = pr_orgin + first_folder + "/" + replacement + "/" + l_PDB_ref[i] + "/" 
            pr_ligand = pr_orgin + first_folder + "/" + replacement + "/" + l_PDB_ref[i] + "/LGD/"
            pr_BS =  pr_orgin + first_folder + "/" + replacement + "/" + l_PDB_ref[i] + "/BS/"
            pr_sust = pr_orgin + first_folder + "/" + replacement + "/" + l_PDB_ref[i] + "/" + str (subst) + "/"
            
            if not path.isdir(pr_final):
                makedirs (pr_final)
            
            if not path.isdir(pr_ligand):
                makedirs (pr_ligand)
            
            if not path.isdir(pr_BS):
                makedirs (pr_BS)
                
            if not path.isdir(pr_sust):
                makedirs (pr_sust)   
            
            # list file
            p_list_smile_queries = pr_sust + "list.smile"
            if not path.exists(p_list_smile_queries) : 
                file_smile_queries = open (p_list_smile_queries, "w")
            else : 
                file_smile_queries = open (p_list_smile_queries, "a")
            file_smile_queries.write (str(smile) + "\n")
            file_smile_queries.close ()
            
            # lig de la query
            writePDBfile.coordinateSection(pr_ligand + "LGD_" + p_lig_query.split ("/")[-1], lig_query_parsed, recorder = "HETATM", header = "LCG_" + p_lig_query.split ("/")[-1], connect_matrix = 1)
            runOtherSoft.babelConvertPDBtoSMILE(pr_ligand + "LGD_" + p_lig_query.split ("/")[-1])
            # lig de reference + smile
            copy2(p_ligand_ref, pr_ligand + "LGD_REF_" + p_ligand_ref.split ("/")[-1])
            runOtherSoft.babelConvertPDBtoSMILE(pr_ligand + "LGD_REF_" + p_ligand_ref.split ("/")[-1])
            # LSR de ref
            copy2(p_frag_ref, pr_sust + "LSR_REF_" + name_ligand + "_" + l_PDB_ref[i] + ".pdb")
            # protein query
            #copy2(p_protein_query, pr_final)
            # LSR query -> p_lig_ref only for the name
            copy2(p_lig_substituate, pr_sust + "LSR_" + p_lig_query.split ("/")[-1])
            # BS query
            copy2(p_BS, pr_BS)   
            
            # BS from reference
            l_atom_BS = parsePDB.computeBS (PDB_ref, p_ligand_ref, thresold = 4.50, option_onlyATOM = 0)
            writePDBfile.coordinateSection(pr_BS + "BS_REF_" + name_ligand + "_" + PDB_ref.split ("/")[-1], l_atom_BS, recorder = "ATOM", header = "BS_REF_" + name_ligand + "_" + PDB_ref, connect_matrix = 0)
            
            i = i + 1
    
    return 1



def controlResult (l_name_ligand):
    
    filout = open(pathManage.result() + "sheap_control.txt", "w")
    
    for name_ligand in l_name_ligand :
        count_sheap = 0
        count_sheap_out = 0 
        count_ribose = 0
        pr_result = pathManage.result(name_ligand)
        
        l_ref = listdir(pr_result)
        for ref_PDB in l_ref : 
            if len(ref_PDB) == 4 : 
                print ref_PDB
                pr_ref = pr_result + ref_PDB
                l_file = listdir(pr_ref)
                for file_ref in l_file : 
                    if search(".hit", file_ref) : 
                        count_sheap = count_sheap + 1
                        
                        if path.getsize(pr_ref +"/" + file_ref ) < 100 : 
                            count_sheap_out = count_sheap_out + 1
                        if search("ribose", file_ref) : 
                            count_ribose = count_ribose + 1
        filout.write (name_ligand + "\n")
        filout.write ("count Shaep:" + str (count_sheap) + "\n")
        filout.write ("count Shaep wrong:" + str (count_sheap_out) + "\n")
        filout.write ("count Shaep ribose:" + str (count_ribose) + "\n")
        filout.write ("******************\n")
        
        

    
def qualityExtraction (l_ligand, name_folder, p_list_ligand, thresold_sheap) : 
    
    pr_result = pathManage.result("final_" + name_folder)
    
    filout = open(pr_result + "quality_extraction.txt", "w")
    
    # number PDB by ligand, without filter
    filout.write ("Number PDB by ligand:\n")
    
    d_dataset =  tool.parseLigandPDBList(p_list_ligand)
    for ligand in l_ligand : 
        filout.write (str (ligand) + ": " + str (len (d_dataset[ligand])) + "\n")
    
    # number references
    filout.write ("\n*************\n\nNumber references by ligands:\n")
    for ligand in l_ligand : 
        pr_result_ligand = pathManage.result(ligand)
        nb_ref = -2
        l_file = listdir(pr_result_ligand)
        for f in l_file : 
            if path.isdir (pr_result_ligand + "/" + f) : 
                nb_ref = nb_ref + 1
        filout.write (ligand + ": " + str (nb_ref) + "\n")
        
    # number of query by ref in means and max and min (after blast)
    filout.write ("\n*************\n\nNumber means queries by references:\n")
    for ligand in l_ligand : 
        d_nb_query = {}
        d_family = {}
        p_filout_family = pathManage.result() + "reference_family_" + ligand + ".txt"
        filout_family = open (p_filout_family, "w")
        pr_result_ligand = pathManage.result(ligand)
        nb_ref = 0
        l_file = listdir(pr_result_ligand)
        for f in l_file : 
            if path.isdir (pr_result_ligand + "/" + f) and len (f) == 4: 
                
                # count by family
                family_ref = analysis.findFamily(f, pathManage.findFamilyFile (ligand))
                if not family_ref in d_family.keys () : 
                    d_family[family_ref] = 0
                d_family[family_ref] = d_family[family_ref] + 1
                
                # count number of references
                nb_ref = nb_ref + 1
                d_nb_query[f] = 0
                l_file_queries = listdir(pr_result_ligand + "/" + f + "/")
                for file_query in l_file_queries : 
                    if search ("CX",file_query) : 
                        d_nb_query[f] = d_nb_query[f] + 1
        filout.write (ligand + ": " + str(np.sum(d_nb_query.values ())) + "\n")
        filout.write (ligand + ": " + str(np.mean(d_nb_query.values ())) + "+/-" + str(np.std (d_nb_query.values ())) + "\n")
        filout.write ("MAX " + str (ligand) + ": " + str (max (d_nb_query.values ())) + " " + str (d_nb_query.keys ()[d_nb_query.values ().index (max (d_nb_query.values ()))]) +"\n")
    
        # family
        filout_family.write ("\t".join(d_family.keys ()) + "\n")
        l_values = [str(x) for x in d_family.values ()]
        filout_family.write ("\t".join(l_values) + "\n")
        filout_family.close ()
        runOtherSoft.piePlot(p_filout_family)
            
    
    # number subref by ligand
    filout.write ("\n*************\n\nNumber of subref considered:\n")
    for ligand in l_ligand :
        d_nb_sub = {}
        d_nb_sub_sheap = {}
        pr_result_ligand = pathManage.result(ligand)
        l_ref = listdir(pr_result_ligand)
        for ref in l_ref : 
            if path.isdir (pr_result_ligand + "/" + ref) and len (ref) == 4: 
                l_file_queries = listdir(pr_result_ligand + "/" + ref + "/")
                for file_query in l_file_queries : 
                    if search ("substituent",file_query) and search (".pdb",file_query): 
                        atom_substituate = file_query.split ("_")[-2]
                        try : value_sheap = float(file_query.split ("_")[-1][:-4])
                        except : continue
                        if not atom_substituate in d_nb_sub.keys () : 
                            d_nb_sub[atom_substituate] = 0
                        d_nb_sub[atom_substituate] = d_nb_sub[atom_substituate] + 1
                        
                        if value_sheap > thresold_sheap : 
                            if not atom_substituate in d_nb_sub_sheap : 
                                d_nb_sub_sheap[atom_substituate] = 0
                            d_nb_sub_sheap[atom_substituate] = d_nb_sub_sheap[atom_substituate] + 1
        filout.write ("\n" + ligand + "\n")
        for atom_substituate in d_nb_sub.keys () : 
            filout.write (atom_substituate + ": " + str (d_nb_sub[atom_substituate]) + "\n")
            try : filout.write (atom_substituate + " ShaEP: " + str (d_nb_sub_sheap[atom_substituate]) + "\n")
            except : filout.write (atom_substituate + " ShaEP: 0\n")
    filout.close()
    
    
    
    
def countingSubstituent (name_final, debug = 1):
    
    pr_final_folder = pathManage.result("final_" + name_final)
    
    d_count = {}
    d_lig = {}
    d_by_ref = {}
    d_count_pr = {}
    l_file_final = listdir(pr_final_folder)
    if debug : print "1", pr_final_folder
    for pr_type_substituent in l_file_final :
        sub_query = pr_type_substituent
        # case where pr type is a file not a folder
        try : l_ref_ID = listdir(pr_final_folder + pr_type_substituent + "/")
        except : continue
        if debug: print "2",pr_final_folder +  pr_type_substituent + "/"
        
        for ref_ID in l_ref_ID : 
            # case where pr_type_substituent is a folder
            try : l_ligand_sub = listdir(pr_final_folder + pr_type_substituent + "/" + ref_ID + "/")
            except : continue
            
            ################
            #  folder LSR  #
            ################
            for ligand_sub in l_ligand_sub :
                if ligand_sub == "BS" or ligand_sub == "LGD" : 
                    continue
                else : 
                    #print l_file_ref
                    if not ligand_sub in d_count.keys () : 
                        d_count[ligand_sub] = {}
                    
                    if not sub_query in d_count[ligand_sub].keys () : 
                        d_count[ligand_sub] [sub_query] = 0
                    
                    d_count[ligand_sub][sub_query] = d_count[ligand_sub][sub_query] + 1
                    
                    ################
                    # complet LSR  #
                    ################
                    l_file_LSR = listdir(pr_final_folder + pr_type_substituent + "/" + ref_ID + "/" + ligand_sub + "/")
                    for file_LSR in l_file_LSR :
                        if search ("LSR", file_LSR):
                            # case LSR reference #
                            ######################
                            if search ("REF_", file_LSR) :
                                lig_ref = file_LSR.split ("_")[2]
                                if not lig_ref in d_by_ref.keys () : 
                                    d_by_ref[lig_ref] = {}
                                break
                            
                    for file_LSR in l_file_LSR : 
                        if not search ("REF_", file_LSR) and search ("LSR", file_LSR) :
                            # count case other LSR that ligand #
                            ####################################      
                            if not ligand_sub in d_by_ref[lig_ref].keys () : 
                                d_by_ref[lig_ref][ligand_sub] = 0
                            d_by_ref[lig_ref][ligand_sub] = d_by_ref[lig_ref][ligand_sub] + 1
            
            
            #################    
            #  folder LGD   #
            #################
            l_file_LGD = listdir(pr_final_folder + pr_type_substituent + "/" + ref_ID + "/LGD/")
            for file_LGD in l_file_LGD : 
                # print file_ref
                if search ("LGD", file_LGD):
                    ligand = file_LGD.split ("_")[1]
                    if ligand == "REF" : 
                        continue
                    if not ligand in d_lig.keys () : 
                        d_lig[ligand] = 1
                    else : 
                        d_lig[ligand] = d_lig[ligand] + 1

            
            ###############
            #  folder BS  #
            ###############
            l_file_BS = listdir(pr_final_folder + pr_type_substituent + "/" + ref_ID + "/BS/")
            for file_BS in l_file_BS : 
                if search ("BS_REF", file_BS):
                    lig_ref = file_BS.split ("_")[2]
                    pr_ref = file_BS.split ("_")[3].split (".")[0]
                    print lig_ref, pr_ref, "*****"
                    if not lig_ref in d_count_pr.keys () : 
                        d_count_pr[lig_ref] = {}
                        d_count_pr[lig_ref]["pr ref"] = []
                        d_count_pr[lig_ref]["pr queries"] = []
                        d_count_pr[lig_ref]["lig queries"] = []
                        d_count_pr[lig_ref]["HSP"] = 0
                        d_count_pr[lig_ref]["kinase"] = 0
                        d_count_pr[lig_ref]["others"] = 0
                                   
                    if not pr_ref in d_count_pr[lig_ref]["pr ref"] : 
                        d_count_pr[lig_ref]["pr ref"].append (pr_ref)
                                
                                
                    try:
                        family = analysis.findFamily (pr_ref, pathManage.dataset (lig_ref) + "family_PDB.txt")
                        d_count_pr[lig_ref][family] = d_count_pr[lig_ref][family] + 1
                    except: pass
                
            for file_BS in l_file_BS : 
                # for not reference BS
                if not search ("BS_REF", file_BS) : 
                    lig_querie = file_BS.split ("_")[1]
                    prot_querie = file_BS.split ("_")[2][0:4]
                    print prot_querie, lig_querie, "*******"
                    # find ligand reference
                    # lig ref define in previous step
                    d_count_pr[lig_ref]["pr queries"].append (prot_querie)
                    d_count_pr[lig_ref]["lig queries"].append (lig_querie)


    # write and plot #
    ##################
    pr_result = pathManage.generatePath(pr_final_folder + "counting/")
    for ligand_sub in d_count.keys () : 
        p_filout = pr_result + ligand_sub
        filout = open (p_filout, "w")
        filout.write ("\t".join(d_count[ligand_sub].keys ()) + "\n")
        l_value = [str(x) for x in d_count[ligand_sub].values ()]
        filout.write ("\t".join(l_value) + "\n")
        filout.close ()
        runOtherSoft.piePlot(p_filout)
    
    filout_lig = open (pr_result + "count_ligand", "w")
    for lig in d_lig.keys () : 
        if d_lig[lig] > 1 : 
            filout_lig.write (str (lig) + "\t" + str (d_lig[lig]) + "\n")
    filout_lig.close ()
    
    filout_LSR_lig = open (pr_result + "CountByLigandRef", "w")
    for lig_ref in d_by_ref.keys () : 
        filout_LSR_lig.write ("====" + str (lig_ref) + "====\n")
        for sub_ref in d_by_ref[lig_ref].keys () : 
            filout_LSR_lig.write (str (sub_ref) + ": " + str (d_by_ref[lig_ref][sub_ref]) + "\n")
    filout_LSR_lig.close ()

    filout_pr_count = open (pr_result + "count_pr", "w")
    for lig in d_count_pr.keys () : 
        filout_pr_count.write ("====" + str (lig) + "====\n")
        filout_pr_count.write ("nb ref pr: " + str (len (d_count_pr[lig]["pr ref"])) + "\n")
        filout_pr_count.write ("nb querie pr: " + str (len (d_count_pr[lig]["pr queries"])) + "\n")
        filout_pr_count.write ("nb ligand queries: " + str (len (d_count_pr[lig]["lig queries"])) + "\n")
        filout_pr_count.write ("Ref kinase: " + str (d_count_pr[lig]["kinase"]) + "\n")
        filout_pr_count.write ("Ref HSP: " + str (d_count_pr[lig]["HSP"]) + "\n")
        filout_pr_count.write ("Ref other: " + str (d_count_pr[lig]["others"]) + "\n")


    filout_pr_count.close ()

    runOtherSoft.barplot(pr_result + "count_ligand")
        
        
        
def enantiomer(l_ligand, name_folder_final) : 
    "to do file output"
    
    pr_final = pathManage.result("final_" + name_folder_final)
    
    pr_enantiomer = pathManage.generatePath(pr_final + "enantiomer/")
    
    d_filout = {}
    for ligand in l_ligand : 
        d_filout[ligand] = {}
        d_filout[ligand]["O3OP"]= open (pr_enantiomer + ligand + "_" + "O3OP.txt" , "w")
        d_filout[ligand]["O4O5"]= open (pr_enantiomer + ligand + "_" + "O4O5.txt" , "w")
        d_filout[ligand]["OPOP"]= open (pr_enantiomer + ligand + "_" + "OPOP.txt" , "w")
        
    l_p_substruct = listdir(pr_final) 
    for pr_substruct in l_p_substruct : 
        print pr_substruct
        # case where pr_substruct is a file not a folder
        try : l_pr_ref = listdir(pr_final + pr_substruct + "/")
        except : continue
        for pr_ref in l_pr_ref : 
            print pr_ref
            # case no folder
            try : l_file = listdir(pr_final + pr_substruct + "/" + pr_ref + "/LGD/")
            except : continue
            for name_file in l_file : 
                if search("REF_A",name_file) and   search(".pdb",name_file): 
                    ligand = name_file.split ("_")[2]
                    l_atom_ligand = parsePDB.loadCoordSectionPDB(pr_final + pr_substruct + "/" + pr_ref + "/LGD/" + name_file, "HETATM")
                    d_minO3OP = 100
                    for atom_ligand in l_atom_ligand : 
                        if atom_ligand["name"] == "O4'" :
                            atom_O4 = atom_ligand
                        elif atom_ligand["name"] == "O5'" :
                            atom_O5 = atom_ligand
                        elif  atom_ligand["name"] == "O3'" :
                            atom_O3 = atom_ligand
                        elif  atom_ligand["name"] == "O1A" :
                            atom_O1A = atom_ligand
                        elif  atom_ligand["name"] == "O2A" :
                            atom_O2A = atom_ligand
                        elif  atom_ligand["name"] == "O1B" :
                            atom_O1B = atom_ligand
                        elif  atom_ligand["name"] == "O2B" :
                            atom_O2B = atom_ligand
                        #elif  atom_ligand["name"] == "O3B" :
                        #    atom_O3B = atom_ligand
                    
                    # d O4 - O5        
                    try : d_O4O5 = parsePDB.distanceTwoatoms(atom_O4, atom_O5)
                    except : continue
                    d_filout[ligand]["O4O5"].write (pr_ref + "_" + pr_substruct  + "\t" + str (d_O4O5) + "\n")

                    # d O3 - OP
                    for atom_ligand in l_atom_ligand : 
                        if ligand == "AMP" : 
                            if atom_ligand["name"] == "O1P" or atom_ligand["name"] == "O2P" or atom_ligand["name"] == "O3P" : 
                                d_tempO3OP = parsePDB.distanceTwoatoms(atom_O3, atom_ligand)
                                if d_tempO3OP < d_minO3OP : 
                                    d_minO3OP = d_tempO3OP
                                    atom_tempO3OP = deepcopy(atom_ligand)
                        else : 
                            if atom_ligand["name"] == "O1A" or atom_ligand["name"] == "O2A" or atom_ligand["name"] == "O3A" : 
                                d_tempO3OP = parsePDB.distanceTwoatoms(atom_O4, atom_ligand)
                                if d_tempO3OP < d_minO3OP : 
                                    d_minO3OP = d_tempO3OP
                                    atom_tempO3OP = deepcopy(atom_ligand)
                    d_filout[ligand]["O3OP"].write (pr_ref + "_" + pr_substruct  +"_" + str(atom_tempO3OP["name"]) + "\t" + str (d_minO3OP) + "\n")
    
                    # d OP OP
                    d_OP = {}
                    if ligand == "ATP" or ligand == "ADP" : 
                        d_OP ["O1AO1B"] = parsePDB.distanceTwoatoms(atom_O1A, atom_O1B)
                        d_OP ["O1AO2B"] = parsePDB.distanceTwoatoms(atom_O1A, atom_O2B)
                        #d_OP ["O1AO3B"] = parsePDB.distanceTwoatoms(atom_O1A, atom_O3B)
                        d_OP ["O2AO1B"] = parsePDB.distanceTwoatoms(atom_O2A, atom_O1B)
                        d_OP ["O2AO2B"] = parsePDB.distanceTwoatoms(atom_O2A, atom_O2B)
                        #d_OP ["O2AO3B"] = parsePDB.distanceTwoatoms(atom_O2A, atom_O3B)
                        
                        d_minOPOP = min (d_OP.values())
                        #print d_minOPOP
                        k_min = [name for name, age in d_OP.items() if age == min (d_OP.values())][0]
                        #print k_min
                        d_filout[ligand]["OPOP"].write (pr_ref + "_" + pr_substruct  +"_" + str(k_min) + "\t" + str (d_minOPOP) + "\n")
                    
                    try :
                        del d_OP 
                        del atom_O1A
                        del atom_O1B
                        del atom_O2A
                        del atom_O2B
                    except : 
                        pass
                    try : 
                        del atom_O3
                        del atom_O4
                        del atom_O5
                    except :
                        pass
            
    # close files
    for lig in l_ligand : 
        for type_dist in d_filout[lig].keys () : 
            p_file = d_filout[lig][type_dist].name
            d_filout[lig][type_dist].close ()
            runOtherSoft.Rhistogram(p_file, type_dist, brk = 20)
    
                                
                                
def superpositionAllRef (l_ligand, name_folder_final):   
    
    pr_final = pathManage.result("final_" + name_folder_final)
    
    pr_align = pathManage.generatePath(pr_final + "refAlignement/")
    
    d_filout_pdb = {}
    d_filout_RMSE = {}
    d_ref = {}
    l_file_RMSE = []
    for ligand in l_ligand : 
        d_filout_pdb[ligand] = open (pr_align + ligand + "_" + "superimposed.pdb" , "w")
        d_filout_RMSE[ligand] = open (pr_align + ligand + "_" + "RMSE.txt" , "w")
        l_file_RMSE.append (pr_align + ligand + "_" + "RMSE.txt") 
    
    l_p_substruct = listdir(pr_final) 
    for pr_substruct in l_p_substruct : 
        try : l_pr_ref = listdir(pr_final + pr_substruct + "/")
        except : continue
        for pr_ref in l_pr_ref : 
            # case where pr_ref is not a folder
            try : l_file = listdir(pr_final + pr_substruct + "/" + pr_ref + "/LGD/")
            except : continue
            for name_file in l_file : 
                if search("REF_A",name_file) and   search(".pdb",name_file): 
                    ligand = name_file.split ("_")[2]
                    l_atom_ligand = parsePDB.loadCoordSectionPDB(pr_final + pr_substruct + "/" + pr_ref + "/LGD/" + name_file, "HETATM", remove_H=1)
                    l_atom_adenine = substructTools.retrieveAdenine(l_atom_ligand)
                    if not ligand in d_ref.keys () : 
                        # stock in tempory dictionary for the reference
                        d_ref[ligand] = []
                        d_ref[ligand].append (l_atom_ligand)
                        d_ref[ligand].append (l_atom_adenine)
                        writePDBfile.coordinateSection(d_filout_pdb[ligand], l_atom_ligand, "HETATM", connect_matrix = 1)
                        continue
                    else : 
                        rotation, translocation =  superimpose.rigid_transform_3D(l_atom_adenine, d_ref[ligand][-1])
                        if rotation == None or translocation == None : 
                            continue
                        # rotation + translation
                        l_atom_lig_rotated = superimpose.applyTranformation(rotation, translocation, l_atom_in=l_atom_ligand)
                        # write PDB file and RMSE
#                         print "============"
#                         print ligand, pr_ref
#                         print len (l_atom_lig_rotated)
#                         print len (d_ref[ligand][0])
#                         print "============"
                        if len (l_atom_lig_rotated) != len (d_ref[ligand][0]) : 
                            continue
                        
                        writePDBfile.coordinateSection(d_filout_pdb[ligand], l_atom_lig_rotated, "HETATM", connect_matrix = 1)
                        RMSE_ligand = superimpose.rmse(d_ref[ligand][0], l_atom_lig_rotated)
                        d_filout_RMSE[ligand].write (str (pr_ref) + pr_substruct  + "\t" + str(RMSE_ligand) + "\n")
            
    # close files
    for lig in d_filout_pdb.keys () : 
        d_filout_pdb[lig].close ()
        d_filout_RMSE[lig].close ()

    for file_RMSE in l_file_RMSE : 
        runOtherSoft.Rhistogram(file_RMSE, "RMSE_Adenine")                                 
            
