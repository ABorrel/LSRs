import smileAnalysis
import pathManage
import parsePDB
import analysis
import superposeStructure
import writePDBfile
import tool
import runOtherSoft

from os import makedirs, listdir, path
from shutil import copy2
from re import search
import numpy as np

def globalArrangement (pr_orgin, p_smile, p_family, name_ligand):
    
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
#         print smile, l_PDB_query, l_PDB_ref, l_ligand
        replacement, metal = smileAnalysis.searchReplacement (smile, l_PDB_query[0], l_PDB_ref[0], name_ligand)
        
        if replacement == "metal" : 
            print metal, l_PDB_query, l_PDB_ref, name_ligand
        
        len_find = len (l_PDB_ref)
        i = 0
        while i < len_find : 
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
            
            pr_final = pr_orgin + replacement + "/" + l_PDB_ref[i] + "/" 
            pr_ligand = pr_orgin + replacement + "/" + l_PDB_ref[i] + "/LGD/"
            pr_BS =  pr_orgin + replacement + "/" + l_PDB_ref[i] + "/BS/"
            pr_sust = pr_orgin + replacement + "/" + l_PDB_ref[i] + "/" + str (subst) + "/"
            
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
        
        

    
def qualityExtraction (l_ligand, p_list_ligand, thresold_sheap) : 
    
    filout = open(pathManage.result() + "quality_extraction.txt", "w")
    
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
    
    
    
    
def countingSubstituent (pr_final_folder, debug = 1):
    
    
    d_count = {}
    d_lig = {}
    l_file_final = listdir(pr_final_folder)
    if debug : print "1", pr_final_folder
    for pr_type_substituent in l_file_final :
        sub_query = pr_type_substituent
        l_file_sub = listdir(pr_final_folder + pr_type_substituent + "/")
        if debug: print "2",pr_final_folder +  pr_type_substituent + "/"
        for ref_ID in l_file_sub : 
            l_ligand_sub = listdir(pr_final_folder + pr_type_substituent + "/" + ref_ID + "/")
            for ligand_sub in l_ligand_sub :
                l_file_ref  = listdir(pr_final_folder + pr_type_substituent + "/" + ref_ID + "/" + ligand_sub + "/")
                if not ligand_sub in d_count.keys () : 
                    d_count[ligand_sub] = {}
                if not sub_query in d_count[ligand_sub].keys () : 
                    d_count[ligand_sub] [sub_query] = 0
		d_count[ligand_sub][sub_query] = d_count[ligand_sub][sub_query] + 1
                
                for file_ref in l_file_ref : 
                    # print file_ref
                    if search ("LGD", file_ref):
                        ligand = file_ref.split ("_")[1]
                        if not ligand in d_lig.keys () : 
                            d_lig[ligand] = 0
                        else : 
                            d_lig[ligand] = d_lig[ligand] + 1
            
        
    # write and plot
    pr_result = pathManage.result("counting")
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
    
    runOtherSoft.barplotQuantity(pr_result + "count_ligand")
        
        
        
def enantiomer(l_ligand, pr_final) : 
    "to do file output"
    
    pr_enantiomer = pathManage.result("enantiomer")
    
    d_filout = {}
    l_filout = []
    for ligand in l_ligand : 
        d_filout[ligand] = open (pr_enantiomer + ligand + "_" + "distOx.txt" , "w")
        l_filout.append (pr_enantiomer + ligand + "_" + "distOx.txt")

    l_p_substruct = listdir(pr_final) 
    for pr_substruct in l_p_substruct : 
        print pr_substruct
        l_pr_ref = listdir(pr_final + pr_substruct + "/")
        for pr_ref in l_pr_ref : 
            print pr_ref
            l_file = listdir(pr_final + pr_substruct + "/" + pr_ref + "/LGD/")
            for name_file in l_file : 
                if search("REF_A",name_file) and   search(".pdb",name_file): 
                    ligand = name_file.split ("_")[2]
                    l_atom_ligand = parsePDB.loadCoordSectionPDB(pr_final + pr_substruct + "/" + pr_ref + "/LGD/" + name_file, "HETATM")
                    for atom_ligand in l_atom_ligand : 
                        if atom_ligand["name"] == "O4'" :
                            atom_O4 = atom_ligand
                        elif atom_ligand["name"] == "O5'" :
                            atom_O5 = atom_ligand
                    d_out = parsePDB.distanceTwoatoms(atom_O4, atom_O5)
                    d_filout[ligand].write (pr_ref + "_" + pr_substruct + "\t" + str (d_out) + "\n")
    
    # close files
    for lig in d_filout.keys () : 
        d_filout[lig].close ()
    
    for file_dist in l_filout : 
        runOtherSoft.Rhistogram(file_dist, "DistanceO4-O5")
                
                                
                                
                                
            
        
        
        
        
        
        
        
    
    
    
    
    
    
    
            
        
        
        
    
