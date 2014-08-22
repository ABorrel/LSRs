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
            
            pr_final = pr_orgin + replacement + "/"+ name_ligand + "_" + subst + "/" + l_PDB_ref[i] + "_" + family + "/" + l_PDB_query[i] + "/" 
            
            if not path.isdir(pr_final):
                makedirs (pr_final)
            
            copy2(PDB_ref, pr_orgin + replacement + "/" + name_ligand + "_" + subst + "/" + l_PDB_ref[i] + "_" + family + "/")
            writePDBfile.coordinateSection(pr_final + p_lig_query.split ("/")[-1], lig_query_parsed, recorder = "HETATM", header = p_lig_query.split ("/")[-1], connect_matrix = 1)
            copy2(p_ligand_ref, pr_orgin + replacement + "/" + name_ligand + "_" + subst + "/" + l_PDB_ref[i] + "_" + family + "/")
            copy2(p_frag_ref, pr_orgin + replacement + "/" + name_ligand + "_" + subst + "/" + l_PDB_ref[i] + "_" + family + "/")
            copy2(p_protein_query, pr_final)
            copy2(p_lig_substituate, pr_final)
            copy2(p_BS, pr_final)   
            
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
    
    
    
    
    
    
    
    
    
    
    
    




