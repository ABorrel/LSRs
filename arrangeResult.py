import smileAnalysis
import pathManage
import parsePDB
import analysis
import superposeStructure
import writePDBfile


from os import path, makedirs, listdir, path
from shutil import copy2
from re import search

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
        
        

    



