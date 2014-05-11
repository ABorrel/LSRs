import smileAnalysis
import pathManage
import parsePDB
import analysis

from os import path, makedirs
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
        l_PDB_ref = line_smile.split ("\t")[-2].split (" ")
        l_ligand = line_smile.strip().split ("\t")[-1].split (" ")
        
        # search replacement
        smile = line_smile.split ("\t")[0]
        replacement = smileAnalysis.searchReplacement (smile)
        
        
        len_find = len (l_PDB_ref)
        i = 0
        while i < len_find : 
            family = analysis.findFamily(l_PDB_ref[i], p_family)
            
            
            pr_final = pr_orgin + replacement + "/" + family + "/" + name_ligand + "/" + subst + "/" + l_PDB_ref[i] + "/" + l_PDB_query[i] + "/" + l_ligand[i] + "/" 
            
            if not path.isdir(pr_final):
                makedirs (pr_final)
            
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
            p_lig_query = pathManage.findligandQuery(pr_dataset , l_ligand[i], l_PDB_query[i])
            p_lig_substituate = pathManage.findSubstructFind(pr_result, l_ligand[i], l_PDB_query[i])
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
            
            copy2(PDB_ref, pr_orgin + replacement + "/" + family + "/" + name_ligand + "/" + subst + "/" + l_PDB_ref[i] + "/")
            copy2(p_ligand_ref, pr_orgin + replacement + "/" + family + "/" + name_ligand + "/" + subst + "/" + l_PDB_ref[i] + "/")
            copy2(p_frag_ref, pr_orgin + replacement + "/" + family + "/" + name_ligand + "/" + subst + "/" + l_PDB_ref[i] + "/")
            
            copy2(p_protein_query, pr_final)
            copy2(p_lig_query, pr_final)
            copy2(p_lig_substituate, pr_final)
            copy2(p_BS, pr_final)   
            
            i = i + 1
    
    return 1


