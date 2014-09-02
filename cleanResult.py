import analysis
import pathManage
import tool
import parseTMalign


from re import search
from os import listdir, remove
from shutil import rmtree


# clean after run (not clean)
def cleanResultFolder (thresold_sheap, l_lig_out, pr_result):
    
    l_pr_lig = listdir(pr_result)
    for pr_lig in l_pr_lig : 
        if len(pr_lig) != 3 : 
            continue
        else : 
            filout_control = open (pr_result + pr_lig + "/control.txt", "w")
            l_pr_ref = listdir(pr_result + pr_lig + "/")
            for pr_ref in l_pr_ref : 
                if len (pr_ref) == 4 : 
                    l_file_query = listdir(pr_result + pr_lig + "/" + pr_ref + "/")
                    for file_query in l_file_query : 
                        if search ("^all", file_query) : 
                            continue
                        elif search ("txt$", file_query) : 
                            remove (pr_result + pr_lig + "/" + pr_ref + "/" + file_query)   
                            continue
                        l_elem_name = file_query.split ("_")
                        lig = l_elem_name[1]
                        pdb_q = l_elem_name[2]
                        sub = l_elem_name[3]
                        if len(l_elem_name) == 5 : 
                            sheap_score = float (l_elem_name[4][0:-4]) 
                        else : 
                            sheap_score = 1.0 #case file CX-BS
                        if lig in l_lig_out : 
                            #print pr_result + pr_lig + "/" + pr_ref + "/" + file_query
                            remove(pr_result + pr_lig + "/" + pr_ref + "/" + file_query)
                            continue
                        else :
                            if len(l_elem_name) == 5 : filout_control.write (str (sub) + "\t" + str (pr_ref) + "\t" + str (pdb_q) + "\t" + str (lig) + "\t" + str (sheap_score) + "\n")
                            if sheap_score < 0.2 : 
                                #print pr_result + pr_lig + "/" + pr_ref + "/" + file_query
                                remove (pr_result + pr_lig + "/" + pr_ref + "/substituent_" + str (lig) + "_" + str(pdb_q) + "_" + str (sub) + ".hit")
                                remove (pr_result + pr_lig + "/" + pr_ref + "/substituent_" + str (lig) + "_" + str (pdb_q) + "_" + str (sub) + ".smi")
                                remove(pr_result + pr_lig + "/" + pr_ref + "/" + file_query) 
                    
                    # del ref folder
                    l_file_query_new = listdir(pr_result + pr_lig + "/" + pr_ref + "/")
                    f = 0
                    for query_new in l_file_query : 
                        if search ("substituent", query_new) : 
                            f = 1
                            break
                    if f==0 : 
                        rmtree(pr_result + pr_lig + "/" + pr_ref + "/")
                
            filout_control.close ()      
                            
    
    
    
    
    




# step 6
# Analysis  
# - smile, filtering
# - shaep on substructure and ligand

def analysisSmile (substruct):

    l_p_smile = pathManage.findListSmileFile(substruct) 
    for p_smile in l_p_smile : 
        analysis.selectSmileCode(p_smile, minimal_length_smile = 4)
    

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
  

####################
###   MAIN     #####
####################
# constante
thresold_RX = 2.7
thresold_BS = 4.5
thresold_blast = 1e-100
thresold_superimposed_ribose = 2.5
thresold_superimposed_pi = 3
thresold_IDseq = 100
thresold_shaep = 0.2
l_ligand_out = ["AMP", "ADP", "ATP", "TTP", "DCP", "DGT", "DTP", "DUP", "ACP", "AD9", "NAD", "AGS", "UDP", "POP", "APC", "CTP", "AOV"]


# main #
########
pr_result = pathManage.result()
cleanResultFolder (thresold_shaep, l_ligand_out, pr_result)


# filter smile + analyse same BS
# l_ligand = ["AMP", "ADP", "ATP", "POP"]
# for lig in l_ligand : 
#     analysisSmile ("AMP")
#     analysisSameBS ("AMP")


