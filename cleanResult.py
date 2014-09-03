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
                            
    

def cleanSmileFile (thresold_shaep, l_ligand_out, pr_result) : 
    
    l_pr_lig = listdir(pr_result)
    for pr_lig in l_pr_lig : 
        if len(pr_lig) != 3 : 
            continue
        else : 
            filin_control = open (pr_result + pr_lig + "/control.txt", "r")
            l_control = filin_control.readlines()
            filin_control.close ()
            d_control = {}
            for control in l_control : 
                l_element_control = control.strip ().split ("\t")
                sub  = l_element_control[0]
                ref = l_element_control[1]
                query = l_element_control[2]
                ligand = l_element_control[3]
                sheap = float(l_element_control[4])
                
                if not sub in d_control.keys () : 
                    d_control[sub] = {}
                if not ref in d_control[sub].keys () : 
                    d_control[sub][ref] = {}
                if not query in d_control[sub][ref].keys () : 
                    d_control[sub][ref][query] = {}
                if not ligand in d_control[sub][ref][query].keys () : 
                    d_control[sub][ref][query][ligand] = sheap
            
            #print d_control
            l_files = listdir(pr_result + pr_lig + "/")
            for files in l_files : 
                if search ("smile.txt", files) : 
                    sub = files.split ("_")[1]
                    filin_simle = open (pr_result + pr_lig + "/" + files, "r")
                    l_smile = filin_simle.readlines()
                    filin_simle.close ()
                    
                    filout_smile = open (pr_result + pr_lig + "/" + files, "w")
                    for smile in l_smile : 
                        #print smile, files, pr_lig, sub
                        l_elem1 = smile.strip().split ("\t")
                        l_ref = l_elem1[3].split (" ")
                        l_queries =  l_elem1[2].split (" ")
                        l_lig = l_elem1[4].split (" ")
                    
                        nb_queries = len (l_queries)
                        i = 0
                        while i < nb_queries :
                            #print "========="
                            #print sub 
                            #print l_ref[i]
                            #print l_queries[i]
                            #print l_lig[i]
                            #print "--------"
                            #print d_control["pi2"]["1RYR"]
                            try : score_shaep = d_control[sub][l_ref[i]][l_queries[i]][l_lig[i]] 
                            except : score_shaep = 0.0
                            if score_shaep < thresold_shaep : 
                                del l_ref[i]
                                del l_queries[i]
                                del l_lig[i]
                                nb_queries = nb_queries - 1
                            elif l_lig[i] in l_ligand_out : 
                                del l_ref[i]
                                del l_queries[i]
                                del l_lig[i]
                                nb_queries = nb_queries - 1
                            else :
                                i = i + 1
                        if len (l_ref) != 0 : 
                            filout_smile.write (str(l_elem1[0]) + "\t" + str (len (l_queries)) + "\t" + " ".join(l_queries) + "\t" + " ".join (l_ref) + "\t" + " ".join (l_lig) + "\n")
                    filout_smile.close ()
                            



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
#cleanResultFolder (thresold_shaep, l_ligand_out, pr_result)
cleanSmileFile (thresold_shaep, l_ligand_out, pr_result)


