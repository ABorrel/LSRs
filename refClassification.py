from os import listdir, remove
from re import search

import pathManage
import runOtherSoft
import parseTMalign



def classifRefProtein (pr_dataset, l_lig):
    
    pr_out = pathManage.result("clasifRef")
    
    l_file_ref = []
    for lig in l_lig : 
        pr_dataset = pathManage.dataset(lig)
        l_file_by_lig = listdir(pr_dataset)
        l_pr_ref_by_lig =[pr_dataset + x for x in l_file_by_lig]
        for pr_ref_by_lig in l_pr_ref_by_lig : 
            PDB_folder = pr_ref_by_lig.split ("/")[-1]
            try : l_file = listdir(pr_ref_by_lig)
            except : continue
            for file_ref in l_file : 
                if search("^" +PDB_folder, file_ref) :
                    l_file_ref.append (pr_ref_by_lig + "/" +file_ref)
                    break
        
    d_outTMalign = applyTMAlignList(l_file_ref, pr_out + "align/")
    
    writeMatrixTMalign (d_outTMalign, pr_out + "matrixRMSD", "RMSD")
    writeMatrixTMalign (d_outTMalign, pr_out + "matrixIDseq", "IDseq")
    writeMatrixTMalign (d_outTMalign, pr_out + "matrixTMscore1", "TMscore1")
    writeMatrixTMalign (d_outTMalign, pr_out + "matrixTMscore2", "TMscore2")


def writeMatrixTMalign (d_in, p_filout, k):
    
    # list PDB
    l_PDB = d_in.keys ()
    for k1 in d_in.keys () : 
        for k2 in d_in[k1].keys () : 
            if not k2 in l_PDB : 
                l_PDB.append (k2)
    
    filout = open (p_filout, "w")
    filout.write ("\t".join (l_PDB) + "\n")
    for PDB1 in l_PDB : 
        filout.write (PDB1)
        for PDB2 in l_PDB : 
            if PDB1 == PDB2 : 
                filout.write ("\t1")
            else : 
                try : filout.write ("\t" + d_in[PDB1][PDB2][k])
                except : filout.write ("\t" + d_in[PDB2][PDB1][k])
        filout.write ("\n")
    filout.close ()
                
                
                
                
                
                 
    

def applyTMAlignList (l_pr_ref, pr_out):
    
    
    pathManage.generatePath(pr_out)
    nb_pr_ref = len (l_pr_ref)
    d_out = {}
    i = 0
    while i < nb_pr_ref : 
        j = i + 1
        PDB1 = l_pr_ref[i].split ("/")[-1][0:4]
        # print PDB1
    
        while j < nb_pr_ref :
            PDB2 =  l_pr_ref[j].split ("/")[-1][0:4]
            # folder TM align
            pr_alignement = pr_out + PDB1 + "__" + PDB2 + "/"
            #print pr_alignement
            print PDB1,i, PDB2,j
            #print l_pr_ref[i]
            # RUN
            out_file = runOtherSoft.runTMalign(l_pr_ref[i], l_pr_ref[j], pr_alignement)
            # clean folders
            CleanResultTMalign (pr_alignement)
            # parse result
            if not PDB1 in d_out.keys () : 
                if not PDB2 in d_out.keys () : 
                    d_out[PDB1] = {}
                    d_out[PDB1][PDB2] = parseTMalign.parseOutputTMalign(out_file[-1])
                else : 
                    d_out[PDB2][PDB1] = {}
                    d_out[PDB2][PDB1] = parseTMalign.parseOutputTMalign(out_file[-1])
            else : 
                d_out[PDB1][PDB2] = {}
                d_out[PDB1][PDB2] = parseTMalign.parseOutputTMalign(out_file[-1])
            j = j + 1
        i = i + 1
    
    return d_out




def CleanResultTMalign (pr_TM_out):

    l_p_filout = listdir (pr_TM_out)
    for p_filout in l_p_filout : 
        if p_filout != "RMSD" and p_filout != "matrix.out" :
            remove (pr_TM_out + p_filout)






