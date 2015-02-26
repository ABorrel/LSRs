from os import listdir, remove
from re import search

import pathManage
import runOtherSoft
import parseTMalign
import downloadFile
import parseEMBOSS



def classifRefProtein (pr_dataset, l_lig):
    
    pr_out = pathManage.result("clasifRef")
    
    # case fasta file
    pr_align_seq = pathManage.generatePath(pr_out + "alignSeq/")
    l_p_fasta = []
    for lig in l_lig : 
        pr_dataset = pathManage.dataset(lig)
        l_file_by_lig = listdir(pr_dataset)
        l_pr_ref_by_lig =[pr_dataset + x for x in l_file_by_lig]
        for pr_ref_by_lig in l_pr_ref_by_lig : 
            PDB_folder = pr_ref_by_lig.split ("/")[-1]
            
            try : l_file = listdir(pr_ref_by_lig)
            except : continue
            for file_ref in l_file : 
                if search("^" + PDB_folder, file_ref) :
                    PDB_ID = file_ref[0:-4]
                    PDB_ID = PDB_ID[0:4].lower () + PDB_ID[4:]
                    # PDB ID with chain associated
                    p_fasta = downloadFile.importFasta(PDB_ID, pr_align_seq, dir_by_PDB = 0, debug = 1, fastaGlobal = "/home/borrel/Yue_project/pdb_seqres.txt")
                    l_p_fasta.append (p_fasta)
                    break
            
                
    d_outNeedle = applyNeedleList (l_p_fasta, pr_align_seq)
    
    
    writeMatrixFromDico (d_outNeedle, pr_out + "matrixSimilarSeq","similarity" )
    writeMatrixFromDico (d_outNeedle, pr_out + "matrixIDSeq","identity" )    
                
    
    # retrieve list of file -> case PDB file
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
                if search("^" + PDB_folder, file_ref) :
                    l_file_ref.append (pr_ref_by_lig + "/" +file_ref)
                    break
    
    
    d_outTMalign = applyTMAlignList(l_file_ref, pr_out + "align/")

    writeMatrixFromDico (d_outTMalign, pr_out + "matrixRMSD", "RMSD")
    writeMatrixFromDico (d_outTMalign, pr_out + "matrixIDseqTMalign", "IDseq")
    writeMatrixFromDico (d_outTMalign, pr_out + "matrixTMscore1", "TMscore1")
    writeMatrixFromDico (d_outTMalign, pr_out + "matrixTMscore2", "TMscore2")


def writeMatrixFromDico (d_in, p_filout, k):
    
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
                if k == "similarity" or k == "identity" : 
                    filout.write ("\t100")
                elif k == "RMSD" : 
                    filout.write ("\t0")
                else : 
                    filout.write ("\t1")
            else : 
                #print PDB1, PDB2
                try : filout.write ("\t" + d_in[PDB1][PDB2][k])
                except : 
                    try : filout.write ("\t" + d_in[PDB2][PDB1][k])
                    except : filout.write ("\tNA")
        filout.write ("\n")
    filout.close ()
                
                
                
def applyNeedleList (l_file_fasta, pr_out): 
    
    
    nb_pr_ref = len (l_file_fasta)
    d_out = {}
    i = 0
    while i < nb_pr_ref : 
        j = i + 1
        PDB1 = l_file_fasta[i].split ("/")[-1][0:4]
        # print PDB1
    
        while j < nb_pr_ref :
            PDB2 =  l_file_fasta[j].split ("/")[-1][0:4]
            # folder TM align
            
            p_outfile = runOtherSoft.needle (l_file_fasta[i], l_file_fasta[j], pr_out + PDB1 + "__" + PDB2 + ".needle")
            
            out = parseEMBOSS.embossFile(p_outfile)
            # list 0-> seq1, 1-> seq2, 2-> similarity, 3->identity
            # print out[2:]
            # parse result
            if not PDB1 in d_out.keys () : 
                if not PDB2 in d_out.keys () : 
                    d_out[PDB1] = {}
                    d_out[PDB1][PDB2] = {}
                    d_out[PDB1][PDB2]["identity"] = out[3].replace ("%", "")
                    d_out[PDB1][PDB2]["similarity"] = out[2].replace ("%", "")
                else : 
                    d_out[PDB2][PDB1] = {}
                    d_out[PDB2][PDB1]["identity"] = out[3].replace ("%", "")
                    d_out[PDB2][PDB1]["similarity"] = out[2].replace ("%", "")
            else : 
                d_out[PDB1][PDB2] = {}
                d_out[PDB1][PDB2]["identity"] = out[3].replace ("%", "")
                d_out[PDB1][PDB2]["similarity"] = out[2].replace ("%", "")
            j = j + 1
        i = i + 1
    
    return d_out          
 
                
                 
    

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
            #print PDB1,i, PDB2,j
            #print l_pr_ref[i]
            # RUN
            out_file = runOtherSoft.runTMalign(l_pr_ref[i], l_pr_ref[j], pr_alignement)
            # clean folders -> pb with several run -> clean too fast -> try / except
            try : CleanResultTMalign (pr_alignement)
            except : pass
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

    try : l_p_filout = listdir (pr_TM_out)
    except : return
    for p_filout in l_p_filout : 
        if p_filout != "RMSD" and p_filout != "matrix.out" :
            remove (pr_TM_out + p_filout)






