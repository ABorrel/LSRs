

def parseOutputTMalign (p_filin):
    
    d_out = {}
    
    filin = open (p_filin, "r")
    f_read = filin.read()
    filin.close ()
    
    d_out["RMSD"] = retrieveRMSD (f_read)
    d_out["IDseq"] = retrieveIDseq (f_read)
    
    return d_out







def retrieveRMSD (filin_read) :
    """
    """
    part_file = filin_read.split ("RMSD=")[1]
    RMSD = part_file.split (",")[0]
    try : RMSD = RMSD.replace (" ", "")
    except : pass
    
    return RMSD
    
    



def retrieveIDseq (filin_read) : 
    """
    Retrieve RMSD in TMalign out file
    args: -> path file RMSD
    return: value of RMSD
    """
    
    IDseq = filin_read.split ("Seq_ID=n_identical/n_aligned=")[1].split ("\n")[0]
    try : IDseq = IDseq.replace (" ", "")
    except : pass
    
    return IDseq



# print parseOutputTMalign("/home/borrel/Yue_project/alignment/12AS_A_2.20__12AS_B_2.20/RMSD")
