from os import path


def parseOutputShaep (p_filin):
    d_out = {}
    
    if path.getsize(p_filin) == 0 : 
        return d_out
    
    filin = open (p_filin, "r")
    
    l_lines = filin.readlines ()
    
    
    l_desc = l_lines[0].strip ().split ("\t")[1:]
    l_val = l_lines[1].strip ().split ("\t")
    
    i = 0
    nb_desc = 3
    
    while i < nb_desc : 
        
        d_out[l_desc[i]] = l_val[i]
        i = i + 1

    return d_out

