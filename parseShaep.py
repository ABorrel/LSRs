from os import path


def parseOutputShaep (p_filin):
    d_out = {}
    
    if path.exists(p_filin) : 
        if path.getsize(p_filin) == 0 : 
            return d_out
    
    filin = open (p_filin, "r")
    
    l_lines = filin.readlines ()
    filin.close ()
    if len (l_lines) < 2 : 
        return d_out
    
    l_desc = l_lines[0].strip ().split ("\t")[1:]
    l_val = l_lines[1].strip ().split ("\t")[1:]
    
    i = 0
    nb_desc = 3
    
    while i < nb_desc : 
        
        d_out[l_desc[i]] = l_val[i]
        i = i + 1

    return d_out




def valueShapeSimilarity (p_filin) : 
    dic = parseOutputShaep (p_filin)
    if not "shape_similarity" in dic.keys () : 
        return 0.0
    return float (dic["shape_similarity"])
    
    
    
