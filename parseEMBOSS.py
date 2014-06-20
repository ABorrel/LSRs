"""
BORREL Alexandre
04-2012
"""
from re import search, sub

def embossFile (path_file_water):
    """Parse water or needle file, retrieve only first alignement
    arg: path of water file
    return: sequences aligned, similariry, identity"""
    
    filin = open(path_file_water, "r")
    lines = filin.read()
    filin.close()
    
    list_in_file = lines.split("#=======================================\n")
    
    carac = list_in_file[1]
    
    list_line_carac = carac.split("\n")

    for line_carac in list_line_carac : 
        # retrieve similarity and identity
        if search ("^# Identity:", line_carac) : 
            line_carac = sub("[ ]{2,}", " ", line_carac)
            list_element = line_carac.split(" ")
            identity = list_element[3]
            identity = identity.replace("(","")
            identity = identity.replace(")","")
        elif search ("^# Similarity:", line_carac) :
            line_carac = sub("[ ]{2,}", " ", line_carac)
            list_element = line_carac.split(" ")
            similarity = list_element[3]
            similarity = similarity.replace("(","")
            similarity = similarity.replace(")","")
    
    # retrieve 2 sequences aligned + initiale position   
    align = list_in_file[2]
    list_line_align = align.split("\n")
    number_lines = len(list_line_align)
    
    seq1 = ""
    seq2 = ""
    begin = int(list_line_align[1][13:21].replace(" ",""))
    for i in range(1,number_lines - 5, 4) : 
        seq1 = seq1 + list_line_align[i] [21:71]
        seq2 = seq2 + list_line_align[i + 2] [21:71]
    seq1 = list(seq1.split()[0])
    seq2 = list(seq2.split()[0])

    return [seq1, seq2, similarity, identity, begin]
    

