from re import search
from os import path


def TSVFiltered (pTSV, lkselect, pfilout = "", debug = 1):
    """Parse file BindingDB"""

    if pfilout != "":
        if path.exists(pfilout):
            return fileFiltered(pfilout)

    l_dout = []
    filin = open (pTSV, "r")
    l_lines = filin.readlines()
    filin.close ()

    #HEADER
    l_header = l_lines[0].strip().split ("\t")

    if debug == 1:
        i = 0
        for h in l_header:
            #print i,h
            i = i + 1

    for lineTSV in l_lines[1:]:
        l_element = lineTSV.strip().split ("\t")
        l_PDB = l_element[27]
        # if PDB null pass
        if l_PDB == "":
            continue
        else:
            l_PDB = l_PDB.split (",")

        #l_identic = list(set(l_PDB).intersection(set(l_PDBsearch)))
        #if l_identic != []:
        d_out = {}
        i = 0
        while i < len (l_header):
            try:
                d_out[l_header[i]] = l_element[i]
                if l_element[i] == "":
                    d_out[l_header[i]] = "-"
            except : d_out[l_header[i]] = "-"
            i = i + 1
        # reduce the output
        for kin in d_out.keys():
            if not kin in lkselect:
                del d_out[kin]
        l_dout.append (d_out)
        #if search("4EOM", lineTSV) or search("3ULI", lineTSV):
        #    print(d_out)
        #    print(lineTSV)
        #    print("================")


    #print(l_dout[0].values())
    #print("\n")
    #print("\n".join(l_dout[0].keys()))
    #print(l_dout)


    # write file reduced
    if pfilout != "":
        filout = open(pfilout, "w")
        header = l_dout[0].keys()
        filout.write("\t".join(header) + "\n")
        for dout in l_dout:
            lw = []
            for h in header:
                lw.append(dout[h])
            filout.write("\t".join(lw) + "\n")
    return l_dout



def fileFiltered(pfilin):

    filin = open(pfilin, "r")
    llinesBinding = filin.readlines()
    filin.close()

    ldout = []

    lheader = llinesBinding[0].strip().split("\t")
    nblines = len(llinesBinding)

    i = 1
    while i < nblines:
        dout = {}
        lelem = llinesBinding[i].strip().split("\t")
        #print len(lheader), len(lelem)
        #print lheader, lelem
        j = 0
        for header in lheader:
            dout[header] = lelem[j]
            j = j + 1
        ldout.append(dout)
        i = i + 1

    return ldout
