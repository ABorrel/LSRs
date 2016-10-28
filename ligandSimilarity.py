from os import listdir
from shutil import copyfile

import pathManage
import generateMCS

def analyseLGDProximity(prclassif):

    print(prclassif)
    nameREF = prclassif.split("/")[-1]
    print(nameREF)

    prout = pathManage.result(nameREF + "_LGDsimilarity")
    print(prout)

    # extract for each reference LGD
    extractLGDfile(prclassif, prout)
    buildMatrixMCS(prout)


def extractLGDfile(prclassif, prresult):
    """Extract from folder classification """

        # test if file in folder result
    if len(listdir(prresult)) > 1:
        return prresult


    lprref = []
    lfoldergroups = listdir(prclassif)
    for foldergroup in lfoldergroups:
        if foldergroup == "cycle":
            lsubtypes = listdir(prclassif + "/cycle/")
            for subtype in lsubtypes:
                lrefprot = listdir(prclassif + "/cycle/" + subtype)
                for refprot in lrefprot:
                    lprref.append(prclassif + "/cycle/" + subtype + "/" + refprot)
        else:
            lrefprot = listdir(prclassif + "/" +  foldergroup + "/")
            for refprot in lrefprot:
                lprref.append(prclassif + "/" + foldergroup + "/" + refprot)


    lout = []
    for prefprot in lprref:
        refprot = prefprot.split("/")[-1]
        if not refprot in lout:
            pathManage.generatePath(prresult + refprot)
            lout.append(refprot)
        # copy file LGD
        lfileLGD = listdir(prefprot + "/LGD/")
        for fileLGD in lfileLGD:
            copyfile(prefprot + "/LGD/" + fileLGD, prresult + refprot + "/" + fileLGD)


    return prresult





def buildMatrixMCS(prin):

    lrefprot = listdir(prin)
    for refprot in lrefprot:
        print(refprot)
        lpsmile = []
        lfileref = listdir(prin + refprot)
        # extract smile
        for fileref in lfileref:
            print(fileref[-3:])
            if fileref[-3:] == "smi":
                lpsmile.append (prin + refprot + "/" + fileref)

        i = 0
        nbsmile = len(lpsmile)
        while i < nbsmile - 1:
            j = i + 1
            while j < nbsmile:
                tanimotoMCS = generateMCS.get_Tanimoto(lpsmile[i], lpsmile[j])
                print(tanimotoMCS)
                # stock
                j = j + 1
            i = i + 1



    return
