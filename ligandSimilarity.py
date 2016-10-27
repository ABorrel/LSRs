from os import listdir
from copy import copy
import pathManage


def analyseLGDProximity(prclassif):

    print(prclassif)
    nameREF = prclassif.split("/")[-1]
    print(nameREF)

    prout = pathManage.result(nameREF + "_LGDsimilarity")
    print(prout)

    # extract for each reference LGD
    extractLGDfile(prclassif, prout)


def extractLGDfile(prclassif, prresult):
    """Extract from folder classification """

    loutref = []
    lfoldergroups = listdir(prclassif)
    for foldergroup in lfoldergroups:
        lrefprots = listdir(prclassif + "/" + foldergroup + "/")
        for refprot in lrefprots:
            if not refprot in loutref:
                pathManage.generatePath(prout + refprot)
            lfileLGD = listdir(prclassif + "/" + foldergroup + "/" + refprot + "/LGD/")
            for fileLGD in lfileLGD:
                copy(prclassif + "/" + foldergroup + "/" + refprot + "/LGD/" + fileLGD, prresult + refprot + "/" + fileLGD)


    return prresult
