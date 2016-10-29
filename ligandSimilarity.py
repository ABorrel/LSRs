from os import listdir
from shutil import copyfile

import pathManage
import generateMCS
import runOtherSoft
import parseShaep

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





def buildMatrixSimilarity(prin, MCS=1, Sheap=1, pfileaffinity=1):

    lrefprot = listdir(prin)
    for refprot in lrefprot:
        print(refprot)
        lpsmile = []
        lppdb = []
        lfileref = listdir(prin + refprot)
        # extract smile
        for fileref in lfileref:
            print(fileref[-3:])
            if fileref[-3:] == "smi":
                lpsmile.append(prin + refprot + "/" + fileref)
            elif fileref[-3:] == "pdb":
                lppdb.append(prin + refprot + "/" + fileref)

        # test same number of pdb files and smi
        if len(lppdb) != len(lpsmile):
            print("Error -> l.80 ligandSimilarity.py")
            return

        dresult={}
        i = 0
        nbsmile = len(lpsmile)
        while i < nbsmile - 1:
            j = i + 1
            while j < nbsmile:
                name1 = lpsmile[i][:-4]
                name2 = lpsmile[j][:-4]
                print(name1, name2)
                if not name1 in dresult.keys():
                    dresult[name1] == {}
                if not name2 in dresult[name1].keys():
                    dresult[name1][name2] = {}
                if MCS == 1:
                    tanimotoMCS = generateMCS.get_Tanimoto(lpsmile[i], lpsmile[j])
                    dresult[name1][name2]["MCS"] = tanimotoMCS
                    print(tanimotoMCS)
                if Sheap == 1:
                    pshaep = prin + refprot + "/outsheap.txt"
                    runOtherSoft.runShaep(lppdb[i], lppdb[j], pshaep)
                    dsheap = parseShaep.parseOutputShaep(pshaep)
                    print(dsheap)

                    runOtherSoft.runShaep(lppdb[j], lppdb[i], pshaep, clean=1)
                    dsheap = parseShaep.parseOutputShaep(pshaep)
                    print(dsheap)
                    dddd
                # stock
                j = j + 1
            i = i + 1

    return
