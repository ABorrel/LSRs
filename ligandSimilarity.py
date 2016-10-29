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
    buildMatrixSimilarity(prout)


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
        lligname = []
        i = 0
        nbsmile = len(lpsmile)
        while i < nbsmile - 1:
            j = i + 1
            while j < nbsmile:
                name1 = lpsmile[i][:-4].split("LGD_")[-1]
                name2 = lpsmile[j][:-4].split("LGD_")[-1]
                if not name1 in lligname:
                    lligname.append(name1)
                if not name2 in lligname:
                    lligname.append(name2)
                print(name1, name2)
                if not name1 in dresult.keys():
                    dresult[name1] = {}
                if not name2 in dresult[name1].keys():
                    dresult[name1][name2] = {}
                if MCS == 1:
                    tanimotoMCS = generateMCS.get_Tanimoto(lpsmile[i], lpsmile[j])
                    dresult[name1][name2]["MCS"] = tanimotoMCS
                    print(tanimotoMCS)
                if Sheap == 1:
                    pshaep = prin + refprot + "/outsheap.txt"
                    runOtherSoft.runShaep(lppdb[i], lppdb[j], pshaep, clean=1)
                    dsheap = parseShaep.parseOutputShaep(pshaep)
                    print(dsheap)
                    dresult[name1][name2]["ESP"] = dsheap["ESP_similarity"]

                    #control same value
                    #runOtherSoft.runShaep(lppdb[j], lppdb[i], pshaep, clean=1)
                    #dsheap = parseShaep.parseOutputShaep(pshaep)
                    #print(dsheap)
                j = j + 1
            i = i + 1

        if MCS == 1:
            filoutMCS = open(prin + refprot + "/matriceMCS", "w")
            filoutMCS.write("\t".join(lligname) + "\n")
        if Sheap == 1:
            filoutSheap = open(prin + refprot + "/matriceSheap", "w")
            filoutSheap.write("\t".join(lligname) + "\n")

        for namelig in lligname:
            if MCS == 1:
                filoutMCS.write(namelig)
            if Sheap == 1:
                filoutSheap.write(namelig)

            for nameligcol in lligname:
                if MCS == 1:
                    if nameligcol == namelig:
                        filoutMCS.write("\t1.0")
                    else:
                        try: filoutMCS.write("\t" + str(dresult[namelig][nameligcol]["MCS"]))
                        except: filoutMCS.write("\t" + str(dresult[nameligcol][namelig]["MCS"]))
                if Sheap == 1:
                    if nameligcol == namelig:
                        filoutSheap.write("\t1.0")
                    else:
                        try: filoutSheap.write("\t" + str(dresult[namelig][nameligcol]["ESP"]))
                        except: filoutSheap.write("\t" + str(dresult[nameligcol][namelig]["ESP"]))
            if MCS == 1:
                filoutMCS.write("\n")
            if Sheap == 1:
                filoutSheap.write("\n")
        if MCS == 1:
            filoutMCS.close()
        if Sheap == 1:
            filoutSheap.close()

        dddd
    return
