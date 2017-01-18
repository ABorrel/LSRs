from os import listdir, remove, path
from shutil import copyfile
from re import search

import pathManage
import generateMCS
import runOtherSoft
import parseShaep
import parseTSV
import tool

PBINDINGDB = "/home/borrel/BindingDB_All.tsv"

def analyseLGDProximity(prclassif):

    print(prclassif)
    nameREF = prclassif.split("/")[-1]
    print(nameREF)

    prout = pathManage.result(nameREF + "_LGDsimilarity")
    print(prout)

    # extract IC550 for PDB and ligand
    pbindingDBfiltered = prout + "bindingDBfiltered.txt"
    lkeep = [ "PDB ID(s) for Ligand-Target Complex", "Ligand HET ID in PDB", "Kd (nM)", "Ki (nM)", "IC50 (nM)"]
    parseTSV.TSVFiltered(PBINDINGDB, lkeep, pfilout=pbindingDBfiltered)

    # extract for each reference LGD
    extractLGDfile(prclassif, prout)
    buildMatrixSimilarity(prout, pfileaffinity=pbindingDBfiltered, MCS=1, Sheap=0)

    # extract MMP
    extractMMP(prout)



def extractLGDfile(prclassif, prresult):
    """Extract from folder classification """

    # test if file in folder result
    #if len(listdir(prresult)) > 1:
    #    return prresult


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
    dLSR = {}
    ltypeLSR = ["pi1", "pi2", "pi3"]
    for prefprot in lprref:#########################################to reduce
        refprot = prefprot.split("/")[-1]
        if not refprot in lout:
            pathManage.generatePath(prresult + refprot)
            lout.append(refprot)
        # copy file LGD
        lfileLGD = listdir(prefprot + "/LGD/")
        for fileLGD in lfileLGD:
            ligid = fileLGD.split("_")[1]
            if ligid == "REF":
                ligid = fileLGD.split("_")[2]
                pdbid = refprot.split("_")[-1]
                LSR = "REF"
            else:
                pdbid = fileLGD.split("_")[2]
                LSR = prefprot.split("/")[-2].replace("_", "")
                if prefprot.split("/")[-3] == "cycle":
                    LSR = "cycle-" + str(LSR)
            nameout = str(LSR) + "_" + str(ligid) + "_" + str(pdbid) + str(fileLGD[-4:])
            copyfile(prefprot + "/LGD/" + fileLGD, prresult + refprot + "/" + nameout)

        # extract SMILES LSR
        folderresult = prresult + refprot + "/"
        if not folderresult in dLSR.keys():
            dLSR[folderresult] = {}

        prLSRin = prefprot + "/LSR/"
        lfileLSR = listdir(prLSRin)
        for fileLSR in lfileLSR:
            #print prefprot + "/LSR/" + fileLSR,"l93===ligandSimilarity"
            if search("^LSR", fileLSR) and search("pdb", fileLSR):
                lelemname = fileLSR.split("_")
                nameLSR = lelemname[1]
                if nameLSR == "REF":
                    continue
                else:
                    lig = lelemname[2]
                    PDBid = lelemname[3]
                    smiles =  runOtherSoft.babelConvertPDBtoSMILE (prLSRin + fileLSR, rm_smi = 1)
                    #print(smiles, "l101 - ligandSimilarity")
                    kin = str(lig) + "-" + PDBid
                    if not kin in dLSR[folderresult].keys():
                        dLSR[folderresult][kin] = {}
                        for typeLSR in ltypeLSR:
                            dLSR[folderresult][kin][typeLSR] = "-"
                    dLSR[folderresult][kin][nameLSR] = smiles
    #print dLSR
    # write filout
    for folderresult in dLSR.keys():
        pfileLSR = folderresult + "listLSRsmiles"

        if path.exists(pfileLSR):
            filoutLSR = open(pfileLSR, "a")
        else:
            filoutLSR = open(pfileLSR, "w")
            filoutLSR.write("\t".join(ltypeLSR) + "\n")

        for kin in dLSR[folderresult].keys():
            lsmiles = [dLSR[folderresult][kin][i] for i in ltypeLSR]
            print lsmiles, "l.122 ligandSimilarity.py"
            filoutLSR.write(kin + "\t" + "\t".join(lsmiles) + "\n")
        filoutLSR.close()
    return prresult



def buildMatrixSimilarity(prin, MCS=1, Sheap=1, pfileaffinity=""):

    lrefprot = listdir(prin)
    for refprot in lrefprot:
        if not path.isdir(prin + refprot):
            continue
        lpsmile = []
        lppdb = []
        lfileref = listdir(prin + refprot)
        # extract smile
        for fileref in lfileref:
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
        while i < nbsmile:
            j = i
            while j < nbsmile:
                name1 = lpsmile[i][:-4].split("/")[-1]
                name2 = lpsmile[j][:-4].split("/")[-1]
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
                    lMCS = generateMCS.get_Tanimoto(lpsmile[i], lpsmile[j])
                    dresult[name1][name2]["MCS"] = lMCS[0]
                    dresult[name1][name2]["MAXdiff"] = lMCS[1]
                    #print(tanimotoMCS)
                if Sheap == 1:
                    pshaep = prin + refprot + "/outsheap.txt"
                    runOtherSoft.runShaep(lppdb[i], lppdb[j], pshaep, clean=1)
                    dsheap = parseShaep.parseOutputShaep(pshaep)
                    #print(dsheap)
                    dresult[name1][name2]["ESP"] = dsheap["ESP_similarity"]

                    #control same value
                    #runOtherSoft.runShaep(lppdb[j], lppdb[i], pshaep, clean=1)
                    #dsheap = parseShaep.parseOutputShaep(pshaep)
                    #print(dsheap)
                j = j + 1
            i = i + 1

        # remove sheap txt
        if Sheap == 1:
            remove(pshaep)

        # load affinity if available
        if pfileaffinity != "":
            laffinity = parseTSV.fileFiltered(pfileaffinity)
            pfiloutaff = prin + refprot + "/affinity"
            filoutaff = open(pfiloutaff, "w")
            filoutaff.write("IC50(nM)\n")
            for namelig in lligname:
                #if search("^REF", namelig):
                #    ligID = namelig.split("_")[1]
                #    PDBid = namelig.split("_")[2]
                #else:
                ligID = namelig.split("_")[1]
                PDBid = namelig.split("_")[2]

                # parse list known
                aff = 0
                for affinity in laffinity:
                    if affinity["Ligand HET ID in PDB"] == ligID:
                        if search(PDBid, affinity["PDB ID(s) for Ligand-Target Complex"]):
                            filoutaff.write(str(namelig) + "\t" + str(affinity["IC50 (nM)"]) + "\n")
                            aff = 1
                            break
                if aff == 0:
                    filoutaff.write(str(namelig) + "\t-\n")
            filoutaff.close()


        # write matrice
        if MCS == 1:
            ptanimoto = prin + refprot + "/matriceMCSTanimoto"
            pnbatomdiff = prin + refprot + "/matriceMCSNbAtomDiff"
            filouttanimoto = open(ptanimoto, "w")
            filoutNbatom = open(pnbatomdiff, "w")
            filouttanimoto.write("\t".join(lligname) + "\n")
            filoutNbatom.write("\t".join(lligname) + "\n")

        if Sheap == 1:
            filoutSheap = open(prin + refprot + "/matriceSheap", "w")
            filoutSheap.write("\t".join(lligname) + "\n")

        for namelig in lligname:
            if MCS == 1:
                filouttanimoto.write(namelig)
                filoutNbatom.write(namelig)
            if Sheap == 1:
                filoutSheap.write(namelig)

            for nameligcol in lligname:
                if MCS == 1:
                    #if nameligcol == namelig:
                    #    filouttanimoto.write("\t1.0")
                    #    filoutNbatom.write("\t-")
                    #else:
                    try: filouttanimoto.write("\t" + str(dresult[namelig][nameligcol]["MCS"]))
                    except: filouttanimoto.write("\t" + str(dresult[nameligcol][namelig]["MCS"]))

                    try: filoutNbatom.write("\t" + str(dresult[namelig][nameligcol]["MAXdiff"]))
                    except: filoutNbatom.write("\t" + str(dresult[nameligcol][namelig]["MAXdiff"]))

                if Sheap == 1:
                    if nameligcol == namelig:
                        filoutSheap.write("\t1.0")
                    else:
                        try: filoutSheap.write("\t" + str(dresult[namelig][nameligcol]["ESP"]))
                        except: filoutSheap.write("\t" + str(dresult[nameligcol][namelig]["ESP"]))
            if MCS == 1:
                filouttanimoto.write("\n")
                filoutNbatom.write("\n")
            if Sheap == 1:
                filoutSheap.write("\n")
        if MCS == 1:
            filouttanimoto.close()
            filoutNbatom.close()
        if Sheap == 1:
            filoutSheap.close()

        # plot matice
        if MCS == 1:
            runOtherSoft.plotMatrice(ptanimoto, pfiloutaff, pnbatomdiff, prin + refprot + "/listLSRsmiles" )
        if Sheap == 1:
            runOtherSoft.plotMatrice(prin + refprot + "/matriceSheap", pfiloutaff)


def extractMMP(prin, maxNbatom = 3, verbose = 0):

    pfilout = prin + "MMP.txt"
    filout = open(pfilout, "w")
    header = "LGD1-PDB1\tLGD2-PDB2\tLGDsmile1\tLGDsmile2\tLSRsLGD1\tLSRLGD2\tIC50(nM) LGD1\tIC50(nM) LGD2\tNb diff atom\tNAMS-Tanimoto\n"
    filout.write(header)

    lrefprot = listdir(prin)
    print len(lrefprot)
    ltemp = []# control repetition
    for refprot in lrefprot:
        print prin + refprot
        if not path.isdir(prin + refprot):
            print "AAAAaA"
            continue
        pfileMSCatomdiff = prin + refprot + "/matriceMCSNbAtomDiff"
        if not path.exists(pfileMSCatomdiff):
            print "Error: " + pfileMSCatomdiff
        else:
            dMSCatomdiff = tool.matrriceFileTODict(pfileMSCatomdiff)
            #print pfileMSCatomdiff
            #print dMSCatomdiff
            #print "###############"

        pfileMSCTanimoto = prin + refprot + "/matriceMCSTanimoto"
        if not path.exists(pfileMSCTanimoto):
            print "Error: " + pfileMSCTanimoto
        else:
            dMCSTanimoto = tool.matrriceFileTODict(pfileMSCTanimoto)
            #print pfileMSCatomdiff
            #print dMSCatomdiff
            #print "###############"




        paffinity = prin + refprot + "/affinity"
        if not path.exists(paffinity):
            print "Error: " + paffinity
        else:
            daff = tool.matrriceFileTODict(paffinity)

        plistSMILES = prin + refprot + "/listLSRsmiles"
        if not path.exists(plistSMILES):
            print "Error: " + plistSMILES
        else:
            dLSR = tool.matrriceFileTODict(plistSMILES)
            #print dLSR

        # extract smile LGD
        lfileref = listdir(prin + refprot + "/")
        dlig = {}
        for fileref in lfileref:
            if search(".smi", fileref):
                fsmile = open(prin + refprot + "/" + fileref, "r")
                smile = fsmile.readlines()[0].strip()
                fsmile.close()
                smile = smile.split("\t")[0]
                dlig[fileref.split(".")[0]] = smile

        #print len(dMSCatomdiff)
        # find MMP
        #print dMSCatomdiff
        for k1 in dMSCatomdiff.keys():
            for k2 in dMSCatomdiff[k1].keys():
                if k1 == k2:
                    continue
                else:
                    #print dMSCatomdiff[k1][k2], "===="
                    nbatomdiff = int(dMSCatomdiff[k1][k2].split("-")[-1])
                    scoreMCS = float(dMCSTanimoto[k1][k2])
                    if nbatomdiff <= maxNbatom and scoreMCS >= 0.9:
                        LGD1 = k1.split("_")[1]
                        PDB1 = k1.split("_")[2]
                        LGD2 = k2.split("_")[1]
                        PDB2 = k2.split("_")[2]
                        # if same ligand continue
                        if LGD1 == LGD2:
                            continue
                        # affinity
                        for kaff in daff.keys():
                            if kaff == k1:
                                aff1 = daff[kaff]["IC50(nM)"]
                            elif kaff == k2:
                                aff2 = daff[kaff]["IC50(nM)"]

                        # LSR
                        LSR1 = ""
                        LSR2 = ""
                        for klsr in dLSR.keys():
                            if klsr == str(LGD1 + "-" + PDB1):
                                for klsr2 in dLSR[klsr].keys():
                                    LSR1 = str(LSR1) + " " + str(dLSR[klsr][klsr2])
                            elif klsr == str(LGD2 + "-" + PDB2):
                                for klsr2 in dLSR[klsr].keys():
                                    LSR2 = str(LSR2) + " " + str(dLSR[klsr][klsr2])
                        if LSR1 == "":
                            LSR1 = "- - - "
                        if LSR2 == "":
                            LSR2 = "- - - "


                        if verbose ==1:
                            print "#########################"
                            print prin + refprot
                            print LGD1, "LGD1"
                            print LGD2, "LGD2"
                            print dlig[k1], "SMI1"
                            print dlig[k2], "SMI2"
                            print PDB1, "PDB1"
                            print PDB2, "PDB2"
                            print LSR1, "LSR1"
                            print LSR2, "LSR2"
                            print aff1, "aff1"
                            print aff2, "aff2"
                            print nbatomdiff, "diff"
                            print "$$$$$$$$$$$$$$$$$$$$$$$$$"
                        ktemp = PDB1 + "-" + LGD1 + "_" + PDB2 + "-" + LGD2
                        if ktemp in ltemp:
                            continue
                        if str(PDB2 + "-" + LGD2 + "_" + PDB1 + "-" + LGD1) in ltemp:
                            continue
                        else:
                            ltemp.append(ktemp)

                        filout.write(LGD1 + "-" + PDB1 + "\t" + LGD2 + "-" + PDB2 + "\t" + dlig[k1] + "\t" + dlig[k2] + "\t" + LSR1 + "\t" + LSR2 + "\t" + str(aff1) + "\t" + str(aff2) + "\t" + str(nbatomdiff) + "\t" + str(scoreMCS) + "\n")
    filout.close()
