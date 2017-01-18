import pathManage
import parseShaep
import runOtherSoft

from os import listdir, path
from re import search

def SheapScoreToClass(prclassif):
    nameREF = prclassif.split("/")[-1]

    prout = pathManage.result(nameREF + "_SheapClassif")
    pfilout = prout + "PiLSRType"
    filout = open(pfilout, "w")
    filout.write("ClassLSR\tESP\tshape\tName\n")

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
            lrefprot = listdir(prclassif + "/" + foldergroup + "/")
            for refprot in lrefprot:
                lprref.append(prclassif + "/" + foldergroup + "/" + refprot)

    for reffolder in lprref:
        #print reffolder
        classcycle = reffolder.split("/")[-3]
        if classcycle == "cycle":
            classif = classcycle + "-" + reffolder.split("/")[-2]
        else:
            classif = reffolder.split("/")[-2]

        # PDB reference
        PDBref = reffolder.split("/")[-1]
        PDBref = PDBref.split("_")[-1]

        lLSR = listdir(reffolder + "/LSR")
        lgdREF = ""
        for fileLSR in lLSR:
            if search("LSR_REF", fileLSR):
                lgdREF = fileLSR.split("_")[2]
                break

        if lgdREF == "":
            print "Error reference l.49 classifResult.py"

        for fileLSR in lLSR:
            if search(".pdb", fileLSR):
                lelemsplit = fileLSR.split("_")
                typeLSR = lelemsplit[1]
                if typeLSR == "REF":
                    continue
                lgd = lelemsplit[2]
                PDBLSR = lelemsplit[3]
                #print classif, PDBref, typeLSR, lgd, PDBLSR
                # file sheap in result folder
                psheap = pathManage.result() + lgdREF + "/" + PDBref + "/substituent_" + lgd + "_" + PDBLSR + "_" + typeLSR + ".hit"
                #print psheap
                if not path.exists(psheap):
                    continue
                dsheap = parseShaep.parseOutputShaep(psheap)
                filout.write(classif + "\t" + str(dsheap["ESP_similarity"]) + "\t" + str(dsheap["shape_similarity"]) + "\t" + lgd + "_" + PDBLSR + "_" + typeLSR + "\n")
    filout.close()

    # plot R to do

    runOtherSoft.plotClassifSheap(pfilout)
