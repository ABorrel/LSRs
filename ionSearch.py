import parsePDB
import pathManage
import runOtherSoft

from os import listdir



# global
l_ions = ['CR', 'AL', 'MG', 'MN', 'ZN', 'CA', 'FE', 'CU', 'CD', 'NI', 'PD', 'CO', 'BA', 'HG', 'BE', 'LI', 'NA','CS', 'K']



def retrieveTwoAtomForAngle (lig_parsed, substruct):
    
    l_phosphate = []
    for atom_lig in lig_parsed : 
        if atom_lig["element"] == "P" : 
            if atom_lig["resName"] == substruct : 
                l_phosphate.append (atom_lig)
    
    
    if substruct == "ADP" or substruct == "POP" : 
        if len (l_phosphate) == 2 : 
            return l_phosphate
        else : 
            return []
    elif substruct == "AMP" : 
        if len (l_phosphate) > 1 : 
            print "ERROR ION l29 ionSearch.py"
        for atom_lig in lig_parsed : 
            if atom_lig["name"] == "O5'" : 
                l_phosphate.append (atom_lig)
    
    # case where ATP not complet
    if substruct == "ATP" and len (l_phosphate) != 3 : 
        return []
    return l_phosphate




def analyseIons (pr_dataset, name_ligand, p_filout, thresold_max_interaction = 4.0) : 

    l_folder_ref = listdir(pr_dataset)

    filout = open (p_filout, "w")
    if name_ligand == "ATP" : 
        filout.write ("PDB\tIon\tD1\tD2\tD3\tAngle1\tAngle2\tAt1\tAt2\tA3\n")
    else : 
        filout.write ("PDB\tIon\tD1\tD2\tAngle\tAt1\tAt2\n")
    
    # dictionnary of counting
    d_count = {}
    d_count["CX"] = 0
    d_count["CX + ions"] = 0
    d_count["BS + ions"] = 0
    d_count["BS + 1-ion"] = 0
    d_count["BS + 2-ions"] = 0
    d_count["BS + more-ions"] = 0
    d_count["Interact-1"] = 0
    d_count["Interact-2"] = 0

    
    # dictionnary by ions
    d_ions = {}
    
    for ref_folder in l_folder_ref  :
        only_one = 0
        if len (ref_folder) != 4 : 
            continue
        d_count["CX"] = d_count["CX"] + 1
        l_temp = []
        # path and complex
        p_lig_ref = pathManage.findligandRef(pr_dataset + ref_folder + "/", name_ligand)
        p_complex = pathManage.findPDBRef(pr_dataset + ref_folder + "/")
    
        # parsing
        lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref, "HETATM")
        l_het_parsed = parsePDB.loadCoordSectionPDB(p_complex, "HETATM")
    
        # retrieve phosphate
        l_pi = retrieveTwoAtomForAngle (lig_ref_parsed, name_ligand)
        if l_pi == [] : # case ligand without phosphate 
            continue 
        flag_interact = 0
        flag_between_1 = 0
        flag_between_2 = 0
        for het_parsed in l_het_parsed : 
            if het_parsed["resName"] in l_ions : 
                d_count["CX + ions"] = d_count["CX + ions"] + 1
                if not het_parsed ["resName"] in d_ions.keys () : 
                    d_ions[het_parsed["resName"]] = 0
                if not het_parsed["resName"] in l_temp :  
                    d_ions[het_parsed["resName"]] = d_ions[het_parsed["resName"]] + 1
                    l_temp.append (het_parsed["resName"])
                PDB_id = ref_folder
                d1 = parsePDB.distanceTwoatoms(l_pi[0], het_parsed)
                d2 = parsePDB.distanceTwoatoms(l_pi[1], het_parsed)
                if name_ligand == "ATP" : 
                    # print len(l_pi), ref_folder, p_lig_ref
                    d3 = parsePDB.distanceTwoatoms(l_pi[2], het_parsed)
                    angle_bis = parsePDB.angleVector(l_pi[1], het_parsed, l_pi[2])
                angle = parsePDB.angleVector(l_pi[0], het_parsed, l_pi[1])
            
                if d1 < 10 and d2 < 10 : 
                    if not het_parsed["resName"] in d_count.keys () : 
                        d_count[het_parsed["resName"]] = 0
                    if only_one == 0 : 
                        d_count[het_parsed["resName"]] = d_count[het_parsed["resName"]] + 1
                        only_one = 1
                    d_count["BS + ions"] = d_count["BS + ions"] + 1
                    flag_interact = flag_interact + 1
                    if d1 < thresold_max_interaction and d2 < thresold_max_interaction : 
                        flag_between_1 = flag_between_1 + 1


                    if name_ligand == "ATP" :
                        if d3 < thresold_max_interaction and d2 < thresold_max_interaction : 
                            flag_between_2 = flag_between_2 + 1
                        filout.write (str (PDB_id) + "\t" + str(het_parsed["resName"]) + "\t" + str(d1) + "\t" + str(d2) + "\t" + str (d3) + "\t" + str(angle) + "\t" + str(angle_bis) + "\t" + str(l_pi[0]["serial"]) + "\t" + str(l_pi[1]["serial"]) + "\t" + str(l_pi[2]["serial"]) + "\n")
                    else : 
                        filout.write (str (PDB_id) + "\t" + str(het_parsed["resName"]) + "\t" + str(d1) + "\t" + str(d2) + "\t" + str(angle) + "\t" + str(l_pi[0]["serial"]) + "\t" + str(l_pi[1]["serial"]) + "\n")
    
        if flag_interact == 1 : 
            d_count["BS + 1-ion"] = d_count["BS + 1-ion"] + 1
        elif flag_interact == 2 : 
            d_count["BS + 2-ions"] = d_count["BS + 2-ions"] + 1
        elif flag_interact > 2 : 
            d_count["BS + more-ions"] = d_count["BS + more-ions"] + 1

        if flag_between_1 >= 1 : 
            d_count["Interact-1"] = d_count["Interact-1"] + flag_between_1
        if flag_between_2 >= 1 : 
            d_count["Interact-2"] = d_count["Interact-2"] + flag_between_2 

    filout.close ()
    
    filout_count = open (p_filout[0:-4] + "count.txt", "w")
    filout_count.write ("CX: " + str (d_count["CX"]) + "\n")
    filout_count.write ("CX + ions: " + str (d_count["CX + ions"]) + "\n")
    filout_count.write ("BS + ions: " + str(d_count["BS + ions"]) + "\n")
    filout_count.write ("BS + 1-ion: " + str(d_count["BS + 1-ion"]) + "\n")
    filout_count.write ("BS + 2-ions: " + str(d_count["BS + 2-ions"]) + "\n")
    filout_count.write ("BS + more-ions: " + str(d_count["BS + more-ions"]) + "\n")
    filout_count.write ("Interact Pi-alpha + Pi-beta: " + str(d_count["Interact-1"]) + "\n")
    filout_count.write ("Interact Pi-beta + Pi-gama: " + str(d_count["Interact-2"]) + "\n")
    filout_count.close ()

    filout_by_ion = open(p_filout[0:-4] + "byIons_" + name_ligand, "w")
    l_k = d_ions.keys ()
    for k in l_k : 
        filout_by_ion.write (str (k.capitalize()) + "\t" + str (d_ions[k]) + "\n")
    filout_by_ion.close ()
   
    runOtherSoft.barplot (p_filout[0:-4] + "byIons_" + name_ligand)

