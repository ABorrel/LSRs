from re import search 
from os import listdir, path
from math import sqrt
from copy import deepcopy

import pathManage
import parseShaep
import parsePDB
import writePDBfile
import runOtherSoft


def selectSmileCode (p_file_smile, minimal_length_smile = 3):
    """
    selection: -> length 
               -> separated by .
    """
    p_filout = p_file_smile[0:-4] + ".filter"
    p_filout_by_ligand = p_file_smile[0:-4] + "_by_ligand.filter"
    filout = open (p_filout, "w")
    filout_by_ligand = open (p_filout_by_ligand, "w")    

    d_smile = {}
    
    filin = open (p_file_smile, "r")
    l_line_smile = filin.readlines ()
    filin.close ()
    
    for line_smile in l_line_smile : 
        
        l_element = line_smile.strip ().split ("\t")
        l_pdb = l_element [2].split (" ")
        l_ligand = l_element[4].split (" ")
        smile = l_element[0]
        
        # length 
        if len (smile) < minimal_length_smile : 
            continue
        elif search("\.", smile) : 
            l_smile = smile.split (".")
            for smile_separated in l_smile : 
                if len (smile_separated) < minimal_length_smile : 
                    continue
                elif not smile_separated in d_smile.keys () : 
                    d_smile[smile_separated] = {}
                    d_smile[smile_separated]["ligand"] = l_ligand
                    d_smile[smile_separated]["PDB"] = l_pdb
                else :
                    for lig in l_ligand : 
                        if not lig in d_smile[smile_separated]["ligand"] : 
                            d_smile[smile_separated]["ligand"].append (lig)
                    for pdb in l_pdb : 
                        if not pdb in d_smile[smile_separated]["PDB"] : 
                            d_smile[smile_separated]["PDB"].append (pdb)
        else : 
            if not smile in d_smile.keys () : 
                d_smile[smile] = {}
                d_smile[smile]["ligand"] = l_ligand
                d_smile[smile]["PDB"] = l_pdb
            else :
                for lig in l_ligand : 
                    if not lig in d_smile[smile]["ligand"] : 
                        d_smile[smile]["ligand"].append (lig)
                for pdb in l_pdb : 
                    if not pdb in d_smile[smile]["PDB"] : 
                        d_smile[smile]["PDB"].append (pdb)
        
    for smile_code in d_smile.keys () : 
        filout.write (str (smile_code) + "\t" + str (len (d_smile[smile_code]["PDB"])) + "\t" + " ".join (d_smile[smile_code]["PDB"]) + "\t" + " ".join(d_smile[smile_code]["ligand"]) + "\n")
            
    filout.close ()
    
    # classe by ligand
    d_l = {}
    for smile_code in d_smile.keys () : 
        #print smile_code
        #print d_smile[smile_code]["ligand"]
        for lig in d_smile[smile_code]["ligand"] : 
            if not lig in d_l.keys () : 
                d_l[lig] = [str (smile_code)]
            else : 	
                d_l[lig].append (str(smile_code))
    
    for ligand in d_l.keys () : 
        filout_by_ligand.write (ligand + "\t" + " ---  ".join(d_l[ligand]) + "\n")
    filout_by_ligand.close ()
  
    return p_filout
    

def globalShaepStat (substruct):
    
    pr_result = pathManage.result(substruct)
    
    p_filout = pr_result + "shaep_global.txt"
    filout = open (p_filout, "w")
    filout.write ("best_similarity\tshape_similarity\tESP_similarity\n")
    
    l_folder = listdir(pr_result)
    
    
    for ref_folder in l_folder  :
        if not path.isdir(pr_result + ref_folder + "/") : continue
        l_file_result = listdir(pr_result + ref_folder + "/")
        for file_result in l_file_result : 
            if search(".hit", file_result) :
                d_shaep_parsed = parseShaep.parseOutputShaep(pr_result + ref_folder + "/" + file_result) 
                if d_shaep_parsed != {} : 
                    filout.write (ref_folder + "_" + file_result[10:-4] + "\t" + str(d_shaep_parsed["best_similarity"]) + "\t" + str(d_shaep_parsed["shape_similarity"]) + "\t" + str(d_shaep_parsed["ESP_similarity"]) + "\n")
    filout.close ()
    runOtherSoft.RhistogramMultiple (p_filout, "Shaep_score")
                
        
        
def computeRMSDBS (p_ref, p_query, p_substruct, pr_result, thresold_BS = 6) :
    
    
    l_atom_query_parsed = parsePDB.loadCoordSectionPDB(p_query, "ATOM")
    l_atom_ref_parsed = parsePDB.loadCoordSectionPDB(p_ref, "ATOM")
    
    l_atom_substruct = parsePDB.loadCoordSectionPDB(p_substruct)
    
    
    
    l_BS_ref = []
    
    for atom_substruct in l_atom_substruct : 
        for atom_ref in l_atom_ref_parsed : 
            d_atom = parsePDB.distanceTwoatoms(atom_substruct, atom_ref)
            if d_atom <= thresold_BS : 
                l_BS_ref.append (atom_ref)
    # retrieve residue full
    l_BS_ref = parsePDB.getResidues(l_BS_ref, l_atom_ref_parsed)
    
#     print len (l_BS_ref)
#     print len (l_atom_query_parsed)
    
    
    l_BS_query = []
    flag_identic_crystal = 1
    for atomBS_ref in l_BS_ref :
#         print  atomBS_parsed 
        d_max = 100.0 
        for atom_query in l_atom_query_parsed :
            if atom_query["resName"] ==  atomBS_ref["resName"] and atom_query["name"] ==  atomBS_ref["name"] : 
                d = parsePDB.distanceTwoatoms(atom_query, atomBS_ref)
                if d < d_max : 
                    d_max = d
                    res_temp = atom_query
                
        
        if d_max < thresold_BS : 
            l_BS_query.append (deepcopy(res_temp))
            # identic check number
            if res_temp["resSeq"] != atomBS_ref["resSeq"] : 
                flag_identic_crystal = 0
        else : 
            # case structure not found
            return []
    
    
#     print len (l_BS_query), len (l_BS_ref)
    l_RMSD = RMSDTwoList (l_BS_query, l_BS_ref)
    
    # write PDB
    p_filout_pdb = pr_result + p_query.split ("/")[-1][0:-4] + "_" + str (flag_identic_crystal) + "_" + p_substruct.split ("_")[-2] + "_" + p_ref.split ("/")[-1]
    filout_pdb = open (p_filout_pdb, "w")
    writePDBfile.coordinateSection(filout_pdb, l_BS_ref, recorder = "ATOM")
    writePDBfile.coordinateSection(filout_pdb, l_BS_query, recorder = "ATOM", header = 0 )
    filout_pdb.close ()
    
    if l_RMSD == [] : 
        return []
    else : 
        return l_RMSD + [flag_identic_crystal]
    
    
    
def RMSDTwoList (l_atom1, l_atom2) : 
    
    nb_ca = 0.0
    d_max = {"value": 0.0}
    diff_position_all = 0.0
    diff_position_ca = 0.0
    
    if len (l_atom1) != len (l_atom2) or len (l_atom2) == 0 : 
        print "ERROR - RMSD: list length different or null"
        return []
    else : 
        i = 0
        while i < len (l_atom1): 
            if l_atom1[i]["name"] != l_atom2[i]["name"] and l_atom1[i]["resName"] != l_atom2[i]["resName"]: 
                print l_atom1[i]["name"] , l_atom2[i]["name"]
                print "ERROR"
                return []
            else : 
                d_atom = parsePDB.distanceTwoatoms(l_atom1[i], l_atom2[i])
                diff_position_all = diff_position_all + d_atom
                
                if l_atom1[i]["name"] == "CA" : 
                    diff_position_ca = diff_position_ca + d_atom
                    nb_ca = nb_ca + 1
                
                if d_atom > d_max["value"] : 
                    d_max["value"] = d_atom
                    d_max["atom"] = l_atom1[i]["name"] + "-" +  l_atom2[i]["name"] + "_" + l_atom1[i]["resName"] + "-" +  l_atom2[i]["resName"]
                    
            i = i + 1
#     print d_max
    return [sqrt(diff_position_all / len (l_atom1)), sqrt (diff_position_ca / nb_ca), d_max["value"], len (l_atom1)]
 
 
 
def familyPDBRef (d_dataset, p_filout) : 
    
    print d_dataset
    filout = open (p_filout, "w")
    filout.write ("PDBID\tUniprot\tName protein\tkwords\n")
    
    for PDB in d_dataset.keys () : 
        p_pdb = d_dataset[PDB]['p_pdb']
        name_prot = parsePDB.nameProtein(p_pdb)
        uniprot_id = parsePDB.UniProtID(p_pdb)
        kwords = parsePDB.keywords(p_pdb)
        filout.write (PDB + "\t" + uniprot_id + "\t" + name_prot + "\t" + kwords + "\n")
        
    filout.close ()
    
    
        
def findFamily (PDB_ID, p_file_family):
    
    filin = open (p_file_family, "r")
    l_line_flin = filin.readlines ()
    filin.close ()
    
    for line_filin in l_line_flin [1:]: 
        name_pr = line_filin.strip ().split ("\t")[-2]
        kwords = line_filin.strip ().split ("\t")[-1]
        PDB = line_filin.strip ().split ("\t")[0]
        if PDB == PDB_ID : 
            # kinase
            if search("KINASE", name_pr.upper ()) or search("KINASE", kwords.upper ()) : 
                return "kinase"
            
            elif search("HEAT SHOCK PROTEIN", name_pr.upper ()) or search("HEAT SHOCK PROTEIN", kwords.upper ()) or search("CHAPERONE", name_pr.upper ()) or search("CHAPERONE", kwords.upper ()) or         search("CHAPERONIN", name_pr.upper ()) or search("CHAPERONIN", kwords.upper ()) or search("HSP", name_pr.upper ()) or search("HSP", kwords.upper ()) : 
                return "HSP"
            else : 
                return "others"
    
    return "others"
    
    

            
