import parsePDB
import pathManage









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
    
    return l_phosphate




def analyseIons (p_data_ref, substruct, filout) : 
    
    # path out file
    
    
    # path and complex
    p_lig_ref = pathManage.findligandRef(p_data_ref, substruct)
    p_complex = pathManage.findPDBRef(p_data_ref)
    
    
    # parsing
    lig_ref_parsed = parsePDB.loadCoordSectionPDB(p_lig_ref, "HETATM")
    l_het_parsed = parsePDB.loadCoordSectionPDB(p_complex, "HETATM")
    
    # retrieve phosphate
    l_pi = retrieveTwoAtomForAngle (lig_ref_parsed, substruct)
    if l_pi == [] : # case ligand without phosphate 
        return 
    
    for het_parsed in l_het_parsed : 
        if het_parsed["resName"] in l_ions : 
            PDB_id = p_data_ref.split ("/")[-2][0:4]
            d1 = parsePDB.distanceTwoatoms(l_pi[0], het_parsed)
            d2 = parsePDB.distanceTwoatoms(l_pi[1], het_parsed)
            angle = parsePDB.angleVector(l_pi[0], het_parsed, l_pi[1])
            
            if d1 < 10 and d2 < 10 : 
                filout.write ( str (PDB_id) + "\t" + str(het_parsed["resName"]) + "\t" + str(d1) + "\t" + str(d2) + "\t" + str(angle) + "\t" + str(l_pi[0]["serial"]) + "\t" + str(l_pi[1]["serial"]) + "\n")
    
    




  