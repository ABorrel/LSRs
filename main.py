import pathManage
import parsePDB
import runOtherSoft
from os import listdir, path
import writePDBfile



# step 1
# first step -> preparation
# - extract the ligand


def datasetPreparation (substruct):
    
    p_dir_dataset = pathManage.dataset(substruct)
    
    l_folder = listdir(p_dir_dataset)
    
    
    for ref_folder in l_folder  :
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        
        for pdbfile in l_pdbfile : 
            p_file_pdb = p_dir_dataset + ref_folder + "/" + pdbfile
            
            # extract ligand in PDB
            l_ligand = parsePDB.retrieveListLigand(p_file_pdb)
            if l_ligand == [] : 
                continue
            else : 
                l_atom_pdb_parsed = parsePDB.loadCoordSectionPDB(p_file_pdb)
                for name_ligand in l_ligand : 
                    l_lig_parsed = parsePDB.retrieveLigand(l_atom_pdb_parsed, name_ligand)
                    p_filout_ligand = p_dir_dataset + ref_folder + "/" + name_ligand + "_" + path.split(p_file_pdb)[1]
                    writePDBfile.coordinateSection(p_filout_ligand , l_lig_parsed[0], "HETATM", name_ligand + "_" + p_file_pdb , connect_matrix = 1)
    return 1


# step 2
# second step -> surperimposed all protein on reference
# - tmaling
# - generated new folder alignement with the TMalign ouput


def applyTMAlign (substruct):
    
    p_dir_dataset = pathManage.dataset(substruct)
    l_folder = listdir(p_dir_dataset)
    
    
    for ref_folder in l_folder  :
        l_pdbfile = listdir(p_dir_dataset + ref_folder + "/")
        p_pdb_ref = pathManage.findPDBRef(p_dir_dataset + ref_folder + "/")
        
        
        for pdbfile in l_pdbfile : 
            # try if PDB not ligand
            if len(pdbfile.split ("_")[0]) != 4 : 
                continue
            else : 
                p_file_pdb = p_dir_dataset + ref_folder + "/" + pdbfile
                p_dir_align = pathManage.alignmentOutput(p_pdb_ref.split ("/")[-1][:-4] + "__" + p_file_pdb.split ("/")[-1][:-4])
                
                # superimpose 
                runOtherSoft.runTMalign(p_file_pdb, p_pdb_ref, p_dir_align)
    return 1


# step 3
# Third step -> generate the SMART code
# - apply rotated matrix on the ligand
# - retrieve substruct at 4a
# - convert in smart
# - generated pdb files for each reference the ligand superimposition
# - generated a list of SMART by global phosphate


def retrieveSubstructSuperimposed ():
    
    
    
    
    
    
    
    return 0


# step 4
# shaep run for every reference
# - input 1 Phosphate with oxygen
# - input 2 substructure 

def applyShaep ():
    
    return 0

# step 5
# manage results  
# - table result
# - folder tree

def manageResult ():
    
    return 0


 




#################
# RUN MAIN !!!! #
#################

# datasetPreparation ("AMP")
applyTMAlign ("AMP")











