'''
RunOtherSoft 10-03-2014
'''
import tool
import superposeStructure
import os





TMalign = "/home/borrel/softwares/TMalign/TMalign"
shaep = "/home/borrel/softwares/shaep/shaep"




def runTMalign(path_pr1, path_pr2, path_dir_out, debug = 1) : 
    
    # delet chain in PDB + remove ligand
    p_pr1 = tool.removeChain (path_pr1, path_dir_out)
    p_pr2 = tool.removeChain (path_pr2, path_dir_out)
    
    cmd_run = TMalign + " " + str (p_pr1) + " " + str (p_pr2) + " -o " + path_dir_out + "align.out -m " + path_dir_out + "matrix.out" +" > " + path_dir_out + "RMSD"
    if debug : 
        print cmd_run
    os.system(cmd_run)
    
    return [path_dir_out + "align.out", path_dir_out + "align.out_all", path_dir_out + "align.out_atm",path_dir_out + "align.out_all_atm", path_dir_out + "RMSD" ]


def  babelConvertPDBtoSMILE (p_file_pdb) : 
    
    path_filout = p_file_pdb[0:-4] + ".smi"
    
    if not os.path.exists(path_filout) : 
        cmd_convert = "babel " + p_file_pdb + " " + path_filout
        print cmd_convert
        os.system (cmd_convert)
    
    filin = open (path_filout, "r")
    l_Fline = filin.readlines ()
    filin.close ()
    smile = l_Fline[0].split ("\t")[0]
    
    return smile



def runShaep (p_struct1, p_struct2, p_out):
    
    cmd = shaep + " --output-file "  + p_out + " " + p_struct1 + " " + p_struct2 
    
    print cmd
    os.system (cmd)
    
    # supp others files
    cmd_rm = "rm " + p_out[0:-4] + "_hits.txt"
    
    os.system (cmd_rm)
    
    
    return p_out
    
    
    
    
