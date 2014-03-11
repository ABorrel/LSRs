'''
RunOtherSoft 10-03-2014
'''
import tool
import superposeStructure
import os





TMalign = "/home/borrel/softwares/TMalign/TMalign"





def runTMalign(path_pr1, path_pr2, path_dir_out, debug = 1) : 
    
    # delet chain in PDB
    path_pr1 = tool.removeChain (path_pr1)
    path_pr2 = tool.removeChain (path_pr2)
    
    path_pr1 = superposeStructure.manageTMalign (path_pr1)
    path_pr2 = superposeStructure.manageTMalign (path_pr2)
    
    cmd_run = TMalign + " " + str (path_pr1) + " " + str (path_pr2) + " -o " + path_dir_out + "align.out -m " + path_dir_out + "matrix.out" +" > " + path_dir_out + "RMSD"
    if debug : 
        print cmd_run
    os.system(cmd_run)
    
    return [path_dir_out + "align.out", path_dir_out + "align.out_all", path_dir_out + "align.out_atm",path_dir_out + "align.out_all_atm", path_dir_out + "RMSD" ]



