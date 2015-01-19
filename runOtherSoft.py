'''
RunOtherSoft 10-03-2014
'''
import tool
import superposeStructure
import os
import subprocess
import time




TMalign = "/home/borrel/softwares/TMalign/TMalign"
shaep = "/home/borrel/softwares/shaep/shaep"


def subprocessTimeControl(command, time_out=60):
    """executing the command with a watchdog"""

    # launching the command
    c = subprocess.Popen(command, shell=True)

    # now waiting for the command to complete
    t = 0
    while t < time_out and c.poll() is None:
        time.sleep(1)  # (comment 1)
        t += 1

    # there are two possibilities for the while to have stopped:
    if c.poll() is None:
        # in the case the process did not complete, we kill it
        c.terminate()
        # and fill the return code with some error value
        returncode = -1  # (comment 2)

    else:                 
        # in the case the process completed normally
        returncode = c.poll()

    return returncode 


def runTMalign(path_pr1, path_pr2, path_dir_out, debug = 1) : 
    
    # if exist doesnt run again
    
    if not os.path.exists(path_dir_out + "align.out") and not  os.path.exists( path_dir_out + "matrix.out")  : 
        p_pr1 = tool.removeChain (path_pr1, path_dir_out)
        p_pr2 = tool.removeChain (path_pr2, path_dir_out)
        cmd_run = TMalign + " " + str (p_pr1) + " " + str (p_pr2) + " -o " + path_dir_out + "align.out -m " + path_dir_out + "matrix.out" +" > " + path_dir_out + "RMSD"
        if debug : 
            print cmd_run
        os.system(cmd_run)
    
    return [path_dir_out + "align.out", path_dir_out + "align.out_all", path_dir_out + "align.out_atm",path_dir_out + "align.out_all_atm", path_dir_out + "RMSD" ]


def babelConvertPDBtoSMILE (p_file_pdb) : 
    
    path_filout = p_file_pdb[0:-4] + ".smi"
    
    if not os.path.exists(path_filout) : 
        cmd_convert = "babel " + p_file_pdb + " " + path_filout
        os.system(cmd_convert)
        print cmd_convert
#         subprocessTimeControl(cmd_convert, time_out=10)
    
    try : filin = open (path_filout, "r")
    except : return "0"
    l_Fline = filin.readlines ()
    filin.close ()
    try : smile = l_Fline[0].split ("\t")[0]
    except : return "0"
    
    return smile



def runShaep (p_struct1, p_struct2, p_out, clean = 0):
    if clean == 1 : 
        if os.path.exists(p_out) : 
            os.remove(p_out)
        else : 
            pass
    elif os.path.exists(p_out) :
        return p_out
    
    
    # run
    cmd = shaep + " --output-file "  + p_out + " " + p_struct1 + " " + p_struct2  + " --noOptimization" 
    print cmd
    os.system (cmd)
#     subprocessTimeControl(cmd, time_out=30) 
            
    # supp others files
    cmd_rm = "rm " + p_out[0:-4] + "_hits.txt"
   
    try : os.system (cmd_rm)
    except : pass
    
    return p_out


def RhistogramMultiple (p_filin, brk = 20) : 
    
    cmd_run = "./histograms.R " + p_filin + " " + str (brk) 
    print cmd_run
    os.system (cmd_run)


def RhistogramRMSD (p_filin, brk = 20) : 
    
    cmd_run = "./histogramsRMSD.R " + p_filin + " " + str (brk) 
    print cmd_run
    os.system (cmd_run)



def barplot (p_filin) :
    cmd_run = "./barplotQuantity.R " + p_filin + " "
    print cmd_run
    os.system (cmd_run)

    
def Rhistogram (p_filin, name_main, brk = 100) : 
    
    cmd_run = "./histogram.R " + p_filin + " " + name_main + " " + str (brk) 
    print cmd_run
    os.system (cmd_run)
    
def water(path_file_fasta1, path_file_fasta2, path_filout, gapopen = 10, gapextend = 0.5, debug = 1):
    """
    Run water from emboss with 2 fasta files and gap open option and gap extend
    args: -> file fasta 1
          -> file fasta 2
          -> gap open value
          -> gap extend value
     return: -> water file
     """
    
    cmd = "water -asequence " + path_file_fasta1 + " -bsequence " + path_file_fasta2 + " -outfile " + path_filout + " -gapopen " + str(gapopen) + " -gapextend " + str(gapextend)
    
    if debug : 
        print cmd
    os.system (cmd)   
    return path_filout    



def needle (path_file_fasta1, path_file_fasta2, path_filout, gapopen = 10, gapextend = 0.5, debug = 1):
    """
    Run water from emboss with 2 fasta files and gap open option and gap extend
    args: -> file fasta 1
          -> file fasta 2
          -> gap open value
          -> gap extend value
     return: -> water file
     """
    
    cmd = "needle -asequence " + path_file_fasta1 + " -bsequence " + path_file_fasta2 + " -outfile " + path_filout + " -gapopen " + str(gapopen) + " -gapextend " + str(gapextend)
    
    if debug : 
        print cmd
    os.system (cmd)   
    return path_filout



def babelPDBtoMOL2 (path_file_pdb) : 
    
    path_filout = path_file_pdb[0:-4] + ".mol2"
    if not os.path.exists(path_filout) : 
        cmd_convert = "babel " + path_file_pdb + " "+ path_filout 
#         os.system(cmd_convert)
        subprocessTimeControl(cmd_convert, time_out=5)
    return path_filout


def piePlot (p_filin):
    
    cmd = "./piePlot.R " + p_filin
    print cmd
    os.system (cmd)
 
