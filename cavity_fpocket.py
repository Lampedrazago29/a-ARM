import shutil
import os
import subprocess
from lib import yes_no

##################################################
#Fpocket cavity generation
##################################################

def fpocket(pdbARMFix):
    """ First, the program receives as input the pdbARMFix file and automatically executes the Fpocket software 
    using the options by default. As output, a new folder with the files of the $N$ protein pockets in PDB format
    is obtained, along with a file which summarizes the information of all the pockets, including the Score and 
    some characteristics of each pocket. Then, N pocket lists (pl) with the residues of each of the N pockets 
    are generated, sorted from higher to lower score and saved as internal storage. To identify the correct 
    chromophore cavity, the different lists are explored and those which contains the PSB residue are pre-selected.
    The pocket list which contains the PSB and has the highest score is selected and named as cavity list. Finally, 
    the main and second counter-ion residues are added to the cavity list and this new list is printed in the 
    cavity file using the ARM format.
    """
    from initial_Setup import pdbARM
    from lib import chainName
    from initial_Setup import counterion_ID, linker_aa_ID, second_counterion_ID

    os.system('fpocket -f '+pdbARMFix)
    shutil.move(pdbARMFix[:-4]+"_out", pdbARM[:-4]+"_out")

    print( "\n The folder", pdbARM[:-4]+"_out", "containing the pockets has been generated")

#LMPG 29-05-2018
#This function identifies the pocket which contains the linker AA

    pocketFile= subprocess.getoutput("grep -lr "+"\""+chainName+'{:>4}'.format(str(linker_aa_ID))+"\""+" "+pdbARM[:-4]+"_out/pockets/*pdb")
# In fpocket the pockets are numbered from highest to lowest score. Here, the pocket of highest score which contains the linker aa is selected.
    pocketFile= pocketFile.split("\n")[0] 

    with open(pocketFile) as cavityFile:
        listcavity = [counterion_ID, linker_aa_ID, second_counterion_ID]
        for line in cavityFile:
            if "ATOM" in line.split()[0]:
                listcavity.append(line.split()[5])                    
                cavityList0 = list(set(listcavity))
                cavityList = sorted(cavityList0, key = int)
                globals().update({"cavityList" : cavityList})

    with open("cavity", "w") as cavity:
        for i in range(0, len(cavityList)):
            cavity.writelines(cavityList[i]+"\n")
            
            # print in screen the cavity
    print( '---> The standard Fpocket cavity file has been generated. The residues are: \n', cavityList)

    question_add = yes_no('\n Do you want to add residues to the standard cavity?')
#    while(1):
    if question_add == True:
        while(1):
            additional_residues = input('\n <-> Introduce the NUMBER ID of the residues, separated by space: \t').split(' ')
            for i in range(0,len(additional_residues)):
                cavityList.append(additional_residues[i])
                cavityListMod = sorted(list(set(cavityList)), key= int)
                globals().update({"cavityList" : cavityListMod})
                   
            with open("cavity", "w") as cavity:
                for i in range(0, len(cavityListMod)):
                    cavity.writelines(cavityListMod[i]+"\n")

            print( '---> The modified cavity file has been generated. The residues are: \n', cavityListMod)
            return

    if question_add == False:
        pass
#    else:
#        continue
                

    
