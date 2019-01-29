import lib
from lib import ChooseNumOption, yes_no, replaceLine, deleteLine
import shutil
import Lists_dictionaries
import initial_Setup

def searchRotamers(newFile = "TempFile"):
    """ This function identify the residues with different rotamers and ask the user for select the correct one """   

    from Lists_dictionaries import rmList
    from initial_Setup import pdbARMTemp

    occupancyList1 = []
    occupancyList = []
    with open(pdbARMTemp, "r") as pdbTemp:
        for line in pdbTemp:
            if 'ATOM' in line and float(line.split()[9]) < 1.0 and line.split()[3] not in rmList:
                occupancyList1.append(line.split()[3]+" "+line.split()[5]+" Ocuppancy= "+line.split()[9])
                occupancyList = sorted(list(set(occupancyList1)))
        ChooseNumOption(occupancyList,"", ":", '\n The following residues have multiple conformations', '', '', False, "")
        pdbTemp.seek(0)

        for line in pdbTemp:
            if 'ATOM' in line and ('CA' in line  or 'C1 ' in line) and float(line.split()[9]) < 1.0 and line.split()[3] not in rmList and line.split()[3] not in aaList:
                conformerName = line.split()[3]
                conformerResNum = line.split()[5]
                question= yes_no('\n <-> Do you want to keep the residue '+conformerName+' '+conformerResNum+"? ("+"Ocuppancy= "+line.split()[9]+")")
                if question == True:
                    replaceLine(pdbARMTemp, conformerName, conformerResNum, { conformerName : " "+str(conformerName)[+-3:]}, newFile, False)
                    shutil.move(newFile, pdbARMTemp)
                    print("---> The "+"\x1b[0;33;49m"+str(conformerName+" "+conformerResNum)+'\x1b[0m'+" will be KEPT.")
                else:
                    deleteLine(pdbARMTemp, conformerName, conformerResNum, { conformerName : " "+str(conformerName)[+-3:]}, newFile)
                    print("---> The "+"\x1b[0;33;49m"+str(conformerName+" "+conformerResNum)+'\x1b[0m'+" will be REMOVED.")
