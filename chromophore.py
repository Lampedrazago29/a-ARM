import initial_Setup
import shutil
from initial_Setup import FormatPDB1


def LYR_to_RET():
    """                                                                                                                                                                                                                                       
    This subroutine converts the LYR chromophore to the standard LYS+RET format                                                                                                                                                               
    """

    from lib import chromophoreName, chainName
    from initial_Setup import pdbARM, pdbARMTemp
    from Lists_dictionaries import LYRtoRETDictionary, LYRtoLysDictionary

#    if chromophoreName == "RET":

#        chromophoreName = "RET"
#        globals().update({"chromophoreName" : chromophoreName})


#    if chromophoreName != "LYR":
    globals().update({"chromophoreName" : chromophoreName})

    if chromophoreName == "LYR":
#        FormatPDB1(pdbARM, pdbARMTemp ) #Temporal file with spaces                                                                                                                                                                            
        pdbARMTemp1 = pdbARMTemp+".1"

        with open(pdbARMTemp, "r") as pdbARMTempOld, open(pdbARMTemp1, "w") as pdbARMTempNew:

            for line in pdbARMTempOld:
                ls = line.split()

                if 'LYR' in line:
                    if ls[2] in LYRtoLysDictionary:
                    #Here, we change the labels according to LYS format                                                                                                                                                                       
                        newLine = ls[0]+"\t"+ls[1]+"\t"+ls[2]+"\t"+"LYS"+"\t"+ls[4]+"\t"+ls[5]+"\t"+ls[6]+"\t"+ls[7]+"\t"+ls[8]+"\t"+ls[9]+"\t"+ls[10]+"\t"+ls[11]

                        pdbARMTempNew.writelines(newLine+"\n")
#                        print(newLine)                                                                                                                                                                                                       

                    if ls[2] in LYRtoRETDictionary:
                    #Here, we change the labels according to RET format                                                                                                                                                                       

                        newLine = ls[0]+"\t"+ls[1]+"\t"+LYRtoRETDictionary.get(ls[2])+"\t"+"RET"+"\t"+ls[4]+"\t"+ls[5]+"\t"+ls[6]+"\t"+ls[7]+"\t"+ls[8]+"\t"+ls[9]+"\t"+ls[10]+"\t"+ls[11]

                        pdbARMTempNew.writelines(newLine+"\n")
 #                       print(newLine)                                                                                                                                                                                                       

                        chromophoreName = "RET"
                        globals().update({"chromophoreName" : chromophoreName})


                else:
                    pdbARMTempNew.writelines(line)
            shutil.move( pdbARMTemp1, pdbARMTemp)

