import os

def countPosNeg(pdbARM,aa1, aa2, aa3, aa4, ION, Charge, List):
    from initial_Setup import counterion_ID, linker_aa_ID, second_counterion_ID
    from chromophore import chromophoreName

    global topCharge, bottomCharge
    topList = []
    topList2 = []
    bottomList = []
    bottomList2 = []
    topCharge = 0
    bottomCharge = 0

#09-05-2018 The count is based on the position of the LA    
#Position of the linker amino acid                                              
    with open(pdbARM) as file:
        for line in file:
            if (aa1 in line or aa2 in line or aa3 in line or aa4 in line and chromophoreName not in line) and "CA" in line: #The linker_aa and main_counterion should be conserved in cases where one of them is protonated 
                if signlinker_aa > 0:
                    resta = (float(line.split()[8]) - float(signlinker_aa))
                    if  (float(line.split()[8]) - float(signlinker_aa)) > 0:
                        topCharge = topCharge + 1
                        topList.append(line.split()[5])
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"TOP"+"\t"+str(float(resta)))
                    else:
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"BOTTOM"+"\t"+str(float(resta)))
                        bottomCharge = bottomCharge + 1
                        bottomList.append(line.split()[5])
                if signlinker_aa < 0:
                    resta = (float(line.split()[8]) - float(signlinker_aa))
                    if  (float(line.split()[8]) - float(signlinker_aa)) < 0:
                        topCharge = topCharge + 1
                        topList.append(line.split()[5])
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"TOP"+"\t"+str(float(resta)))
                    else:
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"BOTTOM"+"\t"+str(float(resta)))
                        bottomCharge = bottomCharge + 1
                        bottomList.append(line.split()[5])
#For CL- ions
    with open(pdbARM) as file:
        for line in file:
            if ION in line and chromophoreName not in line: #The linker_aa and main_counterion should be conserved in cases where one of them is protonated 
                if signlinker_aa > 0:
                    resta = (float(line.split()[8]) - float(signlinker_aa))
                    if  (float(line.split()[8]) - float(signlinker_aa)) > 0:
                        topCharge = topCharge + 1
                        topList.append(line.split()[5])
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"TOP"+"\t"+str(float(resta)))
                    else:
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"BOTTOM"+"\t"+str(float(resta)))
                        bottomCharge = bottomCharge + 1
                        bottomList.append(line.split()[5])
                if signlinker_aa < 0:
                    resta = (float(line.split()[8]) - float(signlinker_aa))
                    if  (float(line.split()[8]) - float(signlinker_aa)) < 0:
                        topCharge = topCharge + 1
                        topList.append(line.split()[5])
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"TOP"+"\t"+str(float(resta)))
                    else:
                        print( line.split()[3]+"\t"+line.split()[5]+"\t"+"BOTTOM"+"\t"+str(float(resta)))
                        bottomCharge = bottomCharge + 1
                        bottomList.append(line.split()[5])

#Remove the LA, MC and SC from the target residues lists

    ExcludedResList = [counterion_ID, linker_aa_ID, second_counterion_ID]

    for i in range(0,len(ExcludedResList)):
        if ExcludedResList[i] in bottomList:
            bottomList.remove(ExcludedResList[i])
    for i in range(0,len(ExcludedResList)):
        if ExcludedResList[i] in topList:
            topList.remove(ExcludedResList[i])
##
    globals().update({ "top"+Charge  : topCharge})
    globals().update({ "top"+List  : topList})
    globals().update({ "bottom"+Charge  : bottomCharge})
    globals().update({ "bottom"+List  : bottomList})

    print( "_"*29, "\n top"+Charge+"=", topCharge)
    print( "bottom"+Charge+"=", bottomCharge, "\n")

def IonTopBottom(TotalPosition, Position):
    if TotalPosition > 0:
        IonPosition = "CL"
    else:
        IonPosition = "NA"
    globals().update({ "Ion"+Position : IonPosition})


##################################################
#numberCounterions
##################################################

def numberCounterions(pdbARM):
    """
This subroutine is designed to compute the number of positively and negatively charged amino acids present in the inner and outer faces of the protein and, based on this information, calculate the total charge in each face.  To this aim, the protein structure -oriented along the Z axis- is “subdivided” in two regions called top and bottom, which are populated with the residues upper and downer the chromophore, respectively. The border between the top and bottom regions is defined by the z component of the coordinates of the center of mass of the chromophore (newcenterXYZ[2], calculated in the step 6 of the code). In other words, one amino acid belongs to the top region if the difference between its ð carbon z coordinate and newcenterXYZ[2] is positive, otherwise the amino acid belongs to the bottom region. Following this criteria, the positively charged (ARG, HIS, LYS) and negatively charged (ASP, GLU) amino acids are classified as top or bottom and are grouped into the categories: top positive, top negative, bottom positive and bottom negative. Finally, the total charge of the top and bottom regions is computed as the difference between the number of its positively and negatively charged residues.
    """
    from lib import chainName
    from initial_Setup import counterion_ID, linker_aa_ID, second_counterion_ID

    TotalTop = 0
    TotalBottom = 0
    print( "-"*29, "\n Summary \n", "-"*29)

    signlinker_aa = '' # Define the outer and inner side of the protein 
    with open(pdbARM) as file:
        for line in file:
            if "ATOM" in line and chainName+'{:>4}'.format(str(linker_aa_ID)) in line and "CA" in line:
                signlinker_aa = float(line.split()[8])
                globals().update({ "signlinker_aa"  : signlinker_aa})
                print(signlinker_aa) 



    countPosNeg(pdbARM,"ARG", "HIS", "LYS", "HIP", " NA ", "Positive", "PosList")
    countPosNeg(pdbARM,"ASP", "GLU", "ASP", "GLU", " CL ", "Negative", "NegList")

    TotalTop = topPositive + (topNegative*-1)
    TotalBottom = bottomPositive + (bottomNegative*-1)
    TotalPos = topPositive + bottomPositive
    TotalNeg = topNegative + bottomNegative
        
    globals().update({ "TotalTop"  : TotalTop})
    globals().update({ "TotalBottom"  : TotalBottom})
    globals().update({ "TotalPos"  : TotalPos})
    globals().update({ "TotalNeg"  : TotalNeg})

    print( "Number of positively charged residues: \t", TotalPos)
    print( "Number of negatively charged residues: \t", TotalNeg)
 
        
    IonTopBottom(TotalTop, "Top")
    IonTopBottom(TotalBottom, "Bottom")
        
    print( "-"*29, "\n", "Total charge top: \t", TotalTop, "\t", "|  Suggestion: Add to the top ", abs(TotalTop), IonTop)
    print( "-"*29, "\n", "Total charge bottom: \t", TotalBottom, "\t",  "|  Suggestion: Add to the bottom", abs(TotalBottom), IonBottom)
    print( "-"*29)
#        print "\n", warning, "The linker amino acid", linker_aa, "and the main counterion", main_counterion, "have been excluded from this analysis"

#Exclude the linker amino acid and the main counterion from the list of target residues

    


##################################################
#This function adds the counterions needed to neutralice the system
##################################################
def addCounterIons(pdbARM):
    """
    This step calculates a set of coordinates for the Cl- and Na+ counterions needed to neutralize the complex protein-chromophore. These coordinates are written in the pdbARM  file using the correct format.
    """
    from protonationStep import pH
    from external_Software import pdb2pqr, putIon

    print("FILE", pdbARM)
    global pH

    if pH:
        pass
    else:
        pH = PickNumber(14.0, '<-> Please write the pH-value (suggested value physiological pH 7.4) in the range ',  0, float)
        globals().update({ "pH"  : str(pH)})

#Uncomment the following two lines and comment the third one to ask the user for the force field to generate the pir file        
#    ForceFieldList = ["amber", "charmm", "parse"]
#    ChooseNumOption(ForceFieldList, "FF", "ForceField", '\n Choose the ', 'to generate the pqr file:', 'is selected as the Force Field to run the PDB2PQR software.', True)
        
    global pqrARM
    pqrARM = pdbARM[:-3]+"pqr"
    ForceFieldName = "AMBER"

    os.system(pdb2pqr+" --ff="+ForceFieldName+" --with-ph="+str(pH)+" --ph-calc-method=propka --summary "+str(pdbARM+" "+pqrARM))

########
    print( "The ","\x1b[0;33;49m"+pqrARM+"\x1b[0m", "file has been generated and is ready to be used in the putIon module.")

    PutIon(pdbARM)
    global ErrorION
    while ErrorION == True:
        PutIon(pdbARM)
#        
    writeIonCoord(pdbARM)

###########################
def PutIon(pdbARM):
    from external_Software import pdb2pqr, putIon
    from cavity_fpocket import cavityList
    from initial_Setup import counterion_ID, linker_aa_ID, second_counterion_ID

    global outIon
    global infoIonFile

    infoIonFile = "infoIon"+pdbARM[:-4]+".sh"

    print("topPosList", topPosList)
    print("topNegList", topNegList)
    print("bottomPosList", bottomPosList)
    print("bottomNegList", bottomNegList) 
# Remove the linker aa, the main counterion and the residues of the cavity from the target residues list

#inner
    with open(infoIonFile, "w") as infoIon:
        TargetInner = [] #This list contains the target residues 
        TargetOuter = []
        infoIon.write("#!/bin/bash \n \n")
        infoIon.write(putIon+" << EOF \n")
        infoIon.write(pqrARM+"\n")
        infoIon.write(str(TotalTop)+"\n")

        if TotalTop == 0:
            pass
        else:
            if TotalTop > 0:
                infoIon.write(str(len(topPosList))+"\n")
                for i in range(0,len(topPosList)):
                #                    infoIon.write(str(topPosList[i])+"\n")
                    infoIon.write(str(topPosList[i])+"\n")
                TargetInner = sorted(list(set(topPosList).difference(cavityList)))
#topPosList
            else:
                infoIon.write(str(len(topNegList))+"\n")
                for i in range(0,len(topNegList)):
                    infoIon.write(str(topNegList[i])+"\n")
                TargetInner = sorted(list(set(topNegList).difference(cavityList)))
#topNegList
#outer
        infoIon.write(str(TotalBottom)+"\n")
        if TotalBottom == 0:
            pass
        else:

            if TotalBottom > 0:
                infoIon.write(str(len(bottomPosList))+"\n")
                for i in range(0,len(bottomPosList)):
                    infoIon.write(str(bottomPosList[i])+"\n")
                TargetOuter = sorted(list(set(bottomPosList).difference(cavityList)))
#bottomPosList
            else:
                infoIon.write(str(len(bottomNegList))+"\n")
                for i in range(0,len(bottomNegList)):
                    infoIon.write(str(bottomNegList[i])+"\n")
                TargetOuter = sorted(list(set(bottomNegList).difference(cavityList)))

            globals().update({"TargetInner" : TargetInner})
            globals().update({"TargetOuter" : TargetOuter})

#bottomNegList
        infoIon.write("EOF")
        outIon = "outputPutIon."+pqrARM[:-4]
        os.system("chmod 777 "+infoIonFile)
            
        TargetInner_pymol = '+'.join(TargetInner)
        TargetOuter_pymol = '+'.join(TargetOuter)

        globals().update({"TargetInner_pymol" : TargetInner_pymol})
        globals().update({"TargetOuter_pymol" : TargetOuter_pymol})


#Execution of the PutIon module
    ExecutePutIon()
    identifyError(pdbARM)

#########################                                                                                                                                                                                                                     
def ExecutePutIon():
    os.system("sh "+infoIonFile+" > "+outIon )

#########################                                                                                                                                                                                                                     
def identifyError(pdbARM):
    global topPositive, topNegative, bottomPositive, bottomNegative, ErrorION
    with open(outIon, "r") as outputIon:
        for line in outputIon:
            if ("Error: Charge too small, exiting") in line:
                ErrorION = True
                print("We will jump...")
                outputIon.seek(0)
                resError = (outputIon.readlines()[-2]).split()[0]
                if resError in topPosList:
                    topPosList.remove(resError)
                    topPositive = int(len(topPosList))
                if resError in topNegList:
                    topNegList.remove(resError)
                    topNegative = int(len(topNegList))
                if resError in bottomPosList:
                    bottomPosList.remove(resError)
                    bottomPositive = int(len(bottomPosList))
                if resError in bottomNegList:
                    bottomNegList.remove(resError)
                    bottomNegative = int(len(bottomNegList))
            else:
                ErrorION = False
    print (("\n ---> Running PutIon analysis over the \x1b[0;33;49m"+str(pqrARM)+"\x1b[0m file for the placement of the counterions").rjust(100, '.'))
                
    globals().update({ "ErrorION"  : ErrorION})
    ExecutePutIon()

# ####################
def writeIonCoord(pdbARM):
    from pdb_format import LastResNum, lN
    from lib import chainName
    print( "---> The following coordinates for the ", IonTop, "and", IonBottom, " counterions have been added to the ", pdbARM, "file: \n")
    with open("ions.pdb", "r") as ions:
        contador = LastResNum+1
        line_num = lN+1
        atomName = []
        resName = []
        x_position = []
        y_position = []
        z_position = []
        for line in ions:
            if ("HETATM") in line:
                atomName.append(line.split()[2])
                resName.append(line.split()[3])
                x_position.append(float(line.split()[6]))
                y_position.append(float(line.split()[7]))
                z_position.append(float(line.split()[8]))

    with open(pdbARM, "a") as file2:
        for i in range(0,abs(len(atomName))): 
            cont = contador+i  # continue line number in pdbARM
            line_Num = line_num+i       # continue residue number in pdbARM
            pdbFormat = "%-6s%5d %s %-4s%3s %s%4d %3s%8.3f%8.3f%8.3f"
            IonCoordinates = pdbFormat % ("HETATM", line_Num, "", atomName[i], resName[i], chainName, cont, "", x_position[i], y_position[i], z_position[i])
            print( pdbFormat % ("HETATM", line_Num, "", atomName[i], resName[i], chainName, cont, "", x_position[i], y_position[i], z_position[i]))
            file2.writelines(IonCoordinates+"\n")
        file2.writelines("END")       
        

    
