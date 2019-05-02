import os
import numpy as np
import lib
import shutil
from lib import PickNumber, ChooseNumOption, yes_no, replaceLine
#from add_heavy_atoms import add_Heavy_atoms

##################################################                                                                                                                                          
# Luca's proposal to calculate the charge                                                                                                                                                  
##################################################                                                                                                                                           
def Charge(Name, NameNum, pKa, pH):
    """ This novel approach, proposed by Luca De Vico, consist of the calculation of the amino acid charges using as reference 
    a pH value provided by the user. """

    global resPkaAna
    global The_Charge

    The_Charge = 0.0
    if pKa == 99.99: # CYS residue bonded in di-sulfide bond, no charge                                                                  
        return The_Charge

    #Calculation of the charge based on the calculated pKa and the given pH                                                            
    exponent = np.power(10, (pKa - pH))
    The_Charge = exponent / (1.0 + exponent)

    if (Name == 'ASP') or (Name == 'GLU') or (Name == 'C-') or (Name == 'CYS') or (Name == 'TYR') or (Name == 'Oco') or (Name == 'OXT'):
        The_Charge = The_Charge-1

    The_Charge = round(The_Charge)
    resPkaAna= "%3s%4d%5d" % (Name, int(NameNum), The_Charge)
    return The_Charge, resPkaAna

##################################################                                                                                                                                          
# pH identification                                                                                                                                                  
##################################################                                                                                                                                            
def identify_pH(pdbName):
    """ This subroutine identifies the experimental crystallization pH in the Template PDB 
    """
    pH_Exp= 'NULL'
    with open(pdbName) as pdbFile:
        for line in pdbFile:
            if 'REMARK' in line and 'PH       ' in line:
                if line.split()[4] == "NULL":
                    pH_Exp = "NULL"
                    print("--> The experimental crystallization pH is not reported in the Template PDB.")
                else:
                    pH_Exp = float(line.split()[4])                                                                                                                                               
                globals().update({ "pH_Exp"  : str(pH_Exp)})                                                                         
                return
            else:
                pH_Exp= 'NULL'
                globals().update({ "pH_Exp"  : str(pH_Exp)})                                                                         
                continue

##################################################                                                                                                                                          
# propka analysis                                                                                                                                                  
##################################################                                                                                                                                            
def propKa(pdbName,pdbARM):
    """ This subroutine is designed to run the PROPKA3.1 software, using as input the pdbARM file.
    From the computed pKa, presented in the output file pdbARM.pka, 
    two different approaches are carried out to identify the residues which ionization state should be changed. """

    from Lists_dictionaries import protList, protAApKa
    from external_Software import pdb2pqr, propkaScript
    from lib import chainName
    from initial_Setup import linker_aa, main_counterion 

    global resProtNumID
    ForceFieldName = "amber"
    resPkaAnaList = []
    resProtNumIDList = []

    identify_pH(pdbName)

    #LMPG 09-07-2018
    print( '\x1b[0;33;31m'+"\n ---> The crystallographic structure was crystallized at pH:"'\x1b[0m', pH_Exp)
    pH = PickNumber(14.0, '<-> Please write the pH-value in the range ',  0, float)
    globals().update({ "pH"  : str(pH)})
    
    global pdbPropkaTemp, pdbPropka, pdbARMFix 
    pdbARMFix = pdbARM[:-3]+'fix.pdb'
    pdbPropka = pdbARMFix[:-3]+'pka'
    pdbPropkaTemp = pdbPropka+'.temp'

    globals().update({"pdbARMFix": pdbARMFix })                

    os.system(pdb2pqr+" --chain --ff="+ForceFieldName+" --with-ph="+str(pH)+"  --ph-calc-method=propka --summary "+str(pdbARM+" "+pdbARMFix))
    print( '\n', str('Running PROPKA3.0 analysis for the \x1b[0;33;49m'+str(pdbARM)+'\x1b[0m input file').rjust(100, '.'))
    
    os.system (propkaScript+" --pH="+str(pH)+" "+pdbARMFix+ ">> /dev/null")
    print( str("Done! The files \x1b[0;33;49m"+pdbARM[:-3]+"propka_input\x1b[0m and \x1b[0;33;49m"+pdbPropka+"\x1b[0m has been generated. ").rjust(100, '.'))

        # A temporal file with the summary of the propka analysis including the Buried values
    with open(pdbPropka, "r") as file, open(pdbPropkaTemp, "w") as file2:
        for line in file:
            for i in range(0, len(protList)):
                if protList[i] in line and chainName in line and "%" in line:
                    line2 = line.split()[0]
                    if protList[i] in line2 and chainName in line and "%" in line:
                        contentFile2 = ("%s \t %s \t %s \t  %s \t %s \t %s \n" % (line.split()[0],line.split()[1],line.split()[3],protAApKa.get(line.split()[0]),line.split()[4],line.split()[5])) 
#Coupled residues detected by propKa are marked with * symbol
                        file2.writelines(contentFile2.replace('*', ' '))
                        
    print( "\n ---> At pH ", pH, "the predicted charge of the residues is: \n", '_'*60 +'\n', '{:^9}'.format('RESIDUE')+'{:^6}'.format('CHARGE')+'{:^6}'.format(' pKa')+'{:^18}'.format(' (pKa - pKa-model)')+'{:^12}'.format('BURIED (%)')+ '\n'+ '_'*60)

        # Analysis based on Luca's proposal
    with open(pdbPropkaTemp) as pkaFile:
        for line in pkaFile:
            for i in range(0, len(protList)):
                if protList[i] in line.split()[0]:# and linker_aa_ID not in line.split()[1] and counterion_ID not in line.split()[1]:
                    shift = 1.5
                    perc_buried_min = 55 
                    resProtNumID = line.split()[1]
                    pKa_calc = line.split()[2]
                    pKa_model = line.split()[3]
                    perc_buried = line.split()[4]
                    diff_pKa = '%.2f' %(abs((float(pKa_calc) - float(pKa_model))))
                    Charge(line.split()[0], line.split()[1], float(pKa_calc), pH)
                    print( resPkaAna, '{:^11}'.format(pKa_calc) , '{:^10}'.format(diff_pKa), '{:^10}'.format(perc_buried))
                    if (protList[i] == "ASP" or protList[i] == "GLU") and The_Charge != -1 and int(perc_buried)>= perc_buried_min:# and float(diff_pKa) > 1:
                        resPkaAnaList.append(resPkaAna)
                        resProtNumIDList.append(resProtNumID)
                        
                    if (protList[i] == "ARG" or protList[i] == "LYS") and The_Charge != 1 and int(perc_buried)>= perc_buried_min:# and float(diff_pKa) > 1:
                        resPkaAnaList.append(resPkaAna)
                        resProtNumIDList.append(resProtNumID)

                    # "Traditional" analysis using a shift paramter and buried percentage > 60% (LMPG-31-05-2018)
                    if protList[i] == "HIS" and float(diff_pKa) > shift and int(perc_buried)>= perc_buried_min:
                        resPkaAnaList.append(resPkaAna)
                        resProtNumIDList.append(resProtNumID)
            
    globals().update({"resPkaAnaList": resPkaAnaList })
    globals().update({"resProtNumIDList": resProtNumIDList})
    ChooseNumOption(resPkaAnaList,"", "ionizable residue", '\n Based on the computed charges, the suggested residues to be protonated are:', '', '', False, "")     
    print("WARNING: The linker amino acid", linker_aa, "and the main counterion", main_counterion, "should be excluded from this analysis")
    print("WARNING: Check carefully if HIS residues must be protonated!" )
    

##################################################                                                                                           
#Protonation side chain amino acids
##################################################                                                                                          
def protonation(pdbARM):
    """ This step ask the user to type the Residue sequence number of the amino acids which ionization state should be changed.
    The label of these amino acids is modified in the pdbARM file (ASN --> ASH, GLU --> GLH, HIS --> HIE, LYS --> LYD).
    """

    from lib import chainName
    from Lists_dictionaries import protAA

    print("FILE:", pdbARM)
    protResDictionary = {}

    question = yes_no('\n <-> Based on the computed charges and the experimental information, do you want to change \n the ionization state of any amino acid? (ASP --> ASH, GLU --> GLH, HIS --> (HID or HIE or HIP), LYS --> LYD)')

    if question == True:
#Use the suggested protonation states
        question1 = yes_no('\n Do you want to protonated the suggested residues '+str(resProtNumIDList)+ '?')

        if question1== True:
            protRes = resProtNumIDList
        if question1== False:
#Use other protonation states different to the suggested ones            
            protRes = input('\n <-> Introduce the NUMBER ID of the residues, separated by space \" \" (i.e. 1 2 3): \t').split(' ')
#        if question1 == False:
        while(1):
            with open(pdbARM) as pdbTemp:
                for line in pdbTemp:
                    for i in range(0, len(protRes)):
                        if "ATOM" in line and chainName+'{:>4}'.format(protRes[i]) in line:
                            protResDictionary.update({ protRes[i] : line.split()[3] } )
   
            question2 = yes_no('---> The residues you selected are: '+str(protResDictionary)+'\n Are you sure about your selection?' )
            print( "\n")
            if question2 == True:

                globals().update({"protResDictionary": protResDictionary })                
                print( "*"*80, "\n", "WARNING: Remember that:", "\n HID: Histidine with hydrogen on the delta nitrogen", "\n HIE: Histidine with hydrogen on the epsilon nitrogen", "\n HIP: Histidine with hydrogens on both nitrogens; this is positively charged. \n", "*"*80, "\n")

                for key in protResDictionary:
                    if protResDictionary[key] == "HIS":
                        print( "\n <-> Please select the correct ionization state of the residue \x1b[0;33;49m"+protResDictionary[key]+" "+key+"\x1b[0m:")
                        HISList = ["HIE", "HID", "HIP"]
                        ChooseNumOption(HISList, "HIS", "HIS", '\n Choose the ', 'protonation state.', 'is the new protonation state.', True)
                        from lib import HISName
                        replaceLine(pdbARM, "ATOM", chainName+'{:>4}'.format(key), {"HIS" : HISName})
                    
                    else:
                        replaceLine(pdbARM, "ATOM", chainName+'{:>4}'.format(key), protAA)
                        print( "\n ---> The protonation state of the residue \x1b[0;33;49m"+protResDictionary[key]+" "+key+"\x1b[0m has been changed to \x1b[0;33;49m"+str(protAA[protResDictionary[key]])+"\x1b[0m \n")
                return
                
            if question2 == False:
                continue
    else:
        print( '---> No amino acid requires to change its ionization state')
