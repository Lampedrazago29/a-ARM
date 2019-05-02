import shutil
import os
import re
import textwrap

from lib import yes_no, FormatPDB1
from pdb_format import *
from protonationStep import *
from cavity_fpocket import fpocket
from counterions_Placement import numberCounterions, addCounterIons
from pymol_Script import pyMoLFig

from Lists_dictionaries import threeToOneDic, oneToThreeDic
globals().update({"threeToOneDic" : threeToOneDic})
globals().update({"oneToThreeDic" : oneToThreeDic})

#############################                                                                                                                                                                                    
#
#############################                                                                                                                                                                                  
def Select_mutation_Software(pdbARM):
    """
    This step ask the user for select the software for the mutations                                                                                                                                               
    """
    from lib import ChooseNumOption
    
    mutationSoftwareList = ["Scwrl4", "MODELLER"]
    ChooseNumOption(mutationSoftwareList, "mutation_Software", "mutation_Software", '\n Choose the ', 'to perform the mutations:', 'will be used to perform the mutations.', True)

#############################                                                                                                                                                                                    
#
#############################                                                                                                                                                                                  
def Mutations_ID(pdbARM):
    """
    This subroutine searches the .seqmut file and make a list with the mutations
    """
    from lib import SearchFileType
    
    SearchFileType('seqmut', '\n The following', 'files with a list of mutations are found: \n', 'file will be used for the mutations..')  #seqmutName                                                          
    
    from lib import seqmutName

    mutation_seqmutDic={}
    with open(seqmutName, "r") as file:
        content = file.read().splitlines()
        for line in content:
            if "Mutation" in line:
                key=line
                mutation_seqmutDic[key]=[]
            else:
                mutation_seqmutDic[key].append(line)
        globals().update({"mutation_seqmutDic" : mutation_seqmutDic})

    global wt_pdbARM
    wt_pdbARM = "wt_"+pdbARM
    shutil.copyfile(pdbARM, wt_pdbARM)
    globals().update({"wt_pdbARM" : wt_pdbARM})

#############################                                                                                                                                                                                    
#
#############################                                                                                                                                                                                  
def Insert_mutations(pdbARM,mutation_seqmut):
    """
    This step ask the user for the list of mutations. Similar to the seqmut file                                                                                                                              
    """

    from lib import workingFolder

    for i in range(0, len(mutation_seqmut)):
        mutation_seqmut[i] = mutation_seqmut[i].upper()

    globals().update({"mutation_seqmut" : mutation_seqmut})
    mutationsFormat(pdbARM,mutation_seqmut)

#Create a new working folder and a new pdb file for the mutation                                                                                                                                                                                                                
    mutFile=''
    for i in range(0,len(mutation_seqmut)):
        mutFile = mutFile+mutation_seqmut[i]+"-"
    globals().update({ "mutFile" : mutFile})

    mut_Folder = mutFile+workingFolder
    os.system("mkdir "+mut_Folder)
    os.system("cp "+pdbARM+" "+mut_Folder)
    os.chdir(mut_Folder)

    global mut_pdbARM, mut_pdbARMTemp, mutation_output
    mut_pdbARM = mutFile+pdbARM
    mut_pdbARMTemp = mut_pdbARM[:-3]+"temp"
    mutation_output = mut_pdbARM[:-3]+"output"
    shutil.copyfile(pdbARM, mut_pdbARM) #working mutation File                                                                                                                                                                                                                       
    NumIDmutations = []
    with open("seqmut"+mut_pdbARM[:-8], "w") as seqmutFile:
        for i in range(0,len(mutation_seqmut)):
            seqmutFile.writelines(mutation_seqmut[i]+"\n")
            NumID = re.findall('\d+', str(mutation_seqmut[i]))[0]
#List with the ResID numbers of the residues to be mutated is stored as NumIDmutationsList                                                                                                                                                                                      
            NumIDmutations.append(NumID)
    globals().update({"NumIDmutationsList" : NumIDmutations})
#    globals().update({"NumIDmutationsList" : })

    print( "\n ---> The following mutation(s) will be performed: ")
    for i in range(0,len(mutation_seqmut)):
        print( str(i+1)+") ",  mutation_seqmut[i])

#############################                                                                                                                                                                                       
#This step recognices non-standard residues in the mutations and unifies the format to 1 letter amino acid                                                                                                          
#############################                                                                                                                                                                                       
def mutationsFormat(pdbARM,mutation_seqmut):
    """#Recognizes non-standard residues and ask the user to insert the mutation again                                                                                                        
    #Unifies the format to 1 letter amino acid                                                                                                                   
    """
    for i in range(0,len(mutation_seqmut)):
        if (mutation_seqmut[i][0:3]) in threeToOneDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][0:3],threeToOneDic[mutation_seqmut[i][0:3]] )
        if (mutation_seqmut[i][-3:]) in threeToOneDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][-3:],threeToOneDic[mutation_seqmut[i][-3:]] )

            globals().update({"mutation_seqmut" : mutation_seqmut})

    #Unifies the format to 3 letter amino acid                                                                                                                                                                      
    for i in range(0,len(mutation_seqmut)):
        if (mutation_seqmut[i][0]) in oneToThreeDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][0],oneToThreeDic[mutation_seqmut[i][0]] )
        if (mutation_seqmut[i][-1:]) in oneToThreeDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][-1:],oneToThreeDic[mutation_seqmut[i][-1:]] )

            globals().update({"mutation_seqmut" : mutation_seqmut})


def NonStandarddResID(pdbARM,mutation_seqmut):
    #Recognizes non-standard residues 3 letter format                                                                                                                                                               
    for i in range(0,len(mutation_seqmut)):
        if (mutation_seqmut[i][0:3]).isalpha() == True:
            if (mutation_seqmut[i][0:3]) not in threeToOneDic:
                print( "\n",  "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][0:3])+'\x1b[0m', "\n Try again!")
                mutation_seqmut[i] = input('Insert correctly the mutation:')

        else:
            if (mutation_seqmut[i][0]).isalpha() == True:
                if (mutation_seqmut[i][0]) not in oneToThreeDic:
                    print( "\n",  "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][0])+'\x1b[0m', "\n Try again!") 
                    mutation_seqmut[i] = input('Insert correctly the mutation:')

        if (mutation_seqmut[i][-3:]).isalpha() == True:
            if (mutation_seqmut[i][-3:]) not in threeToOneDic:
                print( "\n",  "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][0:3])+'\x1b[0m', "\n Try again!")
                mutation_seqmut[i] = input('Insert correctly the mutation:')

        else:
             if (mutation_seqmut[i][-1:]).isalpha() == True:
                 if (mutation_seqmut[i][-1:]) not in oneToThreeDic:
                     print( "\n",  "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][-1:])+'\x1b[0m', "\n Try again!")
                     mutation_seqmut[i] = input('Insert correctly the mutation:')

    globals().update({"mutation_seqmut" : mutation_seqmut})

#############################                                                                                                                                                                                                                                                                                                 
#This step generates the getpir.py script. This script is then executed to obtain the .pir file.                                                                                                                                                                                                                              
#############################                                                                                                                                                                                                                                                                                                 
def get_pir_Script(pdbFile, resNumID, chainName):

    from external_Software import modellerScript

    global sequenceWT
    pdbFileTemp = pdbFile[:-3]+"temp"
    pirFile = pdbFile[:-3]+"pir"
    pirFileTemp = pirFile+".temp"

    getpir = ["import os \n",
              "import sys \n \n",
              "from modeller import * \n \n",
              "env = environ() \n",
              "aln = alignment(env) \n"
              "mdl = model(env, file='"+pdbFile+"') \n",
              "aln.append_model (mdl, align_codes='wt_"+pdbFile+"',atom_files='"+pdbFile+"') \n",
              "aln.write(file='"+pirFile+"', alignment_format='PIR') \n" ]

    getpirScript = "getpyr.py"
    with open(getpirScript, "w") as getpirFile:
        getpirFile.writelines(getpir)
    os.system(modellerScript +" "+getpirScript  )

#Identify missing residues and write the pir file in correct format                                                                                                                                                                                                                                                           
    FormatPDB1(pdbFile, pdbFileTemp)
    realResNumDic = {}
    with open(pdbFileTemp, "r") as file:
        for line in file:
            if "ATOM" in line:
                realResNumDic.update({int(line.split()[5]) : threeToOneDic[str(line.split()[3])]})
        globals().update({ "realResNumDic": realResNumDic})

    sequence = ''
    sequenceWT = ''
    sequenceWTList = []
    with open(pirFile, "r") as pir, open(pirFileTemp, "w") as temp:
        for line in pir:
            if str(pdbFile) in line:
                temp.writelines(line)
        missResList = []
        for i in range(0,resNumID):
            i = i+1
            if i not in realResNumDic:
                missResList.append(i)
                sequence = sequence+"-"
                sequenceWTList.append("-")
            else:
                res = str(realResNumDic[i]).lower()
                sequence = sequence+res
                sequenceWT = sequenceWT+res
                sequenceWTList.append(res)
        temp.writelines('\n'.join(textwrap.wrap(sequence, width=75)))
        temp.writelines('* \n')
    globals().update({ "sequenceWTList": sequenceWTList})
    globals().update({ "missResList": missResList})
    print( "missResList", missResList)

    shutil.move(pirFileTemp, pirFile)

    print( "\n The file "+pirFile+" has been generated using the MODELLER 9.19 software. The missing residues "+str(missResList)+" were considered.")

    os.remove(pdbFileTemp)

#############################                                                                                                                                                                                                               
#SCWRL4 routine                                                                                                                                                                                                                             
#############################                                                                                                                                                                                                                 
def Scwrl4_mutations(mut_pdbARM):
    
    from lib import chainName
    from pdb_format import resNumID
    from download_PDB import TitlePDB
    from external_Software import scrl4Script

    pdbHETATM = mut_pdbARM[:-3]+"HETATM.pdb"
    with open(mut_pdbARM, "r") as pdb, open(pdbHETATM, "w") as hetatm:
        for line in pdb:
            if "HETATM" in line:
                hetatm.writelines(line)

    global sequenceWTList
    seqFileName = mut_pdbARM[:-7]+"seqFileName"

    get_pir_Script(mut_pdbARM, resNumID, chainName)

    for i in range(0,len(mutation_seqmut)):
        NumIDmut = NumIDmutationsList[i]
        globals().update({"NumIDmut" : NumIDmut})
        sequenceWTList[int(NumIDmut)-1] = threeToOneDic[mutation_seqmut[i][-3:]]

        print("seqFileName is:", seqFileName)
        os.system("pwd")
        with open(seqFileName, "w") as seqFile:
            for j in range(0, len(sequenceWTList)):
                if sequenceWTList[j] != "-":
                    seqFile.writelines(sequenceWTList[j])

        print( str("Running SCWRL4 for the mutation number " + str(i+1)).rjust(100, '.'))
        os.system(scrl4Script+" -i "+mut_pdbARM+" -o "+mutation_output+" -h -s "+seqFileName+" > scwrl4_mut."+str(mutation_seqmut[i])+".log" )                                                                                               
        print(scrl4Script+" -i "+mut_pdbARM+" -o "+mutation_output+" -h -s "+seqFileName+" > scwrl4_mut."+str(mutation_seqmut[i])+".log" )                                                                                               

#################
#NO BORRAR
#        os.system(scrl4Script+" -i "+mut_pdbARM+" -o "+mutation_output+" -s "+seqFileName+" -h -t -f "+pdbHETATM+" > scwrl4_mut."+str(mutation_seqmut[i])+".log" )                                                                                               

#        print( (scrl4Script+" -i "+mut_pdbARM+" -o "+mutation_output+" -h -t -f "+pdbHETATM+" -s "+seqFileName+" > scwrl4_mut."+str(mutation_seqmut[i])+".log" ))
##################


        print( "\n The mutation ", mutation_seqmut[i], "has been succesfully generated!")

        sequenceWTList[int(NumIDmut)-1] = sequenceWTList[int(NumIDmut)-1].lower()

        #Fix the format of the mutation_ouput file                                                                                                                                                                                            
        FormatPDB1(mutation_output, "mut_temp")

        with open("mut_temp", "r") as out, open("mutation_output", "w") as temp:
            for line in out:
                if "ATOM" in line:
                    temp.writelines(line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[2]+"\t"+line.split()[3]+"\t"+line.split()[4]+"\t"+line.split()[5]+"\t"+line.split()[6]+"\t"+line.split()[7]+"\t"+line.split()[8]+"\t"+line.split()[9]+"\t"+str("0.0")+"\t"+line.split()[10]+"\n")
        FormatPDB("mutation_output", mutation_output, mutation_output,TitlePDB)
#        os.remove("mut_temp")

        MutatedToARMFormat(mut_pdbARM,NumIDmut, chainName)
        os.system("pwd")
#    os.chdir("../")


#############################                                                                                                                                                                                                               
#MODELLER routine                                                                                                                                                                                                                             
#############################                                                                                                                                                                                                                 
def modeller_mutations(mut_pdbARM):
    """
    This subroutine perfoms the mutation of a single residue, by using MODELLER
    """
    from lib import chainName
    from pdb_format import resNumID
    from download_PDB import TitlePDB
    from external_Software import template

    shutil.copyfile(template+"mutate_model.py", "mutate_model.py")
    script_modeller_mutation = "mutate_model.py"

    print("mutation_seqmut", mutation_seqmut)
    for i in range(0,len(mutation_seqmut)):
        NumIDmut = NumIDmutationsList[i]
        globals().update({"NumIDmut" : NumIDmut})
        ResMutName = mutation_seqmut[i][-3:]


        print("python3.7 "+script_modeller_mutation+" "+str(mut_pdbARM)+" "+str(NumIDmut)+" "+str(ResMutName)+" "+chainName+" > "+str(ResMutName)+str(NumIDmut)+".log")
        os.system("python3.7 "+script_modeller_mutation+" "+str(mut_pdbARM)+" "+str(NumIDmut)+" "+str(ResMutName)+" "+chainName+" > "+str(ResMutName)+str(NumIDmut)+".log")
        print( "\n The mutation ", mutation_seqmut[i], "has been succesfully generated!")

#Move the output pdb of modeller to mut_pdbARMTemp
        shutil.move(mut_pdbARM+"_"+str(ResMutName)+str(NumIDmut)+".pdb",mutation_output)

#Fix the format of the mutation_ouput file                                                                                                                                                                                            
        FormatPDB1(mutation_output, "mut_temp")

        with open("mut_temp", "r") as out, open("mutation_output", "w") as temp:
            for line in out:
                if "ATOM" in line:
                    temp.writelines(line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[2]+"\t"+line.split()[3]+"\t"+line.split()[4]+"\t"+line.split()[5]+"\t"+line.split()[6]+"\t"+line.split()[7]+"\t"+line.split()[8]+"\t"+line.split()[9]+"\t"+str("0.0")+"\t"+line.split()[2][:1]+"\n")
        FormatPDB("mutation_output", mutation_output, mutation_output,TitlePDB)
        os.remove("mut_temp")

        MutatedToARMFormat(mut_pdbARM,NumIDmut, chainName)
        os.system("pwd")

#############################                                                                                                                                                                                    
#
#############################                                                                                                                                                                                  
def MutatedToARMFormat(pdbFile,NumIDmut, chainName):
    """
    Modeller and SCWR4 complete the missing atoms of the amino acids.                                                                                                                                   
    To preserve the structure of the wild type is necessary to include the geometry of the substitution in the wild type file. 
    """
#Obtain the geometry of the new residue mutation                                                                                                                                                                                                                                
    from download_PDB import TitlePDB

    geometry_Mutation=[]
    with open(mutation_output, "r") as oldfile, open(mut_pdbARMTemp, "w") as newfile:
        for line in oldfile:
            if 'ATOM' in line and chainName+'{:>4}'.format(str(NumIDmut)) in line:
                newfile.writelines(line)
                geometry_Mutation.append(line)

    shutil.move(mut_pdbARMTemp, mutation_output)

#Calculate the number of atoms of old residue                                                                                                                                                                                                                                   
    i = 0
    numAtomOldRes = ''
    with open(mut_pdbARM, "r") as oldfile:
        for line in oldfile:
            if 'ATOM' in line and chainName+'{:>4}'.format(str(NumIDmut)) in line:
                i = i+1
        numAtomOldRes = i

#Insert the geometry of the new residue mutation in the wild type geometry                                                                                                                                                                                                      
    i = 0
    with open(mut_pdbARM, "r") as oldfile, open(mut_pdbARMTemp, "w") as newfile:
        for line in oldfile:
            if 'ATOM' in line and chainName+'{:>4}'.format(str(NumIDmut)) in line:
                i = i+1
                if i == numAtomOldRes:
                    newfile.writelines(geometry_Mutation)
            else:
                newfile.writelines(line)

#Write the new pdb using the correct format                                                                                                                                                                                                                                     
    shutil.move(mut_pdbARMTemp, mut_pdbARM)
    FormatPDB1(mut_pdbARM, mut_pdbARMTemp)
    FormatPDB(mut_pdbARMTemp, mut_pdbARM, mut_pdbARM, TitlePDB)


#############################                                                                                                                               
#
#############################                                                                                                                                                                                  
def Mutant_prot_addCounterions(mut_pdbARM):

#Start the protonation and counter-ion placement
    #propKa and charges analysis
    propKa(mut_pdbARM,mut_pdbARM)
    #Protonation states assignment
    protonation(mut_pdbARM)
    #Redefinition of the cavity
    mut_pdbARMFix = mut_pdbARM[:-3]+'fix.pdb'
    fpocket(mut_pdbARMFix)
#Add mutated residue to cavity file                                                                                                                        
    with open("cavity") as cavity:
        listcavity = list(cavity.read().splitlines())
        listcavity.append(NumIDmut)
        cavityListMut0 = list(set(listcavity))
        cavityListMut = sorted(cavityListMut0, key = int)
        globals().update({ "cavityListMut": cavityListMut})
        print( "New cavity including the mutated residues: \n", cavityListMut)

        with open("cavity", "w") as cavity:
            for i in range(0, len(cavityListMut)):
                cavity.writelines(cavityListMut[i]+"\n")

#Calculation of type and number of counter-ions
    numberCounterions(mut_pdbARM)
    addCounterIons(mut_pdbARM)

#Generation of pymol script
    pyMoLFig(mut_pdbARM,cavityListMut)
    pyMoLScript = mut_pdbARM[:-4]+".pml"
    with open(pyMoLScript, "a") as pyMoLScript:
        pyMoLScript.writelines(["sele resi "+NumIDmut+" \n",
                                "create Mutation, sele \n",
                                "show sticks, Mutation \n",
                                "color orange, Mutation \n",
                                ])

#############################
#
#############################
def Mutations_procedure(pdbARM,mutation_seqmut):
#This step ask the user for the list of mutations. Similar to the seqmut file                                                                                                                                       
    NonStandarddResID(pdbARM,mutation_seqmut)
    Insert_mutations(pdbARM,mutation_seqmut)

    from lib import mutation_SoftwareName
    if mutation_SoftwareName == "MODELLER":
        modeller_mutations(mut_pdbARM)
#Start the protonation and counter-ion placement
        Mutant_prot_addCounterions(mut_pdbARM)
        os.chdir("../")

# #SCWRL4 routine                                                                                                                                                                                               
    if mutation_SoftwareName == "Scwrl4":
         Scwrl4_mutations(mut_pdbARM)
#Start the protonation and counter-ion placement
         Mutant_prot_addCounterions(mut_pdbARM)
         os.chdir("../")


#############################
#main Mutations function
#############################
def Mutations():
    """
    By exploiting the backbone-dependent rotamer library implemented in the software SCWRL4, 
    m-ARM has the ability to performing amino acid substitutions on rhodopsin structures to
    generate QM/MM m-ARM models mutants,
    """

    from initial_Setup import pdbARM

    question = yes_no('\n <-> Do you want to perform mutations of the wild type '+pdbARM+' file?')
    if question == True:

        print("Aqui deberia iniciaaaar")

        Select_mutation_Software(pdbARM)
        Mutations_ID(pdbARM)

        for key in mutation_seqmutDic:
            mutation_seqmut = mutation_seqmutDic[key]
            Mutations_procedure(pdbARM,mutation_seqmut)
            

    else:
        print( "---> No mutations requested.")
