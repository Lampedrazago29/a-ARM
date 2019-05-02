import sys
import logging
import glob
import shutil
import os
import tempfile
import re

logger = logging.getLogger("lib")
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger.warning('This library is useful for the execution of the a-ARM protocol')

##################################################
# yes_no
##################################################          
def yes_no(question, default= None):
    """This function asks the user a yes/no question via input(), and returns the answer
    as the variable *user_choice*. """
    valid_option = {"yes": True,
                    "y": True,
                    "no": False,
                    "n": False
                    }
    if default is None:
        option = " [y/n] "
    elif default == "yes":
        option == " [Y/n] "
    elif default == "no":
        option == " [y/N] "
    else:
        raise ValueError("'%s'  is a non valid option. Please select " % default)
    
    while True:
        print(question + option)
        user_choice = input().lower()
        if default is not None and user_choice == '':
            return valid_option[default]
        elif user_choice in valid_option:
            return valid_option[user_choice]
        else:
            print("Please respond with 'yes' or 'no' (or 'y' or 'n'). \t")
        
##################################################
# SearchFileType
##################################################                                      
def SearchFileType(ext, message0 = "",  message1 = "", message2 = ""):
    """ This function identify in the working folder all the files with an specific extention (.ext),
    makes an enumerate list with them and ask the user to choose the correct one """
    extList = glob.glob('*'+ext)
    ChooseNumOption(extList, "file", ext, message0,  message1, message2, True)

##################################################                                                                             
# ChooseNumOption                                                                                                              
##################################################
def ChooseNumOption(nameList, element, type, message0, message1, message2, pick, dictionary = {}):
    """ This function identify all the elements of a list (*nameList*), creates an enumerate list, ask the user to 
    pick the number of the correct one and creates a global variable with that information """
    if nameList:
        print( message0, '\x1b[0;31;43m'+str(type)+'\x1b[0m', message1)
        for i, element in enumerate(nameList, 1):
            if dictionary:
                print( str(i)+str(')'), element, dictionary.get(element))
            else:
                print( str(i)+str(')'), element)

        if pick == True:
            number = PickNumber(len(nameList))
            for i, element in enumerate(nameList, 1):
                if i == number:
                    print( '\n ---> The', '\x1b[0;31;43m'+str(element)+'\x1b[0m', message2) # This file is the pdbFile 
                    globals().update({type+str("Name") : element})
                    return element
        else:
            return
#in case of files                                                                                                                                                                    
    else:
        if element == "file":
            print( '\n No '+'\x1b[6;30;42m'+ '.'+type +'\x1b[0m', 'files found. Please put a ' +'\x1b[6;30;42m'+ '.'+type +'\x1b[0m', 'file in the current folder and start again. \n Good bye! \n')
            sys.exit()
        if element != "file" and element != "chromophore" and element != "chain":
            print( '\n No '+'\x1b[6;30;42m'+ type +'s'+'\x1b[0m', 'found in the ', '\x1b[0;31;43m'+pdbName+'\x1b[0m', "file")
        else:
            print( '\n No '+'\x1b[6;30;42m'+ type +'s'+'\x1b[0m', 'found in the ', '\x1b[0;31;43m'+pdbName+'\x1b[0m', "file. The file is corrupt,  verify your PDB and start again! \n Good bye! \n")
            sys.exit()

##################################################                                                                                        
# PickNumber                                                                                                                                
##################################################                                                                                  
def PickNumber(lenList, message = ' To select the correct option pick a number in range ',min = 1, typeInput = int):
    """ This function ask for a number via input() and return the number """
    while True:
        try:
            input1 = typeInput(input('\n'+message+str(min)+'-'+str(lenList)+': \t'))
        except ValueError:
            print( 'That\'s not a number!')
        else:
            if min <= input1 <= lenList:
                return input1
            else:
                print( 'Number out of range. Try again!')

##################################################
# createWorkingFolder
##################################################
def createWorkingFolder(newpath):
    """ This function creates a new folder called *newpath*, and moves there the PDB and .seqmut files """
    newpath = pdbName[:-4]+'_ARM_input'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        globals().update({ "workingFolder" : newpath+"/"})

    shutil.move( pdbName, newpath+"/"+pdbName) 
    if glob.glob("*.seqmut"):
        MutFile = glob.glob("*.seqmut")[0]
        shutil.copyfile(str(MutFile), newpath+"/"+pdbName[:-4]+".seqmut")
    else:
        pass
    os.chdir(newpath)

##################################################
# FormatPDB1
##################################################
def FormatPDB1(pdb, pdbTemp):
    """ This function generates a temporal formatted pdb file,in which the columns are separated by tabs.
    This separation is fundamental to read the content of the columns in the PDB file """
    with open(pdb, "r") as pdbOld, open(pdbTemp, "w") as pdbNew:
        for line in pdbOld:
            if "TER" not in line and "END" not in line:
                file = line[:+16]+"\t"+line[16:22]+"\t"+line[22:60]+"\t"+line[60:]
                pdbNew.writelines(file)

##################################################
# replaceLine
##################################################
def replaceString(oldString, newString):
    regex = re.compile("(%s)" % "|".join(map(re.escape, newString.keys())))
    return regex.sub(lambda x: str(newString[x.string[x.start() :x.end()]]), oldString)

def replaceLine(oldFile, string1, string2, newString, newFile = "TempFile", mvFile = True):
    """ This function replaces a character or string with information contained in a dictionary """
    with open(oldFile, "r") as oldfile, open(newFile, "w") as newfile:
        oldfile_read = oldfile.readlines()
        for line in oldfile_read:
            line_number = oldfile_read.index(line)
            if string1 in line and string2 in line:
                oldfile_read[line_number] = replaceString(oldfile_read[line_number],newString)
                newfile.writelines(oldfile_read[line_number])
            else:
                newfile.writelines(oldfile_read[line_number])

    if mvFile == True:
        shutil.move(newFile, oldFile)

##################################################
# deleteLine
##################################################
def deleteLine(oldFile, string1, string2, newString,  newFile = "TempFile"):
    """ This function deletes the lines with a determined pattern  """
    with open(oldFile, "r") as oldfile, open(newFile, "w") as newfile:
        oldfile_read = oldfile.readlines()
        for line in oldfile_read:
            line_number = oldfile_read.index(line)
            if 'ATOM' in line and not (string1 in line and string2 in line.split()[5]):
                newfile.writelines(oldfile_read[line_number])
    shutil.move(newFile, oldFile)

##################################################
# DicResNameNum 
##################################################
def DicResNameNum(group, type, typeList, eraseList, typeDictionary, pick):
    """ This function identifies the residue sequence number and creates a dictionary with the residue 
    name and the residue sequence number """
    typeList = sorted(list(set(typeList).difference(eraseList)))
    globals().update({type+str("List") : typeList})
    pdbARMTemp = pdbName[:-3]+chainName+'.temp'

    with open(pdbARMTemp) as pdbTemp:
        for line in pdbTemp:
            for i in range(0, len(typeList)):
                if "ATOM" in line and typeList[i] in line:
                    typeDictionary.update({ typeList[i] : line.split()[5]} )

    globals().update({type+str("Dic") : typeDictionary})

    ChooseNumOption(typeList,type, type, '\n The following are identified as possible', 'in the '+str(pdbName)+ ' file:', 'is selected as the '+type, pick, typeDictionary)


##################################################                                                                                                                                          
#Wild type analysis                                                                                                                                                                         
##################################################                                                                                                                                            
def wT_analysis(pdbARM, pdbName, pdbARMTemp):

    wt_workingFolder = "wt_"+workingFolder
    if not os.path.exists(wt_workingFolder):
        os.makedirs(wt_workingFolder)
        shutil.copyfile(pdbARM, wt_workingFolder+"/"+pdbARM)
        shutil.copyfile(pdbARMTemp, wt_workingFolder+"/"+pdbARMTemp)
        shutil.copyfile(pdbName, wt_workingFolder+"/"+pdbName)

    else:
        pass
    os.chdir(wt_workingFolder)
