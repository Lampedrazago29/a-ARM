import os
import shutil
import Lists_dictionaries
import inspect
import initial_Setup
import download_PDB

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

def FormatPDB(oldFile, newFile, FinalFile, TitlePDB):
    """  This function is used to rewrite the information contained in a PDB file, using the correct 
    format required by the ARM protocol. """

    from Lists_dictionaries import aaList
    from lib import chromophoreName
    from initial_Setup import pdbARMTemp, pdbARM
    from download_PDB import TitlePDB

    global resNumChr # New assigned residue sequence number of the chromophore 

    #Identify gaps

    resNumList = []
    resNumList1 = []
    with open(oldFile, "r") as pdb:
        for line in pdb:
            ls = line.split()
            for i in range(0, len(aaList)):
                if aaList[i] in line:
                    resNumList.append(int(ls[5])) # Residue sequence number                                                                                 
                    resNumList1 = sorted(list(set(resNumList)))                
    
    missResList = [] # This is the list of missing residues
    for i in range(1,resNumList1[-1]+1):
            if i not in resNumList1:
                missResList.append(i)
                
    globals().update({"missResList" : missResList})

    # Write a new file in PDB format
    pdbFormat = "%-6s%5d %s %-4s%3s %s%4d %3s%8.3f%8.3f%8.3f%6.2f%s%6.2f %11s" # Standard PDB format

    with open(oldFile, "r") as oldfile, open(newFile, "w") as newfile:
        newfile.writelines("TEMPLATE: "+ TitlePDB )        


        oldfile_read = oldfile.readlines()
        ln = 0 # Atom serial number

        #For keeping the ACE residues
        for line in oldfile_read:
            ls = line.split()    
            if "ACE" in line:
                ln = ln+1 # Atom serial number
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(ls[5]), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""), float(ls[10]), ls[11])
                newfile.writelines(line+"\n")

        #For writing the a.a.
        j = 0
        for line in oldfile_read:
            ls = line.split()
            j = j+1
            for i in range(0, len(aaList)):
                if aaList[i] in line:
                    resNum = int(ls[5]) # Residue sequence number
                    globals().update({"resNumID" : resNum}) # Last residue sequence number used of amino acids

                    ln = ln+1
                    line1 = pdbFormat  % ("ATOM", int(ln), str(""), ls[2], ls[3],ls[4], int(ls[5]), str(""), float(ls[6]), 
                                         float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                    newfile.writelines(line1+"\n")

### This part writes the TER word between gaps
                    if int(ls[5])+1 in missResList:
                        this_line = oldfile_read[j-1]
                        next_line = oldfile_read[j]
                        this_resNum = this_line.split()[5]
                        next_resNum = next_line.split()[5]

                        if this_resNum != next_resNum:  
                            print("--> The following gaps have been identified:", this_resNum, " to ",  next_resNum)
#Descomentar
                            
                            newfile.writelines("TER \n")
####
        newfile.write("TER \n")

        #For writing the chromophore                                                                                              
        for i in range(1,21):
            for line in oldfile_read:
                ls = line.split()
                # To order the chromophore atoms which come from LYR                                                              
                if chromophoreName in line and ("C"+str(i)+"\t" in line or "C"+str(i)+" \t" in line):
                    resNumChr = resNum+1
                    ln = ln+1
                    line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(resNumChr), str(""), float(ls[6]),
                                        float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                    newfile.writelines(line+"\n")
        newfile.write("TER \n")
#        
        #For writing the waters
        resNumHOH = resNumChr # Residue sequence number        
        for line in oldfile_read:
            ls = line.split()
            if "HOH" in line and "O " in line:
                orgResNum = ls[5]
                
                resNumHOH = resNumHOH+1
                ln = ln+1
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(resNumHOH), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                newfile.writelines(line+"\n")
                
        newfile.write("TER \n")

#LMPG 13-12-18 Include the CL ions which are part of the X-ray structure (i.e. 5B2N)
        #For writing the Cl- ions
        anyCl = False
        for line in oldfile_read:    
            ls = line.split()
            if ('CL' in line) and (('HETATM' in line) or ('ATOM' in line)):
                anyCl = True
                resNumHOH = resNumHOH+1
                ln = ln+1
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(resNumHOH), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                newfile.writelines(line+"\n")
        #        
        if anyCl:
          newfile.write("TER \n")
#
#        os.remove(oldFile)
        shutil.move(newFile, FinalFile)
    
    globals().update({"LastResNum" : resNumHOH}) # Last residue sequence number used. Employed later for the numeration of the counterions
    globals().update({"lN" : ln}) # Last atom serial number used. Employed later for the numeration of the counterions
