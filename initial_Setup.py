from lib import ChooseNumOption, FormatPDB1, DicResNameNum
import vmd_interface
from vmd_interface import *
import Lists_dictionaries
import shutil

from add_heavy_atoms import add_Heavy_atoms

##################################################                                                                              
# SearchChain
##################################################                                                                            
def SearchChain(type = "chain", column = 4, pick = True, group = "ATOM"):
    """ This function identify all the possible chains in the *pdbName* file, ask the user to choose 
    one of them and stores the chain identifier as the global variable *chainName* """

    from lib import pdbName
    typeList = []
    listtype = []
    with open(pdbName) as pdbInitial:
        for line in pdbInitial:
            if group in line.split()[0] and (line.split()[column]).isalpha():
                listtype.append(line.split()[column])
                typeList = sorted(list(set(listtype)))

        ChooseNumOption(typeList,type, type, '\n The following are identified as possible', 'in the '+str(pdbName)+ ' file:', 'is selected a\
s the '+type, pick, "")

##################################################                                                                              
# createpdbARM
##################################################                                                                        
def createpdbARM():
    """ This function generates a clean working PDB file, named *pdbARM*, which contains
    only the selected chain """

    from lib import pdbName, chainName
    global pdbARM # This is the final *ARM.pdb file                                                                                          
    global pdbARMTemp # This is the temporal .pdb file with the columns separated by tabs                                                    
    pdbARM = pdbName[:-3]+'ARM.pdb'
    pdbARMTemp = pdbName[:-3]+chainName+'.temp'
    
    string0 = ['mol load pdb %s \n' % pdbName,
               'set sel [atomselect top \"chain '+chainName+ '\" ]\n',
               '$sel writepdb ' +pdbARM+' \n']
    
    VMDTempFile("centerVMD", string0, "center",  False)
    
    # we remove hydrogen atoms if present
    stringH = ['mol load pdb ' +pdbARM+ ' \n',  
               'set sel [atomselect top \"not hydrogen\" ]\n',
               '$sel writepdb ' +pdbARM+ ' \n']
    
    VMDTempFile("centerVMD", stringH, "center", False)
    
    FormatPDB1(pdbARM, pdbARMTemp) # The *ARM.pdb file is re-written with the columns separated by tabs 
    
##################################################                                                                              
# searchChromophore
##################################################                                                                        
def SearchChromophore(type = "chromophore", pick = True, group = 'ATOM'):
    """ This function identifies the possible chromophores and ask the user to select the correct one """
    from Lists_dictionaries import totalList
    listtype = []
    typeList = []
    with open(pdbARMTemp) as pdbTemp:
        for line in pdbTemp:
            if group in line:
                listtype.append(line.split()[3])
                typeList = sorted(list(set(listtype)))
    globals().update({type+str("List") : typeList})

    DicResNameNum(" ", type, typeList, totalList, {}, pick)


##################################################                                                                              
# SearchLinkerAA
##################################################                                                                        
def SearchLinkerAA(diffPar = 1.5):
    """ This function is designed to simultaneously identify the linker atom of the chromophore and the linker amino acid residue. 
    To this aim, the x, y and z coordinates of the atoms of the chromophoreName are stored in individual lists 
    (coordXChr, coordYChr, coordZChr). Then, the difference between each coordinate of the chromophore and each 
    coordinate of all the amino acids is calculate. This exploration allows to identify the atoms closest to the 
    chromophore. Therefore, the criteria for identifying the linker atom and the linker amino acid is to find the 
    closest pair of atoms considering the x, y and z coordinates. """

    from chromophore import chromophoreName
    from Lists_dictionaries import aaList

    numAtomsChr = 0
    with open(pdbARMTemp) as file:
        coordXChr = []
        coordYChr = []
        coordZChr = []
        labelsChr = []
        for line in file:
            numAtomsChr = numAtomsChr+1
            if chromophoreName in line:
                coordXChr.append(float(line.split()[6]))
                coordYChr.append(float(line.split()[7]))
                coordZChr.append(float(line.split()[8]))
                labelsChr.append(line.split()[2])
                globals().update({"coordXChr" : coordXChr})
                globals().update({"coordYChr" : coordYChr})
                globals().update({"coordZChr" : coordZChr})
                globals().update({"labelsChr" : labelsChr})

    with open(pdbARMTemp) as file:
        for line in file:
            for i in range(0, len(aaList)):
                if aaList[i] in line and ("NZ" in line.split()[2] or "O" in line.split()[2]):
                    for j in range(0, len(coordXChr)):                    
                        diffX = float(abs((float(line.split()[6])) - coordXChr[j]))
                        diffY = float(abs((float(line.split()[7])) - coordYChr[j]))
                        diffZ = float(abs((float(line.split()[8])) - coordZChr[j]))
                        if diffX < diffPar and diffY < diffPar and diffZ < diffPar:
                            linker_aa_ID = line.split()[5]
                            linker_aa = str(line.split()[3]+" "+linker_aa_ID)
                            globals().update({"linker_aa_ID" : linker_aa_ID})
                            globals().update({"linker_aa" : linker_aa})
        
                            coordXYZChr = [coordXChr[j], coordYChr[j], coordZChr[j]]
                            globals().update({"coordXYZChr" : coordXYZChr})
    print("The linker amino acid is:", linker_aa)

##################################################                                                                              
# alignRotAxis 
##################################################                                                                        
def alignRotAxis():
    """ This subroutine, designed by Dr. Luca De Vico, is a function which allows to align the geometry of the complex
    protein+chromophore in the center of mass of the chromophore and to oriente it along the Z axis, using the VMD package """

    from lib import chainName
    from chromophore import chromophoreName
    from Lists_dictionaries import aaList
#    from initial_Setup import pdbARM 

#### This subroutine add the geometry of the heavy atoms for incomplete residues                                                                                                                                                                                           
    add_Heavy_atoms()
###

    string1 = ['mol load pdb %s \n' % pdbARM,  
               'molinfo 0 get center\n']
    print( "\n", str("Centering the protein+"+chromophoreName+" on the center of mass of the complex").rjust(100, "."))
    VMDTempFile("centerVMD", string1, "center", True)    

    from vmd_interface import centerXYZ

    string2 = ['mol load pdb %s \n' % pdbARM,
               'set sel [atomselect 0 \"all\"]\n',
               'atomselect0 moveby {%.6f %.6f %.6f}\n' %(float(centerXYZ[0]), float(centerXYZ[1]), float(centerXYZ[2])),
               'package require Orient\n',
               'namespace import Orient::orient\n',
               'set sel [atomselect top \"all\"]\n',
               'set I [draw principalaxes $sel]\n',
               'set A [orient $sel [lindex $I 2] {0 0 1}]\n',
               '$sel move $A\n',
               'set I [draw principalaxes $sel]\n',
               'set A [orient $sel [lindex $I 1] {0 1 0}]\n',
               '$sel move $A\n',
               'set I [draw principalaxes $sel]\n',
               '$sel writepdb %s \n' % pdbARM]
    VMDTempFile("orientVMD", string2, "orient", False)

    string3 = ['mol load pdb %s \n' % pdbARM,
               'set sel [atomselect top \"resname '+chromophoreName+'\"]\n',
               '$sel writepdb ' +chromophoreName+"."+chainName+'.pdb \n',
               'mol delete top\n',
               'mol load pdb ' +chromophoreName+"."+chainName+'.pdb \n',
               'molinfo top get center\n' ]
# create file with chromophore coordinates
    VMDTempFile("newcenterVMD", string3, "newcenter", True) 
    from vmd_interface import newcenterXYZ
#center the proteinx on the retinal center of mass
    string4 = ['mol load pdb %s \n' % pdbARM,
               'set sel [atomselect 0 \"all\"]\n',
               'atomselect0 moveby {%.6f  %.6f  %.6f}\n' %( float(newcenterXYZ[0]), float(newcenterXYZ[1]), float(newcenterXYZ[2])),
               'molinfo top get center\n',
               '$sel writepdb %s\n' % pdbARM] #VMD saves the new pdb file                                                                                               
    VMDTempFile("orientVMD", string4, "orient", True)
    print( str("Done!").rjust(100, "."))

    os.remove(chromophoreName+"."+chainName+'.pdb')
    shutil.copy(pdbARM, pdbARMTemp)
    from pdb_format import resNumChr
    print( "\n ---> The chromophore label has been changed to: ", '\x1b[0;33;49m'+str(chromophoreName)+" "+str(resNumChr)+'\x1b[0m')
#identify the linker atom and the linker amino acid
    SearchLinkerAA()
    FindCounterIons(linker_aa_ID)
    print( "\n ----> The", '\x1b[0;33;49m'+str(pdbARM)+'\x1b[0m', "ARM input file has been generated.")

##################################################                                                                                       
#This function locates the two nearest possible counterions to the RET                                                                     
#Luca De Vico                                                                                                                    
##################################################                                                                                                 
def FindCounterIons(linker_aa_ID):

    from lib import chainName

    counterionsfilename = "output_counter_ions.dat"
    string0 = ['mol load pdb %s \n' % pdbARMTemp,
               'set negatives [atomselect top \"(resname GLU or resname ASP) and type CG\" ]\n',
               'set linker [[atomselect top \"resid ' + linker_aa_ID + ' and type NZ\"] list]\n',
               'set outfile [open ' + counterionsfilename +  ' w]\n',
               'puts $outfile [$negatives get {resname resid}]\n',
               'foreach coord [atomselect0 list] {\n',
               'puts $outfile [measure bond [list $coord $linker]]\n',
               '}\n',
               'close $outfile\n']

#execute VMD                                                                                                                                                                                                                                  
    VMDTempFile("counterions", string0, "counter", False)

# reads in the output file generated by VMD                                                                                                                                                                                                   
    counterionsfile = open(counterionsfilename, 'r')
    outputlines = counterionsfile.readlines()
    counterionsfile.close()

# first line is the list of possible counterion AA                                                                                                                                                                                            
    counters = outputlines[0].replace("{","").replace("}","").replace("\n","").split()

# now creates the list of possible counterions and their distance from the linker AA                                                                                                                                                          
# format distances[ [AA_name, AA_id, distance], ... ]                                                                                                                                                                                         
    distances = []
    for i in range(1, len(outputlines)):
        distances.append([counters[(i-1)*2], counters[1+(i-1)*2], float(outputlines[i].replace("\n",""))])

#order the distances list by the actual distance, item [2]                                                                                                                                                                                    
    distances.sort(key=lambda distance: distance[2])

# retrieves the best and second best guesses at the counter ion                                                                                                                                                                               
    bestguess = [distances[0][0], distances[0][1]]
    secondbest = [distances[1][0], distances[1][1]]

    counterion_ID = distances[0][1]
    second_counterion_ID = distances[1][1]

    globals().update({"counterion_ID" : counterion_ID})
    globals().update({"second_counterion_ID" : second_counterion_ID})
    with open(pdbARMTemp) as pdbTemp:
        for line in pdbTemp:
            if "ATOM" in line and chainName+'{:>4}'.format(str(counterion_ID)) in line:
                main_counterion = str(line.split()[3]+" "+counterion_ID)
                globals().update({"main_counterion" : main_counterion})

                second_counterion = str(line.split()[3]+" "+second_counterion_ID)
                globals().update({"second_counterion" : second_counterion})

    print("The main and second counter-ion residues are:", main_counterion, "and", second_counterion )
    return [bestguess, secondbest]
