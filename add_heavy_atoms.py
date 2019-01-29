
def add_Heavy_atoms():
    """ This subroutine identifies the geometry of the missing heavy atoms present in the file "missing-atoms.dat". 
    Then, the geometry corresponding to each missing heavy atom is placed in the correct place in the pdbARM file. 
    """
    from initial_Setup import pdbARM
    from lib import chainName


    pdbARM_TER = pdbARM+".temp2"

    pdbFormat = "%-6s%5d %s %-4s%3s %s%4d %3s%8.3f%8.3f%8.3f%6.2f%s%6.2f %11s" # Standard PDB format                                                                                                                                                                                                                          

    #This is the list with the incomplete residues
    missingHeavyAtoms=[]

    print("\n --> The following coordinates of missing heavy atoms have been added to the ", pdbARM, "file. \n")

    with open("missing-atoms.dat", "r") as outputPQR, open("new_coordinates.dat", "w") as newAtoms:
        for line in outputPQR:
            ls = line.split()
            if "Added" and "coordinates" in line:

                AtomID = ls[1]
                line = pdbFormat % ("ATOM", int("1"), str(""), ls[1], ls[3],ls[4], int(ls[5]), str(""), float(ls[8]),
                                    float(ls[9]), float(ls[10]), float("1.00"),str(""), float("1.00"), AtomID[:1])
                print(line)
                newAtoms.writelines(line+"\n")
                missingHeavyAtoms.append(ls[5])


    missingHeavyAtoms = sorted(list(set(missingHeavyAtoms)), key = int)
    globals().update({"missingHeavyAtoms" : missingHeavyAtoms})

###############                                                                                                                                                                                                                                                                                                               
    with open(pdbARM, 'r') as oldTemp, open(pdbARM_TER, 'w') as newTemp:
        for line in oldTemp:
            if 'TER' not in line:
                newTemp.writelines(line)

    with open(pdbARM_TER) as oldfile, open(pdbARM, "w") as newfile:

        j=0
        oldfile_read = oldfile.readlines()
        for line in oldfile_read:

            j=j+1

            newfile.writelines(line)

            for i in missingHeavyAtoms:

                if 'ATOM' in line and str(chainName)+'{:>4}'.format(str(i)) in line:

                    this_line = oldfile_read[j-1]
                    next_line = oldfile_read[j]
                    this_resNum = this_line.split()[5]
                    next_resNum = next_line.split()[5]

                    if this_resNum != next_resNum:

                        with open("new_coordinates.dat", "r") as coord:
                            for line1 in coord:

                                if 'ATOM' in line1 and str(chainName)+'{:>4}'.format(str(this_resNum)) in line1:

                                    newfile.writelines(line1)
                                    
    FormatPDB1( initial_Setup.pdbARM, initial_Setup.pdbARMTemp )
    FormatPDB(initial_Setup.pdbARMTemp, initial_Setup.pdbARM, initial_Setup.pdbARM, download_PDB.TitlePDB)
