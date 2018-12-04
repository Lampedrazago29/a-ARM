import time
from lib import *
import initial_Setup
from initial_Setup import *
from rotamers import *
import download_PDB
from download_PDB import *
from pdb_format import *
import vmd_interface
from vmd_interface import *
import protonationStep
from protonationStep import *
import cavity_fpocket
from cavity_fpocket import *
from counterions_Placement import *
from mutate_PDB import *
from pymol_Script import *

def main():
    t0 = time.time()

    DownloadPDB()
    SearchFileType('pdb', '\n The following', 'files are found: \n', 'file will be used for preparing the ARM input file.')
    createWorkingFolder("workingFolder")         
    SearchChain()
    createpdbARM()
    searchRotamers()
    SearchChromophore()
    FormatPDB(initial_Setup.pdbARMTemp, initial_Setup.pdbARM, initial_Setup.pdbARM, download_PDB.TitlePDB)
    alignRotAxis()
    FormatPDB1( initial_Setup.pdbARM, initial_Setup.pdbARMTemp )
    FormatPDB(initial_Setup.pdbARMTemp, initial_Setup.pdbARM, initial_Setup.pdbARM, download_PDB.TitlePDB)
    wT_analysis(initial_Setup.pdbARM, lib.pdbName, initial_Setup.pdbARMTemp)
    propKa(lib.pdbName,initial_Setup.pdbARM)
    protonation(initial_Setup.pdbARM)
    fpocket(protonationStep.pdbARMFix)
    numberCounterions(initial_Setup.pdbARM)
    addCounterIons(initial_Setup.pdbARM)
    pyMoLFig(initial_Setup.pdbARM, cavity_fpocket.cavityList)
    os.chdir("../")
    Mutations()


    print("Total excecution time (min):", (time.time() - t0)/60)

if __name__ == '__main__':
    main()
