import urllib.request
from urllib.request import urlopen
import shutil
import lib
from lib import yes_no

##################################################
# DownloadPDB
##################################################          
def DownloadPDB():
    """ This function is used to download PDB files directly from the 
    RCSB Protein Data Bank web page """
    url_rcsb = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="
    question = yes_no('\n Do you want to download a PDB file from the RCSB Protein Data Bank webpage? \t')
    global TitlePDB
    TitlePDB = "\n"
    while(1):
        if question == True:
            while(1):
                pdbID = input("<-> Please type the 4-character unique identifier of entry in the Protein Data Bank -PDB ID-. (i.e. 1U19): \t")
#For a valid PDB ID
                try:
                    with urllib.request.urlopen(url_rcsb+str(pdbID)) as pdb_id_content:
                        open( pdbID+".pdb", "w" ).write(pdb_id_content.read().decode('utf-8'))
                        #This part allows the user to re-name the project (new!)
                        renamePDB = input("Please, introduce the name of the project: \t")
                        shutil.move( pdbID+".pdb", renamePDB+".pdb")
                        
                    #Print the Title of the PDB ID in screen        
                        with open(renamePDB+".pdb", "r") as pdbFile:
                            for line in pdbFile:
                                if 'TITLE' in line:
                                    TitlePDB = str(line)
                                    print( "\n <---> The protein you have requested is: \t", TitlePDB)
                                    return TitlePDB
                        
#For an invalid PDB ID
                except:
                    question2 = yes_no("\n The requested "+'\x1b[0;31;43m'+pdbID+'\x1b[0m'+" PDB ID was not found on the server. \n Do you want to download another PDB file?")
                    if question2 == True:
                        continue
                    else:
                        question3 = yes_no("\n ---> No file has been downloaded. \n Do you want to continue to the next step?")
                        if question3 == True:
                            return
                        else:
                            print( "\n OK. Have a nice day and do not try to build the ARM input manually!")
                            sys.exit()

        else:
            question3 = yes_no("\n ---> No file has been downloaded. \n Do you want to continue to the next step?")
            if question3 == True:
                return
            else:
                print( "\n OK. Have a nice day and do not try to build the ARM input manually!")
                sys.exit()
