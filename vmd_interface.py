import os
import tempfile

##################################################                                                                                           
# VMDTempFile                                                                                                                                
################################################## 
def VMDTempFile(tempName, content, position, default=False):
    """ This function executes VMD without needed of graphics """
    with tempfile.NamedTemporaryFile(mode='w+') as tempNameI:
        tempNameI.writelines(content)
        tempNameInp = tempNameI.name
        tempNameI.seek(0)

        with tempfile.NamedTemporaryFile(mode='w+') as tempNameO:
            tempNameOut = tempNameO.name
            os.system('vmd -dispdev text -eofexit <' + tempNameInp + '> ' + tempNameOut )

            if default == True:
                for line in tempNameO:
                    if  '{' in line:
                        line = line.replace("{", "").replace("}", "").replace("\n", "")
                        positionXYZ = line.split(" ")
                        positionXYZ[:] = [float(x)*-1 for x in positionXYZ] ## This step is fundamental to obtain a good alignment 

                        globals().update({position+str("XYZ") : positionXYZ})
                        print("We should compare this step with the oldest version:", position+str("XYZ"), positionXYZ)
                        return positionXYZ
            
        

                        
