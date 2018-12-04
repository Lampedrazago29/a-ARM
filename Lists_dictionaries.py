
##################################################  
# My Lists and Dictionaries                                                                                       
##################################################  

rmList =sorted([ "ACE", "HG", "HOH", "ZN", "HTG", "HTO", "MAN", "NAG", "BMA", "SO4", "BNG", "A", "CL", "NA", "BOG",
                 "PLM", "TWT", "OLA", "PEE", "GLC", "GAL", "L2P", "L3P"])

aaList=["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE",
        "PRO","SER","THR","TRP","TYR","VAL"]

totalList= aaList+rmList

chargedList= ["C-", "Oco", "N+", "ACE"]

protList = ['ASP', 'HIS', 'LYS', 'GLU', 'ARG']+chargedList


protAA = {'ASP': 'ASH', 
          'LYS':'LYD', 
          'GLU':'GLH',
          'ARG':'ARN'}


##################################################
# propka modeled pKa 
##################################################  
protAApKa = {'ASP': '3.80',
             'LYS':'10.5', 
             'GLU':'4.50',
             'ARG':'12.5',
             'HIS': '6.50',
             'C-': '3.20',
             'N+': '8.00'}

##################################################
# AA
##################################################  

threeToOneDic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'UNK': '*'}
oneToThreeDic = {}

for key in threeToOneDic:
    oneToThreeDic.update({threeToOneDic[key] : key})
