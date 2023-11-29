import os
import time

acc = "python /Bioinformatics_Machine_Learning/Machine_Learning/feature_extraction/Pse-in-One-2.0/acc.py"
pse = "python /Bioinformatics_Machine_Learning/Machine_Learning/feature_extraction/Pse-in-One-2.0/pse.py"
sc = "python /Bioinformatics_Machine_Learning/Machine_Learning/feature_extraction/Pse-in-One-2.0/sc.py"

path = '/Bioinformatics_Machine_Learning/Machine_Learning/data/protein/Special_Protein_Identification/enzyme/cd-hit0.3/featurePseinall/'

acc_mode = ['AC', 'CC', 'ACC', 'PDT']
pse_mode = ['PC-PseAAC', 'SC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC-General']

files = ['multipleClass_0.3_.txt', 'noClass_0.3_.txt', 'singleClass_0.3_.txt']
for f in files:
    print time.ctime(), 'name:', f
    for m in acc_mode:
        command = acc + ' ' + path+f + ' Protein '+ m + ' -out ' + path+f[:f.rindex('.')]+m+'.txt'
#         print command
        try:
            os.popen(command).read()
        except:
            continue
            
    for m in pse_mode:
        command = pse + ' ' + path+f + ' Protein '+ m + ' -out ' + path+f[:f.rindex('.')]+m+'.txt'
#         print command
        try:
            os.popen(command).read()
        except:
            continue