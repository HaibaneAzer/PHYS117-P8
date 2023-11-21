from LHCO_reader import LHCO_reader
import os
background = "background/"
events1 = LHCO_reader.Events(f_name=background+"ttbar.lhco") # 100000
events2 = LHCO_reader.Events(f_name=background+"ttbar_largejet.lhco")

# print(len(events1)) # 100000
# print(len(events2)) # 100000

dirpath = os.path.join(os.path.dirname(__file__), "sphaleron")
path_list = []
for _, _, filename in os.walk(dirpath):
    path_list = ["{}\{}".format(dirpath,file) for file in filename]    

for path in path_list:
    print(path + ": ", len(LHCO_reader.Events(f_name=path)))

