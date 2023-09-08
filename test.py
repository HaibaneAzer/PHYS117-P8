import sys
directory_path = "/Users/marcu/anaconda3/envs/phys117_py2/Lib/site-packages"
sys.path.append(directory_path)
from LHCO_reader import LHCO_reader

events = LHCO_reader.Events(f_name="example.lhco")
print (events)