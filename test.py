
from LHCO_reader import LHCO_reader

events = LHCO_reader.Events(f_name="example.lhco")

events.plot("electron", "PT")

""" for event in events: """
    
    