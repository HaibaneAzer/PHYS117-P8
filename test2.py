from LHCO_reader import LHCO_reader
import pandas as pd

events = LHCO_reader.Events(f_name="BH_n4_M8.lhco")

# Determine the number of events
num_events = len(events)

print(events.column('electron', 'PT'))