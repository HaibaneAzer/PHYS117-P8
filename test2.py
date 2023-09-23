from LHCO_reader import LHCO_reader
import pandas as pd

events = LHCO_reader.Events(f_name="example.lhco")

# Determine the number of events
num_events = len(events)

for n in range(num_events):
    length = len(events[n]['electron'])
    for k in range(length):
        print(events[n]['electron'][k]['PT'])