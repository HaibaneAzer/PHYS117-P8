from LHCO_reader import LHCO_reader
import numpy as np

events = LHCO_reader.Events(f_name="example.lhco")

# Determine the number of events
num_events = len(events)

particles = [
    'photon',
    'electron',
    'muon',
    'tau',
    'jet',
    'MET',
]
tot_obj = 0
obj_per_event = []
for event in events:
    objects_event = event.number()
    for particle_name in particles:
        tot_obj += objects_event[particle_name]
    obj_per_event.append(tot_obj)
    tot_obj = 0
    
print(obj_per_event[49])
""" i = 0
for n in events:
    if i < 50:
        print(n.number())
    i += 1 """
