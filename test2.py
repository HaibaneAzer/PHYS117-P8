from LHCO_reader import LHCO_reader
import numpy as np

events = LHCO_reader.Events(f_name="example.lhco")

# Determine the number of events

print(len(events))

particles = [
    'photon',
    'electron',
    'muon',
    'tau',
    'jet',
    ]
es = events
max_length = max(len(es.column(particles[0], 'PT')),
                 len(es.column(particles[1], 'PT')),
                 len(es.column(particles[2], 'PT')),
                 len(es.column(particles[3], 'PT')),
                 len(es.column(particles[4], 'PT')))
print(max_length)
HT_sum_list = []
for idx in range(max_length):
    print(idx)
    HT = 0
    for part in particles:
        if idx < len(events.column(part, "PT")):
            """ print("{}: ".format(part), len(events.column(part, 'PT'))) """
            HT += events.column(part, 'PT')[idx]
    HT_sum_list.append(HT)
print(HT_sum_list)
    