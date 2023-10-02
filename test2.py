from LHCO_reader import LHCO_reader

events = LHCO_reader.Events(f_name="example.lhco")

objects_1 = ['electron', 'muon', 'tau']

print(events.count('electron'))
print(events[0].HT())
print(events[0].ET())
print(events[0]) # event 1
print(events[0]['jet']) # jets at event 1
print(events[0]['jet'][0]) # first jet at event 1
print(events[0]['jet'][0]['PT']) # value of PT at first jet at event 1

# number of electrons/muons/taus in an event
electrons_per_event = []
muons_per_event = []
taus_per_event = []
for event_idx in range(len(events)):
    electrons_per_event.append(len(events[event_idx][objects_1[0]]))
    muons_per_event.append(len(events[event_idx][objects_1[1]]))
    taus_per_event.append(len(events[event_idx][objects_1[2]]))
        
print(electrons_per_event)

# PT of objects with largest PT in event (electron, muon, tau, jet, MET)
# also extracting phi values
objects_2 = ['electron', 'muon', 'tau', 'jet']

Largest_PT_per_event = []
phi_per_event_L_PT = []
phi_met_per_event = []
for event_idx in range(len(events)):
    current_PT = 0 # reset PT
    current_phi_highest_PT = 0 # reset phi
    if events[event_idx]['MET']:
        phi_met = events[event_idx]['MET'][0]['phi']
        phi_met_per_event.append(phi_met)
    # compare every object in event
    for particle in objects_2: 
        if events[event_idx][particle]: # check if the object is in event
            for object_idx in range(len(events[event_idx][particle])):
                if current_PT < events[event_idx][particle][object_idx]['PT']:
                    current_PT = events[event_idx][particle][object_idx]['PT']
                    current_phi_highest_PT = events[event_idx][particle][object_idx]['phi']
    
    Largest_PT_per_event.append(current_PT)
    phi_per_event_L_PT.append(current_phi_highest_PT)


print(Largest_PT_per_event[1])
print(phi_per_event_L_PT[1])
print(events[1]) # check if correct

# difference in phi between such objects and met

delta_phi_per_event = []

for i in range(len(phi_met_per_event)):
    delta_phi_value = abs(phi_per_event_L_PT[i] - phi_met_per_event[i])
    delta_phi_per_event.append(delta_phi_value)

print(delta_phi_per_event[1])