from LHCO_reader import LHCO_reader

events = LHCO_reader.Events(f_name="example.lhco") # 10000

objects_1 = ['electron', 'muon', 'tau']

new_events_odd = events[1::2]
new_events_even = events[0::2]
print(len(new_events_even)) # 5000 
print(len(new_events_odd)) # 5000

