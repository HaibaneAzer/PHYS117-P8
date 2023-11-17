from LHCO_reader import LHCO_reader

events = LHCO_reader.Events(f_name="example.lhco")

objects_1 = ['electron', 'muon', 'tau']

""" print(events.count('electron'))
print(events[0].HT())
print(events[0].ET()) """
print(events[0]) # event 1
""" print(events[0]['jet']) # jets at event 1
print(events[0]['jet'][0]) # first jet at event 1
print(events[0]['jet'][0]['PT']) # value of PT at first jet at event 1 """

