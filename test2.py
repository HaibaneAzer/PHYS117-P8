from LHCO_reader import LHCO_reader

events = LHCO_reader.Events(f_name="example.lhco")

objects = ['photon', 'electron', 'muon', 'tau', 'jet']

HT_sum_list = []
meff_sum_list = []
for event in events:
    HT = 0
    for particle in objects:
        if particle in event:
            HT += sum(p["PT"] for p in event[particle])
    HT_sum_list.append(HT)
    missing_ET = sum(p["PT"] for p in event['MET'])
    meff_sum_list.append(missing_ET + HT)

print(HT_sum_list)
print(meff_sum_list)