import numpy as np

### main computational functions ###

def calculate_HT_and_meff(events):
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
    
    return HT_sum_list, meff_sum_list

def objects_per_event(events):
    particles = [
    'photon',
    'electron',
    'muon',
    'tau',
    'jet',
    'MET',
    ]
    obj_per_event = []

    for event in events:
        tot_obj = 0
        for particle_name in particles:
            tot_obj += event.number()[particle_name]
        obj_per_event.append(tot_obj)
    
    return obj_per_event

def electrons_muons_taus_per_event(events):
    particles = [
    'electron',
    'muon',
    'tau',
    ]
    # number of electrons/muons/taus in an event
    electrons_per_event = []
    muons_per_event = []
    taus_per_event = []
    for event_idx in range(len(events)):
        electrons_per_event.append(len(events[event_idx][particles[0]]))
        muons_per_event.append(len(events[event_idx][particles[1]]))
        taus_per_event.append(len(events[event_idx][particles[2]]))

    return electrons_per_event, muons_per_event, taus_per_event

def largest_PT_in_event(events):
    # PT of objects with largest PT in event (electron, muon, tau, jet, MET)
    # also extracting phi values
    objects_2 = ['electron', 'muon', 'tau', 'jet']

    Largest_PT_per_event = [] # Largest PT variables
    phi_per_event_L_PT = [] # phi variables
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
    
    delta_phi = delta_phi_per_event(phi_met_per_event, phi_per_event_L_PT)
    
    return Largest_PT_per_event, delta_phi

def delta_phi_per_event(phi_met_per_event, phi_per_event_L_PT):
    # difference in phi between such objects and met

    delta_phi_per_event = []
    for i in range(len(phi_met_per_event)):
        phi_L = phi_per_event_L_PT[i] # angles phi varies from -pi to pi
        phi_met = phi_met_per_event[i]
        delta_phi_value = abs(phi_L - phi_met)
        if delta_phi_value > np.pi: # handle wrap around the circle
            delta_phi_value = 2*np.pi - delta_phi_value
        delta_phi_per_event.append(delta_phi_value)

    return delta_phi_per_event

def delta_eta_per_event(eta_met_per_event, eta_per_event_L_PT):
    
    delta_eta_per_event = []
    for i in range(len(eta_met_per_event)):
        eta_L = eta_per_event_L_PT[i] # angles eta varies from -pi to pi
        eta_met = eta_met_per_event[i]
        delta_eta_value = abs(eta_L - eta_met)
        if delta_eta_value > np.pi: # handle wrap around the circle
            delta_eta_value = 2*np.pi - delta_eta_value
        delta_eta_per_event.append(delta_eta_value)

    return delta_eta_per_event

def delta_R_per_event(events, particle):

    # PT of objects with largest PT in event (electron, muon, tau, jet, MET)
    # also extracting phi values

    phi_per_event_L_PT = [] # phi variables
    phi_met_per_event = []
    eta_per_event_L_PT = [] # Eta variables
    eta_met_per_event = []
    for event_idx in range(len(events)):
        current_PT = 0 # reset PT
        current_phi_highest_PT = 0 # reset phi
        current_eta_highest_PT = 0 # reset eta
        if events[event_idx]['MET']:
            phi_met = events[event_idx]['MET'][0]['phi']
            eta_met = events[event_idx]['MET'][0]['eta']
            phi_met_per_event.append(phi_met)
            eta_met_per_event.append(eta_met)
        # compare every object in event
    
        if events[event_idx][particle]: # check if the object is in event
            for object_idx in range(len(events[event_idx][particle])):
                if current_PT < events[event_idx][particle][object_idx]['PT']:
                    current_PT = events[event_idx][particle][object_idx]['PT']
                    current_phi_highest_PT = events[event_idx][particle][object_idx]['phi']
                    current_eta_highest_PT = events[event_idx][particle][object_idx]['eta']

    
        phi_per_event_L_PT.append(current_phi_highest_PT)
        eta_per_event_L_PT.append(current_eta_highest_PT)

    delta_eta = delta_eta_per_event(eta_met_per_event, eta_per_event_L_PT)
    delta_phi = delta_phi_per_event(phi_met_per_event, phi_per_event_L_PT)
    

    # delta_R = sqrt(delta_phi^2 + delta_eta^2)
    delta_R_per_event = []
    for i in range(len(delta_eta)):
        delta_R_value = np.sqrt(delta_phi[i]**2 + delta_eta[i]**2)
        delta_R_per_event.append(delta_R_value)

    return delta_R_per_event

def calculate_epsilon(x_data_1, y_data_1, x_data_2, y_data_2, t_cut, sum_direction):
    if sum_direction == "1":
        boxes_b = y_data_1[t_cut:]
        epsilon_b = sum(boxes_b)
        boxes_s = y_data_2[:t_cut]
        epsilon_s = sum(boxes_s)
        
    else:
        boxes_b = y_data_2[:t_cut]
        epsilon_b = sum(boxes_b)
        boxes_s = y_data_1[t_cut:]
        epsilon_s = sum(boxes_s)

    return epsilon_b, epsilon_s

def signal_efficiency(x_data_file1, y_data_file1, x_data_file2, y_data_file2, num_events): # NB!: bh = file1, sph = file2
    # f(x|H_0)_0 and f(x|H_1)_1
    # where x represents the value from x-axis (currently either meff or HT)
    # H_0 is rejected signal and H_1 is wanted signal (either sphaleron or BH, or opposite).
    # t_cut is the cut point on the graph chosen and represents where we either start or end the riemann sum.
    # to find epsilon_b = area of f_0, and epsilon_s = area of f_1. W determines where the interval stops/starts.
    # b = epsilon_b * N_b and s = epsilon_s * N_s (N is total number of events from each compared data)
    # signal efficiency formula: optimal_selection = s / sqrt(s + b)

    # don't need to iterate to the edges since t_cut will most likely be somewhere in the middle.
    # also prevents index error.

    ### NB: s / (s + b) zero_to_tcut + b / (b + s) t_cut_to_end = 1 doesn't happen. s of file 1 =/= s of file 2?

    t = max(len(x_data_file1), len(x_data_file2))
    t_cut1 = 0
    t_cut2 = 0
    N_b = num_events[0]
    N_s = num_events[1]
    signal_eff = 0
    signal_eff_list = []
    optimal_t_cut = 0
    optimal_t_list = []
    # lists for plotting significance
    signal_efficiencies_list = []
    y_value_s_b_list = []
    x_values_list = []
    x_value_s_b_list = []
    print("y_data1: ", y_data_file1)
    print("y_data2: ", y_data_file2)
    # turn y_data into pdfs
    y_data_1_pdf = [float(y) / np.sum(y_data_file1) for y in y_data_file1]
    y_data_2_pdf = [float(y) / np.sum(y_data_file2) for y in y_data_file2] 
    print(y_data_1_pdf)
    print(y_data_2_pdf)
    # choose summing direction
    sum_direction = None
    while sum_direction not in ['1', '2']:
        print("sum from:")
        print("1: zero to t-cut")
        print("2: t-cut to end")
        sum_direction = raw_input("Selected integer corrosponding to your choice: ")
        if sum_direction not in ['1', '2']:
            print("Wrong value. Try again")
    # pick t_cut
    
    for i, x in enumerate(x_data_file1):
        print(i, x)
    t_cut = int(raw_input("select t_cut: "))
    t_cut1 = t_cut
    # sum all events from selected direction to t_cut
    ### NB: Choose blackhole file as first!!!
    if sum_direction == "1":
        blackhole_events1_a = sum(y_data_file1[:t_cut])
        sphaleron_events1_a = sum(y_data_file2[t_cut:])
        blackhole_events1_r = sum(y_data_file1[t_cut:])
        sphaleron_events1_r = sum(y_data_file2[:t_cut])
    else:
        blackhole_events1_a = sum(y_data_file1[t_cut:])
        sphaleron_events1_a = sum(y_data_file2[:t_cut])
        blackhole_events1_r = sum(y_data_file1[:t_cut])
        sphaleron_events1_r = sum(y_data_file2[t_cut:])

    #
    epsilon_b, epsilon_s = calculate_epsilon(x_data_file1, y_data_1_pdf, x_data_file2, y_data_2_pdf, t_cut, sum_direction)
    b = epsilon_b * N_b
    s = epsilon_s * N_s
    
    current_signal_eff = (s / np.sqrt(s + b))
    # get highest signal_eff method
    if signal_eff < current_signal_eff:
        signal_eff = current_signal_eff
        optimal_t_cut = x_data_file1[t_cut]

    signal_efficiencies_list.append(current_signal_eff)
    x_values_list.append(x_data_file1[t_cut])

    y_value_s_b_list.append(signal_efficiencies_list)
    x_value_s_b_list.append(x_values_list)
    signal_eff_list.append(signal_eff)
    optimal_t_list.append(optimal_t_cut)


    ###### sum for opposite s and b ######
    # invert direction
    if sum_direction == "1":
        sum_direction = "2"
    else: 
        sum_direction = "1"

    signal_efficiencies_list = []
    x_values_list = []
    signal_eff = 0
    optimal_t_cut = 0
    
    for i, x in enumerate(x_data_file1):
        print(i, x)
    t_cut = int(raw_input("select t_cut: "))
    t_cut2 = t_cut
    # sum all events from selected direction to t_cut
    if sum_direction == "2":
        blackhole_events2_a = sum(y_data_file1[t_cut:])
        sphaleron_events2_a = sum(y_data_file2[:t_cut])
        blackhole_events2_r = sum(y_data_file1[:t_cut])
        sphaleron_events2_r = sum(y_data_file2[t_cut:])
    else:
        blackhole_events2_a = sum(y_data_file1[:t_cut])
        sphaleron_events2_a = sum(y_data_file2[t_cut:])
        blackhole_events2_r = sum(y_data_file1[t_cut:])
        sphaleron_events2_r = sum(y_data_file2[:t_cut])
    #
    epsilon_b, epsilon_s = calculate_epsilon(x_data_file1, y_data_1_pdf, x_data_file2, y_data_2_pdf, t_cut, sum_direction)
    b = epsilon_b * N_b
    s = epsilon_s * N_s
    
    current_signal_eff = (s / np.sqrt(s + b))
    # get highest signal_eff method
    if signal_eff < current_signal_eff:
        signal_eff = current_signal_eff
        optimal_t_cut = x_data_file1[t_cut]

    signal_efficiencies_list.append(current_signal_eff)
    x_values_list.append(x_data_file1[t_cut])
    #######

    y_value_s_b_list.append(signal_efficiencies_list)
    x_value_s_b_list.append(x_values_list)
    signal_eff_list.append(signal_eff)
    optimal_t_list.append(optimal_t_cut)
    print(signal_eff_list)
    print(optimal_t_list)

    print("Signal efficiency given t_cut: {}".format(optimal_t_cut))
    
    print(" t_cut 1:")
    print(" S/B accept: {} / {}".format(blackhole_events1_a, sphaleron_events1_r))
    print(" S/B reject: {} / {}".format(sphaleron_events1_a, blackhole_events1_r))
    print(" S/B total: {} / {}".format(sphaleron_events1_a + sphaleron_events1_r, blackhole_events1_a + blackhole_events1_r))
    print(" uncertainty of t_cut 1 based on 3 bins from chosen cut: {} +- {}".format(x_data_file1[t_cut1], abs(x_data_file1[t_cut1] - x_data_file1[t_cut1 + 3])))
    print(" t_cut 2:")
    print(" S/B accept: {} / {}".format(blackhole_events2_a, sphaleron_events2_r))
    print(" S/B reject: {} / {}".format(sphaleron_events2_a, blackhole_events2_r))
    print(" S/B total: {} / {}".format(sphaleron_events2_a + sphaleron_events2_r, blackhole_events2_a + blackhole_events2_r))
    print(" uncertainty of t_cut2 based on 3 bins from chosen cut: {} +- {}".format(x_data_file1[t_cut2], abs(x_data_file1[t_cut2] - x_data_file1[t_cut2 + 3])))
    print("Signal_eff =", signal_eff)
    return signal_eff_list, optimal_t_list, y_value_s_b_list, x_value_s_b_list
