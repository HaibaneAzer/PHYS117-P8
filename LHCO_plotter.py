from LHCO_reader import LHCO_reader
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from collections import OrderedDict

### NB: this code is only compatible with python 2.x (solely used to read LHCO files, which its
### known library "LHCO_reader" has stopped its support since 2015).

####### path #########
def read_events_from_file(file_path):
    """ 
    Reads .lhco files using :mod:`LHCO_reader` and
    converts them to :func:`Events` objects which can
    be used to produce tables, plots and calculate
    useful values. 
    """
    try:
        events = LHCO_reader.Events(f_name=file_path)
        return events
    except Exception as e:
        print("Error reading file {}: {}".format(file_path, e))
        return None

def validate_events(events):
    """
    Checks if object :class:`Events` has correct type. Mean't mostly to
    further validate file format.

    Args:
        events (LHCO_reader.Events): The LHCO events object.
    """
    if not isinstance(events, list):
        raise ValueError("Events must be a list")
    for event in events:
        if not isinstance(event, dict):
            raise ValueError("Each event must be a dictionary")
        # Additional validation checks for event data

################ MATH ###############

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
        for particle in objects_2: 
            if events[event_idx][particle]: # check if the object is in event
                for object_idx in range(len(events[event_idx][particle])):
                    if current_PT < events[event_idx][particle][object_idx]['PT']:
                        current_PT = events[event_idx][particle][object_idx]['PT']
                        current_phi_highest_PT = events[event_idx][particle][object_idx]['phi']
                        current_eta_highest_PT = events[event_idx][particle][object_idx]['eta']

    
        Largest_PT_per_event.append(current_PT)
        phi_per_event_L_PT.append(current_phi_highest_PT)
        eta_per_event_L_PT.append(current_eta_highest_PT)
    
    delta_eta = delta_eta_per_event(eta_met_per_event, eta_per_event_L_PT)
    delta_phi = delta_phi_per_event(phi_met_per_event, phi_per_event_L_PT)
    delta_R = delta_R_per_event(delta_phi, delta_eta)
    return Largest_PT_per_event, delta_phi, delta_R

def delta_phi_per_event(phi_met_per_event, phi_per_event_L_PT):
    # difference in phi between such objects and met

    delta_phi_per_event = []

    for i in range(len(phi_met_per_event)):
        delta_phi_value = abs(phi_per_event_L_PT[i] - phi_met_per_event[i])
        delta_phi_per_event.append(delta_phi_value)

    return delta_phi_per_event

def delta_eta_per_event(eta_met_per_event, eta_per_event_L_PT):
    delta_eta_per_event = []

    for i in range(len(eta_met_per_event)):
        delta_eta_value = abs(eta_per_event_L_PT[i] - eta_met_per_event[i])
        delta_eta_per_event.append(delta_eta_value)

    return delta_eta_per_event

def delta_R_per_event(delta_phi, delta_eta):
    # delta_R = sqrt(delta_phi^2 + delta_eta^2)
    delta_R_per_event = []
    for i in range(len(delta_eta)):
        delta_R_value = np.sqrt(delta_phi[i]**2 + delta_eta[i]**2)
        delta_R_per_event.append(delta_R_value)

    return delta_R_per_event

def calculate_F1_score(signal_eff, background_rejection):
    # recommended method to find optimal t_cut by chatgpt
    return 2 * (signal_eff * background_rejection) / (signal_eff + background_rejection)

def signal_efficiency(x_data_file1, y_data_file1, x_data_file2, y_data_file2, num_events): # NB!: b = file1, s = file2
    # f(x|H_0)_0 and f(x|H_1)_1
    # where x represents the value from x-axis (currently either meff or HT)
    # H_0 is rejected signal and H_1 is wanted signal (either sphaleron or BH, or opposite).
    # t_cut is the cut point on the graph chosen and represents where we either start or end the riemann sum.
    # to find epsilon_b = area of f_0, and epsilon_s = area of f_1. W determines where the interval stops/starts.
    # b = epsilon_b * N_b and s = epsilon_s * N_s (N is total number of events from each compared data)
    # signal efficiency formula: optimal_selection = s / sqrt(s + b)

    # don't need to iterate to the edges since t_cut will most likely be somewhere in the middle.
    # also prevents index error.
    t = min(len(x_data_file1), len(x_data_file2))
    N_b = num_events[0] # (might change for actual event number)
    N_s = num_events[1]
    signal_eff = 0
    optimal_t_cut = None
    max_F1_score = 0
    attempts = []
    print("x_data_file1:", x_data_file1)
    print("x_data_file2:", x_data_file2)
    # pick t_cut
    for t_cut in range(t-1):
        epsilon_b = sum((x_data_file1[t_cut+1] - x_data_file1[t_cut]) * y_data_file1[:t_cut])
        epsilon_s = sum((x_data_file2[t_cut+1] - x_data_file2[t_cut]) * y_data_file2[:t_cut])

        b = epsilon_b * N_b
        s = epsilon_s * N_s
        current_signal_eff = (s / np.sqrt(s + b))
        background_rejection = 1 - b / N_b

        F1_score = calculate_F1_score(current_signal_eff, background_rejection)
        # f1 score method
        """ if F1_score > max_F1_score:
            max_F1_score = F1_score
            optimal_t_cut = x_data_file1[t_cut]
            signal_eff = current_signal_eff
            attempts.append(optimal_t_cut) """
        
        # get highest signal_eff method
        if signal_eff < current_signal_eff:
            signal_eff = current_signal_eff
            optimal_t_cut = x_data_file1[t_cut]
            attempts.append(optimal_t_cut) # find when signal efficiency is last defined (aka. when max)
        print("current signal efficiencies:", current_signal_eff)

    print("Signal efficiency given t_cut: {}".format(optimal_t_cut))
    print("all attempts:", attempts)
    print("Signal_eff =", signal_eff)
    return signal_eff, optimal_t_cut

    
#####################################

#### prompts #####

def plot_selection_prompt():
    """
    Prompt user to select the type of plots they want to generate.

    Returns:
        str: The selected plot type.
    """
    print("Select the type of plots:")
    print("1. Particle Properties (PT, phi, eta)")
    print("2. HT and signal_efficiency")
    print("3. meff and signal_efficiency")
    print("4. Objects per Event (histogram)")
    print("5. Electrons, Muons, Taus per Event (histrogram)")
    print("6. largest PT and delta_phi per Event")
    print("7: delta_R per event")
    
    choice = raw_input("Enter the number corresponding to your choice: ")

    if choice == '1':
        return 'particle_properties'
    elif choice == '2':
        return 'HT'
    elif choice == "3":
        return 'meff'
    elif choice == '4':
        return 'objects_per_event'
    elif choice == '5':
        return 'electrons_muons_taus_per_event'
    elif choice == '6':
        return 'largest_PT_in_event'
    elif choice == '7':
        return 'delta_R_per_event'
    else:
        print("Invalid choice. Please enter a valid option.")
        return plot_selection_prompt()

def locate_file_prompter():
    """
    Uses path and directory information to locate necessary files
    NB: Important that the files you want to use is located in the 
    same directory as this code.

    Returns:
        str: Path of subdirectory.
        list(str): a list of files names from path.
        list(str): a list of indicies for corrosponding file in file list
    """
    # define path and file name
    data_path = os.getcwd()

    # get list of all subdirectories in the data directory
    subdirectories = next(os.walk(data_path))[1]
    subdirectories.remove('.git') # ignore git folder

    subdir_dict = OrderedDict() # necessary for 3.6 and lower python versions
    
    # do new prompt
    file_list = []
    for idx in range(len(subdirectories)):
        selected_subdir = subdirectories[idx]
        # List files from all subdirectory
        subdir_path = os.path.join(data_path, selected_subdir)
        subdir_file_list = os.listdir(subdir_path)

        # store in dict
        subdir_dict[subdir_path] = subdir_file_list

        # selected_subdir is all
        selected_subdir = "all files"
        for files in subdir_file_list:
            file_list.append(files)
    
    # Prompt user to select files
    print("\nSelect files from {}:".format(selected_subdir))
    for idx, file_name in enumerate(file_list):
        print("{}. {}".format(idx + 1, file_name))

    # create file_list including both path and file name
    file_path_list = create_file_path_list(subdir_dict)

    # handle input exceptions
    try:
        # compatibility check
        if (sys.version_info[0] == 2):
            # python 2.x
            selected_files = raw_input("Enter the indices of the files (comma-separated): ")
        else:
            # python 3.x
            selected_files = input("Enter the indices of the files (comma-separated): ")
    except ValueError:
        print("Invalid input. Please enter a valid line of comma-separated integers.")
        return locate_file_prompter() # prompt retry
    
    selected_files = selected_files.split(',')
    selected_files = [int(idx.strip()) - 1 for idx in selected_files]
    
    return selected_files, file_path_list

def create_file_path_list(subdir_dict):
    file_path_list = []

    for subdir_path, file_names in subdir_dict.items():
        for file_name in file_names:
            file_path = os.path.join(subdir_path, file_name)
            file_path_list.append(file_path)

    return file_path_list

#################

def select_particle_and_property():
    particles = ['electron', 'jet', 'muon', 'tau', 'MET']
    props = ['PT', 'phi', 'eta']

    print("\nSelect particle to look at")
    for idx, particle in enumerate(particles):
        print("{}. {}".format(idx + 1, particle))

    particle_idx = int(input("Enter valid index of particle: ")) - 1
    if not (0 <= particle_idx < len(particles)):
        print("Invalid index. Please try again.")
        return select_particle_and_property()

    print("\nSelect property to view for {}".format(particles[particle_idx]))
    for idx, prop in enumerate(props):
        print("{}. {}".format(idx + 1, prop))

    prop_idx = int(input("Enter valid index of property: ")) - 1
    if not (0 <= prop_idx < len(props)):
        print("Invalid index. Please try again.")
        return select_particle_and_property()

    return particles[particle_idx], props[prop_idx]

def process_selected_file(file_path):
    events = read_events_from_file(file_path)

    if events is not None:
        validate_events(events)
        return events
    else:
        return None

def processed_file_to_data(file_path_list, selected_files, plot_type):
    legend_labels = []
    x_data = []
    y_data = []
    num_events_list = []
    signal_eff = None
    t_cut_optimal = None

    particle_name, prop_name = select_particle_and_property()
    for idx in selected_files:
        if 0 <= idx < len(file_path_list):
            file_name = file_path_list[idx]
            events = process_selected_file(file_name)
            num_events_list.append(len(events))

            if events is not None:
                x_data, y_data = process_data_for_plot(events, x_data, y_data, plot_type, particle_name, prop_name)
                legend_labels.append(file_name)
    # get signal efficiency if only 2 files and HT, meff data is analyzed
    if (len(selected_files) == 2) and (plot_type in ['HT', 'meff']):
        signal_eff, t_cut_optimal = signal_efficiency(x_data[0], y_data[0], x_data[1], y_data[1], num_events_list) # file 1, file 2

    return x_data, y_data, legend_labels, particle_name, prop_name, signal_eff, t_cut_optimal

def data_to_bincenter_y_norm(data, bins):
    y, binEdges = np.histogram(data, bins=bins)
    bincenters = 0.5*(binEdges[1:] + binEdges[:-1])
    y_norm = y / np.sum(data)
    
    return bincenters, y_norm

def data_to_bincenter_histogram(data, bins):
    """
    Calculates histogram and returns bin centers and histogram values.

    Args:
        data (list): List of data points.
        num_bins (int): Number of bins for the histogram.

    Returns:
        tuple: Tuple containing bin centers and histogram values.
    """
    bin_edges = np.linspace(min(data), max(data), bins+1)
    hist, _ = np.histogram(data, bins=bin_edges, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, hist

def process_data_for_plot(events, x_data, y_data, plot_type, particle_name, prop_name):
    if plot_type == 'HT':
        HT_list, _ = calculate_HT_and_meff(events)
        bincenters_HT, normalized_y_HT = data_to_bincenter_y_norm(HT_list, 50)

        x_data.append(bincenters_HT)
        y_data.append(normalized_y_HT)
    elif plot_type == 'meff':
        _, meff_list = calculate_HT_and_meff(events)
        bincenters_meff, normalized_y_meff = data_to_bincenter_y_norm(meff_list, 50)

        x_data.append(bincenters_meff)
        y_data.append(normalized_y_meff)
    elif plot_type == 'objects_per_event':
        obj_per_event_data = objects_per_event(events)

        x_data.append(list(range(1, len(obj_per_event_data) + 1)))  # Event numbers
        y_data.append(obj_per_event_data)
    elif plot_type == 'particle_properties':
        data = events.column(particle_name, prop_name)  
        bincenters, normalized_y = data_to_bincenter_y_norm(data, 50)

        x_data.append(bincenters)
        y_data.append(normalized_y)
    elif plot_type == 'electrons_muons_taus_per_event':
        electrons, muons, taus = electrons_muons_taus_per_event(events)

        x_data.append(list(range(1, len(electrons) + 1)))
        y_data.append(electrons)
        x_data.append(list(range(1, len(muons) + 1)))
        y_data.append(muons)
        x_data.append(list(range(1, len(taus) + 1)))
        y_data.append(taus)
    elif plot_type == 'largest_PT_in_event':
        largest_PT, delta_phi, _ = largest_PT_in_event(events)
        
        # Divide the events into bins
        bin_size = len(largest_PT) // 100
        binned_events = np.arange(0, len(largest_PT), bin_size)
        binned_PT = [np.mean(largest_PT[i:i+bin_size]) for i in binned_events]
        binned_phi = [np.mean(delta_phi[i:i+bin_size]) for i in binned_events]

        x_data.append([binned_events, binned_events])
        y_data.append([binned_PT, binned_phi])
    elif plot_type == 'delta_R_per_event':
        largest_PT, _, delta_R = largest_PT_in_event(events)
         # Divide the events into bins
        bincenters_R, y_normalized = data_to_bincenter_y_norm(delta_R, 50)
    
        x_data.append(bincenters_R)
        y_data.append(y_normalized)
    return x_data, y_data

def plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name, plot_type, signal_eff, t_cut_optimal):
    """
    plot all data into one line plot 

    Args:
        x_data (list(float)): List containing x-coordinates.
        y_data (list(float)): List containing y-coordinates.
        legend_labels (list(str)): List containing file names
        p_name (str): The particle name.
        p_prop (str): The particle property.
    """
    if plot_type == 'electrons_muons_taus_per_event':
        fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

        particles = ['Electrons', 'Muons', 'Taus']
        for i, ax in enumerate(axes):
            # Count the occurrences of each value in y_data[i]
            counts = {}
            for value in y_data[i]:
                if value in counts:
                    counts[value] += 1
                else:
                    counts[value] = 1

            # Separate the counts and values
            values = list(counts.keys())
            frequencies = list(counts.values())

            # Create a bar plot
            ax.bar(values, frequencies, align='center')
            ax.set_title('{} per Event'.format(particles[i]))
            ax.set_ylabel('Frequency')

        plt.xlabel('Total Particles')
        plt.tight_layout()
        plt.show()
    elif plot_type == 'largest_PT_in_event':
        # New code for binned events line plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))

        # Extract file names from legend labels
        legend_labels = [os.path.basename(label) for label in legend_labels]

        # Make plot for each data file
        for label, (x_PT, x_phi), (y_PT, y_phi) in zip(legend_labels, x_data, y_data):
            ax1.plot(x_PT, y_PT, '-', label=label)
            ax2.plot(x_phi, y_phi, '-', label=label)
        
        ax1.grid()
        ax1.set_xlabel('Events')
        ax1.set_ylabel('Largest PT (GeV)')
        ax1.legend()

        ax2.grid()
        ax2.set_xlabel('Event')
        ax2.set_ylabel('Delta Phi (angle)')
        ax2.legend()

        plt.tight_layout()
        plt.show()
    elif plot_type == 'delta_R_per_event':
        # New code for binned events line plot
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        # Extract file names from legend labels
        legend_labels = [os.path.basename(label) for label in legend_labels]

        # Make plot for each data file
        for label, x_R, y_R in zip(legend_labels, x_data, y_data):
            ax1.plot(x_R, y_R, '-', label=label)
        
        ax1.grid()
        ax1.set_xlabel('Events')
        ax1.set_ylabel('Delta R (angle)')
        ax1.legend()

        plt.tight_layout()
        plt.show()
    elif plot_type == 'objects_per_event':
        # Count the occurrences of each value in y_data
        counts = {}
        for value in y_data[0]:
            if value in counts:
                counts[value] += 1
            else:
                counts[value] = 1

        # Separate the counts and values
        values = list(counts.keys())
        frequencies = list(counts.values())

        # Create a bar plot
        plt.bar(values, frequencies, align='center')

        plt.xlabel('Total Objects')
        plt.ylabel('Frequency')
        plt.title('Objects per Event')

        plt.show()
    else: # all other plot types as line graph
        fig = plt.figure()
        ax = fig.add_subplot(111)
        unit = ""

        # Extract file names from legend labels
        legend_labels = [os.path.basename(label) for label in legend_labels]
        # make plot for each data file
        for label, x, y in zip(legend_labels, x_data, y_data):
            ax.plot(x, y, '-', label=label)
            if signal_eff and t_cut_optimal != None:
                ax.axvline(x=t_cut_optimal, color='r', linestyle='--', label='T_cut Optimal ({}), with efficiency ({})'.format(t_cut_optimal, signal_eff))
        
        ax.grid()

        if plot_type == 'particle_properties':
            ax.set_title('{} {}'.format(particle_name, prop_name))
        elif plot_type in ['HT', 'meff']:
            ax.set_title('{}'.format(particle_name))

        # give appropiate unit scale
        if prop_name == 'PT':
            unit = "(GeV)"
        elif prop_name == 'phi':
            unit = "(angle)"

        if plot_type == 'HT':
            unit = "(GeV)"
        elif plot_type == 'meff':
            unit = "(GeV)"

        # set y and x labels
        if plot_type in ['HT', 'meff', 'particle_properties']:
            ax.set_xlabel('{} {}'.format(prop_name, unit))
            ax.set_ylabel('Relative Frequency')
        ax.legend()

        plt.show()

def main():
    selected_files, file_list = locate_file_prompter()

    plot_type = plot_selection_prompt()

    if plot_type == 'particle_properties':
        x_data, y_data, legend_labels, particle_name, prop_name = processed_file_to_data(file_list, selected_files, plot_type)
        plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name, plot_type, None, None)
    elif plot_type in ['HT', 'meff', 'objects_per_event', 'electrons_muons_taus_per_event', 'largest_PT_in_event', 'delta_R_per_event']:
        x_data, y_data, legend_labels, particle_name, prop_name, signal_eff, t_cut_optimal = processed_file_to_data(file_list, selected_files, plot_type)
        plot_line_graph(x_data, y_data, legend_labels, plot_type, "", plot_type, signal_eff, t_cut_optimal)

# run
if __name__ == "__main__":
    main()
