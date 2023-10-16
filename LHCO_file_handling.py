from collections import OrderedDict as OD
from LHCO_reader import LHCO_reader
import os
import sys
import numpy as np

from LHCO_comp_functions import (
    calculate_HT_and_meff,
    objects_per_event,
    electrons_muons_taus_per_event,
    largest_PT_in_event,
    delta_R_per_event,
    signal_efficiency
)

print("path: ", os.getcwd())

### main file handling functions ### 
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

    subdir_dict = OD() # necessary for 3.6 and lower python versions
    
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
    
def processed_file_to_data(file_path_list, selected_files, plot_type):
    legend_labels = []
    x_data = []
    y_data = []
    num_events_list = []
    signal_eff = None
    t_cut_optimal = None
    signal_list = None
    t_list = None

    particle_name, prop_name = select_particle_and_property(plot_type)
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
        signal_eff, t_cut_optimal, signal_list, t_list = signal_efficiency(x_data[0], y_data[0], x_data[1], y_data[1], num_events_list) # file 1, file 2

    return x_data, y_data, legend_labels, particle_name, prop_name, signal_eff, t_cut_optimal, signal_list, t_list


### main prompt selection functions ###

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

def select_particle_and_property(plot_type):
    # NB! update this function for other parts of code aswell!!!
    particles = ['electron', 'jet', 'muon', 'tau', 'MET']
    props = ['PT', 'phi', 'eta']
    if plot_type == 'delta_R_per_event':
        print("\nSelect particle to look at")
        for idx, particle in enumerate(particles):
            print("{}. {}".format(idx + 1, particle))

        particle_idx = int(input("Enter valid index of particle: ")) - 1
        if not (0 <= particle_idx < len(particles)):
            print("Invalid index. Please try again.")
            return select_particle_and_property()
        # we look at highest PT
        prop_idx = 0     
    else:
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


### helper functions ###
### used in module only ###

def process_selected_file(file_path):
    events = read_events_from_file(file_path)
    if events is not None:
        validate_events(events)
        return events
    else:
        return None
    
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

        x_data.append([list(range(1, len(electrons) + 1)), list(range(1, len(muons) + 1)), list(range(1, len(taus) + 1))])
        y_data.append([electrons, muons, taus])
    elif plot_type == 'largest_PT_in_event':
        largest_PT, delta_phi = largest_PT_in_event(events)
        
        # Divide the events into bins
        bin_size = len(largest_PT) // 100
        binned_events = np.arange(0, len(largest_PT), bin_size)
        binned_PT = [np.mean(largest_PT[i:i+bin_size]) for i in binned_events]
        binned_phi = [np.mean(delta_phi[i:i+bin_size]) for i in binned_events]

        x_data.append([binned_events, binned_events])
        y_data.append([binned_PT, binned_phi])
    elif plot_type == 'delta_R_per_event':
        delta_R = delta_R_per_event(events, particle_name)
         # Divide the events into bins
        bincenters_R, y_normalized = data_to_bincenter_y_norm(delta_R, 50)
    
        x_data.append(bincenters_R)
        y_data.append(y_normalized)
    return x_data, y_data

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

def create_file_path_list(subdir_dict):
    file_path_list = []

    for subdir_path, file_names in subdir_dict.items():
        for file_name in file_names:
            file_path = os.path.join(subdir_path, file_name)
            file_path_list.append(file_path)

    return file_path_list

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
