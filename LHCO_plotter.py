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
    tot_obj = 0
    obj_per_event = []
    for event in events:
        objects_event = event.number()
        for particle_name in particles:
            tot_obj += objects_event[particle_name]
        obj_per_event.append(tot_obj)
        tot_obj = 0
    
    return obj_per_event

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
    print("2. HT")
    print("3. meff")
    print("4. Objects per Event")
    choice = raw_input("Enter the number corresponding to your choice: ")

    if choice == '1':
        return 'particle_properties'
    elif choice == '2':
        return 'HT'
    elif choice == "3":
        return 'meff'
    elif choice == '4':
        return 'objects_per_event'
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
    x_data = []
    y_data = []
    legend_labels = []

    particle_name, prop_name = select_particle_and_property()
    for idx in selected_files:
        if 0 <= idx < len(file_path_list):
            file_name = file_path_list[idx]
            events = process_selected_file(file_name)

            if events is not None:
                process_data_for_plot(events, plot_type, x_data, y_data, particle_name, prop_name)
                legend_labels.append(file_name)

    return x_data, y_data, legend_labels, particle_name, prop_name

def process_data_for_plot(events, plot_type, x_data, y_data, particle_name, prop_name):

    if plot_type == 'HT':
        HT_list, _ = calculate_HT_and_meff(events)
        y_HT, binEdges_HT = np.histogram(HT_list, bins=50)
        bincenters_HT = 0.5*(binEdges_HT[1:] + binEdges_HT[:-1])
        normalized_y_HT = y_HT / np.sum(HT_list)
        x_data.append(bincenters_HT)
        y_data.append(normalized_y_HT)

    elif plot_type == 'meff':
        _, meff_list = calculate_HT_and_meff(events)
        y_meff, binEdges_meff = np.histogram(meff_list, bins=50)
        bincenters_meff = 0.5*(binEdges_meff[1:] + binEdges_meff[:-1])
        normalized_y_meff = y_meff / np.sum(meff_list)
        x_data.append(bincenters_meff)
        y_data.append(normalized_y_meff)

    elif plot_type == 'objects_per_event':
        obj_per_event_data = objects_per_event(events)
        x_data.append(list(range(1, len(obj_per_event_data) + 1)))  # Event numbers
        y_data.append(obj_per_event_data)
    
    elif plot_type == 'particle_properties':
        data = events.column(particle_name, prop_name)  
        y, binEdges = np.histogram(data, bins=50)
        bincenters = 0.5*(binEdges[1:] + binEdges[:-1])
        normalized_y = y / np.sum(data)
        x_data.append(bincenters)
        y_data.append(normalized_y)
    

    return x_data, y_data

def plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name, plot_type):
    """
    plot all data into one line plot 

    Args:
        x_data (list(float)): List containing x-coordinates.
        y_data (list(float)): List containing y-coordinates.
        legend_labels (list(str)): List containing file names
        p_name (str): The particle name.
        p_prop (str): The particle property.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    unit = ""

    # Extract file names from legend labels
    legend_labels = [os.path.basename(label) for label in legend_labels]

    # make plot for each data file
    for label, x, y in zip(legend_labels, x_data, y_data):
        ax.plot(x, y, '-', label=label)
        
    ax.grid()

    if plot_type == 'particle_properties':
        ax.set_title('{} {}'.format(particle_name, prop_name))
    elif plot_type in ['HT', 'meff']:
        ax.set_title('{}'.format(particle_name))
    elif plot_type == 'objects_per_event':
        ax.set_title('{}'.format(plot_type))

    # give appropiate unit scale
    if prop_name == 'PT':
        unit = "(GeV)"
    elif prop_name == 'phi':
        unit = "(angle)"

    if plot_type == 'HT':
        unit = "(GeV)"
    elif plot_type == 'meff':
        unit = "(GeV)"
    elif plot_type == 'objects_per_event':
        unit = "(amount)"

    # set y and x labels
    if plot_type in ['HT', 'meff', 'particle_properties']:
        ax.set_xlabel('{} {}'.format(prop_name, unit))
        ax.set_ylabel('Relative Frequency')
    elif plot_type == 'objects_per_event':
        ax.set_xlabel('{} {}'.format('Events', unit))
        ax.set_ylabel('Total Objects')
    ax.legend()

    plt.show()

def main():
    selected_files, file_list = locate_file_prompter()

    plot_type = plot_selection_prompt()

    if plot_type == 'particle_properties':
        x_data, y_data, legend_labels, particle_name, prop_name = processed_file_to_data(file_list, selected_files, plot_type)
        plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name, plot_type)
    elif plot_type in ['HT', 'meff', 'objects_per_event']:
        x_data, y_data, legend_labels, particle_name, prop_name = processed_file_to_data(file_list, selected_files, plot_type)
        plot_line_graph(x_data, y_data, legend_labels, plot_type, "", plot_type)

# run
if __name__ == "__main__":
    main()
