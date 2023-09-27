from LHCO_reader import LHCO_reader
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

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

def LHCO_plot(events, file, x_data, y_data, p_name, p_prop):
    """
    Calculate and append data points for plotting.

    Args:
        events (LHCO_reader.Events): The LHCO events object.
        file (str): The name of the file being processed.
        x_data (list): List to store x-coordinates.
        y_data (list): List to store y-coordinates.
        p_name (str): The particle name.
        p_prop (str): The particle property.

    Returns:
        str: The name of the file.
    """
    data = events.column(p_name, p_prop)  # Make data into a column

    # create x/y-coordinates using data
    y, binEdges = np.histogram(data, bins=50) # number of x-axis-segments 
    bincenters = 0.5*(binEdges[1:] + binEdges[:-1])

    # normalize y
    normalized_y = y / np.sum(data)

    x_data.append(bincenters)
    y_data.append(normalized_y)
    return file

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

    # Prompt user to select a subdirectory
    print("Select a subdirectory:")
    for idx, subdir in enumerate(subdirectories):
        print("{}. {}".format(idx + 1 , subdir))

    # handle input exceptions
    try:
        subdir_idx = int(input("Enter the index of the subdirectory: ")) - 1
        if subdir_idx < 0 or subdir_idx >= len(subdirectories):
            raise ValueError("Invalid subdirectory index")
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return locate_file_prompter() # prompt retry
    
    selected_subdir = subdirectories[subdir_idx]

    # List files in the selected subdirectory
    subdir_path = os.path.join(data_path, selected_subdir)
    file_list = os.listdir(subdir_path)

    # Prompt user to select files
    print("\nSelect files from {}:".format(selected_subdir))
    for idx, file_name in enumerate(file_list):
        print("{}. {}".format(idx + 1, file_name))

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
    
    return subdir_path, file_list, selected_files

def file_processer(subdir_path, file_list, selected_files):
    """
    Turns file data into :mod:`matplotlib` variables that can be used
    in :func:`plot`.

    Args:
        subdir_path (str): path of subdirectory.
        file_list (list(str)): list of file names found from path.
        selected_files (list(int)): list of indicies for corrosponding file from file list.
    
    Return:
        x_data (list(float)): List containing x-coordinates.
        y_data (list(float)): List containing store y-coordinates.
        legends_labels (list(str)): List containing file names
        particles[particle_idx] (str): Name of particle
        props[prop_idx] (str): name of property
    """
    # select particle and prop, then add to list
    particles = [
        'electron',
        'jet',
        'muon',
        'tau',
        'MET',
    ]
    props = [
        'PT',
        'phi',
        'eta',
    ]

    # print list of valid particles
    print("\nSelect particle to look at")
    for idx, particle in enumerate(particles):
        print("{}. {}".format(idx + 1, particle))
    
    # handle prompt exceptions
    try:
        particle_idx = int(input("Enter valid index of particle: ")) - 1
        if particle_idx < 0 or particle_idx >= len(particles):
            raise ValueError("Invalid particle-list index")
    except ValueError:
        print("Invalid input. Please enter a valid index")
        return file_processer(subdir_path, file_list, selected_files)
    # print list of valid properties
    print("\nSelect properties to {}".format(particles[particle_idx]))
    for idx, prop in enumerate(props):
        print("{}. {}".format(idx + 1, prop))
    
    # handle prompt exceptins
    try:
        prop_idx = int(input("Enter valid prop name of particle: ")) - 1
        if prop_idx < 0 or prop_idx >= len(props):
            raise ValueError("Invalid props-list index")
    except ValueError:
        print("Invalid input. Please enter a valid index")
        return file_processer(subdir_path, file_list, selected_files)

    # initialize lists to store data
    x_data = []
    y_data = []

    legend_labels = []
    # Process selected files
    for idx in selected_files:
        if 0 <= idx < len(file_list):
            file_name = file_list[idx]
            file_path = os.path.join(subdir_path, file_name)

            # Read and convert file to event object
            events = read_events_from_file(file_path)
            if events is not None:
                # Validate events
                validate_events(events)

                # Plot
                legend_label = LHCO_plot(events, file_name, x_data, y_data, particles[particle_idx], props[prop_idx])
                legend_labels.append(legend_label)
    return x_data, y_data, legend_labels, particles[particle_idx], props[prop_idx]

def plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name):
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

    # make plot for each data file
    for x, y, label in zip(x_data, y_data, legend_labels):
        ax.plot(x, y, '-', label=label)
        
    ax.grid()
    ax.set_title('{} {}'.format(particle_name, prop_name))

    # give appropiate unit scale
    if prop_name == 'PT':
        unit = "(GeV)"
    elif prop_name == 'phi':
        unit = "(angle)"

    ax.set_xlabel('{} {}'.format(prop_name, unit))
    ax.set_ylabel('Relative Frequency')
    ax.legend()

    plt.show()

def main():
    subdir_path, file_list, selected_files = locate_file_prompter()
    x_data, y_data, legend_labels, particle_name, prop_name = file_processer(subdir_path, file_list, selected_files)
    plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name)

# run
if __name__ == "__main__":
    main()
