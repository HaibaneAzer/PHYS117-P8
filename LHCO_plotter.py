from LHCO_reader import LHCO_reader
import numpy as np
import collections
import matplotlib.pyplot as plt
import os

def read_events_from_file(file_path):
    try:
        events = LHCO_reader.Events(f_name=file_path)
        return events
    except Exception as e:
        print("Error reading file: {}".format(e))
        return None
    

    

def validate_events(events):
    if not isinstance(events, list):
        raise ValueError("Events must be a list")
    for event in events:
        if not isinstance(event, dict):
            raise ValueError("Each event must be a dictionary")
        # Additional validation checks for event data

def LHCO_plot(events):
    data = events.column('jet', 'PT')  # Make data into a column

    fig = plt.figure()
    ax = fig.add_subplot(111)
    y, binEdges = np.histogram(data, bins=100)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    ax.plot(bincenters, y, '-')
    """ ax.hist(data, 50, normed=1, facecolor='Crimson', alpha=0.9) """
    ax.grid()  # Add grid lines

    ax.set_title('jet PT')  # Titles etc
    ax.set_xlabel('PT')
    ax.set_ylabel("Frequency")
    plt.rc('font', size=20)  # Fonts - LaTeX possible, but slow

    plt.show()  # Show the plot on screen

def main(idx, jdx):
    # define path and file name
    data_path = os.getcwd()

    subdirectories = next(os.walk(data_path))[1]
    subdirectories.remove('.git') # ignore git folder

    # error check
    if idx >= len(subdirectories):
        raise ValueError("Invalid idx value. There is only {} subdirectories".format(len(subdirectories)))
    
    folder_name = subdirectories[idx]
    folder_path = os.path.join(data_path, folder_name)

    file_list = os.listdir(folder_path)

    # error check
    if jdx >= len(file_list):
        raise ValueError("Invalid jdx value. There are only {} files in {}".format(len(file_list), folder_name))

    file_path = os.path.join(data_path, folder_name, file_list[jdx])
    """ print("filepath:", file_path) # current path """

    print("file:", str(file_list[jdx])) # current file

    # read and convert file to event object
    print()
    events = read_events_from_file(file_path)
    print()
    
    if events is not None:
        # validate events
        validate_events(events)

        # plot
        LHCO_plot(events)
        

# run
if __name__ == "__main__":
    print()
    folder = int(input("folder (int): "))
    file = int(input("file (int): "))
    print()
    main(folder, file)

    

    
