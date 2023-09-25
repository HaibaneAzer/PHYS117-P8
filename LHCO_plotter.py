from LHCO_reader import LHCO_reader
import numpy as np
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

def LHCO_plot(events, file, x_data, y_data):
    data = events.column('jet', 'PT')  # Make data into a column

    # create x/y-coordinates using data
    y, binEdges = np.histogram(data, bins=50) # number of x-axis-segments 
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

    # normalize y
    normalized_y = y / np.sum(data)

    x_data.append(bincenters)
    y_data.append(normalized_y)

    return file

def main():
    # define path and file name
    data_path = os.getcwd()

    # get list of all subdirectories in the data directory
    subdirectories = next(os.walk(data_path))[1]
    subdirectories.remove('.git') # ignore git folder

    # Prompt user to select a subdirectory
    print("Select a subdirectory:")
    for idx, subdir in enumerate(subdirectories):
        print("{}. {}".format(idx+1 , subdir))

    subdir_idx = int(input("Enter the index of the subdirectory: ")) - 1
    selected_subdir = subdirectories[subdir_idx]

    # List files in the selected subdirectory
    subdir_path = os.path.join(data_path, selected_subdir)
    file_list = os.listdir(subdir_path)

    # Prompt user to select files
    print("\nSelect files from {}:".format(selected_subdir))
    for idx, file_name in enumerate(file_list):
        print("{}. {}".format(idx+1, file_name))

    selected_files = raw_input("Enter the indices of the files (comma-separated): ") # 2.x python only
    selected_files = selected_files.split(',')
    selected_files = [int(idx.strip()) - 1 for idx in selected_files]

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
                legend_label = LHCO_plot(events, file_name, x_data, y_data)
                legend_labels.append(legend_label)

    # plot all data into one graph
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for x, y, label in zip(x_data, y_data, legend_labels):
        ax.plot(x, y, '-', label=label)
        
    ax.grid()
    ax.set_title('jet PT')
    ax.set_xlabel('PT (GeV)')
    ax.set_ylabel('Relative Frequency')
    ax.legend()

    plt.show()

# run
if __name__ == "__main__":
    main()
