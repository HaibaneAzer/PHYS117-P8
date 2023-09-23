from LHCO_reader import LHCO_reader
import os
import pandas as pd
# a

def read_events_from_file(file_path):
    try:
        events = LHCO_reader.Events(f_name=file_path)
        return events
    except Exception as e:
        print("Error reading file: {}".format(e))
        return None

def process_events(events):
    # Determine the number of events
    num_events = len(events)

    # Define the initial dictionary with empty lists
    d = {
        'Events': [],
        'jet PT': [],
        'jet Phi': [],
        'jet eta': [],
    }

    # Loop through the events and update the values
    for n in range(num_events):    
        event_label = 'e{}'.format(n+1) 
        d['Events'].append(event_label)
        if events[n]["electron"]:
            d['jet PT'].append(events[n]["electron"][0]['PT'])
            d['jet Phi'].append(events[n]["electron"][0]['phi'])
            d['jet eta'].append(events[n]["electron"][0]['eta'])
        else:
            d['jet PT'].append('-')
            d['jet Phi'].append('-')
            d['jet eta'].append('-')

    # Create the DataFrame
    table = pd.DataFrame(data=d)
    print(table)

def validate_events(events):
    if not isinstance(events, list):
        raise ValueError("Events must be a list")
    for event in events:
        if not isinstance(event, dict):
            raise ValueError("Each event must be a dictionary")
        # Additional validation checks for event data


def main(idx, jdx):
    # define folder and file

    # locate file with index (idx = folder, jdx = file)
    idx = 0
    jdx = 0

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
    print("filepath:", file_path) # current path

    print("file:", str(file_list[jdx])) # current file

    # read and convert file to event object
    events = read_events_from_file(file_path)

    if events is not None:
        # validate events
        validate_events(events)
        
        # run processes
        process_events(events)

# run
if __name__ == "__main__":
    folder = int(input("folder (int): "))
    file = int(input("file (int): "))
    main(folder, file)