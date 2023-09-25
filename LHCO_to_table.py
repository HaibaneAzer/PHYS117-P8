from LHCO_reader import LHCO_reader
import matplotlib as plt
import os
import pandas as pd

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

    # Define a list of particle/unit names and their corresponding labels
    particles = [
        {'name': 'electron', 'label': 'electron'},
        {'name': 'jet', 'label': 'jet'},
        {'name': 'MET', 'label': 'MET'},
        {'name': 'muon', 'label': 'muon'},
        {'name': 'tau', 'label': 'tau'},
        # Add more particles/units as needed
    ]

    # Define the attributes (PT, phi, eta) you want to extract for each particle/unit
    attributes = ['PT', 'phi', 'eta']

    # Initialize the dictionary to store the data
    d = {'Events': []}

    # Add columns for each particle/unit and attribute combination
    for particle_info in particles:
        for attribute in attributes:
            column_name = '{} {}'.format(particle_info['label'], attribute)
            d[column_name] = []

    # Loop through the events and update the values
    # Loop through the events and update the values
    for n in range(num_events):    
        event_label = 'e{}'.format(n+1) 
        d['Events'].append(event_label)
    
        for particle_info in particles:
            particle_name = particle_info['name']
            particle_label = particle_info['label']
        
            if events[n][particle_name]:
                for attribute in attributes:
                    value = events[n][particle_name][0][attribute]
                    d['{} {}'.format(particle_label, attribute)].append(value)
            else:
                for attribute in attributes:
                    d['{} {}'.format(particle_label, attribute)].append('-')  # Append None instead of '-' (optional)

    # Print the counts
    counters = event_counter(events, particles)

    print("Counts of Events:")
    for particle_label, count in counters.items():
        print("{} events occurred with {}".format(count, particle_label))
    print()

    # Create the DataFrame
    table = pd.DataFrame(data=d)
    print(table)

def event_counter(events, particles):
    # Initialize counters for each particle/unit
    counters = {particle_info['label']: 0 for particle_info in particles}

    for n in range(len(events)):
        empty_event = True

        for particle_info in particles:
            particle_name = particle_info['name']
            particle_label = particle_info['label']

            if events[n][particle_name]:
                counters[particle_label] += 1
                empty_event = False

        # If all particles are missing, mark the event as empty
        if empty_event:
            counters['Empty Event'] += 1

    return counters

def validate_events(events):
    if not isinstance(events, list):
        raise ValueError("Events must be a list")
    for event in events:
        if not isinstance(event, dict):
            raise ValueError("Each event must be a dictionary")
        # Additional validation checks for event data


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
        
        # run processes
        process_events(events)

# run
if __name__ == "__main__":
    print()
    folder = int(input("folder (int): "))
    file = int(input("file (int): "))
    print()
    main(folder, file)