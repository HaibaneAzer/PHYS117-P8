from LHCO_reader import LHCO_reader
import os
import sys
from sklearn.metrics import confusion_matrix
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
from LHCO_comp_functions import signal_efficiency, filter_events_by_HT, filter_events_by_Largest_PT
from LHCO_file_handling import create_file_path_list, process_selected_file

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

   
    selected_files = "3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,21,22,1,2"
    
    selected_files = selected_files.split(',')
    selected_files = [int(idx.strip()) - 1 for idx in selected_files]
    
    return selected_files, file_path_list

def select_value_prompt(type):
    print("Select GeV value of {} selection".format(type))
    selected_value = float(raw_input())
    return selected_value

def filter_events_by_selection(events, HT_cut, Largest_PT_cut_L, Largest_PT_cut_R):
    
    events_left_HT, events_right_HT = filter_events_by_HT(events, HT_cut)
    events_L_L, events_L_R = filter_events_by_Largest_PT(events_left_HT, Largest_PT_cut_L)
    events_R_L, events_R_R = filter_events_by_Largest_PT(events_right_HT, Largest_PT_cut_R)

    sphaleron_selection = (events_L_L, events_R_L)
    black_hole_selection = (events_L_R, events_R_R)

    return events_L_L, events_R_L, events_L_R, events_R_R

def process_files_to_selection_table(file_path_list, selected_files):
    typeList = ["HT", "Largest PT"]
    legend_labels = []
    num_events_list = []
    events_list = []
    sphaleron_event_list_L = []
    sphaleron_event_list_R = []
    blackhole_event_list_L = []
    blackhole_event_list_R = []
    
    for idx in selected_files:
        file_name = file_path_list[idx]
        events_list.append(process_selected_file(file_name))
        legend_labels.append(file_name)

    if events_list is not None:
        print("HT cut")
        HT_cut = select_value_prompt(typeList[0])
        print("L PT left Cut")
        Largest_PT_cut_L = select_value_prompt(typeList[1])
        print("L PT right cut")
        Largest_PT_cut_R = select_value_prompt(typeList[1])
        for event in events_list:
            events_L_sph, events_R_sph, events_L_BH, events_R_BH = filter_events_by_selection(event, HT_cut, Largest_PT_cut_L, Largest_PT_cut_R)
            sphaleron_event_list_L.append(len(events_L_sph))
            sphaleron_event_list_R.append(len(events_R_sph))
            blackhole_event_list_L.append(len(events_L_BH))
            blackhole_event_list_R.append(len(events_R_BH))
    

    return sphaleron_event_list_L, sphaleron_event_list_R, blackhole_event_list_L, blackhole_event_list_R, legend_labels


def main():
    # prompts
    selected_files, file_path = locate_file_prompter()

    sphaleron_selection_L, sphaleron_selection_R, black_hole_selection_L, black_hole_selection_R, legend_labels = process_files_to_selection_table(file_path, selected_files)

    legend_labels = [os.path.basename(label) for label in legend_labels]

    events_table = [
        black_hole_selection_L,
        black_hole_selection_R,
        sphaleron_selection_L,
        sphaleron_selection_R
    ]

    max_label_length = max(len(label) for label in legend_labels)

    transposed_matrix = list(map(list, zip(*events_table)))
    idx = 0
    print(" "*(max_label_length + 1)+"BH Left | BH Right | Sph Left | Sph Right")
    for idx, label in enumerate(legend_labels):
        row_values = " ".join(str(val).ljust(9) for val in transposed_matrix[idx])
        label_spaces = " " * (max_label_length - len(label))  # Calculate spaces for label alignment
        print("{}{} {}".format(label, label_spaces, row_values))

# run
if __name__ == "__main__":
    main()