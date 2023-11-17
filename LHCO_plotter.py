# Necessary modules for plotting
from LHCO_file_handling import (
    locate_file_prompter, 
    processed_file_to_data, 
    plot_selection_prompt
)
from LHCO_data_to_plots import (
    plot_line_graph
)

### NB: this code is only compatible with python 2.x (solely used to read LHCO files, which its
### known library "LHCO_reader" has stopped its support since 2015).

def main():
    # prompts
    selected_files, file_list = locate_file_prompter()
    plot_type = plot_selection_prompt()

    # file to data
    x_data, y_data, legend_labels, particle_name, prop_name, signal_eff, t_cut_optimal, signal_eff_list = processed_file_to_data(file_list, selected_files, plot_type)
    
    # plotting
    plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name, plot_type, signal_eff, t_cut_optimal, signal_eff_list)

# run
if __name__ == "__main__":
    main()
