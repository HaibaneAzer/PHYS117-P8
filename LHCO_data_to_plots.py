import matplotlib.pyplot as plt
import os


### main plotting functions ###
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
    # shorten labels for all plots
    legend_labels = [os.path.basename(label) for label in legend_labels]

    if plot_type == 'electrons_muons_taus_per_event':
        fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
        particles = ['Electrons', 'Muons', 'Taus']

        for i, ax in enumerate(axes):
            # Create a list to store data for each file
            data_per_file = []

            for j, (x, y) in enumerate(zip(x_data, y_data)):
                # Count the occurrences of each value in y_data[i]
                counts = {}
                for value in y[i]:
                    if value in counts:
                        counts[value] += 1
                    else:
                        counts[value] = 1

                # Separate the counts and values
                values = list(counts.keys())
                frequencies = list(counts.values())

                # Store the frequencies for this file
                data_per_file.append((values, frequencies))

            num_files = float(len(legend_labels))

            # Calculate the width based on the number of files
            base_width = 0.7  # Set a base width
            width_scaling_factor = 1 / num_files  # Adjust this factor to control the scaling
            total_width = num_files * base_width * width_scaling_factor

            for j, (values, frequencies) in enumerate(data_per_file):
                x_values_with_offset = [value + (j - (num_files - 1) / 2) * (total_width / num_files) - total_width / (2 * num_files) for value in values]
                adjusted_width = total_width / num_files
                ax.bar(x_values_with_offset, frequencies, width=adjusted_width, align='edge', alpha=0.7, label=legend_labels[j])

            ax.set_title('{} per Event'.format(particles[i]))
            ax.set_ylabel('Frequency')
            ax.legend()

        plt.xlabel('Total Particles')
        plt.tight_layout()
    elif plot_type == 'largest_PT_in_event':
        # New code for binned events line plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))

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
    elif plot_type == 'delta_R_per_event':
        # New code for binned events line plot
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        # Make plot for each data file
        for label, x_R, y_R in zip(legend_labels, x_data, y_data):
            ax1.plot(x_R, y_R, '-', label=label)
        
        ax1.grid()
        ax1.set_title("Delta_R(Largest {}, MET)".format(particle_name))
        ax1.set_xlabel('Delta R (Angle)')
        ax1.set_ylabel('Relative Frequency')
        ax1.legend()

        plt.tight_layout()
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

    elif plot_type in ['HT', 'meff']:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        unit = ""

        # make plot for each data file
        chosen_file = ""
        colorlist = ["red", "blue"]
        if signal_eff and t_cut_optimal != None: 
            t_cut_idx = 0 
            for x, y in zip(x_data, y_data):
             
                ax1.axvline(x=t_cut_optimal[t_cut_idx], linestyle='--', color=colorlist[t_cut_idx], label='T_cut Optimal ({}), with Significance ({})'.format(t_cut_optimal[t_cut_idx], signal_eff[t_cut_idx]))
                t_cut_idx += 1

        loop = 0
        # make plot of Ht/meff
        for label, x, y in zip(legend_labels, x_data, y_data):
            if len(x_data) == 2 and loop == 0:
                chosen_file = "first choice: "
            elif len(x_data) == 2 and loop == 1:
                chosen_file = "second choice: "
            ax1.plot(x, y, '-', color=colorlist[loop], label="{}{}".format(chosen_file, label))
            loop += 1
            
        ax1.grid()
        ax1.set_title('{}'.format(plot_type))
        ax1.set_ylabel('Relative Frequency')
        ax1.legend()

    else: # all other plot types as line graph
        fig = plt.figure()
        ax = fig.add_subplot(111)
        unit = ""

        # make plot for each data file
        for label, x, y in zip(legend_labels, x_data, y_data):
            ax.plot(x, y, '-', label=label)
        
        ax.grid()
        ax.set_title('{} {}'.format(particle_name, prop_name))
        
        # give appropiate unit scale
        if prop_name == 'PT':
            unit = "(GeV)"
        elif prop_name == 'phi':
            unit = "(angle)"

        # set y and x labels
        
        ax.set_xlabel('{} {}'.format(prop_name, unit))
        ax.set_ylabel('Relative Frequency')
        ax.legend()
    # for all plots
    plt.show()