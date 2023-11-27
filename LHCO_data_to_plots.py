import matplotlib.pyplot as plt
import os


### main plotting functions ###
def plot_line_graph(x_data, y_data, legend_labels, particle_name, prop_name, plot_type, signal_eff, t_cut_optimal, signal_eff_list):
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

            markers=[">", ".", "*", "o"]
            # Calculate the width based on the number of files
            base_width = 0.7  # Set a base width
            width_scaling_factor = 1 / num_files  # Adjust this factor to control the scaling
            total_width = num_files * base_width * width_scaling_factor

            for j, (values, frequencies) in enumerate(data_per_file):
                x_values_with_offset = [value + (j - (num_files - 1) / 2) * (total_width / num_files) - total_width / (2 * num_files) for value in values]
                adjusted_width = total_width / num_files
                ax.scatter(values, frequencies, alpha=0.7, label=legend_labels[j], marker=markers[j])

            ax.set_title('{} per Event'.format(particles[i]))
            ax.set_ylabel('Frequency')

            # adjust legends
            box = ax.get_position()
            """ ax.set_position([box.x0, box.y0, box.width*0.8, box.height]) """
            """ ax.legend(bbox_to_anchor=(1, 0.5)) """
            ax.legend()

        plt.xlabel('Total Particles')
        plt.tight_layout()
    elif plot_type == 'largest_PT_in_event':
        # New code for binned events line plot
        #####
        fig = plt.figure()
        sig_vs_tcut_plot = ""
        while sig_vs_tcut_plot not in ["y", "n"]:
            print("Add signal vs tcut plot (y/n)?")
            sig_vs_tcut_plot = raw_input()
            if sig_vs_tcut_plot not in ["y", "n"]:
                print("Not valid choice. Please choose 'y' for yes or 'n' for no")
        if sig_vs_tcut_plot == "y":
            ax1 = fig.add_subplot(311)
            ax1_a = fig.add_subplot(312)
            ax1_b = fig.add_subplot(313)
            """ ax2 = fig.add_subplot(223) """
            """ ax2_a = fig.add_subplot(322)
            ax2_b = fig.add_subplot(323) """
        else:
            ax1 = fig.add_subplot(111)
            """ ax2 = fig.add_subplot(212) """
        unit = ""

        # make plot for each data file
        chosen_file = ""
        colorlist = ["red", "blue", "purple"]
        if signal_eff and t_cut_optimal != None: 
            t_cut_idx = 0 
            for _ in colorlist:
                label_ax1 = 'T_cut ({}), with Significance ({})'.format(t_cut_optimal[t_cut_idx], signal_eff[t_cut_idx]) # removed for now...
                ax1.axvline(x=t_cut_optimal[t_cut_idx], linestyle='--', color=colorlist[t_cut_idx])
                t_cut_idx += 1

            loop = 0
            # make plot of Ht/meff
            for label, (x_PT, x_phi), (y_PT, y_phi) in zip(legend_labels, x_data, y_data):
                if loop == 0:
                    chosen_file = "blackhole: "
                elif loop == 1:
                    chosen_file = "sphaleron: "
                ax1.plot(x_PT, y_PT, '-', color=colorlist[loop], label="{}{}".format(chosen_file, label))
                """ ax2.plot(x_phi, y_phi, '-', color=colorlist[loop], label="{}{}".format(chosen_file, label)) """
                if sig_vs_tcut_plot == "y":
                    ax1_a.plot(x_PT, signal_eff_list[loop][0], color=colorlist[loop])
                    ax1_b.plot(x_PT, signal_eff_list[loop][1], color=colorlist[loop])
                loop += 1
            if sig_vs_tcut_plot == "y":
                ax1_a.grid()
                ax1_a.set_ylabel("Significance")
                
                ax1_b.grid()
                ax1_b.set_ylabel("Efficiency")
                
            
        else: 
        #####
        # Make plot for each data file
            for label, (x_PT, x_phi), (y_PT, y_phi) in zip(legend_labels, x_data, y_data):
                ax1.plot(x_PT, y_PT, '-', label=label)
                """ ax2.plot(x_phi, y_phi, '-', label=label) """
        
        ax1.grid()
        ax1_b.set_xlabel('Largest PT [GeV]')
        ax1.set_ylabel("Frequency")
        
        """ box = ax1.get_position( )
        ax1.set_position([box.x0, box.y0, box.width*0.6, box.height])
        ax1.legend(bbox_to_anchor=(1, 0.5)) """

        """ ax2.grid()
        ax2.set_xlabel('Delta Phi []')
        ax2.set_ylabel("Relative Frequency") """

        """ box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax2.legend(bbox_to_anchor=(1, 0.5)) """

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
        
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax1.legend(bbox_to_anchor=(1, 0.5))

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
        sig_vs_tcut_plot = ""
        while sig_vs_tcut_plot not in ["y", "n"]:
            print("Add signal vs tcut plot (y/n)?")
            sig_vs_tcut_plot = raw_input()
            if sig_vs_tcut_plot not in ["y", "n"]:
                print("Not valid choice. Please choose 'y' for yes or 'n' for no")
        if sig_vs_tcut_plot == "y":
            ax1 = fig.add_subplot(311)
            ax2 = fig.add_subplot(312)
            ax3 = fig.add_subplot(313)
        else:
            ax1 = fig.add_subplot(111)
        unit = ""

        # make plot for each data file
        chosen_file = ""
        colorlist = ["red", "blue", "purple"]
        if signal_eff and t_cut_optimal != None: 
            t_cut_idx = 0 
            for _ in colorlist:
                label_ax1 = 'T_cut ({}), with Significance ({})'.format(t_cut_optimal[t_cut_idx], signal_eff[t_cut_idx]) # removed for now...
                ax1.axvline(x=t_cut_optimal[t_cut_idx], linestyle='--', color=colorlist[t_cut_idx])
                t_cut_idx += 1

            loop = 0
            # make plot of Ht/meff
            for label, x, y in zip(legend_labels, x_data, y_data):
                if loop == 0:
                    chosen_file = "blackhole: "
                elif loop == 1:
                    chosen_file = "sphaleron: "
                ax1.plot(x, y, '-', color=colorlist[loop], label="{}{}".format(chosen_file, label))
                if sig_vs_tcut_plot == "y":
                    ax2.plot(x, signal_eff_list[loop][0], color=colorlist[loop])
                    ax3.plot(x, signal_eff_list[loop][1], color=colorlist[loop])
                loop += 1
            if sig_vs_tcut_plot == "y":
                ax2.grid()
                ax2.set_ylabel("Significance")
                
                ax3.grid()
                ax3.set_ylabel("Efficiency")
                ax3.set_xlabel("{} [GeV]".format(plot_type))
            
        else: 
            # make plot of Ht/meff
            for label, x, y in zip(legend_labels, x_data, y_data):
                ax1.plot(x, y, '-', label="{}".format(label))
            ax1.set_xlabel("{} [GeV]".format(plot_type))
            
        ax1.grid()
        ax1.set_title('{}'.format(plot_type))
        ax1.set_ylabel('Frequency')
        
        box = ax1.get_position( )
        ax1.set_position([box.x0, box.y0, box.width*0.6, box.height])
        ax1.legend(bbox_to_anchor=(1, 0.5))

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

        box = ax.get_position( )
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(bbox_to_anchor=(1, 0.5))
    # for all plots
    plt.show()