import argparse
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import os
import sys

def plot_ss_data(filename, dpi=300, tu='us', dt=1e-04, H_col='darkviolet', E_col='yellow', C_col='aqua', out='ssmtx.png', width=10, height=8, x_label_fontsize=14, y_label_fontsize=14, y_num_major_ticks=12,
                 x_num_major_ticks=10, x_tick_size=12, y_tick_size=12, title=None, title_size=16):
    """
    Plot secondary structure data from an input file.

    Args:
        filename (str): The input file name.
        dpi (int): Dots per inch for saving the figure (default: 300).
        tu (str): Time unit, either 'us' (default) or 'ns'.
        dt (float): Time between consecutive frames (default: 1e-04).
        H_col (str): Alpha-helix color (default: 'darkviolet').
        E_col (str): Beta-sheet color (default: 'yellow').
        C_col (str): Coil color (default: 'aqua').
        out (str): Image output name (default: 'ssmtx.png').
        width (float): Image width in inches (default: 10).
        height (float): Image height in inches (default: 8).
        x_label_fontsize (int): Fontsize for x-axis labels (default: 14).
        y_label_fontsize (int): Fontsize for y-axis labels (default: 14).
        y_num_major_ticks (int): Number of major tick locators on the y-axis (default: 10).
        x_num_major_ticks (int): Number of major tick locators on the x-axis (default: 10).
        x_tick_size (int): Size of ticks on the x-axis (default: 12).
        y_tick_size (int): Size of ticks on the y-axis (default: 12).
        title_size (int): Size of the plot title (default: 16).
    """
    
    # Edit the font, font size, and axes width
    mpl.rcParams['font.family'] = 'Nimbus Sans'

    # Configure the plot size
    plt.figure(figsize=(width, height))

    # Color scheme
    colors = [H_col, E_col, C_col]
    cmap = plt.matplotlib.colors.ListedColormap(colors)

    # Define the time unit (tu) and adjust dt accordingly
    if tu == 'ns':
        dt = 1e-01
    else:
        dt = 1e-04

    if filename == "ss_by_frame.xvg":

        data = pd.read_csv(filename, delim_whitespace=True, skiprows=2,
                           names=['frame', 'H(%)', 'E(%)', 'C(%)'])

        # Convert frames to microseconds
        data['frame'] = data['frame'] * dt

        # Plot
        line_h, = plt.plot(data['frame'], data['H(%)'], label='H', linestyle='-', lw=0.5, color=colors[0])
        line_e, = plt.plot(data['frame'], data['E(%)'], label='E', linestyle='-', lw=0.5, color=colors[1])
        line_c, = plt.plot(data['frame'], data['C(%)'], label='C', linestyle='-', lw=0.5, color=colors[2])

        # Automatically adjust axis limits
        plt.xlim(data['frame'].min(), data['frame'].max())
        plt.ylim(data[['H(%)', 'E(%)', 'C(%)']].min().min(), data[['H(%)', 'E(%)', 'C(%)']].max().max())

        # Add legend with squares and matching colors
        legend_elements = [
            Line2D([0], [0], marker='s', color='w', label='H', markersize=12, markerfacecolor=colors[0], markeredgecolor=colors[0]),
            Line2D([0], [0], marker='s', color='w', label='E', markersize=12, markerfacecolor=colors[1], markeredgecolor=colors[1]),
            Line2D([0], [0], marker='s', color='w', label='C', markersize=12, markerfacecolor=colors[2], markeredgecolor=colors[2])
        ]

        plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.99, 0.5), frameon=False, fontsize=14)

    elif filename == "ss_by_res.xvg":

        # If the file contains "res" in its name, create a bar plot

        # Read data from the file and create a DataFrame
        data = pd.read_csv(filename, delim_whitespace=True, skiprows=2,
                           names=['ResidueNumber', 'H(%)', 'E(%)', 'C(%)'])

        bars_h = plt.bar(data['ResidueNumber'], data['H(%)'], label='H', color=colors[0])
        bars_e = plt.bar(data['ResidueNumber'], data['E(%)'], label='E', bottom=data['H(%)'], color=colors[1])
        bars_c = plt.bar(data['ResidueNumber'], data['C(%)'], label='C', bottom=data['H(%)'] + data['E(%)'], color=colors[2])

        # Automatically adjust axis limits
        plt.xlim(data['ResidueNumber'].min(), data['ResidueNumber'].max())
        plt.ylim(data[['H(%)', 'E(%)', 'C(%)']].min().min(), data[['H(%)', 'E(%)', 'C(%)']].max().max())

        # Add legend with squares and matching colors
        legend_elements = [
            Line2D([0], [0], marker='s', color='w', label='H', markersize=12, markerfacecolor=colors[0], markeredgecolor=colors[0]),
            Line2D([0], [0], marker='s', color='w', label='E', markersize=12, markerfacecolor=colors[1], markeredgecolor=colors[1]),
            Line2D([0], [0], marker='s', color='w', label='C', markersize=12, markerfacecolor=colors[2], markeredgecolor=colors[2])
        ]

        plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.99, 0.5), frameon=False, fontsize=14)

    elif filename == "ss.mtx":

        data = pd.read_csv(filename, delim_whitespace=True, skiprows=1, header=None)

        # Replace 'H' with 1, 'E' with 2, and 'C' with 3 in all columns except the first (frame)
        data.iloc[:, 1:] = data.iloc[:, 1:].replace({'H': -100, 'E': 0, 'C': 100})

        # Transpose the DataFrame so that residues are on the x-axis and time on the y-axis
        data_transposed = data.transpose()

       
        # Replace Frames by time in xticks (columns)
        xticks_labels = data_transposed.iloc[0]
        xticks_labels *= dt
        xticks_labels = round(xticks_labels, 1)

        data_transposed.columns = xticks_labels

       
        sns.heatmap(data_transposed.iloc[1:, :], cmap=cmap, cbar=False)

        # Define colors for H, E, and C
        legend_elements = [Line2D([0], [0], marker='s', color='w', label='H', markersize=12, markerfacecolor=colors[0]),
                           Line2D([0], [0], marker='s', color='w', label='E', markersize=12, markerfacecolor=colors[1]),
                           Line2D([0], [0], marker='s', color='w', label='C', markersize=12, markerfacecolor=colors[2])]

        # Add legend with squares and matching colors
        plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.99, 0.5), frameon=False, fontsize=14)

        y_max = data_transposed.shape[0]

        plt.ylim(0, y_max)

        plt.xticks(rotation=0)
        plt.yticks(rotation=0)


    # Configure axis labels
    if filename.lower() in ['ss.mtx', 'ss_by_frame.xvg']:

        if tu == 'us':
            plt.xlabel('Time ($\mu$s)', fontsize=x_label_fontsize)
        else:
            plt.xlabel('Time (ns)', fontsize=x_label_fontsize)

        if filename.lower() != 'ss.mtx':
            plt.ylabel('Percentage', fontsize=y_label_fontsize)
        else:
            plt.ylabel('Residue', fontsize=y_label_fontsize)

    else:
        plt.xlabel('Residue Number', fontsize=y_label_fontsize)
        plt.ylabel('Percentage', fontsize=y_label_fontsize)

    # Configure tick sizes on x and y axes
    plt.xticks(fontsize=x_tick_size)
    plt.yticks(fontsize=y_tick_size)

    # Configure title and title size
    if args.title is not None:
        plt.title(args.title, fontsize=args.title_size)

    # Configure the number of major tick locators on the y-axis
    plt.locator_params(axis='y', nbins=y_num_major_ticks)

    # Configure the number of major tick locators on the x-axis
    plt.locator_params(axis='x', nbins=x_num_major_ticks)

    # Save the plot as a PNG image
    plt.savefig(out, dpi=dpi, bbox_inches='tight')

    # Show the plot
    plt.show()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a PNG image from an input file (ss.mtx, ss_by_frame.xvg, or ss_by_res.xvg) generated by SIRAH.')
    parser.add_argument('-i', dest='filename', metavar='[input]', type=str, default='ss.mtx', help='Input file name (ss.mtx, ss_by_frame.xvg, ss_by_res.xvg)')
    parser.add_argument('-d', dest='dpi', metavar='[dpi]', type=int, default=300, help='DPI (dots per inch) for saving the figure (default: 300)')
    parser.add_argument('-tu', dest='tu', metavar='[tu]', type=str, default='us', help='Time unit (us or ns) (default: us)')
    parser.add_argument('-dt', dest='dt', metavar='[dt]', type=float, default=1e-04, help='Time between consecutive frames (default: 1e-04)')
    parser.add_argument('-H', dest='H_col', metavar='[helix-color]', type=str, default='darkviolet', help='Alpha-helix color (default: darkviolet)')
    parser.add_argument('-E', dest='E_col', metavar='[beta-sheet color]', type=str, default='yellow', help='Beta-sheet color (default: yellow)')
    parser.add_argument('-C', dest='C_col', metavar='[coil color]', type=str, default='aqua', help='Coil color (default: aqua)')
    parser.add_argument('-o', dest='out', metavar='[out name]', type=str, default=None, help='Image output name (default: input name)')
    parser.add_argument('-wt', dest='width', metavar='[width]', type=float, default=10, help='Image width (inches) (default: 10)')
    parser.add_argument('-ht', dest='height', metavar='[height]', type=float, default=8, help='Image height (inches) (default: 8)')
    parser.add_argument('-xfs', dest='x_label_fontsize', metavar='[xlab fontsize]', type=int, default=14, help='Fontsize for x-axis labels (default: 14)')
    parser.add_argument('-yfs', dest='y_label_fontsize', metavar='[ylab fontsize]', type=int, default=14, help='Fontsize for y-axis labels (default: 14)')
    parser.add_argument('-yticks', dest='y_num_major_ticks', metavar='[y # ticks]', type=int, default=10, help='Number of major tick locators on the y-axis (default: 10)')
    parser.add_argument('-xticks', dest='x_num_major_ticks', metavar='[x # ticks]', type=int, default=10, help='Number of major tick locators on the x-axis (default: 10)')
    parser.add_argument('-xtsize', dest='x_tick_size', metavar='[xtsize]', type=int, default=12, help='Size of ticks on the x-axis (default: 12)')
    parser.add_argument('-ytsize', dest='y_tick_size', metavar='[ytsize]', type=int, default=12, help='Size of ticks on the y-axis (default: 12)')
    parser.add_argument('-title', dest='title', metavar='[title]', type=str, default=None, help='Title of the plot (default: None)')
    parser.add_argument('-ttsize', dest='title_size', metavar='[title_size]', type=int, default=16, help='Size of the plot title (default: 16)')
    parser.add_argument('--version', action='store_true', help='Print version and exit')

    args = parser.parse_args()
    
    if args.version:
        print("Version 1.0 [September 2023]")
        sys.exit(0)

    if args.out is None:
        input_filename = os.path.splitext(args.filename)[0]
        args.out = input_filename


    # Check if the input file is valid
    if args.filename not in ['ss.mtx', 'ss_by_frame.xvg', 'ss_by_res.xvg']:
        print(f"Input file '{args.filename}' is not a valid input file.")
        sys.exit(1)

    # Call the function to plot secondary structure data
    plot_ss_data(args.filename, args.dpi, args.tu, args.dt, args.H_col, args.E_col, args.C_col, args.out, args.width, args.height, args.x_label_fontsize, args.y_label_fontsize, args.y_num_major_ticks,
                 args.x_num_major_ticks, args.x_tick_size, args.y_tick_size, args.title_size, args.title)
