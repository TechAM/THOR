# THOR
Data Analysis Tools for the THOR Catalogue

PRE-REQUISITES:
- Python 3.6

Modules:
- Tkinter
- matplotlib
- pandas
- pandastable
- numpy
- scipy
- astropy
- mpl_toolkits

1) Download both Python files, the THOR catalogue CSV file and store in the same directory.
2) Open a terminal and cd into the directory.
3) Run data_analysis_gui5 using the command "python data_analysis_gui5.py"

Here is a brief description for each button
- Plot Coordinates: select either galactic or J2000 coordinates, and can select a third variable with a colourmap to colour the points
- Intensity Histogram: select the frequency of the spectral window from the options shown, and change some parameters of the histogram
- Intensity Histogram: choose two spectral windows and the intensities will be plotted from them

When your own data is loaded, the above buttons get disabled as they are only applicable for THOR.

- Clear Figure: clears all plots and the figure including the title
- Open figure settings: currently only allows you to set the title of the figure

- Custom Plot: can select up to 4 variables and adjust graph parameters as needed (note that buttons get disabled/enabled as you click other buttons so you can't e.g. only select a colour variable as that doesn't make sense)
- Show table: shows the table
- Load own table: opens file dialogue allowing you to load in a different catalogue - must be CSV format so the dialogue only shows CSV files
- Reload THOR: loads the THOR catalogue back in and enables the THOR buttons
- Filter sources: select a column, specifying the minimum and maximum value
- Remove outliers: select a column and enter the number of standard deviations from the mean beyond which outliers should be removed
- Remove blanks: select a column and any rows which have a blank in that column will be dropped
- Remove filters: removes all filters

Any filters you apply will be shown in the little box under "List of currently applied filters". This box will be emptied when filters are removed explicitly or a new file is loaded.

You can pan around and zoom in/out of the plot using the buttons under the canvas in the top

Currently, there isn't much input handling but entering non-numeric letters won't crash the program.
