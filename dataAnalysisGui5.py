#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 14:05:39 2018

@author: avi
"""

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasAgg, FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
from tkinter import filedialog
from tkinter.colorchooser import *
from analysis_tools_update2 import *
from pandastable import Table

import pandas
import scipy
import astropy

root = Tk()
catalogue, sources, numeric_cols = loadThor()

#draw figure into canvas
def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = PhotoImage(master=canvas, width=figure_w, height=figure_h)

    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    #return handle containing reference to photo - must be kept live else picture disappears
    return photo


class CoordPlot():
    def __init__(self, master):
        self.plotting=False
        self.coord_scatter_win = Toplevel()
        self.coord_scatter_win.title("Plot coordinates")
        
        Label(self.coord_scatter_win, text="Select type of coordinates to plot: ").grid(row=0, column=0, sticky=N)
        self.coord_type = IntVar(self.coord_scatter_win)
        Radiobutton(self.coord_scatter_win,text="Galactic", variable=self.coord_type, value=0).grid(row=1, column=0, sticky=N)
        Radiobutton(self.coord_scatter_win, text="J2000", variable=self.coord_type, value=1).grid(row=2, column=0, sticky=N)
                
        Label(self.coord_scatter_win, text="Select the variable to colour the points by: ").grid(row=0, column=1)
        self.color_var = StringVar(self.coord_scatter_win)
        self.color_var.set(list(numeric_cols[0]).insert(0, ""))
        OptionMenu(self.coord_scatter_win, self.color_var, *numeric_cols).grid(row=1, column=1)
        
        
        Label(self.coord_scatter_win, text="Select colour map: ").grid(row=0, column=2)
        self.color_map = StringVar(self.coord_scatter_win)
        self.color_map.set(color_map_options[0])
        OptionMenu(self.coord_scatter_win, self.color_map, *color_map_options).grid(row=1, column=2)
        
        Button(self.coord_scatter_win, text="Plot", command=self.plot).grid(row=2, column=4)
        Button(self.coord_scatter_win, text="Cancel", command=self.coord_scatter_win.destroy).grid(row=2, column=3)
    def plot(self):
        print("plot pressed")
        self.plotting=True
        self.coord_scatter_win.destroy()
        

class IntensityHistPlot():
    def __init__(self, master):
        self.plotting = False
        
        self.intensity_hist_win = Toplevel()
        self.intensity_hist_win.title("Plot intensity histogram")
        
        Label(self.intensity_hist_win, text="Select frequency in MHz: ").grid(row=0, column=0)
        self.freq = IntVar(self.intensity_hist_win)
        self.freq.set(valid_frequencies[0])
        freq_option_menu = OptionMenu(self.intensity_hist_win, self.freq, *valid_frequencies)
        freq_option_menu.grid(row=1, column=0)
        
        
        self.fit_var = IntVar(self.intensity_hist_win)
        fit_line_check = Checkbutton(self.intensity_hist_win, text="Show best fit line", variable=self.fit_var)
        fit_line_check.grid(row=0, column=1)
        
        
        Label(self.intensity_hist_win, text="Number of bins").grid(row=0, column=2)
        self.num_bins = StringVar(self.intensity_hist_win)
        self.num_bins.set("50")
        num_bins_spin = Spinbox(self.intensity_hist_win, from_=1, to=200, textvariable=self.num_bins)
        num_bins_spin.grid(row=1, column=2)
        
        Label(self.intensity_hist_win, text="Label for histogram").grid(row=0, column=3)
        self.label_var = StringVar(self.intensity_hist_win)
        label_entry = Entry(self.intensity_hist_win, textvariable=self.label_var)
        label_entry.grid(row=1, column=3)
                   
        self.color_var = StringVar(self.intensity_hist_win)
        self.color_var.set("#000000")
        def getColor():
            color = askcolor() 
            self.color_var.set(str(color[1]))
        Button(self.intensity_hist_win, text="Select colour", command=getColor).grid(row=0, column=4)
        Label(self.intensity_hist_win, textvariable=self.color_var).grid(row=1, column=4)
        
        Button(self.intensity_hist_win, text="Plot", command=self.plot).grid(row=2, column=6)
        Button(self.intensity_hist_win, text="Cancel", command=self.intensity_hist_win.destroy).grid(row=2, column=5)
        
    def plot(self):
        print("plot pressed")
        self.plotting=True
        self.intensity_hist_win.destroy()
    

class IntensityScatter():
    def __init__(self, master):
        self.plotting = False
        self.intensity_scat_win = Toplevel()
        self.intensity_scat_win.title("Plot intensity scatter")
          
        Label(self.intensity_scat_win, text="Select frequencies in MHz: ").grid(row=0, column=0)
        self.freq_1 = IntVar(self.intensity_scat_win)
        self.freq_1.set(valid_frequencies[0])
        freq_1option_menu = OptionMenu(self.intensity_scat_win, self.freq_1, *valid_frequencies)
        freq_1option_menu.grid(row=1, column=0)
        
        Label(self.intensity_scat_win, text="Select frequencies in MHz: ").grid(row=0, column=1)
        self.freq_2 = IntVar(self.intensity_scat_win)
        self.freq_2.set(valid_frequencies[0])
        self.freq_2option_menu = OptionMenu(self.intensity_scat_win, self.freq_2, *valid_frequencies)
        self.freq_2option_menu.grid(row=1, column=1)
        
        Button(self.intensity_scat_win, text="Plot", command=self.plot).grid(row=3, column=3)
        Button(self.intensity_scat_win, text="Cancel", command=self.intensity_scat_win.destroy).grid(row=3, column=2)
    
    def plot(self):
        print("plot pressed")
        self.plotting=True
        self.intensity_scat_win.destroy()
        

#THOR BUTTONS
btn_group_thor = Frame(root)
btn_group_thor.grid(row=0, column=0, sticky=W)
convenience_plots = Label(btn_group_thor, text="THOR buttons:")
plot_coord_btn = Button(btn_group_thor, text="Plot coordinates")
plot_intensity_hist_btn = Button(btn_group_thor, text="Intensity histogram")
plot_intensity_scat_btn = Button(btn_group_thor, text="Intensity scatter")
convenience_plots.pack(side=LEFT)
plot_coord_btn.pack(side=LEFT)
plot_intensity_hist_btn.pack(side=LEFT)
plot_intensity_scat_btn.pack(side=LEFT)

btn_group_figure = Frame(root)
btn_group_figure.grid(row=1, column=0, sticky=W)
figure_settings_label = Label(btn_group_figure, text="Figure settings:")
clear_figure_btn = Button(btn_group_figure, text="Clear figure")
figure_settings_btn = Button(btn_group_figure, text="Open figure settings")
figure_settings_label.pack(side=LEFT)
clear_figure_btn.pack(side=LEFT)
figure_settings_btn.pack(side=LEFT)

#SECOND ROW OF BUTTONS
btn_group_2 = Frame(root)
btn_group_2.grid(row=2, column=0, sticky=W)
other_plots = Label(btn_group_2, text="Other tools")
remove_filter_btn = Button(btn_group_2, text="Remove filters")
plot_custom = Button(btn_group_2, text="Custom plot")
show_table_btn = Button(btn_group_2, text="Show table")
load_data_btn = Button(btn_group_2, text="Load own table")
load_thor_btn = Button(btn_group_2, text="Reload THOR")
filter_general_btn = Button(btn_group_2, text="Filter sources")
outlier_btn = Button(btn_group_2, text="Remove outliers")
drop_blanks_btn = Button(btn_group_2, text="Remove blanks")
plot_custom.pack(side=LEFT)
show_table_btn.pack(side=LEFT)
load_data_btn.pack(side=LEFT)
load_thor_btn.pack(side=LEFT)
filter_general_btn.pack(side=LEFT)
outlier_btn.pack(side=LEFT)
drop_blanks_btn.pack(side=LEFT)
remove_filter_btn.pack(side=LEFT)

def open_file(event):
    ftypes = [('Commas Separated Variable', '*csv')]
    dlg = filedialog.Open(root, filetypes = ftypes)
    file_name = dlg.show()
    if file_name != "":
        global catalogue
        global sources
        global numeric_cols
        catalogue, sources, numeric_cols = loadCustom(file_name)
        update_cols_list()
        plot_coord_btn.configure(state="disabled")
        plot_intensity_hist_btn.configure(state="disabled")
        plot_intensity_scat_btn.configure(state="disabled")
        plot_coord_btn.unbind("<Button-1>")
        plot_intensity_hist_btn.unbind("<Button-1>")
        plot_intensity_scat_btn.unbind("<Button-1>")
load_data_btn.bind("<Button-1>", open_file)

def open_thor(event):
    global catalogue
    global sources
    global numeric_cols
    catalogue, sources, numeric_cols = loadThor()
    update_cols_list()

    plot_coord_btn.configure(state="normal")
    plot_intensity_hist_btn.configure(state="normal")
    plot_intensity_scat_btn.configure(state="normal")

    plot_coord_btn.bind("<Button-1>", coord_plotter)
    plot_intensity_hist_btn.bind("<Button-1>", intensity_hist_plotter)
    plot_intensity_scat_btn.bind("<Button-1>", intensity_scatter_plotter)
load_thor_btn.bind("<Button-1>", open_thor)

#THE FIGURE AND THE CANVAS IN WHICH IT WILL BE DISPLAYED
fig = Figure(figsize=(8,5), dpi=80)
ax = fig.subplots()

canvas_frame = Frame(root)
canvas = FigureCanvasTkAgg(figure=fig, master=canvas_frame)
canvas.show()
canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)
toolbar = NavigationToolbar2TkAgg(canvas, canvas_frame)
toolbar.update()
canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
canvas_frame.grid(row=1, column=2)

#PLOT COORDINATES FROM THOR CATALOGUE
def coord_plotter(event):
    coord_plot_window = CoordPlot(root)
    coord_plot_window.coord_scatter_win.wait_window()
    if coord_plot_window.plotting:
        coord_type = coord_plot_window.coord_type.get()
        color_var = coord_plot_window.color_var.get()
        color_map = coord_plot_window.color_map.get()
        
        scatterPlotCoords(ax, sources, coord_type==1, color_var,str(color_map))
        canvas.draw()
plot_coord_btn.bind("<Button-1>", coord_plotter)

#PLOT INTENSITY HISTOGRAM FROM THOR CATALOGUE
def intensity_hist_plotter(event):
    intensity_hist_window = IntensityHistPlot(root)
    intensity_hist_window.intensity_hist_win.wait_window()
    
    if intensity_hist_window.plotting:
        freq = intensity_hist_window.freq.get()
        num_bins = intensity_hist_window.num_bins.get()
        label = intensity_hist_window.label_var.get()
        color_var = intensity_hist_window.color_var.get()
        fit = intensity_hist_window.fit_var.get()
        plotIntensityHistogram(ax, sources, freq, num_bins, label, color_var, fit=(fit==1))
        canvas.draw()
plot_intensity_hist_btn.bind("<Button-1>", intensity_hist_plotter)

#INTENSITY SCATTER PLOT FROM THOR CATALOGUE
def intensity_scatter_plotter(event):
    intensity_scat_window = IntensityScatter(root)
    intensity_scat_window.intensity_scat_win.wait_window()
    
    if intensity_scat_window.plotting:
        freq_1 = intensity_scat_window.freq_1.get()
        freq_2 = intensity_scat_window.freq_2.get()
        scatterIntensity(ax, sources, freq_1, freq_2)
        canvas.draw()
plot_intensity_scat_btn.bind("<Button-1>", intensity_scatter_plotter)

class FilterGeneral():
    def __init__(self, master):
        self.filtering = False
    
        self.filter_win = Toplevel()
        self.filter_win.title("Filter sources") 
        
        Label(self.filter_win, text="Filter column: ").grid(row=0, column=0)
        self.filter_var = StringVar(self.filter_win)
        self.filter_var.set("")
        filter_col_menu = OptionMenu(self.filter_win, self.filter_var, *numeric_cols)
        filter_col_menu.grid(row=0, column=1)
        
        self.min_val = DoubleVar(self.filter_win)
        self.max_val = DoubleVar(self.filter_win)
        Label(self.filter_win, text="Min value:").grid(row=1, column=0)
        Label(self.filter_win, text="Max value:").grid(row=2, column=0)
        Entry(self.filter_win, textvariable=self.min_val).grid(row=1, column=1)
        Entry(self.filter_win, textvariable=self.max_val).grid(row=2, column=1)
        
        Label(self.filter_win, text="In your sources...").grid(row=0, column=2)
        self.min_sources = DoubleVar(self.filter_win)
        Label(self.filter_win, textvariable=self.min_sources).grid(row=1, column=3)
        self.max_sources = DoubleVar(self.filter_win)
        Label(self.filter_win, textvariable=self.max_sources).grid(row=2, column=3)

        
        self.filter_var.trace("w", self.updateMinMax)
        
        Button(self.filter_win, text="Confirm filters", command=self.filter_).grid(row=5, column=4)
        Button(self.filter_win, text="Cancel", command=self.filter_win.destroy).grid(row=5, column=3)  
        
    def updateMinMax(self, *args):
        self.min_sources.set(sources[self.filter_var.get()].min())
        self.max_sources.set(sources[self.filter_var.get()].max())

    def filter_(self):
        print("filter pressed")
        self.filtering=True
        self.filter_win.destroy()

#FILTER BY SPECIFYING COLUMN, MIN AND MAX
def filter_general(event):
    filter_window = FilterGeneral(root)
    filter_window.filter_win.wait_window()
    
    if filter_window.filtering:
        global sources
        col = filter_window.filter_var.get()
        min_val = filter_window.min_val.get()
        max_val = filter_window.max_val.get()
        sources = filterByRange(sources, col, min_val, max_val)
        filters_list.append("Restricted " + str(col) + " to between " + str(min_val) + " and " + str(max_val))
        update_filters_list()
filter_general_btn.bind("<Button-1>", filter_general)

class Outliers:
    def __init__(self, master):
        self.removing = False
        self.outlier_win = Toplevel()
        self.outlier_win.title("Filter sources") 


        Label(self.outlier_win, text="Name of column").grid(row=0, column=0)
        Label(self.outlier_win, text="Number of standard deviations").grid(row=1, column=0)

        self.col_var = StringVar(self.outlier_win)
        self.col_var.set(numeric_cols[0])
        OptionMenu(self.outlier_win, self.col_var, *numeric_cols).grid(row=0, column=1)
        self.std_devs = DoubleVar(self.outlier_win)
        self.std_devs.set("2")
        Spinbox(self.outlier_win, from_=0, to=10, textvariable=self.std_devs).grid(row=1, column=1)
        
        self.log = IntVar(self.outlier_win)
        Checkbutton(self.outlier_win, text="Logarithmic", variable=self.log).grid(row=2, column=1)      
        
        Button(self.outlier_win, text="Remove outliers", command=self.remove).grid(row=3, column=4)
        Button(self.outlier_win, text="Cancel", command=self.outlier_win.destroy).grid(row=3, column=3)  

    def remove(self):
        print("removing outliers")
        self.removing = True
        self.outlier_win.destroy()

#REMOVE OUTLIERS BY SPECIFYING COLUMN, NUMBER OF STANDARD DEVIATIONS
def remove_outliers(event):
    outlier_window = Outliers(root)
    outlier_window.outlier_win.wait_window()
    
    if outlier_window.removing:
        global sources
        col = outlier_window.col_var.get()
        std_devs = outlier_window.std_devs.get()
        log = outlier_window.log.get()
        sources = removeOutliers(sources, col, std_devs, log)
        filters_list.append("Removed outliers beyond " + str(std_devs) + " from the mean for " + str(col))
        update_filters_list()
outlier_btn.bind("<Button-1>", remove_outliers)

#RESET BY SETTING SOURCES EQUAL TO ORIGINAL CATALOGUE TO REMOVE ALL FILTERS
def remove_filters(event):
    global sources
    sources = catalogue
    filters_list.clear()
    update_filters_list()
    print(len(sources))
remove_filter_btn.bind("<Button-1>", remove_filters)

#CLEAR THE FIGURE AND ANY PLOTS AND UPDATE THE CANVAS
def clear_fig(event):
    global ax
    fig.gca().cla()
    fig.clf()
    ax = fig.subplots()
    canvas.draw()
clear_figure_btn.bind("<Button-1>", clear_fig)


class FigureSettings():
    def __init__(self, master):
        self.saving = False
        self.fig_settings_win= Toplevel()
        self.fig_settings_win.title("Figure settings")
        
        Label(self.fig_settings_win, text="Title of figure: ").grid(row=0, column=0)
        self.title_var = StringVar(self.fig_settings_win)
        Entry(self.fig_settings_win, textvariable=self.title_var).grid(row=0, column=1)
        
        Button(self.fig_settings_win, text="Save", command=self.save).grid(row=5, column=5)
        Button(self.fig_settings_win, text="Cancel", command=self.fig_settings_win.destroy).grid(row=5, column=4)
    def save(self):
        self.saving=True
        self.fig_settings_win.destroy()

#SETTING FOR CHANGING TITLE OF FIGURE
def add_title(event):
    fig_settings_window = FigureSettings(root)
    fig_settings_window.fig_settings_win.wait_window()
    
    if fig_settings_window.saving:
        fig.suptitle(fig_settings_window.title_var.get())
        canvas.draw()
figure_settings_btn.bind("<Button-1>", add_title)

class CustomPlot():
    def __init__(self, master):
        self.plotting = False
        self.custom_win = Toplevel()
        self.custom_win.title("Create a custom plot")
        
        Label(self.custom_win, text="x column: ").grid(row=0, column=0)
        self.col_1 = StringVar(self.custom_win)
        self.col_1.set("")
        col_1_menu = OptionMenu(self.custom_win, self.col_1, *numeric_cols)
        col_1_menu.grid(row=0, column=1)
        self.fit_var = IntVar(self.custom_win)
        fit_var_check = Checkbutton(self.custom_win, text="Show best fit line", variable=self.fit_var)
        fit_var_check.grid(row=0, column=3)
        Label(self.custom_win, text="Order of fitted polynomial: ").grid(row=1, column=3)
        self.order = IntVar(self.custom_win)
        self.order.set("1")
        order_spin = Spinbox(self.custom_win, from_=1, to=5, textvariable=self.order)
        order_spin.grid(row=1, column=4)
        self.color_var = StringVar(self.custom_win)
        self.color_var.set("#000000")
        def getColor():
            color = askcolor() 
            self.color_var.set(str(color[1]))
        color_btn = Button(self.custom_win, text="Select colour", command=getColor)
        color_btn.grid(row=0, column=4)
        Label(self.custom_win, textvariable=self.color_var).grid(row=0, column=5)
        self.stepped = IntVar(self.custom_win)
        stepped_chk = Checkbutton(self.custom_win, text="Stepped histogram?", variable=self.stepped)
        stepped_chk.grid(row=0, column=6)
        Label(self.custom_win, text="Label: ").grid(row=0, column=7)
        self.label_var = StringVar(self.custom_win)
        label_entry = Entry(self.custom_win, textvariable=self.label_var).grid(row=0, column=8)
        self.log_freq = IntVar(self.custom_win)
        log_freq_check = Checkbutton(self.custom_win, text="Logarithmic frequency?", variable=self.log_freq)
        log_freq_check.grid(row=0, column=9)
        self.cumulative = IntVar(self.custom_win)
        cumulative_chk = Checkbutton(self.custom_win, text="Cumulative?", variable=self.cumulative)
        cumulative_chk.grid(row=0, column=10)
        
        Label(self.custom_win, text="y column: ").grid(row=1, column=0)
        self.col_2 = StringVar(self.custom_win)
        self.col_2.set("")
        col_2_menu = OptionMenu(self.custom_win, self.col_2, *numeric_cols)
        col_2_menu.grid(row=1, column=1)
        col_2_menu.configure(state="disabled")
        self.spearman = IntVar(self.custom_win)
        spearman_chk = Checkbutton(self.custom_win, text="Spearman's rank coefficient", variable=self.spearman, state=DISABLED)
        spearman_chk.grid(row=1, column=5)
        Label(self.custom_win, text="color column: ").grid(row=2, column=0)
        self.col_3 = StringVar(self.custom_win)
        self.col_3.set("")
        col_3_menu = OptionMenu(self.custom_win, self.col_3, *numeric_cols)
        col_3_menu.grid(row=2, column=1)
        col_3_menu.configure(state="disabled")
    
        Label(self.custom_win, text="size column: ").grid(row=3, column=0)
        self.col_4 = StringVar(self.custom_win)
        self.col_4.set("")
        col_4_menu = OptionMenu(self.custom_win, self.col_4, *numeric_cols)
        col_4_menu.grid(row=3, column=1)
        col_4_menu.configure(state="disabled")
            
        self.log_x = IntVar(self.custom_win)
        Checkbutton(self.custom_win, text="Logarithmic", variable=self.log_x,).grid(row=0, column=2)
        self.log_y = IntVar(self.custom_win)
        chk_y = Checkbutton(self.custom_win, text="Logarithmic", variable=self.log_y, state=DISABLED)
        chk_y.grid(row=1, column=2)
        self.log_color = IntVar(self.custom_win)
        chk_color = Checkbutton(self.custom_win, text="Logarithmic", variable=self.log_color, state=DISABLED)
        chk_color.grid(row=2, column=2)
        self.log_size = IntVar(self.custom_win)
        chk_size = Checkbutton(self.custom_win, text="Logarithmic", variable=self.log_size, state=DISABLED)
        chk_size.grid(row=3, column=2)
        
        Label(self.custom_win, text="Select colour map: ").grid(row=2, column=3)
        self.color_map = StringVar(self.custom_win)
        self.color_map.set(color_map_options[0])
        opt_color = OptionMenu(self.custom_win, self.color_map, *color_map_options)
        opt_color.grid(row=2, column=4)
        opt_color.configure(state="disabled")
        def enableOption2(*args):
            col_2_menu.configure(state="active")
            chk_y.config(state=ACTIVE)
            spearman_chk.config(state=ACTIVE)
        def enableOption3(*args):
            col_3_menu.configure(state="active")
            chk_color.config(state=ACTIVE)
            opt_color.config(state=ACTIVE)
            color_btn.config(state=DISABLED)
            stepped_chk.config(state=DISABLED)
            log_freq_check.config(state=DISABLED)
            cumulative_chk.config(state=DISABLED)
        def enableOption4(*args):
            col_4_menu.configure(state="active")
            chk_size.config(state=ACTIVE)
        self.col_1.trace("w", enableOption2)
        self.col_2.trace("w", enableOption3)
        self.col_3.trace("w", enableOption4)
        
        Button(self.custom_win, text="Plot", command=self.plot).grid(row=5, column=5)
        Button(self.custom_win, text="Cancel", command=self.custom_win.destroy).grid(row=5, column=4)
    
    
    def plot(self):
        print("plot pressed")
        self.plotting = True
        self.custom_win.destroy()
#CUSTOM PLOT WITH 1 - 4 VARIABLES
def custom_plot(event):
    custom_window = CustomPlot(root)
    custom_window.custom_win.wait_window()
    
    if custom_window.plotting:
        scatter = customPlot(ax, sources, custom_window.col_1.get(), custom_window.col_2.get(), custom_window.col_3.get(), custom_window.col_4.get(), custom_window.log_x.get(), custom_window.log_y.get(), custom_window.log_color.get(), custom_window.log_size.get(), custom_window.color_map.get(), custom_window.fit_var.get(), custom_window.order.get(), custom_window.color_var.get(), custom_window.color_var.get(), custom_window.label_var.get(), custom_window.log_freq.get(), custom_window.cumulative.get(), custom_window.spearman.get())
        canvas.draw()
plot_custom.bind("<Button-1>", custom_plot)

#REMOVE ROWS WHICH HAVE BLANK VALUE IN SELECTED COLUMN
class DropBlanks():
    def __init__(self, master):
        self.dropping = False
        self.drop_win = Toplevel()
        self.drop_win.title("Remove blanks")
        
        Label(self.drop_win, text="Column: ").grid(row=0, column=0)
        self.drop_var = StringVar(self.drop_win)
        self.drop_var.set("")
        drop_col_menu = OptionMenu(self.drop_win, self.drop_var, *catalogue.columns)
        drop_col_menu.grid(row=0, column=1)
    
        Button(self.drop_win, text="Drop", command=self.drop).grid(row=2, column=3)
        Button(self.drop_win, text="Cancel", command=self.drop_win.destroy).grid(row=2, column=2)
    
    def drop(self):
        self.dropping = True
        self.drop_win.destroy()

def drop_blanks(event):
    drop_window = DropBlanks(root)
    drop_window.drop_win.wait_window()
    
    if drop_window.dropping:
        global sources
        drop_var = drop_window.drop_var.get()
        sources = sources.dropna(subset = [drop_var])
        filters_list.append("Dropped blanks in column: " + str(drop_var))
        update_filters_list()
drop_blanks_btn.bind("<Button-1>", drop_blanks)

#SUMMARY SECTION FOR COLUMNS
col_frame = Frame(root)
col_frame.grid(row=3, column=0, sticky=NW)
scrollbar_cols = Scrollbar(col_frame)
listbox_cols = Listbox(col_frame, yscrollcommand=scrollbar_cols.set)

#UPDATE LIST OF COLUMNS EACH TIME A NEW CATALOGUE IS LOADED IN
def update_cols_list():
    if listbox_cols.size() > 0:
        listbox_cols.delete("0", "end")
    for col in sources.columns:
        listbox_cols.insert(END, str(col))
update_cols_list()
scrollbar_cols.config(command = listbox_cols.yview)
listbox_cols.pack(side=LEFT, fill=BOTH)
scrollbar_cols.pack(side=LEFT, fill=Y)

summary_frame = Frame(root)
summary_frame.grid(row=4, column=0, sticky=NW)
scrollbar_summary = Scrollbar(summary_frame)
listbox_summary = Listbox(summary_frame, yscrollcommand=scrollbar_summary.set)
scrollbar_summary.config(command = listbox_summary.yview)
listbox_summary.pack(side=LEFT, fill=BOTH)
scrollbar_summary.pack(side=LEFT, fill=Y)

#SHOW LIST OF FILTERS THAT HAVE BEEN APPLIED
filters_list = []
filter_frame = Frame(root)
Label(filter_frame, text="List of currently applied filters").pack(side=TOP)
filter_frame.grid(row=3, column=1, sticky=NW, columnspan=3)
scrollbar_filters = Scrollbar(filter_frame)
listbox_filters = Listbox(filter_frame, yscrollcommand=scrollbar_filters.set)
def update_filters_list():
    if listbox_filters.size()>0:
        listbox_filters.delete("0", "end")
    if len(filters_list) != 0:
        for filter_ in filters_list:
            listbox_filters.insert(END, str(filter_))
    else:
        if listbox_filters.size()>0:
            listbox_filters.delete("0", "end")
scrollbar_filters.config(command=listbox_filters.yview)
listbox_filters.pack(side=LEFT, fill="both", expand=True)
scrollbar_filters.pack(side=LEFT)

#UPDATE SUMMARY BOX TO SHOW INFO FOR SELECTED COLUMN
def display_summary(event):
    listbox_summary.delete('0','end')
    col_name = str(listbox_cols.get(listbox_cols.curselection())) 
    listbox_summary.insert(END, "Column " + col_name)
    rows = str(sources[col_name].describe()).split("\n")    
    for row in rows:
        listbox_summary.insert(END, row)
listbox_cols.bind("<<ListboxSelect>>", display_summary)

class TableView:
    def __init__(self, master):
        toplevel = Toplevel()
        toplevel.title("Table")
        frame=Frame(toplevel)
        frame.grid(row=0, column=0)
        table = Table(frame, dataframe=sources)
        table.show()
        
#OPEN WINDOW SHOWING TABLE
def show_table(event):
    table_view = TableView(root)
show_table_btn.bind("<Button-1>", show_table)

root.title("THOR Data Analysis Tools")
root.mainloop()
