import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterExponent
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u
from astropy.coordinates import SkyCoord

fig, ax = plt.subplots()

ERR_INVALID = "Invalid choice. Try again."
ERR_INVALID_NUM ="Not a number. Try again."
ERR_MAX_MIN = "Maximum {0} must be greater than minimum {0}"
valid_frequencies = [1060, 1310, 1440, 1690, 1820, 1950] 
color_map_options = ["Greys_r", "viridis","plasma","inferno", "magma", "gist_heat", "Wistia", "summer", "bwr", "RdGy", "Set1", "tab10", "Pastel1","CMRmap", "brg", "nipy_spectral"]

catalogue = None
sources = None
numeric_cols = None
def loadThor():
    global catalogue
    global sources
    global numeric_cols
    catalogue = pd.read_csv("thor_continuum_catalog.csv")
    sources = catalogue
    numeric_cols = list(catalogue.select_dtypes(include=["float64", "int64"]).columns)
    return catalogue, sources, numeric_cols

def loadCustom(file_name):
    global catalogue
    global sources
    global numeric_cols
    catalogue = pd.read_csv(file_name)
    sources = catalogue
    numeric_cols = list(catalogue.select_dtypes(include=["float64", "int64"]).columns)
    return catalogue, sources, numeric_cols

plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=20)  # fontsize of the figure title


#remove outliers by checking which data points lie beyond n standard deviations from the mean
def removeOutliers(sources, col_name, n_std_devs, log):
    n_std_devs = float(n_std_devs)
    #if user wishes, take log to base 10 before removing outliers and create a temporary column for the logs
    if log:
        sources = sources[sources[col_name] > 0]
        col_log = "col_log"
        #remove all rows for which the data in specified column isn't available
        sources = sources[np.isfinite(sources[col_name])]
        sources[col_log] = np.log10(sources[col_name])
        col_mean = sources[col_log].mean()
        sources =sources[np.abs(sources[col_log]-col_mean) <= (n_std_devs*sources[col_log]).std()]
    else:
    #else remove outliers normally
        sources = sources[np.abs(sources[col_name]-sources[col_name].mean()) <= (n_std_devs*sources[col_name].std())] 
    return sources

#generic function for filtering by specifying min and max values for a column
def filterByRange(sources, col_filter, min_val, max_val):
    min_condition = sources[col_filter] >= min_val
    max_condition = sources[col_filter] <= max_val
    return sources[min_condition & max_condition]

#convert epoch J2000 coordinates into galactic latitude and longitude - used for plotting galactic coordinates
def getGalacticCoords(sources):
    coords = SkyCoord(ra = sources["RA"].values * u.degree, dec = sources["Dec"].values * u.degree)
    return {"l" : coords.galactic.l.deg, "b": coords.galactic.b.deg}

def scatterIntensity(ax, sources, freq_1, freq_2):
    if freq_1 == freq_2:
        return "Both intensities are the same"
    
    #convert frequencies into the column names that correspond to the frequencies and set the axes labels
    intensity_attrib_1 = intensityAttrib(freq_1)
    intensity_attrib_2 = intensityAttrib(freq_2)
    
    ax.set_xlabel(createLabel(intensity_attrib_1))
    ax.set_ylabel(createLabel(intensity_attrib_2))
        
    #get the columns and plot
    intensities_1 = sources[intensity_attrib_1]
    intensities_2 = sources[intensity_attrib_2]
    ax.errorbar(intensities_1, intensities_2, fmt="k,")
    
    print("Plotted the intensity scatter graph")
    ax.set_xscale("log")
    ax.set_yscale("log")

def plotIntensityHistogram(ax, sources, freq_mhz, num_bins, label, color ="#338768", fit=False):
    intensity_attrib = intensityAttrib(freq_mhz)
    
    #if the frequency has not been recorded, then no such column will exist so return
    try:
        peak_intensity = sources[intensity_attrib]
    except:
        print("Intensities for frequency " + str(freq_mhz) + "MHz have not been recorded")
        return 
    #filtering out negative intensities
    peak_intensity = peak_intensity[peak_intensity>0]
    
    #creating the log-scaled bins
    bin_max = np.log10(peak_intensity.max())
    bin_min = np.log10(peak_intensity.min())
    bins = 10**np.linspace(bin_min, bin_max, num_bins)
    
    #set axis label and x scale
    ax.set_xlabel(createLabel(intensity_attrib))
    ax.set_ylabel("Frequency")
    ax.set_xscale("log")
    
    counts, bin_edges, ignored = ax.hist(peak_intensity, bins = bins, rwidth=1.0, color=color, label=str(freq_mhz)+"MHz " + label)
    #draw best fit logarithmic Gaussian curve
    if fit:
        shape, loc, scale = lognorm.fit(peak_intensity, floc=0)
        bins_log_len = np.r_[bin_edges[1:] - bin_edges[:-1], 0]
        # get pdf-values for same intervals as histogram
        samples_fit_log = lognorm.pdf(bins, shape, loc=loc, scale=scale)
        # plot the fit line
        ax.plot(bins, samples_fit_log * bins_log_len * counts.sum(), 'k-', label="Fit line " + label,  linewidth=2)   
        
        #display mean and std dev in a textbox
        mean = round(scale, 2)
        std = round(np.log10(shape), 2)        
        ax.legend((dummyObj(), dummyObj()), ("Mean = " + str(mean), "SD = " + str(std)))

    ax.legend(loc = "lower left", bbox_to_anchor = (0.1, 1.01))
    print("Plotted the histogram")


def scatterPlotCoords(ax, sources, ra, colourVar = "None", colourMap = "Greys_r"):
    #user can plot galactic coordinates or J2000
    if not ra:
        if colourVar != "resolved_source":
            galactic_coords = getGalacticCoords(sources)
            coords = {"x" : galactic_coords["l"], "y" :galactic_coords["b"]}
            
            #doesn't work because galactic longitude and latitude columns aren't completely numeric so their data types are "object" on which numeric operations can't be performed
            #coords = {"x": np.array(sources["GAL_LON"]), "y": np.array(sources["GAL_LAT"]}
        #set axes labels and flip x axis
        ax.set_xlabel(createLabel("GAL_LON"))
        ax.set_ylabel(createLabel("GAL_LAT"))
        ax.invert_xaxis()
    else:       
        coords = {"x": sources["RA"], "y": sources["Dec"]}
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("Dec (deg)")
            
    fig = plt.gcf()
    marker_size = calcMarkerSize(fig)
    
    #if 3rd variable is specified, use it to colour the sources and display a colourbar
    if colourVar == "None":
        scatter = ax.scatter(coords["x"], coords["y"], c="k", s=marker_size)
    elif isLogarithmic(colourVar):
        #if the 3rd variable is better represented on a logarithmic scale then make the colourbar scale logarithmic
        scatter = ax.scatter(coords["x"], coords["y"], c=sources[colourVar],  s=marker_size, cmap=colourMap, norm=LogNorm())
        createColorbar(ax, scatter, log10Label(createLabel(colourVar)), True)
    elif colourVar == "resolved_source":
        #if the 3rd variable is a Boolean (resolved or not), filter the sources into two categories and plot them as separate colours on same set of axes        
        extreme_colour_1 = getColorFromCMAP(colourMap, 0.0)
        extreme_colour_2 = getColorFromCMAP(colourMap, 1.0)

        resolvedSources = filterByRange(sources, colourVar, 1, 1)
        unresolvedSources = filterByRange(sources, colourVar, 0, 0)
        
        #get coordinates of the resolved sources and of the unresolved sources and store in dictionary
        if not ra:
            coords_1, coords_2 = getGalacticCoords(resolvedSources), getGalacticCoords(unresolvedSources)    
            coords_resolved = {"x": coords_1["l"], "y": coords_1["b"]}
            coords_unresolved = {"x": coords_2["l"], "y": coords_2["b"]}
        else:
            coords_resolved = {"x": resolvedSources["RA"], "y": resolvedSources["Dec"]}
            coords_unresolved = {"x": unresolvedSources["RA"], "y": unresolvedSources["Dec"]}
        
        #one plot for resolved, one for unresolved
        scatter = ax.scatter(coords_resolved["x"], coords_resolved["y"], s=marker_size, c=extreme_colour_1, label="Resolved sources")
        scatter = ax.scatter(coords_unresolved["x"], coords_unresolved["y"], s=marker_size, c=extreme_colour_2, label="Unresolved sources")
        
        ax.legend(loc = "lower left", bbox_to_anchor = (0.1, 1.01))
        makeLegendActuallyAppear(ax)
    else:
        scatter = ax.scatter(coords["x"], coords["y"], c=sources[colourVar],  s=marker_size, cmap=colourMap)
        createColorbar(ax, scatter, createLabel(colourVar))
    print("Plotted the coordinates")


def customHistogram(ax, sources, x_col, log_x, fit, color, step_filled, histlabel, log_freq, cumulative):
        #logarithmic y axis?
        if log_freq:
            ax.set_yscale("log")
        ax.set_ylabel("Frequency")
        
        #stepfilled means histogram is coloured, step just shows the outline
        if step_filled:
            histtype = "stepfilled"
        else:
            histtype = "step"
            
        #get column to plot and removing blanks
        values = sources[x_col]
        values = values.dropna()
        
        #create bins depending on whether logarithmic or linear
        if log_x:
            values = values[values>0]
            bins_min = np.log10(values.min())
            bins_max = np.log10(values.max())
            bins = 10**np.linspace(bins_min, bins_max, 50)
            ax.set_xscale("log")
        else:        
            bins = np.linspace(sources[x_col].min(), sources[x_col].max(), 50)
        
        if fit:
            #best fit Normal or LogNormal curve
            if log_x:
                counts, bin_edges, ignored = ax.hist(values, bins=bins, color=color, histtype=histtype, label=histlabel, linewidth=1, cumulative=cumulative)
                shape, loc, scale = lognorm.fit(values, floc=0)
                bins_log_len = np.r_[bin_edges[1:] - bin_edges[:-1], 0]
                # get pdf-values for same intervals as histogram
                samples_fit_log = lognorm.pdf(bins, shape, loc=loc, scale=scale)
                # plot the fit line
                ax.plot(bins, samples_fit_log * bins_log_len * counts.sum(), 'k-',  linewidth=2)                   
                mean = round(scale, 2)
                std = round(np.log10(shape), 2)
            else:
                counts, bin_edges, ignored = ax.hist(values, bins=bins, color=color, normed=True, histtype=histtype, label=histlabel, linewidth=1, cumulative=cumulative)
                mean, std = norm.fit(values)
 
                xmin, xmax = ax.get_xlim()
                x = np.linspace(xmin, xmax, 100)
                y = norm.pdf(x, mean, std)
                print(y)
                ax.plot(x, y, "k-")
     
                mean = round(mean, 2)
                std = round(std, 2)
            
            #display mean and standard deviation in the legend
            ax.legend((dummyObj(), dummyObj()), ("Mean = " + str(mean), "SD = " + str(std)))
        else:
            #draw histogram without best fit line
            ax.hist(values, bins=bins, color=color, histtype=histtype, label=histlabel, linewidth=2, edgecolor=color, cumulative=cumulative)
        
        #display legend if user wants a label - this 
        if histlabel!="":
            ax.legend()
            
            
def customPlot(ax, sources, x_col, y_col, color_col, size_col, log_x = False, log_y = False, log_color = False, log_size=False, color_map="Greys_r", fit = False,order=1, color="r", stepped=False, histlabel="", log_freq=False, cumulative=False, spearman=False):
    fig = plt.gcf()
    marker_size = calcMarkerSize(fig)
            
    ax.set_xlabel(createLabel(x_col))
    
#    try:
    #if only one variable chosen then plot histogram and return
    if y_col == "":
        customHistogram(ax, sources, x_col, log_x, fit, color, stepped, histlabel, log_freq, cumulative)
        return
    #if two variables, plot simple scatter
    elif color_col == "":
        scatter = ax.scatter(sources[x_col], sources[y_col], c="k", s=marker_size, label="")
    #if three variables, colour the points and display colourbar
    elif size_col == "":
        label = createLabel(color_col)
        if log_color:
            label = log10Label(label)
        
        scatter = ax.scatter(sources[x_col], sources[y_col], c=sources[color_col], cmap=color_map, s=marker_size, label="")
        createColorbar(ax, scatter, label, log_color)
    #if four variables, adjust size of points
    else:
        #calculate the quartiles of the size column  
        sizes = sources[size_col]
        color = getColorFromCMAP(color_map, 0.5)
        if log_size:
            #can only log positive
            sizes = sizes.fillna(0.00001)
            sizes = sizes[sizes>0]
            sizes = np.log10(sizes)
            temp_sizes = getTempSizes(sizes)
            #get quartiles of logged sizes to show 3 sizes in the legend
            q1, q2, q3 = getQuartiles(temp_sizes)
            q1_plot, q2_plot, q3_plot = plotsForLegend(q1, q2, q3, color, True)
            q1, q2, q3 = getQuartiles(sizes)
            #create label for legend
            size_col = log10Label(size_col)
            #the actual sizes the points will be
            sizes = 10*(10**temp_sizes)
        else:
            #for sizes that aren't available, set equal to 0 (so they won't get plotted)
            sizes = sizes.fillna(0)
            temp_sizes = getTempSizes(sizes)
            #get quartiles to show 3 sizes in the legend
            q1, q2, q3 = getQuartiles(temp_sizes)
            q1_plot, q2_plot, q3_plot = plotsForLegend(q1, q2, q3, color, False)
            q1, q2, q3 = getQuartiles(sizes)
            #the actual sizes the points will be
            sizes = 100*sizes
        
        #round the quartiles to 3 decimal places for displaying in the legend and create the legend
        dp = 3
        q1 = round(q1, dp)
        q2 = round(q2, dp)
        q3 = round(q3, dp)
        ax.legend((dummyObj(), q1_plot,q2_plot,q3_plot), (size_col, str(q1), str(q2), str(q3)), scatterpoints=1, loc = "lower left", bbox_to_anchor = (0.01, 1.01))
        makeLegendActuallyAppear(ax)
        
        #create colourbar and appropriate label for it
        label = createLabel(color_col)
        if log_color:
            scatter = ax.scatter(sources[x_col], sources[y_col], c=sources[color_col], cmap=color_map, s=sizes, label="")
            label = log10Label(label)
        else:
            scatter = ax.scatter(sources[x_col], sources[y_col], c=sources[color_col], cmap=color_map, s=sizes, label="")
        createColorbar(ax, scatter, label, log_color)
    
    ax.set_ylabel(createLabel(y_col))
    
    #adjust scale on axes if logged so all data points are visible
    if log_x:
        ax.set_xscale("log")
        xs = sources[x_col]
        xs = xs[xs>0]
        ax.set_xlim(xmin=xs.min(), xmax=xs.max())
    if log_y:
        ax.set_yscale("log")
        ys = sources[y_col]
        ys = ys[ys>0]
        ax.set_ylim(ymin=ys.min(), ymax=ys.max())
        
    #if needed, get and draw polynomial line of best fit
    if fit:#
        #drop blank values in the two columns and any values <= 0 if taking log
        sources = sources.dropna(subset=[x_col, y_col])
        if log_x:
            sources = sources[sources[x_col]>0]
        if log_y:
            sources = sources[sources[y_col]>0]
        xs, ys, coeff = getPolynomialFunc(sources[x_col], sources[y_col], log_x, log_y, order)    
        ax.plot(xs, ys(xs), "r--")
        ax.text(0.1, 1.02, "Fitted line coefficients\n" + str(np.round(coeff, 3)), transform=ax.transAxes)
    
    #if needed, get Spearman's Rank Correlation Coefficient and display in a textbox
    if spearman:
        ax.text(0.5, 1.02, r"r$_s$ = " +str(np.round(spearmanr(sources[x_col], sources[y_col])[0], 4)), transform=ax.transAxes)
    return scatter

#    except:
#        print("exception")
#        pass
    

#create label for putting on axis by combining name of quantity with the unit - mainly for THOR catalogue
def createLabel(quantity):
    quantity_return = quantity
    
    #make the labels S_p... and S_int friendlier to read
    if "S_p" in quantity:
        if "delta" in quantity:
            quantity_return = "Delta intensity at " + quantity[14:-1] + "MHz"
        elif len(quantity) > 0:
            quantity_return = "Intensity at " + quantity[8:-1] + "MHz"
        else:
            quantity_return = "Intensity"
    elif quantity == "S_int":
        quantity_return = "Flux density"
    elif quantity == "GAL_LON":
        quantity_return = "Galactic longitude"
    elif quantity == "GAL_LAT":
        quantity_return ="Galactic latitude"
        
    #attach unit if applicable
    unit = getUnit(quantity)
    if unit is None:
        label = quantity_return
    else:
        label = str(quantity_return) + " (" + str(getUnit(quantity)) + ")"
    return label

#return the unit of quantities using Astropy units package
def getUnit(quantity):
    if "S_p" in quantity:
        return u.Jy / u.beam
    elif quantity == "RA" or quantity == "Dec" or quantity == "BPA" or quantity=="GAL_LON" or quantity=="GAL_LAT":
        return u.degree
    elif quantity == "S_int":
        return u.Jy
    elif quantity == "BMAJ" or quantity == "BMIN":
        return u.arcsec
    else:
        return None

def createColorbar(ax, scatter, label, log=False):
    #append the colour bar axes onto the main set of axes
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="10%", pad=0.05)
    cbar = plt.colorbar(scatter, cax=cax, shrink = 0.8, pad = 0.05)    
    cbar.set_label(label)
    
    #make scale of colourbar logarithmic if needed
    if log:
        cbar.formatter = LogFormatterExponent()
        cbar.update_ticks()
    return cbar

def makeLegendActuallyAppear(ax):
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height*0.8])    

def getPolynomialFunc(x, y,log_x, log_y, order):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    #create and return a polynomial function into which x values will be passed to
    if log_x and log_y:
        logx = np.log10(x)
        logy = np.log10(y)
        coefficients = np.polyfit(logx, logy, order)
        polynomial = np.poly1d(coefficients)
        ys = lambda x: np.exp(polynomial(np.log(x)))
    #TODO: fix middle two elifs
    elif log_x:
        logx = np.log10(x)
        coefficients= np.polyfit(logx, y, order)
        polynomial = np.poly1d(coefficients)
        ys = lambda x: polynomial(x)
    elif log_y:
        logy = np.log10(y)
        coefficients = np.polyfit(x, logy, order)
        polynomial = np.poly1d(coefficients)
        ys = lambda x: np.exp(polynomial(np.log(x)))
    else:
        coefficients = np.polyfit(x, y, order)
        polynomial = np.poly1d(coefficients)
        ys = lambda x: polynomial(x)
    return x, ys, coefficients

#get a specific colour from a colourmap
def getColorFromCMAP(color_map, fraction):
    cmap = plt.cm.get_cmap(str(color_map))
    return cmap(fraction)

def log10Label(label) -> str:
    return "log" + r"$_1$$_0$(" + label + ")"

def getQuartiles(arr):
    return np.percentile(arr, 25), np.percentile(arr, 50), np.percentile(arr, 75)
#used for creating legend for sizes
def getTempSizes(sizes):
    col_range = sizes.max() - sizes.min()
    return (sizes-sizes.min())/col_range

#used for adding extra information to a legend
def dummyObj():
    return Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

#used for creating the sizes legend
def plotsForLegend(q1, q2, q3, color, log=False):
    marker = "o"
    print(color)
    if log:
        return plt.scatter([], [], s=10*(10**q1), marker=marker, color=color), plt.scatter([], [], s=10*(10**q2), marker=marker, color=color), plt.scatter([], [], s=10*(10**q3), marker=marker, color=color)
    else:
        return plt.scatter([], [], s=100*q1, marker=marker, color=color), plt.scatter([], [], s=100*q2, marker=marker, color=color), plt.scatter([], [], s=100*q3, marker=marker, color=color)

def calcMarkerSize(fig):
    return (200./fig.dpi)**2

#convert a frequency into the corresponding column attribute in the THOR catalogue
#e.g. 1060 returns S_p(spw-1060)
def intensityAttrib(freq_mhz):
    return "S_p(spw-" + str(freq_mhz) + ")"

def deltaIntensityAttrib(freq_mhz):
    return "delta_" + intensityAttrib(freq_mhz)

#used for the THOR catalogue to determine whether scale should be logarithmic
def isLogarithmic(col_name):
    if col_name == "S_int" or "S_p" in col_name or col_name.upper() == "INTENSITY" or col_name=="alpha" or col_name=="SNR" or col_name=="n_pix":
        return True
    
#not quite working but not being used at the moment
#def filterByGalacticCoords(sources, min_lon, max_lon, min_lat, max_lat):
#    min_coord = SkyCoord(frame="galactic", l=min_lon*u.degree, b=min_lat*u.degree)
#    max_coord = SkyCoord(frame="galactic", l=max_lon*u.degree, b=max_lat*u.degree)
#    min_ra = float(min_coord.icrs.ra.deg)
#    max_ra = float(max_coord.icrs.ra.deg)
#    min_dec = float(min_coord.icrs.dec.deg)
#    max_dec = float(max_coord.icrs.dec.deg)
#    return filterByRADec(sources, min_ra, max_ra, min_dec, max_dec)

#convenience function for returning sources within specified region of the sky
#def filterByRADec(sources, min_ra, max_ra, min_dec, max_dec):
#    restricted = filterByRange(sources, "RA", min_ra, max_ra)
#    return filterByRange(restricted, "Dec", min_dec, max_dec)
