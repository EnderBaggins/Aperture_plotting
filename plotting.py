# %% [markdown]
# # plotting_classes.py

# %%
import matplotlib
import matplotlib.pyplot as plt #cause I use both plt and random mpl functions
import numpy as np



# matplotlib.rc("text", usetex=True)
from IPython.display import display # for displaying figs without clearing them

import inspect # for comparing arguments and parameters

# data handling and file processing
import h5py
import os
import sys
sys.path.append('.')

from datalib_logsph import DataSph, flag_to_species
from datalib import Data

from typing import Any #used in __setattr__ to allow any attribute type to be set

#from scipy.ndimage import gaussian_filter1d
#from matplotlib import rcParams
#import matplotlib.colors as colors

# making things look nice

from mpl_toolkits.axes_grid1 import make_axes_locatable # for colorbar
from matplotlib.colors import LinearSegmentedColormap
matplotlib.rc("font", family="serif")

cdata = {
    "red": [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (0.55, 1.0, 1.0), (1.0, 1.0, 1.0),],
    "green": [(0.0, 1.0, 1.0), (0.45, 0.0, 0.0), (0.55, 0.0, 0.0), (1.0, 1.0, 1.0),],
    "blue": [(0.0, 1.0, 1.0), (0.45, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0),],
}

hot_cold_cmap = LinearSegmentedColormap("hot_and_cold", cdata, N=1024, gamma=1.0)
if "hot_and_cold" not in matplotlib.colormaps:
    matplotlib.colormaps.register(hot_cold_cmap)


# %%
################################################################################

#A debug class object
class debugging:
    def __init__(self, enabled=False):
        self.enabled = enabled
        self.level = 0 #allows for different debug messages
        # 0 is for all messages
        # 1 is for important messages
        # 2 is for only the most important messages

debug = debugging(enabled=False)

# %%
aperture_figure_objects = {}
# this is all the apt_fig objects that have been created
# I will need to fix this

# %%
################################################################################
'''
Defines the function that will be used to calculate the quantity
This class will also store the arguments necessary for plotting the function
'''
class apt_plot:
    def __init__(self,  func, name, plot_function,**kwargs):
        self.func = func
        self.plot_object = None
        self.plot_function = plot_function # e.g colorplot using pcolormesh
        self.parameters = {} # all possible parameters for the plot
        self.position = None # a (row,col) tuple (starts at 0->N-1)
        self.name = name
        self.ax = None

        self.set_default_parameters()
        self.parameters.update(kwargs) #override the defaults
    
    
    def set_default_parameters(self):
        self.parameters['cmap'] = 'hot_and_cold'
        self.parameters['aspect'] = 'equal'
        self.parameters['fontsize'] = 24

    # overrides the parameters of the class object without changing them
    def override_params(self, **kwargs):
        parameters = self.parameters.copy()
        parameters.update(kwargs)

        if debug.enabled and debug.level <= 0:
            overridden = {k: (self.parameters[k], v) for k, v in kwargs.items() if k in self.parameters and self.parameters[k] != v}
    
            for key, (original_value, new_value) in overridden.items():
                print(f"    {self.name}'s Parameter '{key}' was overridden: Original = {original_value}, New = {new_value}")
    
        return parameters

    def make_plot(self,data, **kwargs):
        parameters = self.override_params(**kwargs)
        # creates the plot_object desired
        self.plot_object = self.plot_function(self,data, **parameters)
        # sets the parameters of the axis
        self.set_plot_attr(**parameters)
    
    def set_plot_attr(self, **kwargs):
        parameters = self.parameters.copy()
        parameters.update(kwargs)
        # This requires that fld_val has the following attributes:
        # xlim, ylim, aspect, title
        attrs = ['xlim', 'ylim', 'aspect', 'title']
        for param in parameters:
            if param in attrs:
                try:
                    getattr(self.ax, f"set_{param}")(parameters[param])
                except Exception as e:
                    print(f"Could not set {param}: {e}")

    
    '''
    Takes a new axis and copies over all important aspects of self.ax, 
    then sets self.ax to the new axis
    '''
    def copy_ax_attr(self, new_ax):
        if self.ax is None:
            return #nothing to copy
        # this is mainly so I can create a new axis with the same properties
        # for use of changing the shape of the figure
        ax = self.ax
        properties = {}
        methods = inspect.getmembers(ax, predicate=inspect.ismethod)
        
        #specific attributes I want to pass between axes
        attrs = [
            "title",
            "xlabel",
            "ylabel",
            "xlim",
            "ylim",
            "xscale",
            "yscale",
           # "xticks", #these would copy old ones, if you change the shape it looks weird
            #"yticks",
            #"xticklabels",
            #"yticklabels"

        ]
        #filtering for get_ with set_ by the attrs of axis
        for name, method in methods:
            if (name[4:] in attrs and name.startswith('get_')):
                try:
                    properties[name] = method()
                except Exception as e:
                    print(f"Could not execute {name}: {e}")
        #setting the properties
        for name, value in properties.items():
            try:
                set_method =getattr(new_ax, 'set_'+name[4:])
                set_method(value)
            except Exception as e:
                print(f"Could not set {name}: {e}")
        self.ax = new_ax

    '''
    Sets the fontsizes of the plot
    has each type be set to default of fontsize if not specified
    options are: label_fontsize, title_fontsize, tick_fontsize, legend_fontsize
    '''
    def set_fontsize(self, **kwargs):
        parameters = self.override_params(**kwargs)
        print(f"{self.name} parameters: {parameters}")
        fontsize = parameters.get('fontsize', None)
        
        #each possible fontsize type with default as the fontsize
        label_fontsize = parameters.get('label_fontsize', fontsize)
        title_fontsize = parameters.get('title_fontsize', fontsize)
        tick_fontsize = parameters.get('tick_fontsize', fontsize)
        legend_fontsize = parameters.get('legend_fontsize', fontsize)

        if debug.enabled and debug.level <= 0:
            print(f"    Setting fontsize for {self.name}: {fontsize}")
            print(f"    Label Fontsize: {label_fontsize}")
            print(f"    Title Fontsize: {title_fontsize}")
            print(f"    Tick Fontsize: {tick_fontsize}")
            print(f"    Legend Fontsize: {legend_fontsize}")
    

        #title font size
        self.ax.title.set_fontsize(title_fontsize) if title_fontsize is not None else None
        
        #x and y axis label font size
        if label_fontsize is not None:
            self.ax.xaxis.label.set_fontsize(label_fontsize) 
            self.ax.yaxis.label.set_fontsize(label_fontsize)

        #legend font size
        legend = self.ax.get_legend()
        if legend:
            for text in legend.get_texts():
                text.set_fontsize(legend_fontsize) if legend_fontsize is not None else None
        
        # All tick label fontsizes
        for tick_label in (self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            tick_label.set_fontsize(tick_fontsize) if tick_fontsize is not None else None

            # set font size for colorbar
            cbar = self.cbar
            if cbar:
                cbar.ax.tick_params(labelsize=tick_fontsize) if tick_fontsize is not None else None
        
    def __str__(self):
        # A overinformative string representation
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"fld_quantity with function {self.func.__name__} and kwargs: {kwargs_str}"
    

# %%
################################################################################
'''
apt_fig is a one dataset class object storing all the information necessary to plot a figure of subplots
this will have children objects that will store the information necessary to plot the data of one subplot
'''

class apt_fig:
    def __init__(self,data, unique_ident= "default_identifier", **kwargs):

        #ensures uniqueness of the object
        # i.e you can recreate the same object and it overrides the old one
        global aperture_figure_objects
        self.unique_ident = unique_ident
        if aperture_figure_objects.get(unique_ident) is not None:
            del aperture_figure_objects[unique_ident]
            if debug.enabled and debug.level <= 2:
                print(f"Overriding apt_fig \'{unique_ident}\' with new object")
        aperture_figure_objects[unique_ident] = self
            
        self.fig = None
        self.plots = {}        # has axis information inside the class objects
        self.parameters = {}
        self.data = data       #Only one dataset in apt_fig
        self.step = 0
        self.post_process = {} # post_processing functions. Things that act on all axes or figures

        # storing the information about the shape of the figure
        # and it's subplots. even accounts for empty subplots
        self.columns = 1      # number of columns
        self.rows = 1         # number of rows

        self.parameters.update(kwargs) #override the defaults
    
        #self._init_complete = True #This is to allow set_attr to run normally for initialization
    
    # overrides the parameters of the class object without changing them
    def override_params(self, **kwargs):
        parameters = self.parameters.copy()
        parameters.update(kwargs)

        if debug.enabled and debug.level <= 0:
            overridden = {k: (self.parameters[k], v) for k, v in kwargs.items() if k in self.parameters and self.parameters[k] != v}
    
            for key, (original_value, new_value) in overridden.items():
                print(f"    Parameter '{key}' was overridden: Original={original_value}, New={new_value}")
    
        return parameters
    
    #####################################################
    # changing the default setattr and delattr to do other things
    def __setattr__(self, name: str, value: Any) -> None:
        # loading the data at the timestep
        if name == 'step':
            self.data.load(value)
        '''
        # only runs once the object is fully initialized
        if '_init_complete' in self.__dict__:
            pass
        '''
        #This is necessary to avoid infinite recursion
        super().__setattr__(name, value)
    
    def __delattr__(self, name: str) -> None:
        super().__delattr__(name)
    ####################################################
    
    # checks if this position is taken by another subplot
    def check_position_taken(self,pos):
        for plot in self.plots.values():
            if plot.position == pos:
                return True
        return False
    
    # resizes the row and column to be the larger
    # of the current size and the new position
    def resize_row_col(self,pos):
        #updates the rows and columns to be large enough
        old_rows = self.rows
        old_columns = self.columns
        self.rows = max(self.rows,pos[0]+1)
        self.columns = max(self.columns,pos[1]+1)
        return old_rows, old_columns # in case you need them

    # reshapes the figure to have appropriate shape of subplots
    def set_fig_shape(self,num_rows,num_columns):
        new_fig = plt.figure()
        #copies over the old axes
        for plot in self.plots.values():
            pos = plot.position

            new_ax = plt.subplot2grid((num_rows,num_columns),pos,fig= new_fig)
            # setting plot.ax to new axis with old properties
            plot.copy_ax_attr(new_ax) 

        plt.close(self.fig) #closes the old figure
        self.fig = new_fig
        self.fig.set_label(self.unique_ident)

        if debug.enabled and debug.level <= 1:
            print(f"  Reloaded figure to {num_rows}x{num_columns}, with {len(list(self.plots))+1} subplots")
    
    def add_plot(self,apt_plot_object,pos=None, **kwargs):
        #convert string to function
        ap = apt_plot_object
        #enforces that apt_plot_object is an apt_plot object
        assert isinstance(ap, apt_plot), "apt_plot_object must be a apt_plot object, try e.g EdotB() not EdotB"
        name = ap.name
        #if the plot already exists, raise an error
        if name in self.plots:
            raise ValueError(f"{name} already exists as plot")

        # if no position is given, add to a new column on top
        # consider making it fill empty spots first
        if pos is None:
            #check if there are no plots
            if len(self.plots) == 0:
                pos = (0,0)
            else:
                self.columns += 1
                if self.check_position_taken((0,self.columns-1)):
                    pos = (0,self.columns)
                else:
                    pos = (0,self.columns-1)
        
        #resize the shape if necessary
        self.resize_row_col(pos)
        
        # check if the position is already taken
        if self.check_position_taken(pos):
            # consider making this move the old plot to a new position
            # or this plot to a nearby one
            raise ValueError(f"Position {pos} is already taken")

        # sets the self.fig to have correct shape
        # I think redundent if the number of rows and columns are unchanged
        # but I think its cleaner to make sure everything updates
        self.set_fig_shape(self.rows,self.columns)
            
        # add the plot to the dictionary
        self.plots[name] = ap
        ap.position = pos
        #connecting the ax to the position
        ap.ax = plt.subplot2grid((self.rows,self.columns),pos,fig = self.fig)
        

        if debug.enabled and debug.level <= 2:
            print(f"Added plot {name} to position {pos}")

    def del_plot(self, name):
        # first asserts that the plot exists
        if name not in self.plots:
            raise ValueError(f"{name} does not exist as plot to delete")
        
        del_plot = self.plots[name]

        #now we need to see if we can shrink the figure
        # I want the largest row and column from plot positions
        for plot in self.plots.values():
            if plot == del_plot:
                continue
            pos = plot.position
            self.resize_row_col(pos)

        #deletes the plot
        del self.plots[name]
        if debug.enabled and debug.level <= 2:
            print(f"Deleted plot {name}")

        # reshapes the figure
        # again I think redundant if the number of rows and columns are unchanged
        self.set_fig_shape(self.rows,self.columns)
    
    def move_plot(self, name, pos):
        #first asserts that the plot exists
        if name not in self.plots:
            raise ValueError(f"{name} does not exist as plot")
        
        #check if the position is already taken
        if self.check_position_taken(pos):
            raise ValueError(f"Position {pos} is already taken, implement swapping instead?")
        
        #resize the shape if necessary
        self.resize_row_col(pos)

        #moves the plot to the new position
        self.plots[name].position = pos
        self.set_fig_shape(self.rows,self.columns)
        if debug.enabled and debug.level <= 2:
            print(f"Moved plot {name} to position {pos}")


       
    # this just updates the figure with whatever is in the plots
    def make_fig(self, **kwargs): 
        # makes all parameters overriding the defaults with kwargs
        parameters = self.override_params(**kwargs)

        #first make all the plots
        for plot in self.plots.values():
            plot.make_plot(self.data, **parameters)
            if debug.enabled and debug.level <= 0:
                print(f"Made plot {plot.name} \n")

        #then post process
        self.draw_post(**parameters)
        self.set_fontsize(**parameters)

        self.fig.tight_layout()
        #display(self.fig) #shows the fig without clearing it
        return self.fig
        #plt.show(self.fig)

 
   #scales the figure to be more close to target_size while keeping aspect ratio
    def rescale_figure(self,target_size):
        width = self.fig.get_figwidth()
        height = self.fig.get_figheight()
        aspect = width/height

        if aspect > 1:
            new_width = target_size
            new_height = target_size/aspect
        else:
            new_height = target_size
            new_width = target_size*aspect
        
        self.fig.set_figwidth(new_width)
        self.fig.set_figheight(new_height)
        if debug.enabled and debug.level <= 1:
            print(f"  Rescaled figure from {width}x{height} to {new_width}x{new_height}")
   
    def add_post(self,func,**kwargs):
        #adding the function to the post_processing 
        name = func.__name__
        self.post_process[name] = func

        if debug.enabled and debug.level <= 2:
            print(f"Added post processing function {name}")

    
    def draw_post(self,**kwargs):
        for func in self.post_process.values():
            func(self,**kwargs)
            if debug.enabled and debug.level <= 0:
                print(f"    Post processed with {func.__name__}")

    def set_fontsize(self, **kwargs):
        parameters = self.override_params(**kwargs)
        fontsize = parameters.get('fontsize', None)

    # Safely adjust suptitle if it exists
        if hasattr(self.fig, '_suptitle') and self.fig._suptitle is not None and fontsize is not None:
            self.fig._suptitle.set_fontsize(fontsize)
        
        #adjusts each plot's fontsizes
        for plot in self.plots.values():
            plot.set_fontsize(fontsize=fontsize)

    def add_parameters(self, plots, **kwargs):
        if isinstance(plots, str):
            plots = [plots]
        for name,value in kwargs.items():
            for plot in plots:
                self.plots[plot].parameters[name] = value

    def __str__(self):
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"figure_plotting with kwargs: {kwargs_str}"


# %% [markdown]
# # plotting.py

# %%

#from plotting_files.plotting_classes import *

# %%
###################################################################################################
''' 
This function takes a dictionary of parameters and a function object
and returns a dictionary of parameters that match the parameters of the function

'''
def match_param(parameters, obj):
    assert isinstance(parameters, dict), "parameters must be a dictionary"
    
    signature = inspect.signature(obj)
    obj_params = set(signature.parameters.keys())
    # access the callable parameters of the object
    # sadly does not include kwargs possibilities

    accepts_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in signature.parameters.values())
    if accepts_kwargs:
        return parameters

    matching_params = {}

    for attr in parameters.keys():
        if attr in obj_params:
            matching_params[attr] = parameters[attr]

    if debug.enabled and debug.level <= 0:
        print(f"Matched parameters: {matching_params}")
    
    
    return matching_params


    

# %%
#########################################################
#########################################################
# Here there be the main functions for plotting
'''
The Colorplot is a plot_type function
so it needs as input an object of the apt_plot class
as well as the data object
This function plots a pcolormesh plot on the axis
and includes a colorbar nicely placed on the right
'''
def colorplot(apt_plot_object,data,**kwargs): 
    ax = apt_plot_object.ax
    #plotting the colorplot
    params = match_param(kwargs, ax.pcolormesh)
    if debug.enabled and debug.level <= 0:
        print(f"{apt_plot_object.name} is plotting colorplot with parameters {params}")
    
    c = ax.pcolormesh(data.x1, data.x2, apt_plot_object.func(data), **params)
    
    #include the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05) ###############################consider making this a parameter inside fld_val
    cbar = plt.colorbar(c, cax=cax)
    params = match_param(kwargs, cbar.ax.tick_params)
    cbar.ax.tick_params(**params)

    #sets the colorbar object to the plot object
    apt_plot_object.cbar = cbar


    return c

# %%
####################################################################

def line_plot(apt_plot_object,data,**kwargs):
    ap = apt_plot_object
    ax = ap.ax
    params = match_param(kwargs, ax.plot)
    if debug.enabled:
        print(f"{apt_plot_object.name} is plotting equator_Bfield with parameters {params}")
    
    c = ax.plot(data.x1, apt_plot_object.func(data), **params)
    return c

# %%
#########################################################################

'''
default_figure is a function that creates a default figure
this requires an input of the figure_plotting object
it will merely create and return a fig,ax of a certain type

ax is limited to be a column or a row vector for simplicity
'''
#This default figure only allows for horizontal or vertical orientation
##i.e no 2x2 or 2x3 etc
def default_figure(figure_plot, **kwargs):
    fp = figure_plot


    orient = kwargs.get('orientation', 'horizontal')

    data = fp.data
    
    #Loops through the field quantities to find the appropriate size of the figure
    total_xsize = 0
    largest_xsize = 0
    total_ysize = 0
    largest_ysize = 0
    for fld_val in fp.fld_quantities.values():

        #finds the xsize of the subplot
        if fld_val.xlim is not None:
            xsize = fld_val.xlim[1] - fld_val.xlim[0]
        else:
            xsize = np.ceil(data.x1.max()) - np.floor(data.x1.min())

        #finds the ysize of the subplot
        if fld_val.ylim is not None:
            ysize = fld_val.ylim[-1] - fld_val.ylim[0]
        else:
            ysize = np.ceil(data.x2.max()) - np.floor(data.x2.min())
        
        #updates the largest,total sizes
        total_xsize += xsize
        largest_xsize = xsize if xsize > largest_xsize else largest_xsize
        total_ysize = ysize if ysize > total_ysize else total_ysize
        largest_ysize += ysize

    if orient == 'horizontal':
        figsize = total_xsize, largest_ysize
        fig, ax  = plt.subplots(ncols=fp.num_flds)
    elif orient == 'vertical':
        figsize = largest_xsize, total_ysize
        fig, ax  = plt.subplots(nrows=fp.num_flds)
    else:
        raise ValueError("orient must be either 'horizontal' or 'vertical'")
    
    #figsize = rescale_figure(figsize, target_size=kwargs.get('target_figsize', 13))
    if debug.enabled: print(f"Default figure size: {figsize}")

    fig.set_size_inches(figsize)
    
    if fp.num_flds == 1:
        ax = [ax]
    return fig, ax
     


# %%
#########################################################################
#########################################################################
# Here there be the axis functions

'''
draw_field_lines is a post_processing function
so it requires the apt_fig object as input

Current problems is it required Bp inside the config of the dataset
'''
def draw_field_lines_sph(apt_fig,**kwargs):
    
    dataset = apt_fig.data

    Bp      = kwargs.get('Bp',dataset.conf["Bp"])
    flux    = np.cumsum(dataset.B1 * dataset._rv * dataset._rv * np.sin(dataset._thetav) * dataset._dtheta, axis=0)
    clevels = np.linspace(0.0, np.sqrt(Bp), 10)**2
    clevels = np.linspace(0.0, Bp, 10)
    
    for plot in apt_fig.plots.values():    
        plot.ax.contour(dataset.x1, dataset.x2, flux, clevels, colors='green', linewidths=1)
    
    

# %%
#########################################################################

'''
draw_NS is an  post_processing function
so it requires the apt_fig object as input

This function will merely draw the neutron star on every plot
'''
def draw_NS(apt_fig,**kwargs): #just need dataset for convention
    r = kwargs.get("Radius",1)
    for plot in apt_fig.plots.values():
        plot.ax.add_patch(plt.Circle((0,0),r,fill=True, color="black", alpha=0.5))


# %%
#########################################################################
#########################################################################
# Here there be the figure functions

'''
this function just returns the top right subplot
this is used for the time function
'''
def find_top_right_subplot(fig):
    max_x1 = -1
    max_y1 = -1
    top_right_subplot = None
    
    for ax in fig.axes:
        bbox = ax.get_position()
        # Check if this subplot is further right and up than the current max
        if bbox.x1 >= max_x1 and bbox.y1 >= max_y1:
            max_x1, max_y1 = bbox.x1, bbox.y1
            top_right_subplot = ax
            
    return top_right_subplot


# %%
#########################################################################

'''
show_time is a figure function
so it requires as input  the figure to plot on
as well as the dataset being plotted

This draws the time on the top right of the figure
'''
def show_time(fig,dataset,**kwargs):
    time = round(dataset.time,2)
    width,height = fig.get_size_inches()
    add_height = 2 #adjustment for the text
    fig.set_size_inches(width, height+add_height, forward=True)
    # Adjust the top margin to make space for the text
    top_right = find_top_right_subplot(fig)
    #fig.subplots_adjust(top=0.85) 
    top_right.text(0.0, 1.1, r"$t = %a units $" % time, fontsize=14, transform=top_right.transAxes)
    #fig.text(x,y,f"t = {time:.2f}", fontsize=24,color = 'blue')
    return fig

# %%
#########################################################
def EdotB(name='EdotB'):
    return apt_plot(
                     lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3,
                     name = name,
                     plot_function = colorplot,
                     #optional
                     vmin = -1, #default vmin/vmax forthis quantity
                     vmax = 1,
                     )

