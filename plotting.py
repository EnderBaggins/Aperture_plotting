import matplotlib
# matplotlib.rc("text", usetex=True)
matplotlib.rc("font", family="serif")
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
import inspect
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
import os
import sys
from typing import Any
sys.path.append('.')
from scipy.ndimage import gaussian_filter1d

from datalib_logsph import DataSph, flag_to_species
from datalib import Data

from matplotlib.colors import LinearSegmentedColormap

cdata = {
    "red": [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (0.55, 1.0, 1.0), (1.0, 1.0, 1.0),],
    "green": [(0.0, 1.0, 1.0), (0.45, 0.0, 0.0), (0.55, 0.0, 0.0), (1.0, 1.0, 1.0),],
    "blue": [(0.0, 1.0, 1.0), (0.45, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0),],
}

hot_cold_cmap = LinearSegmentedColormap("hot_and_cold", cdata, N=1024, gamma=1.0)
if "hot_and_cold" not in matplotlib.colormaps:
    matplotlib.colormaps.register(hot_cold_cmap)


################################################################################

#A debug class object
class debugging:
    def __init__(self, enabled=False):
        self.enabled = enabled
        self.level = 0 #allows for different debug messages

debug = debugging(enabled=False)


################################################################################

#Defines the function that will be used to calculate the quantity
#This class will also store the arguments necessary for plotting the function
class apt_plot:
    def __init__(self,  func,name,plot_function,**kwargs):
        self.func = func
        self.plot_object = None
        self.plot_function = plot_function # e.g colorplot using pcolormesh
        self.parameters = {}
        #self.axis = None     # the axis object where the plot will be drawn
        self.position = None # a (row,col) tuple
        self.name = name
        self.ax = None

        self.set_default_parameters()
        self.parameters.update(kwargs) #override the defaults
    '''
    # updating the set_attr to call
    # functions every time I set a value
    def __setattr__(self, name, value):
        #This is necessary to avoid infinite recursion
        super().__setattr__(name, value) 
    
    def __delattr__(self, name):
        super().__delattr__(name)
    '''
    def set_default_parameters(self):
        self.parameters['cmap'] = 'hot_and_cold'
        self.parameters['aspect'] = 'equal'
        self.parameters['fontsize'] = 24

    def make_plot(self,data):
        # creates the plot_object desired
        self.plot_object = self.plot_function(self,data, **self.parameters)
        # sets the parameters of the axis
        self.set_plot_attr()
    
    def set_plot_attr(self):
        # This requires that fld_val has the following attributes:
        # xlim, ylim, aspect, title, fontsize
        attrs = ['xlim', 'ylim', 'aspect', 'title']#, 'fontsize']
        for param in self.parameters:
            if param in attrs:
                try:
                    getattr(self.ax, f"set_{param}")(self.parameters[param])
                except Exception as e:
                    print(f"Could not set {param}: {e}")

    def copy_ax_attrfail(self,new_ax,figure):
        #This is a failed attempt because not all artists can be copied
        ax = self.ax

        ax.remove() # separate the axis from the figure
        for artist in ax.get_children():
            new_ax.add_artist(artist.clone())
        
        if not ax.figure.axes:
            plt.close(ax.figure)

        self.ax = new_ax

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
           # "xticks",
            #"yticks",
            #"xticklabels",
            #"yticklabels"

        ]
        #filtering for get_ with set_ by the attrs
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

    def __str__(self):
        # A overinformative string representation
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"fld_quantity with function {self.func.__name__} and kwargs: {kwargs_str}"
    
    ################################################################################
'''
apt_fig is a one dataset class object storing all the information necessary to plot a figure of subplots
this will have children objects that will store the information necessary to plot the data of one subplot
'''
class apt_fig:
    def __init__(self,data, **kwargs):
        self.fig = plt.figure()
        #self.axes = [] # tightly connected with plots # of the form {axis,position}
        self.plots = {}# has axis information inside the class objects
        self.parameters = {}
        self.data = data #Only one dataset in apt_fig
        self.step = 0
        self.draw_on_plot = {}
        self.draw_on_fig = {}

        # storing the information about the shape of the figure
        # and it's subplots. even accounts for empty subplots
        self.columns = 1
        self.rows = 1

        self.parameters.update(kwargs) #override the defaults
    
        self._init_complete = True #This is to allow set_attr to run normally for initialization
    
    
    
    #####################################################
    # changing the default setattr and delattr to do other things
    def __setattr__(self, name: str, value: Any) -> None:
        # loading the data at the timestep
        if name == 'step':
            self.data.load(value)

        # only runs once the object is fully initialized
        if '_init_complete' in self.__dict__:
            pass
        #This is necessary to avoid infinite recursion
        super().__setattr__(name, value)
    
    def __delattr__(self, name: str) -> None:
        #updates fld_quantity if the thing added is a field val
        super().__delattr__(name)
    ####################################################
    
    # checks if this position is taken
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
            plot.copy_ax_attr(new_ax) #this sets plot.ax to this new axis with the old properties
            plot.make_plot(self.data)
        plt.close(self.fig)
        self.fig = new_fig
        self.fig.tight_layout()
        if debug.enabled:
            print(f"Reloaded figure to {num_rows}x{num_columns}, with {len(list(self.plots))+1} subplots")

                
        
    def add_plot(self,apt_plot_object,pos=None, **kwargs):
        ap = apt_plot_object
        #enforces that apt_plot_object is an apt_plot object
        assert isinstance(ap, apt_plot), "apt_plot_object must be a apt_plot object, try e.g EdotB() not EdotB"
        name = ap.name
        #if the plot already exists, raise an error
        if name in self.plots:
            raise ValueError(f"{name} already exists as plot")

        
        # if no position is given, add to a new column on top
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
        
        self.resize_row_col(pos)
        
        # check if the position is already taken
        if self.check_position_taken(pos):
            # consider making this move the old plot to a new position
            # or this plot to a nearby one
            raise ValueError(f"Position {pos} is already taken")

        
        self.set_fig_shape(self.rows,self.columns)
            
        # add the plot to the dictionary
        self.plots[name] = ap
        ap.position = pos
        #setting the new ax
        ap.ax = plt.subplot2grid((self.rows,self.columns),pos,fig =self.fig)
        print(type(ap))
        ap.make_plot(self.data)


        if debug.enabled:
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
        if debug.enabled:
            print(f"Deleted plot {name}")
        #reshapes the figure
        self.set_fig_shape(self.rows,self.columns)
        

        


    def move_plot(self, name, pos):
        #first asserts that the plot exists
        if name not in self.plots:
            raise ValueError(f"{name} does not exist as plot")
        
        #check if the position is already taken
        if self.check_position_taken(pos):
            raise ValueError(f"Position {pos} is already taken, implement swapping instead?")
        
        self.resize_row_col(pos)
        #moves the plot to the new position
        self.plots[name].position = pos
        self.set_fig_shape(self.rows,self.columns)
        if debug.enabled:
            print(f"Moved plot {name} to position {pos}")

    
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
   



    def __str__(self):
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"figure_plotting with kwargs: {kwargs_str}"

###################################################################################################
''' 
This function takes a dictionary of parameters and a function object
and returns a dictionary of parameters that match the parameters of the function

'''
def match_param(parameters, obj):
    assert isinstance(parameters, dict), "parameters must be a dictionary"
    
    obj_params = set(inspect.signature(obj).parameters.keys())
    # access the callable parameters of the object
    # sadly does not include kwargs possibilities
    matching_params = {}

    for attr in parameters.keys():
        if attr in obj_params:
            matching_params[attr] = parameters[attr]

    return matching_params


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
    if debug.enabled:
        print(f"{apt_plot_object.name} is plotting colorplot with parameters {params}")
    
    c = ax.pcolormesh(data.x1, data.x2, apt_plot_object.func(data), **params)
    
    #include the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05) ###############################consider making this a parameter inside fld_val
    cbar = plt.colorbar(c, cax=cax)
    params = match_param(kwargs, cbar.ax.tick_params)
    cbar.ax.tick_params(**params)
    return c

#########################################################
# here are the default apt_plot object generators
def EdotB(name='EdotB'):
    return apt_plot(
                     lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3,
                     name = name,
                     plot_function = colorplot,
                     #optional
                     vmin = -1, #default vmin/vmax forthis quantity
                     vmax = 1,
                     )



    