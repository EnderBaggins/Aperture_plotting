# %% [markdown]
# # plotting_classes.py

# %%
import matplotlib
import matplotlib.pyplot as plt #cause I use both plt and random mpl functions
import numpy as np


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
matplotlib.rc("text", usetex=True)

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
################################################################################
'''
Defines the function that will be used to calculate the quantity
This class will also store the arguments necessary for plotting the function
This will hold a reference to it's data object
'''
class apt_plot:
    def __init__(self,data, fld_val, plot_function,name=None,**kwargs):
        ###self.func = func
        self.plot_object = None
        self.plot_function = plot_function # e.g colorplot using pcolormesh
        self.plot_type = None # plot type for use in updating i.e pcolormesh, lineplot, etc
        self.parameters = {} # all possible parameters for the plot
        self.position = None # a (row,col) tuple (starts at 0->N-1)
        self.made = False

        self.name = name
        self.ax = None
        self.data=data
        # if hasattr(data,'dataname'):
        #     self.dataname = data.dataname

        self._construct_plot(fld_val,**kwargs)
        self.set_default_parameters()
        self.parameters.update(kwargs) #override the defaults
    
    #construct_plot will be called with __init__ to create the plot
    # it will account for both data.keys apt_plot_types and lambda functions
    # to appropriately set the fld_val and name
    def _construct_plot(self,fld_vals,**kwargs):
        # key can be a list of multiple keys
        # in which case there should also have been a list of fld_vals and labels

        global apt_plot_types
        #first we tackle if there is a single of each
        if not isinstance(fld_vals,list):
            if callable(fld_vals):
                self.fld_val = fld_vals #as its a lambda function
                if self.name is None:
                    raise ValueError(f"Name must be specified for lambda function fld_val")
            
            elif fld_vals in apt_plot_types:
                if self.name is None:
                    self.fld_val,self.name = apt_plot_types[fld_vals]()#does not override name
                else:
                    self.fld_val,self.name = apt_plot_types[fld_vals](self.name)#overrides name
            
            elif fld_vals in self.data.keys:
                self.fld_val = lambda data: getattr(data,fld_vals)
                if self.name is None:
                    self.name = fld_vals #since fld_vals is a string it can identify
            elif fld_vals is None:
                pass #This is for plotting functions that have no fld_val (i.e particle stuff)
            else:
                raise ValueError(f"Key {fld_vals} not found in apt_plot_types, data.keys or is not a lambda function")

        # having fld_vals be a list means labels exists and is a list and data is also a list
        # this is for lineplots where you want multiple in one subplot
        # in general the plot_funciton you input needs to be able to deal with the list of items
        else:
            fld_val_list = []
            labels = kwargs.get('labels',None)
            if labels is None:
                raise ValueError(f"Labels must be specified for multiple fld_vals")
            for fld_val,label,datum in zip(fld_vals,labels,self.data):

                if callable(fld_val):
                    
                    print(type(fld_val))
                    fld_val_list.append(fld_val) #as its a lambda function
                elif fld_val in apt_plot_types:
                    print(type(fld_val))
                    fv,name = apt_plot_types[fld_val](label)
                    fld_val_list.append(fv)
                    #name will be default as the label it was given
                elif isinstance(fld_val,str) and fld_val in datum.keys:
                    
                    fv = lambda data: getattr(data,fld_val)
                    fld_val_list.append(fv)
                    #self.name = 
                elif fld_val is None:
                    pass #This is for plotting functions that have no fld_val (i.e particle stuff)
                else:
                    raise ValueError(f"Key {fld_val} not found in apt_plot_types, data.keys or is not a lambda function")
                
            self.fld_val = fld_val_list
            if self.name is None:
                raise ValueError(f"Name must be specified for multiple fld_vals")

        
        
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

    def update_plot(self, data, **kwargs):
        parameters = self.override_params(**kwargs)
        if self.plot_object is None:
            raise ValueError("No plot object to update, make_fig first")
        
        if self.plot_type == "colorplot":
            self.plot_object.set_array(self.func(data).flatten())
            
        elif self.plot_type == "lineplot":
            # any lineplot function must specify xdata and ydata as attributes
            # of the apt_plot object and only run the plotting when linemade is false
            
            self.plot_function(self, data, **parameters)
            if isinstance(self.plot_object,list):
                #then there is a list of different lines stores to update
                #where xdata and ydata are also lists
                for i, line in enumerate(self.plot_object):
                    line.set_data(self.xdata[i], self.ydata[i])
            else:
                self.plot_object.set_data(self.xdata, self.ydata)

        else:
            raise ValueError("Update not implemented for this plot function")


        #self.set_plot_attr(**parameters)
        
    def set_plot_attr(self, **kwargs):
        parameters = self.parameters.copy()
        parameters.update(kwargs)
        # This requires that fld_val has the following attributes:
        # xlim, ylim, aspect, title
        attrs = ['xlim', 'ylim', 'aspect', 'title',"xlabel",'ylabel']
        for param in parameters:
            if param in attrs:
                try:
                    getattr(self.ax, f"set_{param}")(parameters[param])
                except Exception as e:
                    print(f"Could not set {param}: {e}")

        if parameters.get("legend",False):
            self.ax.legend()
        

    
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
        fontsize = parameters.get('fontsize', None)
        
        #each possible fontsize type with default as the fontsize
        label_fontsize = parameters.get('label_fontsize', fontsize)
        title_fontsize = parameters.get('title_fontsize', fontsize)
        tick_fontsize = parameters.get('tick_fontsize', fontsize)
        legend_fontsize = parameters.get('legend_fontsize', fontsize)

        ctick_fontsize = parameters.get('ctick_fontsize',tick_fontsize)

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
       
        self.ax.tick_params(labelsize=tick_fontsize)

        # set font size for colorbar
        if hasattr(self, 'cbar'):
                self.cbar.ax.tick_params(labelsize=ctick_fontsize) if tick_fontsize is not None else None


    
    def __str__(self):
        # A overinformative string representation
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"fld_quantity with function {self.func.__name__} and kwargs: {kwargs_str}"
    

# %%
################################################################################
'''
apt_post is a class that will be used to make single post processing functions
allowing for updating of the plot as well
'''
class apt_post:
    def __init__(self, name, post_func=None, update_func=None, **kwargs):

        self.post_func = post_func
        self.update_func = update_func
        self.post_plots = None # the list of plots to run post on

        self.parameters = {}
        self.name = name

        #self.set_default_parameters() # idk what default for drawing post process
        self.parameters.update(kwargs) #override the defaults

    def override_params(self, **kwargs):
        parameters = self.parameters.copy()
        parameters.update(kwargs)

        if debug.enabled and debug.level <= 0:
            overridden = {k: (self.parameters[k], v) for k, v in kwargs.items() if k in self.parameters and self.parameters[k] != v}
    
            for key, (original_value, new_value) in overridden.items():
                print(f"    {self.name}'s Parameter '{key}' was overridden: Original = {original_value}, New = {new_value}")
    
        return parameters

    def make_post(self,apt_fig_obj, **kwargs):
        parameters = self.override_params(**kwargs)
        # need to pass self into these functions
        # to allow the post_func to access self parameters
        self.post_func(self,apt_fig_obj, **parameters)

    def update_post(self, apt_fig_obj, **kwargs):
        parameters = self.override_params(**kwargs)
        if self.update_func is None:
            if debug.enabled and debug.level <= 0:
                print(f"    {self.name} not updated")
        else:
            self.update_func(self,apt_fig_obj, **parameters)

    def set_fontsize(self, **kwargs):
        parameters = self.override_params(**kwargs)
        fontsize = parameters.get('fontsize', None)
        
        #each possible fontsize type with default as the fontsize
        label_fontsize = parameters.get('label_fontsize', fontsize)
        title_fontsize = parameters.get('title_fontsize', fontsize)
        tick_fontsize = parameters.get('tick_fontsize', fontsize)
        legend_fontsize = parameters.get('legend_fontsize', fontsize)

# %%
################################################################################
aperture_figure_objects = {}
# this is all the apt_fig objects that have been created
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
            plt.close(aperture_figure_objects[unique_ident].fig)
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
        self.made = False      # marks whether make_fig has been run
        # storing the information about the shape of the figure
        # and it's subplots. even accounts for empty subplots
        self.columns = 1      # number of columns
        self.rows = 1         # number of rows

        self.parameters.update(kwargs) #override the defaults
    
    
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
    def set_fig_grid(self,num_rows,num_columns):
        new_fig = plt.figure()
        #copies over the old axes
        for plot in self.plots.values():
            pos = plot.position

            new_ax = plt.subplot2grid((num_rows,num_columns),pos,fig= new_fig)
            # setting plot.ax to new axis with old properties
            plot.copy_ax_attr(new_ax) 

        plt.close(self.fig) #closes the old figure
        self.fig = new_fig
        #self.fig.set_label(self.unique_ident)

        if debug.enabled and debug.level <= 1:
            print(f"  Reloaded figure to {num_rows}x{num_columns}, with {len(list(self.plots))+1} subplots")
    
    '''def construct_plot_obj(self,key,plot_function, **kwargs):
        # default print type for fld obj is colorplot
        global colorplot
        if debug.enabled and debug.level <= 1:
            print(f"Constructing plot object for {key}")
            
        name = kwargs.get('name', key)
        title = kwargs.get('title', key)
        # remove from kwargs
        kwargs.pop('name', None)
        kwargs.pop('title', None)

        return apt_plot(lambda data: getattr(data, key)
                      , name = name
                      , plot_function = plot_function
                      , title = title
                      , **kwargs)
    '''
    def add_plot(self,key,pos=None, plot_function = None, data=None, **kwargs):
        if plot_function is None:
            plot_function = colorplot # default plot function

        if data is None:# allows overriding of fig data for subplot
            data = self.data
    

        if not isinstance(key,apt_plot):
            ap = apt_plot(data,key,plot_function, **kwargs)
        else:
            ap = key
        assert isinstance(ap, apt_plot), "plot must be a apt_plot object, try e.g EdotB() not EdotB or add datakey argument"
        

        name = ap.name
        #if the plot already exists, raise an error
        if name in self.plots:
            raise ValueError(f"{name} already exists as plot (if you are trying to plot the same data with key but with different plot_function), try adding name to change the name")

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
        self.set_fig_grid(self.rows,self.columns)
            
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
        self.set_fig_grid(self.rows,self.columns)
    
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
        self.set_fig_grid(self.rows,self.columns)
        if debug.enabled and debug.level <= 2:
            print(f"Moved plot {name} to position {pos}")

    
    def add_post(self,apt_post_obj,add_to = "all",**kwargs):

        if add_to == "all":
            post_plots = list(self.plots.values())
        elif isinstance(add_to, str):
            post_plots = [self.plots[add_to]]
        else:
            post_plots = [self.plots[plot] for plot in add_to]
        
        global apt_post_types
        # allows for a list of the objects
        if not isinstance(apt_post_obj, list):
            apt_post_obj = [apt_post_obj]
        
        
        for post in apt_post_obj:
            if isinstance(post, str):
                if post in apt_post_types:
                    ap = apt_post_types[post](**kwargs)
            else:
                ap = post
            assert isinstance(ap, apt_post), "apt_post_obj must be a apt_post object"

            name = ap.name
            if name in self.post_process:
                raise ValueError(f"{name} already exists as post process")
            
            self.post_process[name] = ap
            ap.post_plots = post_plots
            
            if debug.enabled and debug.level <= 2:
                print(f"Added post processing function {name}")

    def del_post(self, name):
        if name not in self.post_process:
            raise ValueError(f"{name} does not exist as post process")
        del self.post_process[name]
        if debug.enabled and debug.level <= 2:
            print(f"Deleted post processing function {name}")
       
    # this just updates the figure with whatever is in the plots
    def make_fig(self, **kwargs): 
        # first if it is already made, close the figure and remake it
        # this is to stop drawing over the old figure
        if self.made:
            plt.close(self.fig)
            # then redo the figure
            self.set_fig_grid(self.rows,self.columns)

        # makes all parameters overriding the defaults with kwargs
        parameters = self.override_params(**kwargs)
        
        if kwargs.get("step",None) is not None:
            self.step = kwargs["step"]
        #first make all the plots
        for plot in self.plots.values():
            plot.make_plot(self.data, **parameters)
            if debug.enabled and debug.level <= 0:
                print(f"Made plot {plot.name} \n")

        #then post process
        for post in self.post_process.values():
            post.make_post(self, **parameters)
            if debug.enabled and debug.level <= 0:
                print(f"  Post processed with {post.name}")
        
        # fontsize is rough, so it has its own methods
        self.set_fontsize(**parameters)
        self.set_size()
        self.fig.tight_layout()
        self.made = True # marks that the figure has been made
        return self.fig
    

    def set_size(self, xsize = None,ysize = None,**kwargs):

        if xsize is not None and ysize is not None:
            shape = [xsize, ysize]
            self.shape = shape
        elif hasattr(self, 'shape'):
            shape = self.shape
        else:
            shape = None
        # This allows for make_fig to call set_size even without having a preset shape
        if shape is not None:
            self.fig.set_size_inches(shape,**kwargs)
        

    def update_fig(self, step = None,set_plot_attr=False, **kwargs):

        if step is not None:
            self.step = step

        parameters = self.override_params(**kwargs)
        # updates all the basic plots
        for plot in self.plots.values():
            plot.update_plot(self.data, **parameters)
            if set_plot_attr:
                plot.set_plot_attr(**parameters)
            if debug.enabled and debug.level <= 1:
                print(f"  Updated plot {plot.name} \n")

        # then reruns post processing (these functions need to override)
        for post in self.post_process.values():
            if post.update_post is not None:
                post.update_post(self, **parameters)
                if debug.enabled and debug.level <= 1:
                    print(f"  Updated post process with {post.name}")
            else:
                if debug.enabled and debug.level <= 0:
                    print(f"    {post.name} cannot be updated")

        self.fig.canvas.draw_idle()
        return self.fig
 
    def make_movie(self,save_name="Untitled", start=0, end=None,increment=1,  **kwargs):
        assert all(isinstance(var, (int, type(None))) for var in [start, end, increment]), "start, end, increment must be integers or None"
        if end is None:
            end = len(self.data.fld_steps)
        #checks if untitled already exists in movies
        if os.path.exists(f"movies/{save_name}.mp4"):
            if debug.enabled and debug.level <= 2:
                print(f"Overwriting {save_name}.mp4")
    
        # first ensures that the figure has been made once
        if not self.made:
            self.make_fig(**kwargs)
            # to set all the right things, this specific plot
            # is not saved, because it does not know which step it is

        if not os.path.exists('movie_plots'):
            os.makedirs('movie_plots') #where temp pngs go
        if not os.path.exists('movies'):
            os.makedirs('movies') #where movies go
        #clears out old plots
        if len(os.listdir("movie_plots"))>0:
            os.system("rm movie_plots/*.png")

        # saves each step as a temp png
        for step in range(start, end, increment):
            save_step = step // increment # this is to make sure the steps are always incrementing by 1
            self.update_fig(step, **kwargs)
            self.fig.savefig(f"movie_plots/{save_step:05d}.png")

        os.system(f"ffmpeg -y -loglevel error -r 10 -i movie_plots/%05d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p -threads 0 movies/{save_name}.mp4")
                        

    
    def set_fontsize(self, **kwargs):
        parameters = self.override_params(**kwargs)
        fontsize = parameters.get('fontsize', None)

    # Safely adjust suptitle if it exists
        if hasattr(self.fig, '_suptitle') and self.fig._suptitle is not None and fontsize is not None:
            self.fig._suptitle.set_fontsize(fontsize)
        
        #adjusts each plot's fontsizes
        for plot in self.plots.values():
            plot.set_fontsize(**parameters)

    def add_parameters(self, plots, **kwargs):
        if plots == "all":
            plots = list(self.plots)
        elif isinstance(plots, str):
            plots = [plots]
        for name,value in kwargs.items():
            for plot in plots:
                self.plots[plot].parameters[name] = value

    # other adding plots that wrap add_plot
    def add_colorplot(self, name = None, fld_val= None,**kwargs):
        global colorplot
        self.add_plot(fld_val, name, colorplot, **kwargs)

    def add_lineout_plot(self, name, fld_vals,  restrictions,labels=None, second_restriction=None,data=None, **kwargs):
        global lineout_plot # this is the lineout function for plotting a lineplot

        #restrictions must be either a list of tuples or a tuple
        #second_restriction is entirely for a 3D dataset to cut down to 1D line
        #makes everything a list for purposes of iterating multiple lineouts on one plot
        if not isinstance(fld_vals, list):
            fld_vals = [fld_vals]

        if not isinstance(labels, list):
            labels = [labels]
        if isinstance(restrictions, tuple):
            restrictions = [restrictions]
        if len(restrictions) != len(fld_vals):
            #use the same restrictions for all lineouts
            restrictions = [restrictions] * len(fld_vals)
        else:
            assert all(isinstance(restriction, tuple) for restriction in restrictions), "restrictions must be a list of tuples"
        # similar for second_restriction
        if second_restriction is not None:
            if isinstance(second_restriction, tuple):
                second_restriction = [second_restriction]
            else:
                assert all(isinstance(restriction, tuple) for restriction in second_restriction), "second_restriction must be a list of tuples"
        
        
        # assert that the lengths of the lists are the same
        assert len(fld_vals) == len(labels) == len(restrictions), "fld_vals, labels, and restrictions must be the same length (they are 1 to 1)"

        #asserting that fld_val is either a lambda function of data or a key in data/plot_types
        for fld_val in fld_vals:
            if callable(fld_val):
                pass
            elif fld_val in self.data.keys:#stored data values
                pass
            elif fld_val in apt_plot_types:#stored apt_plot functions
                pass
            else:
                raise ValueError(f"{fld_val} is not a key in data/apt_plot_types or a lambda function of data")

        if data is None:
            data = [self.data] * len(fld_vals) #if data is not specified, it is the same for all lineouts
        elif not isinstance(data, list):
            data = [data] * len(fld_vals) #if data is not a list, it is the same for all lineouts
        elif len(data) != len(fld_vals):
            raise ValueError("If multiple data input then must be same length as fld_vals")
        if labels is None:
            for fld_val in fld_vals:
                if not isinstance(fld_val, str):
                    raise ValueError("If labels is not specified, fld_vals must be strings, if you do a lambda please supply labels as well")
            labels = fld_vals
        #creating the plot object with the lists (apt_plot account for lists seperately)
        ap = apt_plot(data, fld_vals, name=name,plot_function=lineout_plot,labels=labels, **kwargs)
        
        #saving restriction information
        ap.restrictions = restrictions
        ap.second_restriction = second_restriction
        

        self.add_plot(ap, **kwargs)
        

    def add_particle_hist(self, species, x_key, y_key,data=None,**kwargs):
        # first checking species is valid
        if data is None:
            data = self.data

        if species not in ["p", "e", "electron", "positron"] and species not in [0, 1]:
            raise ValueError(f"Species {species} is not valid, must be 'e' or 'p' or 0 or 1, to add ion change add_particle_plot in apt_fig")
        global particle_hist

        name = kwargs.get('name', f"{y_key}_v_{x_key}")
        kwargs.update({'name': name})
        kwargs.setdefault("xlabel", x_key)
        kwargs.setdefault("ylabel", y_key)   

        aplot = apt_plot(data,None, plot_function = particle_hist, **kwargs)

        #Giving correct parameters for the particle_plot function to access
        aplot.y_key = y_key
        aplot.x_key = x_key
        aplot.species = species

        self.add_plot(aplot,**kwargs)

    def add_spectrum(self,name, species,logscale= False,data=None, **kwargs):
        labels = kwargs.get('labels',None)

        # first convert to list of species
        if isinstance(species, str):
            species = [species]
        elif isinstance(species, int):
            species = [species]
        elif isinstance(species, list):
            pass
        else:
            raise ValueError(f"Species {species} is not valid, must be 'e' or 'p' or 0 or 1, to add ion change add_spectrum in apt_fig")
        
        global spectrum_plot
        
        kwargs.setdefault("xlabel", "$E$")
        if logscale:
            kwargs.setdefault("ylabel", "$EdN/dE$")
        else:
            kwargs.setdefault("ylabel", "$dN/dE$")
        
        if data is None:
            data = [self.data] * len(species) #if data is not specified, it is the same for all lineouts
        elif not isinstance(data, list):
            data = [data] * len(species) #if data is not a list, it is the same for all lineouts
        elif len(data) != len(species):
            raise ValueError("If multiple data input then must be same length as species")
        if labels is None:
            labels = []
            for specie,datum in zip(species,data):
                if specie == "e" or specie == 0:
                    labels.append("Electron")
                elif specie == "p" or specie == 1:
                    labels.append("Positron")
                else:
                    labels.append(f"Species {specie}")

        #name = kwargs.get('name', f"{species}_spectrum")
        fld_vals = [None] * len(species)
        kwargs.update({'name': name}) #need to update so kwargs can pass to add_plot nicely
        aplot = apt_plot(data, fld_vals, plot_function = spectrum_plot,labels = labels, **kwargs)

        #Giving correct parameters for the spectrum_plot function to access
        aplot.species = species
        aplot.logscale = logscale

        self.add_plot(aplot,**kwargs)

    def print_info(self):
        afig = self
        from IPython.display import display, Markdown

        def extract_parameters(obj):
            if hasattr(obj, 'parameters') and isinstance(obj.parameters, dict):
                return obj.parameters
            return {}

        def format_parameters(parameters):
            return "\n".join([f"     {key}: {value}" for key, value in parameters.items()])

        output = []
        output.append(f"**{afig.unique_ident}:**")
        output.append(f"   \n ***Overriding Parameters:***")
        output.append(f"\n{format_parameters(afig.parameters)}\n")
        
        output.append(f"**{afig.unique_ident} has the following plots:**\n")
        for name, plot in afig.plots.items():
            output.append(f"  **{name}**\n")
            if hasattr(plot.data, 'name'):
                output.append(f"    Data: {plot.data.name}")
            if isinstance(plot.data,list):
                for i,data in enumerate(plot.data):
                    output.append(f"    {i+1}: Data: {data.name}")
            #output.append(f"  data: {plot.data}")
            output.append(f"    Position: {plot.position}")
            output.append(f"    Parameters:")
            output.append(f"{format_parameters(extract_parameters(plot))}\n")
            if hasattr(plot, 'restrictions'):
                output.append(f"    **Restrictions:**")
                for i, restriction in enumerate(plot.restrictions):
                    output.append(f"      {i+1}: {restriction}")
                if plot.second_restriction is not None:
                    output.append(f"    **Second Restriction:**")
                    for i, restriction in enumerate(plot.second_restriction):
                        output.append(f"      {i+1}: {restriction}")
            
        output.append(f"**Post-processing functions:**\n")
        for name, post in afig.post_process.items():
            output.append(f"  **{name}**\n")
            output.append(f"    Parameters:")
            output.append(f"{format_parameters(extract_parameters(post))}\n")
            
        
        display(Markdown("\n".join(output)))


    def __str__(self):
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"figure_plotting with kwargs: {kwargs_str}"


# %% [markdown]
# # plotting.py

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


# %%
########################################################
def extract_keyword_from_error(error):
    from re import search
    # Use a regular expression to extract the keyword argument name from the error message
    match = search(r"unexpected keyword argument '(\w+)'", str(error))
    if match:
        # If a match is found, the first group contains the keyword argument name
        keyword_error = match.group(1)
        return keyword_error
    else:
        # If no match is found, return None or handle as appropriate
        return None

def run_function_safely(func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except (AttributeError,TypeError) as e:
        # Attempt to extract the problematic attribute name from the error message
        # This is fragile and depends on the error message format
        key_error = extract_keyword_from_error(e)

        # remove this key from the input arguments the recursionally run the function
        if key_error in kwargs:
            #if debug.enabled and debug.level <= 0:
                #print(f"Removing key {key_error} from kwargs")
            del kwargs[key_error]
            return run_function_safely(func, *args, **kwargs)
        else:
            raise ValueError(f"could not remove key {key_error}: {e}")


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
    ap = apt_plot_object
    import warnings
    ax = ap.ax
    ap.plot_type = 'colorplot'
    params = match_param(kwargs, ax.pcolormesh)
    if debug.enabled and debug.level <= 0:
        print(f"{ap.name} is plotting colorplot with parameters {params}")
    
    warnings.filterwarnings("ignore", category=UserWarning, message="The input coordinates to pcolormesh are interpreted as cell centers.*")
    # passing params into pcolormesh crashes, so we need to remove problematic ones with run_function_safely
    c = run_function_safely(ax.pcolormesh, data.x1, data.x2, ap.fld_val(data), **params)
    
    #include the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(c, cax=cax)
    params = match_param(kwargs, cbar.ax.tick_params)
    #cbar.ax.tick_params(**params)

    #sets the colorbar object to the plot object
    apt_plot_object.cbar = cbar


    return c

# %%
####################################################################

def lineout_plot(apt_plot_object,data,**kwargs):
    ap = apt_plot_object
    ax = ap.ax
    labels = kwargs.get('labels', None)
    #checks if restriced exists
    if hasattr(ap, 'restrictions'):
        restrictions = ap.restrictions
    else:
        raise ValueError("Lineout plot must have restrictions")
    
    #checks if second_restriction exists
    if getattr(ap, 'second_restriction') is not None:
        second_restriction = ap.second_restriction
        raise ValueError("3D with second_restrictions not implemented yet")
    
    
    fld_vals = ap.fld_val
    datas = ap.data

    #assert same length
    if len(fld_vals) != len(restrictions) != len(datas):
        raise ValueError("fld_vals and restrictions must be the same length")

    xdata=[]
    ydata=[]

    for restriction,fld_val,data in zip(restrictions,fld_vals,datas):
        
        res_ax = restriction[0]# the axis to restrict
        val = restriction[1]# the restricting value

        #look at each case separately
        if res_ax == "theta":
            restricted_index = np.argmin(np.abs(data._theta-val))
            ys = fld_val(data)[restricted_index,:]
            xs = data._r
            
        elif res_ax == "r":
            restricted_index = np.argmin(np.abs(data._r-val))
            ys = fld_val(data)[:,restricted_index]
            xs = data._theta
        else:
            raise ValueError("Lineout plot must have restrictions on theta or r, if you did ._theta just do .r")
        '''
        axdata = getattr(data,res_ax)
        restricted_index = np.argmin(np.abs(axdata-val))
        #main issue is restricted index finds a single value
        # i need the whole row/column
        
        fld_val = fld_val(data) # the fld function of the data
        line = fld_val[restricted_index] #isolated to restricted index

        if res_ax is 'x1':
            line = fld_val[restricted_index:]
            rem_ax = 'x2'
            xs = getattr(data,rem_ax)[restricted_index,:]
        elif res_ax is 'x2':
            line = fld_val[:,restricted_index]
            rem_ax = 'x1'
            xs = getattr(data,rem_ax)[:,restricted_index]
        else:
            raise ValueError("Lineout plot must have restrictions on x1 or x2 (3D not implemented yet), if you did ._theta just do .x2")
        '''
        xdata.append(xs)
        ydata.append(ys)

    # save values for use in updating
    
    ap.plot_type = 'lineplot' #to allow for lineplot updating
    ap.xdata = xdata
    ap.ydata = ydata
    # now we plot the line
    # this allows us to call the function again to update the line data
    # without redrawing the whole plot as xdata,ydata is the pertinent data
    if hasattr(ap, 'linemade'):
        pass
    else:
        plot_object = []
        params = match_param(kwargs, ax.plot)
        for i in range(len(xdata)):
            label = labels[i] if labels is not None else None
            print(label)
            line, = run_function_safely(ax.plot,xdata[i], ydata[i],label=label, **params)
            plot_object.append(line)
        
        ap.linemade = True
        return plot_object

# %%
def particle_hist(apt_plot_object,data,**kwargs):
    ap = apt_plot_object
    ax = ap.ax

    # First ensures that the plot object has the necessary keys
    if not all(hasattr(ap, key) for key in ['x_key', 'y_key', 'species']):
        raise AttributeError("One or more required attributes ('x_key', 'y_key', 'species') are missing in the object 'ap'")
    x_key = getattr(ap, 'x_key')
    y_key = getattr(ap, 'y_key')
    species = getattr(ap, 'species')
    # The flag requires species to be 0 or 1 (or the ion option)
    if species != 0 or 1:
        if species == "p"or "positron":
            species = 1
        elif species == "e" or "electron":
            species = 0
        else:
            raise ValueError(f"Species {species} is not valid, must be 'e' or 'p', or 0 or 1")

    # Now we can plot the particle plot with the species flag
    x = getattr(data, f"tracked_ptc_{x_key}")[flag_to_species(data.tracked_ptc_flag) == species]
    y = getattr(data, f"tracked_ptc_{y_key}")[flag_to_species(data.tracked_ptc_flag) == species]

    plot = run_function_safely(ax.hist2d, x, y, bins=200,**kwargs)

    return plot
    #ap.plot_type = 'particle_plot' # need to understand how to update first

# %%


# %%
# Spectra looks at the particle data and plots a histogram of the energy
# will need to consider log scale as well
def spectrum_plot(apt_plot_object,data,**kwargs):
    ap = apt_plot_object
    ax = ap.ax
    ap.plot_type = 'lineplot' #to allow for lineplot updating

    #ap.xdata = rs
    #ap.ydata = line

    # First ensures that the plot object has the necessary keys
    if not all(hasattr(ap, key) for key in ['species']):
        raise AttributeError("One or more required attributes ('species') are missing in the object 'ap'")
    species = getattr(ap, 'species')
    if isinstance(species, str):
        species = [species]
        # this is to enable multiple species to be plotted
    for specie in species:
        if specie == "p"or "positron" or 1:
            specie = 1
        elif specie == "e" or "electron" or 0:
            specie = 0
        else:
            raise ValueError(f"Species {specie} is not valid, must be 'e' or 'p', or 0 or 1")
    # The flag requires species to be 0 or 1 (or the ion option)
    
        # Now we can plot the particle plot with the species flag
        energy = getattr(data, f"tracked_ptc_E")[flag_to_species(data.tracked_ptc_flag) == specie]
        logscale = getattr(ap, 'logscale', False)
        kwargs.setdefault('bins', 300)
        if logscale:
            #override the bins which is defaulted to a number to a logarithmic space
            #bins = run_function_safely(np.logspace,np.log10(energy.min()), np.log10(energy.max()), num = kwargs.get("bins"))#,**kwargs)
            bins = np.logspace(np.log10(energy.min()), np.log10(energy.max()), num = kwargs.get("bins"))
            kwargs['bins'] = bins
        
        #making the histogram
        params = match_param(kwargs, ax.hist)
        # saves the histogram data 
        counts, bin_edges = run_function_safely(np.histogram,energy, **params)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        #doesn't work with the multiple species yet
        ap.xdata = bin_centers
        ap.ydata = counts

        if hasattr(ap, 'linemade'):
                pass
        else:
            
            params = match_param(kwargs, ax.plot)
            line, = run_function_safely(ax.step,bin_centers, counts, **params)
            if logscale:
                ax.set_xscale('log')

            ap.linemade = True # to only run the step function once future times will only update the data
            return line



# %%
apt_post_types = {}

# %%
#########################################################################
#########################################################################
# Here there be the axis functions

'''
draw_field_lines is a generator for post processing class object

its func requires apt_fig as input
This stores the contour lines in each plot object

the update is optional and also requires apt_fig
this deletes the old contour lines and recalls func


Current problems is it required Bp inside the config of the dataset
'''
def draw_field_lines1(name='draw_field_lines1',**kwargs):
    
    def func(self,apt_fig,**kwargs):
        data = apt_fig.data
        Bmax = kwargs.get('Bp',data.conf["Bp"])
        flux    = np.cumsum(data.B1 * data._rv * data._rv * np.sin(data._thetav) * data._dtheta, axis=0)
        clevels = np.linspace(0.0, np.sqrt(Bmax), 10)**2
        #clevels = np.linspace(0.0, Bmax, 10)
        
        for plot in self.post_plots:
            contours = plot.ax.contour(data.x1, data.x2, flux, clevels, colors='green', linewidths=1)
            setattr(plot, name, contours)
            # to access for update later
        return None
    
    def update_func(self,apt_fig,**kwargs):
        
        for plot in self.post_plots:
            if hasattr(plot, name):
                for c in getattr(plot, name).collections:
                    c.remove()
        func(self,apt_fig,**kwargs)
            
        if debug.enabled and debug.level <= 1:
            print(f"Updated {name} at timestep {apt_fig.step}")
        return None
    
    return apt_post(name, func, update_func, **kwargs)
apt_post_types['draw_field_lines1'] = draw_field_lines1
        
    

# %%
#########################################################################

'''
draw_NS is a generator for post_processing class object
so it requires the apt_fig object as input

This function will merely draw the neutron star on every plot
'''
def draw_NS(name='draw_NS',**kwargs):
    def func(self,apt_fig,**kwargs): 
        r = kwargs.get("Radius",1)
        for plot in self.post_plots:
            plot.ax.add_patch(plt.Circle((0,0),r,fill=True, color="black", alpha=0.5))

    return apt_post(name, func, **kwargs)
apt_post_types['draw_NS'] = draw_NS

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


# %%


# %%
#########################################################
# Here there be the possible plotting types
apt_plot_types = {} # a dict of possible
# If you want to add one, you need to add it to this dictionary
# or else call add_plot(EdotB) as oppossed to add_plot("EdotB")

def EdotB(name="EdotB"):
    fld_val = lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3
    return (fld_val,name)


def EdotB22(name='EdotB',**kwargs):
    vmin = kwargs.get('vmin', -1) #default values
    vmax = kwargs.get('vmax', 1)
    return apt_plot(
                     lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3,
                     name = name,
                     plot_function = colorplot,
                     #optional
                     title = r"$\vec{E} \cdot \vec{B}$",
                     vmin = vmin,
                     vmax = vmax,
                     **kwargs
                     )
apt_plot_types['EdotB'] = EdotB

def Epar(name='Epar',**kwargs):
    vmin = kwargs.get('vmin', -1)
    vmax = kwargs.get('vmax', 1)
    return apt_plot(
                     lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3,
                     name = name,
                     plot_function = colorplot,
                     title = r"$(\vec{E} \cdot \vec{B})/(\mathbf{B})$",
                     #optional
                     vmin =vmin, #default vmin/vmax forthis quantity
                     vmax = vmax,
                     **kwargs
                     )
apt_plot_types['Epar'] = Epar





# %%
def EdotB_eq(name='EdotB_eq',**kwargs):

    return apt_plot(
                     lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3,
                     name = name,
                     plot_function = lineout_plot,
                     xlabel = "r",
                     ylabel = r"$\vec{E} \cdot \vec{B}$",
                     **kwargs
                     )
apt_plot_types['EdotB_eq'] = EdotB_eq
