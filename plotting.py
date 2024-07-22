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
        self.plot_type = None # plot type for use in updating i.e pcolormesh, lineplot, etc
        self.parameters = {} # all possible parameters for the plot
        self.position = None # a (row,col) tuple (starts at 0->N-1)
        self.made = False
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

    def update_plot(self, data, **kwargs):
        parameters = self.override_params(**kwargs)
        if self.plot_object is None:
            raise ValueError("No plot object to update, make_fig first")
        global colorplot
        global equator_plot
        if self.plot_type == "colorplot":
            self.plot_object.set_array(self.func(data).flatten())
            '''# Ensure the colorbar is updated to reflect the new data range
            # This assumes self.cbar exists and is the colorbar associated with self.plot_object
            if hasattr(self, 'cbar'):
                # Update the colorbar's limits based on the new data
                self.cbar.update_normal(self.plot_object)
                # Redraw the colorbar if necessary
                #self.cbar.draw_all()'''
        elif self.plot_function == "lineplot":
            self.plot_function(self, data, **parameters)
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
    def __init__(self, name, post_func, update_func=None, **kwargs):

        self.post_func = post_func
        self.update_func = update_func

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
        self.post_func(apt_fig_obj, **parameters)

    def update_post(self, apt_fig_obj, **kwargs):
        parameters = self.override_params(**kwargs)
        if self.update_func is None:
            if debug.enabled and debug.level <= 0:
                print(f"    {self.name} not updated")
        else:
            self.update_func(apt_fig_obj, **parameters)

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
        #self.fig.set_label(self.unique_ident)

        if debug.enabled and debug.level <= 1:
            print(f"  Reloaded figure to {num_rows}x{num_columns}, with {len(list(self.plots))+1} subplots")
    
    def construct_plot_obj(self,key,plot_function, **kwargs):
        # default print type for fld obj is colorplot
        global colorplot
        if debug.enabled and debug.level <= 1:
            print(f"Constructing plot object for {key}")
            
        name = kwargs.get('name', key)
        # remove from kwargs
        kwargs.pop('name', None)

        return apt_plot(lambda data: getattr(data, key)
                      , name = name
                      , plot_function = plot_function
                      , title = key
                      , **kwargs)
    
    def add_plot(self,name,pos=None,plot_function = colorplot, datakey=None, **kwargs):
        
        global apt_plot_types
        plot = name
        if isinstance(plot, str):
            # first we check if the str is a data.keys
            if plot in self.data.keys:
                ap = self.construct_plot_obj(plot,plot_function, **kwargs)
            elif datakey is not None:
                ap = self.construct_plot_obj(datakey,plot_function,name =name, **kwargs)
            elif plot in apt_plot_types:
                ap = apt_plot_types[plot](**kwargs)
            else:
                raise ValueError(f"{plot} is not in apt_plot_types, try adding it to the dictionary")
        else:
            ap = plot
            #enforces that apt_plot_object is an apt_plot object
        assert isinstance(ap, apt_plot), "plot must be a apt_plot object, try e.g EdotB() not EdotB or add datakey argument"
        

        name = ap.name
        #if the plot already exists, raise an error
        if name in self.plots:
            raise ValueError(f"{name} already exists as plot (if you are trying to plot the same data but with different plot_function), try adding datakey to override the name")

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

    
    def add_post(self,apt_post_obj,**kwargs):
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
            self.set_fig_shape(self.rows,self.columns)

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

        self.fig.tight_layout()
        self.made = True # marks that the figure has been made
        return self.fig
    
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
            self.update_fig(step, **kwargs)
            self.fig.savefig(f"movie_plots/{step:05d}.png")

        os.system(f"ffmpeg -y -loglevel error -r 10 -i movie_plots/%05d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p movies/{save_name}.mp4")
                        

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
    except AttributeError as e:
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
    c = run_function_safely(ax.pcolormesh, data.x1, data.x2, ap.func(data), **params)
    
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

def equator_plot(apt_plot_object,data,**kwargs):
    ap = apt_plot_object
    ax = ap.ax
    # hardcode B3 at equator for B3 vs r plotting
    # mainly to test line_plots
    # I can figure out automated later
    
    #first finding the equator index
    ths = data._theta # just the vector of theta values
    equator_index = np.argmin(np.abs(ths-np.pi/2))

    # from the func we need to isolate the equator
    func = ap.func(data) # the fld function of the data
    equator = func[equator_index] #isolated to equatorial plane

    rs = data._rv[equator_index]

    # save values for use in updating
    
    ap.plot_type = 'lineplot' #to allow for lineplot updating
    ap.xdata = rs
    ap.ydata = equator

    # now we plot the line
    # this allows us to call the function again to update the line data
    # without redrawing the whole plot as xdata,ydata is the pertinent data
    if hasattr(ap, 'linemade'):
        pass
    else:
        
        params = match_param(kwargs, ax.plot)
        line, = run_function_safely(ax.plot,rs, equator, **params)
        
        ap.linemade = True
        return line

# %%
afig.plots["EdotB"].ax.get_children()

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
    
    def func(apt_fig,**kwargs):
        data = apt_fig.data
        Bp = kwargs.get('Bp',data.conf["Bp"])
        flux    = np.cumsum(data.B1 * data._rv * data._rv * np.sin(data._thetav) * data._dtheta, axis=0)
        clevels = np.linspace(0.0, np.sqrt(Bp), 10)**2
        
        for plot in apt_fig.plots.values():
            contours = plot.ax.contour(data.x1, data.x2, flux, clevels, colors='green', linewidths=1)
            setattr(plot, name, contours)
            # to access for update later
        return None
    
    def update_func(apt_fig,**kwargs):
        
        for plot in apt_fig.plots.values():
            if hasattr(plot, name):
                for c in getattr(plot, name).collections:
                    c.remove()
        func(apt_fig,**kwargs)
            
        if debug.enabled and debug.level <= 1:
            print(f"Updated {name}")
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
    def func(apt_fig,**kwargs): 
        r = kwargs.get("Radius",1)
        for plot in apt_fig.plots.values():
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

def EdotB(name='EdotB',**kwargs):
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
                     plot_function = equator_plot,
                     xlabel = "r",
                     ylabel = r"$\vec{E} \cdot \vec{B}$",
                     **kwargs
                     )
apt_plot_types['EdotB_eq'] = EdotB_eq
