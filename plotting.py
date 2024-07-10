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
class fld_quantity:
    def __init__(self, func,**kwargs):
        self.func = func
        self.kwargs = kwargs
        self.name = kwargs.get('name', None)

        self.set_plot_defaults()

    # updating the set_attr to call
    # functions every time I set a value
    def __setattr__(self, name, value):
        #updates parameter lists
        self.update_colormeshlist() if name != 'cplot_params' else None

        #This is necessary to avoid infinite recursion
        super().__setattr__(name, value) 
    
    def __delattr__(self, name):
        #updates parameter lists
        self.update_colormeshlist()
        
        super().__delattr__(name)

    def set_plot_defaults(self):
        self.cplot_params = {}
        #colorbar defaults
        #don't technically need these
        self.vmin = self.kwargs.get('vmin', None)
        self.vmax = self.kwargs.get('vmax', None)
        self.cmap = self.kwargs.get('cmap', "hot_and_cold")
        self.shading = self.kwargs.get('shading', 'gouraud')
        #self.norm = kwargs.get('norm', colors.Normalize(vmin=self.vmin, vmax=self.vmax)) #This is already default
        #If you change norm you would need to ensure it has the current vmin and vmax
        '''        #I don't think we need these cbar, you can make your own plotting function
        #colorbar placement and scaling #shouldn't need to change these much
        self.cbar_position = kwargs.get('cbar_position', 'right')
        self.cbar_size = kwargs.get('cbar_size', '5%')
        self.cbar_pad = kwargs.get('cbar_pad', '0.05')
        '''
        #plot limits
        self.xlim = self.kwargs.get('xlim', None)
        self.ylim = self.kwargs.get('ylim', None)
        self.aspect = self.kwargs.get('aspect', 'equal')

        #plot labels
        self.title  = self.name
        self.fontsize = self.kwargs.get('fontsize', 24) ########################Consider making this a title_fontsize ##or consider a general method
        
    #updating the list of colorplot parameters
    def update_colormeshlist(self):
        self.cplot_params = {} #clears the list in case something got deleted

        #Compares the attributes of self to colorplot_params
        pcolormesh_sig = inspect.signature(plt.pcolormesh)
        for attr in dir(self):
            if not callable(getattr(self, attr)) and not attr.startswith("__") and attr in pcolormesh_sig.parameters and not attr == 'kwargs':
                if attr in pcolormesh_sig.parameters:
                    self.cplot_params[attr] = getattr(self, attr)


    def __str__(self):
        # A overinformative string representation
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"fld_quantity with function {self.func.__name__} and kwargs: {kwargs_str}"
    

################################################################################
class figure_plotting:
    def __init__(self, **kwargs):
        self.fig = None
        self.ax = None
        self.kwargs = kwargs
        
        #dataset implementation
        self.datasets={} # a dict of DataSph or Data objects
        self.data = None #Current dataset
        self.num_datasets = 0

        #fld_quantity implementation
        self.fld_quantities = {} # a dict of fld_quantity class objects
        self.fld_val = None # THe current fld_quantity
        assert isinstance(self.fld_val, fld_quantity) or self.fld_val == None, "fld_val must be a fld_quantity object"
        self.num_flds = 0

        self.draw_on_axes = {} #A dict of functions that draw on the axes needs func(ax,...)
        self.fig_funcs = {} #A dict of functions that adjust the figure needs func(fig,...)
        self.step = 0
        assert isinstance(self.step, int), "Step must be an integer"

        self._init_complete = True #This is to allow set_attr to run normally for initialization
    
    #####################################################
    # changing the default setattr and delattr to do other things
    def __setattr__(self, name: str, value: Any) -> None:
        # only runs once the object is fully initialized
        if '_init_complete' in self.__dict__:
            #updates the lists of datasets and fld_quantities
            self.update_field_quantities(name, value)
            self.update_datasets(name, value)
            
            #updates the lists of functions that act on axes or the figure
            self.update_draw_on_axes(name, value)
            self.update_fig_funcs(name, value)
        #This is necessary to avoid infinite recursion
        super().__setattr__(name, value)
    
    def __delattr__(self, name: str) -> None:
        #updates fld_quantity if the thing added is a field val
        value = None #just to call the function

        self.update_field_quantities(name, value,delete=True)
        self.update_datasets(name, value,delete=True)
        #updates the lists of functions that act on axes or the figure
        self.update_draw_on_axes(name, value,delete=True)
        self.update_fig_funcs(name, value,delete=True)

        #This is necessary to avoid infinite recursion
        super().__delattr__(name)
    #####################################################

    def add_plot(self,fld_val, **kwargs):
        #enforces that fld_val is a fld_quantity object
        assert isinstance(fld_val, fld_quantity), "fld_val must be a fld_quantity object"
        
        #Currently I have it such that the name is inherited from the fld_quantity object
        #This is to allow for easy identification of the plot
        name = kwargs.get("name",fld_val.name)
        if name == None:
            error = "fld_quantity object must have a name attribute or a name must be passed as a keyword argument"
            raise ValueError(error)
        
        setattr(self,name,fld_val)
        #By setting the attribute it will call __setattr_ which does all pertinent updates
    def del_plot(self, name):
        #first asserts that the plot exists
        if name not in self.fld_quantities:
            raise ValueError(f"{name} does not exist as plot")
        else:
            delattr(self,name)
            #__Delattr__ will handle the rest
    
    def add_dataset(self, data, name="data", **kwargs):
        #enforces that data is a DataSph object or a Data
        assert isinstance(data, (DataSph, Data)), "data must be a DataSph or Data object"

        setattr(self,name,data)
        #By setting the attribute it will call __setattr_ which does all pertinent updates
    def del_dataset(self, name):
        #first asserts that the dataset exists
        if name not in self.datasets:
            raise ValueError(f"{name} does not exist as dataset")
        else:
            delattr(self,name)
            #__Delattr__ will handle the rest

    #def add_draw_on_plot(self,func,**kwargs):

    def update_datasets(self,name,value,delete=False):
        if isinstance(value, (DataSph, Data)) and name != 'data':
            self.data = value
            self.datasets[name] = value
            self.num_datasets = len(self.datasets)
            if debug.enabled:
                print(f'Added {name} to datasets')
        #deletes if delete=True where Delet is from the kwargs
        if name in self.datasets and delete:
            del self.datasets[name]
            self.num_datasets = len(self.datasets)
            if debug.enabled:
                print(f'Removed {name} from datasets')
        
        

    def update_field_quantities(self,name,value,delete=False):
        if isinstance(value, fld_quantity) and name != 'fld_val':
            self.fld_val = value
            self.fld_quantities[name] = value
            self.num_flds = len(self.fld_quantities)
            if debug.enabled:
                print(f'Added {name} to fld_quantities')
        #deletes if delete=True where Delet is from the kwargs
        if name in self.fld_quantities and delete:
            del self.fld_quantities[name]
            self.num_flds = len(self.fld_quantities)
            if debug.enabled:
                print(f'Removed {name} from fld_quantities')
    
        
            
    def update_draw_on_axes(self,name, value,delete=False):
        sig = inspect.signature(value) if callable(value) else None
        #chech if its a valid function for drawing on axes
        if callable(value) and 'ax' in sig.parameters and name not in self.draw_on_axes:
            self.draw_on_axes[name] = value
            if debug.enabled:
                print(f'Added {name} to draw_on_axes')
        #deletes if delete=True where Delet is from the kwargs
        if name in self.draw_on_axes and delete:
            del self.draw_on_axes[name]
            if debug.enabled:
                print(f'Removed {name} from draw_on_axes')

    def update_fig_funcs(self,name,value,delete=False):
        sig = inspect.signature(value) if callable(value) else None
        #chech if its a valid function for adjusting the figure
        if callable(value) and 'fig' in sig.parameters and name not in self.fig_funcs:
            self.fig_funcs[name] = value
            if debug.enabled:
                print(f'Added {name} to fig_funcs')  
        #deletes if delete=True where Delet is from the kwargs
        if name in self.fig_funcs and delete:
            del self.fig_funcs[name]
            if debug.enabled:
                print(f'Removed {name} from fig_funcs')


    def __str__(self):
        kwargs_str = ', '.join(f'{key}={value}' for key, value in self.kwargs.items())
        return f"figure_plotting with kwargs: {kwargs_str}"


#########################################################
#########################################################
# Here there be the main functions for plotting
'''
The Colorplot function is an axis function
so it requires as input  the axis to plot on
as well as an object of the class figure_plotting
This function plots a pcolormesh plot on the axis
and includes a colorbar nicely placed on the right
'''
def colorplot(ax,figure_plot,**kwargs): 
    fp = figure_plot

    #printing parameters if debugging is enabled
    if debug.enabled: print(f"{fp.fld_val.name} colorplot attributes", fp.fld_val.cplot_params)

    #plotting the colorplot
    c = ax.pcolormesh(fp.data.x1, fp.data.x2, fp.fld_val.func(fp.data), **fp.fld_val.cplot_params)
    
    #include the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05) ###############################consider making this a parameter inside fld_val
    cbar = plt.colorbar(c, cax=cax)
    cbar.ax.tick_params(labelsize=fp.fld_val.fontsize)
    return c
    #I think I should make colorbar be its own function, other types of plotting might need cbar


#########################################################################
'''
set_plot_characteristics is an axis function
so it requires as input  the axis to plot on
as well as an object of the class figure_plotting
This sets the plot characteristics of the axis
'''
def set_plot_characteristics(ax,figure_plot,**kwargs):
    fp = figure_plot
    #There might be a way to do this with accessing signatures
    #however, it could get complicated with different naming conventions
    #So I will just hardcode the possible parameters
    
    
    # This requires that fld_val has the following attributes:
    # xlim, ylim, aspect, title, fontsize
    attrs = ['xlim', 'ylim', 'aspect', 'title', 'fontsize']
    for attr in attrs:
        #just checking if there else make some error 
        if not hasattr(fp.fld_val, attr):
            raise ValueError(f"fld_val must have attribute {attr}")
    
    if debug.enabled: 
        #where fld_val attr match attrs
        for attr in attrs:#makes a dict of the attributes
            fld_attrs = {attr: getattr(fp.fld_val, attr) for attr in attrs}
        print(f"{fp.fld_val.name} plot attributes", fld_attrs)


    
    ax.set_xlim(kwargs.get('xlim', fp.fld_val.xlim))
    ax.set_ylim(kwargs.get('ylim', fp.fld_val.ylim))
    ax.set_aspect(kwargs.get('aspect', fp.fld_val.aspect))
    ax.tick_params(axis='both', which='major', labelsize=kwargs.get('fontsize', fp.fld_val.fontsize))

    ax.set_title(kwargs.get('title', fp.fld_val.title), fontsize=kwargs.get('fontsize', fp.fld_val.fontsize))


#########################################################################

'''
rescale_figure is a basic function that rescales a figsize 
to be closer to a target_size while retaining the aspect ratio
'''
def rescale_figure(figsize,target_size=13):
    width, height = figsize
    aspect = width/height

    if aspect > 1:
        figsize = (target_size, target_size/aspect)
    else:
        figsize = (target_size*aspect, target_size)

    return figsize


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
     
#########################################################################

'''
plot_fig is the main plotting function
should be as high level as possible calling functions where necessatry
'''
#This will plot the figure of a single dataset with various field quantities
# Currently it will allow for either horizontal  vertical display,
# no 2x2 or 3x3 etc. display
def plot_fig(figure_plot, **kwargs):
    fp = figure_plot
    
    ## setting the current data to what is wanted to be plotted
    kwdata = kwargs.get('data', None)
    if kwdata not in fp.datasets:
        if fp.num_datasets > 1:
            raise ValueError("Multiple datasets, must specify which one to plot: data=datasetname")
        elif fp.num_datasets == 0:
            if fp.data == None:
                raise ValueError("No data set to plot")
            else:
                pass
                #this already is data
        elif fp.num_datasets == 1:
            fp.data = list(fp.datasets.values())[0]
        else:
            raise ValueError("Something went wrong, shouldn't of happened")
    else:
        fp.data = fp.datasets[kwdata]

    if debug.enabled: print(f"Number of fld_quantities: {fp.num_flds}")
    if fp.num_flds == 0:
        raise ValueError("There must be at least one field quantity to plot: object.add_plot(fld_quantity_object)")
    ##setting the subplots depending on the get_figure function
    get_fig_func = kwargs.get("get_figure", default_figure)
    fig,ax = get_fig_func(fp, **kwargs)

    #loading the dataset
    fp.data.load_fld(fp.step)

    for i, fld_val in enumerate(list(fp.fld_quantities.values())):
        fp.fld_val = fld_val
        
        #plotting the colorplot
        colorplot(ax[i],fp)
        #setting the plot characteristics
        set_plot_characteristics(ax[i],fp,**kwargs)
        
        # Doing any post processing on the axis
        # Such as drawing field lines, NS, etc
        #Will update to be object instances
        
        #runs all draw_on_axes functions
        for ax_func in fp.draw_on_axes.values():
            ax_func(ax[i],fp.data,**kwargs)
        #runs all fig_funcs functions
        for fig_func in fp.fig_funcs.values():
            fig = fig_func(fig,fp.data,**kwargs)

    fig.tight_layout(pad=1.0)

    return fig


#########################################################################
#########################################################################
# Here there be the axis functions

'''
draw_field_lines is an axis function
so it requires as input  the axis to plot on
as well as an object of the class figure_plotting

Current problems is it required Bp inside the config of the dataset
'''
def draw_field_lines_sph(ax,dataset,**kwargs):
    Bp      = dataset.conf["Bp"]
    flux    = np.cumsum(dataset.B1 * dataset._rv * dataset._rv * np.sin(dataset._thetav) * dataset._dtheta, axis=0)
    clevels = np.linspace(0.0, np.sqrt(Bp), 10)**2
    ax.contour(dataset.x1, dataset.x2, flux, clevels
             , colors='green', linewidths=1)


#########################################################################

'''
draw_NS is an axis function
so it requires as input  the axis to plot on
as well as an object of the class figure_plotting

This function will merely draw the neutron star on the axis
'''
def draw_NS(ax,dataset,**kwargs): #just need dataset for convention
    r = kwargs.get("Radius",1)
    ax.add_patch(plt.Circle((0,0),r,fill=True, color="black", alpha=0.5))


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


#########################################################
def EdotB(name='EdotB'):
    return fld_quantity(
                     lambda data: data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3,
                     name = name,
                     #ob
                     vmin = -1, #default vmin/vmax forthis quantity
                     vmax = 1,
                     )

def JdotB(name='JdotB'):
    return fld_quantity(
                    lambda data: data.J1*data.B1 + data.J2*data.B2 + data.J3*data.B3,
                    name = name,
                    vmin = -1,
                    vmax = 1,
                    )
