# plotting2.py - Simplified Aperture Plotting Framework
# A streamlined version with unified plotting registration

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import display
import inspect
import os
import functools
from typing import Any, Callable, Dict, List, Optional, Union

# Data imports
import sys
sys.path.append('.')
from datalib_logsph import DataSph, flag_to_species
from datalib import Data

# Matplotlib setup
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
matplotlib.rc("font", family="serif")
matplotlib.rc("text", usetex=True)

# Load colormap
import pickle
try:
    with open("hot_cold_cmap.pkl", "rb") as f:
        cdata = pickle.load(f)
    try:
        matplotlib.cm.get_cmap("hot_and_cold")
    except:
        matplotlib.colormaps.register(cdata, name="hot_and_cold")
except FileNotFoundError:
    pass  # Colormap file not found, continue without it

################################################################################
# Registry for plotting functions
################################################################################

_plotters: Dict[str, Callable] = {}
_post_processors: Dict[str, Callable] = {}


def register(func: Callable, name: str = None) -> Callable:
    """
    Register a plotting function.
    
    The function should have signature: func(data, ax, **kwargs) -> matplotlib_object
    
    Usage:
        @apl.register
        def my_plotter(data, ax, **kwargs):
            return ax.plot(...)
        
        # or
        apl.register(my_plotter)
        apl.register(my_plotter, name="custom_name")
    """
    if name is None:
        name = func.__name__
    _plotters[name] = func
    return func


def register_post(func: Callable, name: str = None) -> Callable:
    """
    Register a post-processing function.
    
    The function should have signature: func(data, ax, **kwargs) -> matplotlib_object
    Post-processors run after the main plot is made.
    """
    if name is None:
        name = func.__name__
    _post_processors[name] = func
    return func


def list_plotters() -> List[str]:
    """List all registered plotter names."""
    return list(_plotters.keys())


def list_post_processors() -> List[str]:
    """List all registered post-processor names."""
    return list(_post_processors.keys())


################################################################################
# apt_plot: Individual subplot
################################################################################

class apt_plot:
    """
    Represents a single subplot in an apt_fig.
    
    Attributes:
        name: Identifier for this plot
        data: Data object for this plot
        ax: Matplotlib axes object
        position: Grid position (row, col)
        step: Current data step
        parameters: Plot parameters
    """
    
    def __init__(self, name: str, data=None, pos: tuple = None, **kwargs):
        self.name = name
        self.data = data
        self.ax = None
        self.position = pos
        self.step = kwargs.get('step', 0)
        self.colspan = kwargs.get('colspan', 1)
        self.rowspan = kwargs.get('rowspan', 1)
        
        # Plot state
        self.parameters = {'fontsize': 24}
        self.parameters.update(kwargs)
        
        self._plotters: List[tuple] = []  # List of (plotter_name, kwargs)
        self._post_processors: List[tuple] = []  # List of (post_name, kwargs)
        self._plot_objects = {}  # Store matplotlib objects by plotter name
        self.made = False
    
    def set_step(self, step: int = None):
        """Set the data step (data is loaded lazily in make_plot)."""
        if step is not None:
            self.step = step
    
    def _load_data(self):
        """Load data for current step, skipping if already loaded."""
        if self.data is None:
            return
        if isinstance(self.data, list):
            for d in self.data:
                if getattr(d, '_current_fld_step', None) != self.step:
                    d.load(self.step)
        else:
            if getattr(self.data, '_current_fld_step', None) != self.step:
                self.data.load(self.step)
    
    def override_params(self, **kwargs) -> dict:
        """Return parameters with kwargs overrides."""
        params = self.parameters.copy()
        params.update(kwargs)
        return params
    
    def add_plotter(self, plotter_name: str, **kwargs):
        """Queue a plotter to run on this subplot."""
        if plotter_name not in _plotters:
            raise ValueError(f"Plotter '{plotter_name}' not registered. Available: {list_plotters()}")
        # Keep data in kwargs so each plotter can use its own data if specified
        self._plotters.append((plotter_name, kwargs))
    
    def add_post(self, post_name: str, **kwargs):
        """Queue a post-processor for this subplot."""
        if post_name not in _post_processors:
            raise ValueError(f"Post-processor '{post_name}' not registered. Available: {list_post_processors()}")
        # Keep data in kwargs so each post-processor can use its own data if specified
        self._post_processors.append((post_name, kwargs))
    
    def make_plot(self, **kwargs):
        """Execute all queued plotters."""
        self.set_step(kwargs.get('step'))
        self._load_data()  # Load data just before plotting (ensures plots can have different steps)
        params = self.override_params(**kwargs)
        
        for plotter_name, plotter_kwargs in self._plotters:
            merged = {**params, **plotter_kwargs}
            # Determine which data to use for this specific plotter
            data_to_use = merged.pop('data', self.data)  # Pop from merged, default to self.data
            # Load the specific data object if needed
            if data_to_use is not None:
                if isinstance(data_to_use, list):
                    for d in data_to_use:
                        if getattr(d, '_current_fld_step', None) != self.step:
                            d.load(self.step)
                else:
                    if getattr(data_to_use, '_current_fld_step', None) != self.step:
                        data_to_use.load(self.step)
            func = _plotters[plotter_name]
            result = func(data_to_use, self.ax, **merged)
            # Store results as lists to support multiple calls with the same plotter name
            if plotter_name not in self._plot_objects:
                self._plot_objects[plotter_name] = []
            self._plot_objects[plotter_name].append(result)
        
        self.set_plot_attr(**params)
        self.made = True
    
    def make_post(self, **kwargs):
        """Execute all queued post-processors."""
        params = self.override_params(**kwargs)
        
        for post_name, post_kwargs in self._post_processors:
            merged = {**params, **post_kwargs}
            # Determine which data to use for this specific post-processor
            data_to_use = merged.pop('data', self.data)  # Pop from merged, default to self.data
            # Load the specific data object if needed
            if data_to_use is not None:
                if isinstance(data_to_use, list):
                    for d in data_to_use:
                        if getattr(d, '_current_fld_step', None) != self.step:
                            d.load(self.step)
                else:
                    if getattr(data_to_use, '_current_fld_step', None) != self.step:
                        data_to_use.load(self.step)
            func = _post_processors[post_name]
            result = func(data_to_use, self.ax, **merged)
            # Store results as lists to support multiple calls with the same post-processor name
            if post_name not in self._plot_objects:
                self._plot_objects[post_name] = []
            self._plot_objects[post_name].append(result)

    def update_plot(self, **kwargs):
        """Update the plot with new data (selective removal, fallback to cla if needed)."""
        self.set_step(kwargs.get('step'))
        
        # Try selective removal first (preserves axis labels/ticks if possible)
        try:
            for artist in list(self.ax.get_children()):
                # Keep axis decorations
                if isinstance(artist, (matplotlib.spines.Spine, matplotlib.axis.Axis)):
                    continue
                if artist is self.ax.title:
                    continue
                if isinstance(artist, matplotlib.legend.Legend):
                    continue
                # Try to remove
                remove = getattr(artist, "remove", None)
                if callable(remove):
                    remove()
        except NotImplementedError as e:
            # If any artist can't be removed, fall back to full cla()
            # print(f"Warning: Could not remove artist {type(artist).__name__}. Falling back to ax.cla()")
            self.ax.cla()
        
        self._plot_objects.clear()
        
        self.make_plot(**kwargs)
        self.make_post(**kwargs)
    
    def set_plot_attr(self, **kwargs):
        """Set axis attributes from parameters."""
        params = self.parameters.copy()
        params.update(kwargs)
        
        attrs = ['xlim', 'ylim', 'aspect', 'title', 'xlabel', 'ylabel', 'xscale', 'yscale']
        for attr in attrs:
            if attr in params and params[attr] is not None:
                try:
                    getattr(self.ax, f'set_{attr}')(params[attr])
                except Exception as e:
                    print(f"Warning: Could not set {attr} for plot '{self.name}': {e}")
        
        if params.get('legend', False):
            self.ax.legend(loc=params.get('legend_loc', 'best'))
        if params.get('hide_yaxis', False):
            self.ax.yaxis.set_visible(False)
        if params.get('hide_xaxis', False):
            self.ax.xaxis.set_visible(False)
    
    def set_fontsize(self, **kwargs):
        """Set fontsizes for the plot."""
        params = self.override_params(**kwargs)
        fontsize = params.get('fontsize')
        
        label_fontsize = params.get('label_fontsize', fontsize)
        title_fontsize = params.get('title_fontsize', fontsize)
        tick_fontsize = params.get('tick_fontsize', fontsize)
        legend_fontsize = params.get('legend_fontsize', fontsize)
        
        if title_fontsize:
            self.ax.title.set_fontsize(title_fontsize)
        if label_fontsize:
            self.ax.xaxis.label.set_fontsize(label_fontsize)
            self.ax.yaxis.label.set_fontsize(label_fontsize)
        if tick_fontsize:
            self.ax.tick_params(labelsize=tick_fontsize)
        
        legend = self.ax.get_legend()
        if legend and legend_fontsize:
            for text in legend.get_texts():
                text.set_fontsize(legend_fontsize)
    
    def copy_ax_attr(self, new_ax):
        """Copy axis attributes to a new axis."""
        if self.ax is None:
            self.ax = new_ax
            return
        
        attrs = ['title', 'xlabel', 'ylabel', 'xlim', 'ylim', 'xscale', 'yscale','rasterized']
        for attr in attrs:
            try:
                val = getattr(self.ax, f'get_{attr}')()
                getattr(new_ax, f'set_{attr}')(val)
            except:
                pass
        self.ax = new_ax


################################################################################
# apt_fig: Figure container
################################################################################

_figure_objects: Dict[str, 'apt_fig'] = {}


class apt_fig:
    """
    Main figure class for Aperture plotting.
    
    Usage:
        afig = apt_fig(data, "my_figure")
        afig.add_plot("plot1", pos=(0, 0))
        afig.add_plot("plot2", pos=(0, 1))
        afig.colorplot("plot1", fld_val="rho", cmap="hot")
        afig.streamlines(["plot1", "plot2"], density=2)
        afig.make_fig()
    """
    
    def __init__(self, data, unique_ident: str = "default", **kwargs):
        global _figure_objects
        plt.ioff()  # Start with interactive mode off to prevent auto-display during setup
        # Handle existing figure with same identifier
        self.unique_ident = unique_ident
        if unique_ident in _figure_objects:
            old_fig = _figure_objects[unique_ident]
            if old_fig.fig is not None:
                plt.close(old_fig.fig)
            del _figure_objects[unique_ident]
        # register this figure
        _figure_objects[unique_ident] = self
        
        self.data = data
        self.step = kwargs.get('step', 0)
        self.fig = None
        self.plots: Dict[str, apt_plot] = {}
        self.parameters = {'fontsize': 24}
        self.parameters.update(kwargs)
        
        self.rows = 1
        self.columns = 1
        self.made = False
        self.plot_folder = kwargs.get('plot_folder', os.getcwd())
    
    def __getattr__(self, name: str):
        """
        Allow calling registered plotters as methods.
        e.g., afig.colorplot("plot1", ...) calls _plotters["colorplot"]
        """
        if name.startswith('_') or name in ('plots', 'parameters', 'data'):
            raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")
        
        if name in _plotters:
            original_func = _plotters[name]
            @functools.wraps(original_func)
            def plotter_method(plots: Union[str, List[str]], **kwargs):
                return self._apply_plotter(name, plots, **kwargs)
            return plotter_method
        
        if name in _post_processors:
            original_func = _post_processors[name]
            @functools.wraps(original_func)
            def post_method(plots: Union[str, List[str]], **kwargs):
                return self._apply_post(name, plots, **kwargs)
            return post_method
        
        raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}' and it's not a registered plotter")
    
    def _apply_plotter(self, plotter_name: str, plots: Union[str, List[str]], **kwargs):
        """
        Apply a plotter to specified plots. Auto-creates plots if they don't exist.
        
        When auto-creating, all kwargs are passed to add_plot (stored in parameters).
        The plotter also receives all kwargs - it uses what it needs.
        """
        if plots == "all":
            plots = list(self.plots.keys())
        elif isinstance(plots, str):
            plots = [plots]
        
        # Add afig object to kwargs so plotters can access apt_fig properties
        kwargs['afig'] = self
        
        for plot_name in plots:
            if plot_name not in self.plots:
                # Auto-create the plot with all kwargs as parameters
                self.add_plot(plot_name, **kwargs)
            
            self.plots[plot_name].add_plotter(plotter_name, **kwargs)
    
    def _apply_post(self, post_name: str, plots: Union[str, List[str]], **kwargs):
        """Apply a post-processor to specified plots."""
        if plots == "all":
            plots = list(self.plots.keys())
        elif isinstance(plots, str):
            plots = [plots]
        
        for plot_name in plots:
            if plot_name not in self.plots:
                raise ValueError(f"Plot '{plot_name}' not found")
            self.plots[plot_name].add_post(post_name, **kwargs)
    
    def override_params(self, **kwargs) -> dict:
        """Return parameters with kwargs overrides."""
        params = self.parameters.copy()
        params.update(kwargs)
        return params
    
    def set_step(self, step: int, plots: Union[str, List[str]] = "all"):
        """Set step for specified plots (data loaded lazily during make_plot)."""
        self.step = step
        
        if plots == "all":
            plots = list(self.plots.keys())
        elif isinstance(plots, str):
            plots = [plots]
        
        for plot_name in plots:
            self.plots[plot_name].set_step(step)
    
    def _check_position_taken(self, pos: tuple) -> bool:
        """Check if grid position is occupied."""
        for plot in self.plots.values():
            if plot.position == pos:
                return True
        return False
    
    def _resize_grid(self, pos: tuple):
        """Expand grid if needed for position."""
        self.rows = max(self.rows, pos[0] + 1)
        self.columns = max(self.columns, pos[1] + 1)
    
    def _set_fig_grid(self, num_rows: int, num_columns: int):
        """Recreate figure with new grid dimensions."""
        new_fig = plt.figure()
        
        for plot in self.plots.values():
            new_ax = plt.subplot2grid(
                (num_rows, num_columns),
                plot.position,
                fig=new_fig,
                colspan=plot.colspan,
                rowspan=plot.rowspan
            )
            plot.copy_ax_attr(new_ax)
        
        if self.fig is not None:
            plt.close(self.fig)
        self.fig = new_fig
    
    def add_plot(self, name: str, pos: tuple = None, data=None, **kwargs):
        """
        Add a subplot to the figure.
        
        Args:
            name: Identifier for this subplot
            pos: Grid position (row, col), auto-assigned if None
            data: Data object for this plot (defaults to figure data)
            **kwargs: Parameters for the subplot
        """
        if name in self.plots:
            raise ValueError(f"Plot '{name}' already exists")
        
        if data is None:
            data = self.data
        
        step = kwargs.pop('step', self.step)
        
        # Auto-assign position if not specified
        if pos is None:
            if len(self.plots) == 0:
                pos = (0, 0) #defaults to single plot
            else:
                found = False
                # Find first available position
                # Scan row-wise
                for col in range(self.columns):
                    for row in range(self.rows):
                        if not self._check_position_taken((row, col)):
                            pos = (row, col)
                            found = True
                            break
                    if found:
                        break
                if not found:
                    self.columns += 1 # default to adding a new column
                    pos = (0, self.columns - 1)
        else:
            if self._check_position_taken(pos):
                raise ValueError(f"Position {pos} is already taken")
        
        self._resize_grid(pos)
        self._set_fig_grid(self.rows, self.columns)
        
        # Create the plot object
        plot = apt_plot(name, data=data, pos=pos, step=step, **kwargs)
        self.plots[name] = plot
        plot.position = pos
        
        # Create matplotlib axis
        plot.ax = plt.subplot2grid(
            (self.rows, self.columns),
            pos,
            fig=self.fig,
            colspan=plot.colspan,
            rowspan=plot.rowspan
        )
        
    
    def del_plot(self, name: str):
        """Remove a subplot from the figure."""
        if name not in self.plots:
            raise ValueError(f"Plot '{name}' does not exist")
        
        del self.plots[name]
        
        # Recalculate grid size
        self.rows = 1
        self.columns = 1
        for plot in self.plots.values():
            self._resize_grid(plot.position)
        
        self._set_fig_grid(self.rows, self.columns)
    
    def make_fig(self, **kwargs):
        """
        Create the figure by executing all queued plotters and post-processors.
        """
        if self.made:
            # Reset for re-making (don't close fig yet - _set_fig_grid needs old axes)
            for plot in self.plots.values():
                plot.made = False
                plot._plot_objects.clear()
            self._set_fig_grid(self.rows, self.columns)  # closes old fig after copying attrs
        
        params = self.override_params(**kwargs)
        
        if 'step' in kwargs: # Note: will override previous set_step calls
            self.set_step(kwargs['step'])
        
        # Make all plots
        for plot in self.plots.values():
            plot.make_plot(**params)
        
        # Run post-processors
        for plot in self.plots.values():
            plot.make_post(**params)
        
        # Set fontsizes
        self.set_fontsize(**params)
        self.set_size()
        
        # Adjust layout
        layout_keys = ['wspace', 'hspace', 'top', 'bottom', 'left', 'right']
        if any(k in params for k in layout_keys):
            self.fig.subplots_adjust(**{k: params.get(k) for k in layout_keys if k in params})
        else:
            self.fig.tight_layout()
        
        self.made = True # mark as made to prevent re-making unless requested
        
        if kwargs.get('save_fig', False):
            save_name = kwargs.get('save_name', 'Untitled')
            self.fig.savefig(f"{self.plot_folder}/{save_name}.png")
        
        return self.fig
    
    def update_fig(self, step: int = None, **kwargs):
        """Update the figure with new data."""
        if step is not None:
            self.set_step(step)# updates ALL plots to this EXACT step
        
        params = self.override_params(**kwargs)
        
        for plot in self.plots.values():
            plot.update_plot(**params)
        
        self.fig.canvas.draw_idle()
        return self.fig
    
    def set_fontsize(self, **kwargs):
        """Set fontsizes for all plots."""
        params = self.override_params(**kwargs)
        
        for plot in self.plots.values():
            plot.set_fontsize(**params)
    
    def set_size(self, xsize: float = None, ysize: float = None, **kwargs):
        """Set figure size in inches."""
        if xsize is not None and ysize is not None:
            self.shape = [xsize, ysize]
        
        if hasattr(self, 'shape'):
            self.fig.set_size_inches(self.shape, **kwargs)
    
    def add_parameters(self, plots: Union[str, List[str]], **kwargs):
        """Add/override parameters for specified plots."""
        if plots == "all":
            plots = list(self.plots.keys())
        elif isinstance(plots, str):
            plots = [plots]
        
        for plot_name in plots:
            self.plots[plot_name].parameters.update(kwargs)
    
    def make_movie(self, save_name: str = "Untitled", start: int = 0, end: int = None,
                   increment: int = 1, **kwargs):
        """Create a movie from the figure."""
        # Strip .mp4 extension if already present to avoid .mp4.mp4
        if save_name.endswith('.mp4'):
            save_name = save_name[:-4]
        
        fps = kwargs.get('fps', 25)
        pic_folder = kwargs.get('pic_folder', 'movie_plots')
        save_folder = kwargs.get('save_folder', f'{self.plot_folder}/movies')
        
        if end is None:
            end = len(self.data.fld_steps)
        
        if not self.made:
            self.make_fig(**kwargs)
        
        os.makedirs(pic_folder, exist_ok=True)
        os.makedirs(save_folder, exist_ok=True)
        
        # Clear old frames
        for f in os.listdir(pic_folder):
            if f.endswith('.png'): #just to not delete things that may be useful
                os.remove(os.path.join(pic_folder, f))
        
        # Generate frames
        for i, step in enumerate(range(start, end, increment)):
            self.update_fig(step=step, **kwargs)
            self.fig.savefig(f"{pic_folder}/{i:05d}.png")
            print(f"\rMovie progress: {i+1}/{(end-start)//increment}", end="")
        
        print()
        
        # Create movie with ffmpeg
        cmd = f"ffmpeg -y -loglevel error -r 10 -i {pic_folder}/%05d.png -c:v libx264 -vf fps={fps} -pix_fmt yuv420p -threads 0 {save_folder}/{save_name}.mp4"
        os.system(cmd)
        print(f"Movie saved to {save_folder}/{save_name}.mp4")
    
    def print_info(self):
        """Print figure information."""
        print(f"apt_fig: {self.unique_ident}")
        print(f"  Grid: {self.rows}x{self.columns}")
        print(f"  Step: {self.step}")
        print(f"  Plots:")
        for name, plot in self.plots.items():
            print(f"    - {name} at {plot.position}")
            print(f"      Plotters: {[p[0] for p in plot._plotters]}")
            print(f"      Post-processors: {[p[0] for p in plot._post_processors]}")


################################################################################

@register
def colorplot(data, ax, fld_val=None, cmap="hot_and_cold", include_colorbar=True, **kwargs):
    """
    Standard 2D colorplot using pcolormesh.
    
    Args:
        fld_val: Field to plot - can be string key, lambda, or callable
        cmap: Colormap name
        include_colorbar: Whether to add a colorbar
        vmin, vmax: Color limits
        logscale: Use log color scale
        aspect: Axis aspect ratio
    """
    from warnings import filterwarnings
    filterwarnings("ignore", category=UserWarning, message="The input coordinates to pcolormesh.*")
    
    # Resolve fld_val
    if callable(fld_val):
        values = fld_val(data)
    elif isinstance(fld_val, str) and hasattr(data, fld_val):
        values = getattr(data, fld_val)
    elif fld_val is None:
        raise ValueError("fld_val must be specified for colorplot")
    else:
        raise ValueError(f"Cannot resolve fld_val: {fld_val}")
    
    # Get grid
    x_grid = kwargs.get('x_grid', data.x1)
    y_grid = kwargs.get('y_grid', data.x2)
    ax.set_aspect('equal')
    
    # Handle log / symlog scaling
    logscale = kwargs.get('logscale', False)
    symlogscale = kwargs.get('symlogscale', False)
    vmin = kwargs.get('vmin')
    vmax = kwargs.get('vmax')
    
    norm = None
    if symlogscale:
        from matplotlib.colors import SymLogNorm
        linthresh = kwargs.get('linthresh', 1e-3)
        linscale = kwargs.get('linscale', 1)
        norm = SymLogNorm(linthresh=linthresh, linscale=linscale, vmin=vmin, vmax=vmax, base=10)
        pcolor_kwargs = {'cmap': cmap, 'norm': norm, 'rasterized': True}
    elif logscale:
        if vmin is None:
            vmin = values[values > 0].min() if np.any(values > 0) else 1e-10
        if vmax is None:
            vmax = values.max()
        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        pcolor_kwargs = {'cmap': cmap, 'norm': norm, 'rasterized': True}
    else:
        pcolor_kwargs = {'cmap': cmap, 'vmin': vmin, 'vmax': vmax, 'rasterized': True}
    
    c = ax.pcolormesh(x_grid, y_grid, values, **pcolor_kwargs)
    
    # if kwargs.get('aspect'):
    #     ax.set_aspect(kwargs['aspect'])
    
    # Skip colorbar creation if one already exists on this axis (for movies)
    if include_colorbar and not getattr(ax, '_colorbar_exists', False):
        cbar_pos = kwargs.get('cbar_pos', 'right')
        cbar_size = kwargs.get('cbar_size', '5%')
        cbar_pad = kwargs.get('cbar_pad', 0.05)
        cbar_shrink = kwargs.get('cbar_shrink', 1.0)
        cbar_label = kwargs.get('cbar_label')
        cbar_label_size = kwargs.get('cbar_label_size', 12)
        cbar_no_edge_label = kwargs.get('cbar_no_edge_label', False)
        ticks = kwargs.get('ticks')
        linthresh = kwargs.get('linthresh', 1e-3)
        linscale = kwargs.get('linscale', 1)
        
        divider = make_axes_locatable(ax)
        orientation = "horizontal" if cbar_pos in ["top", "bottom"] else "vertical"
        cax = divider.append_axes(cbar_pos, size=cbar_size, pad=cbar_pad)
        cbar = plt.colorbar(c, cax=cax, orientation=orientation, shrink=cbar_shrink)
        
        if cbar_label:
            if orientation == "vertical":
                cbar.ax.set_title(cbar_label, pad=10, fontsize=cbar_label_size)
            else:
                cbar.set_label(cbar_label, rotation=0, labelpad=10, fontsize=cbar_label_size)
        
        if logscale:
            from matplotlib.ticker import LogFormatter, LogLocator
            
            # Custom formatter for 10^n format
            class CustomLogFormatter(LogFormatter):
                def __call__(self, value, tickdir=None):
                    return f"$10^{{{int(np.log10(value))}}}$"
            
            formatter = CustomLogFormatter(labelOnlyBase=False)
            if ticks is None and vmin is not None and vmax is not None and vmin > 0:
                ticks = [10 ** i for i in range(int(np.floor(np.log10(vmin))), int(np.ceil(np.log10(vmax))) + 1)]
            if ticks is not None:
                cbar.set_ticks(ticks)
            axis = cbar.ax.xaxis if orientation == "horizontal" else cbar.ax.yaxis
            axis.set_major_formatter(formatter)
            axis.set_major_locator(LogLocator(base=10))
            if orientation == "horizontal":
                if cbar_pos == "top":
                    axis.set_ticks_position("top")
                    axis.set_label_position("top")
                else:
                    axis.set_ticks_position("bottom")
                    axis.set_label_position("bottom")
            else:
                if cbar_pos == "left":
                    axis.set_ticks_position("left")
                    axis.set_label_position("left")
                else:
                    axis.set_ticks_position("right")
                    axis.set_label_position("right")
        elif symlogscale:
            from matplotlib.ticker import SymmetricalLogLocator, FormatStrFormatter
            locator = SymmetricalLogLocator(linthresh=linthresh, base=10)
            formatter = FormatStrFormatter("%.1g")
            axis = cbar.ax.xaxis if orientation == "horizontal" else cbar.ax.yaxis
            axis.set_major_locator(locator)
            axis.set_major_formatter(formatter)
            if orientation == "horizontal":
                if cbar_pos == "top":
                    axis.set_ticks_position("top")
                    axis.set_label_position("top")
                else:
                    axis.set_ticks_position("bottom")
                    axis.set_label_position("bottom")
            else:
                if cbar_pos == "left":
                    axis.set_ticks_position("left")
                    axis.set_label_position("left")
                else:
                    axis.set_ticks_position("right")
                    axis.set_label_position("right")
        else:
            from matplotlib.ticker import ScalarFormatter
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-np.inf, np.inf))  # Force scientific notation always
            axis = cbar.ax.xaxis if orientation == "horizontal" else cbar.ax.yaxis
            axis.set_major_formatter(formatter)
            if orientation == "horizontal":
                if cbar_pos == "top":
                    axis.set_ticks_position("top")
                    axis.set_label_position("top")
                else:
                    axis.set_ticks_position("bottom")
                    axis.set_label_position("bottom")
            else:
                if cbar_pos == "left":
                    axis.set_ticks_position("left")
                    axis.set_label_position("left")
                else:
                    axis.set_ticks_position("right")
                    axis.set_label_position("right")
        
        if cbar_no_edge_label:
            ticklabels = cbar.ax.get_xticklabels() if orientation == "horizontal" else cbar.ax.get_yticklabels()
            if ticklabels:
                ticklabels[-1].set_visible(False)
        
        # Mark that colorbar has been created on this axis
        ax._colorbar_exists = True
    
    return c


@register
def lineplot(data, ax, x_val=None, y_val=None, label=None, **kwargs):
    """
    Simple line plot.
    
    Args:
        x_val: X data - string key, lambda, or array
        y_val: Y data - string key, lambda, or array
        label: Line label for legend
    """
    def resolve_val(val):
        if callable(val):
            return val(data)
        elif isinstance(val, str) and hasattr(data, val):
            return getattr(data, val)
        elif isinstance(val, np.ndarray):
            return val
        else:
            raise ValueError(f"Cannot resolve value: {val}")
    
    x = resolve_val(x_val) if x_val else np.arange(len(resolve_val(y_val)))
    y = resolve_val(y_val)
    
    plot_kwargs = {k: v for k, v in kwargs.items() if k in ['color', 'linestyle', 'linewidth', 'marker', 'markersize', 'alpha']}
    
    line = ax.plot(x, y, label=label, **plot_kwargs)
    return line


@register_post
def draw_time(data, ax, fmt=r"$t = {:.2f}$", pos=(0.05, 0.95), **kwargs):
    """
    Draw current time on the plot.
    
    Args:
        fmt: Format string for time
        pos: Position (x, y) in axes coordinates
    """
    fontsize = kwargs.get('fontsize', 12)
    color = kwargs.get('color', 'black')
    
    text = ax.text(
        pos[0], pos[1],
        fmt.format(data.time),
        transform=ax.transAxes,
        fontsize=fontsize,
        color=color,
        verticalalignment='top'
    )
    return text

def convert_to_index(ax_index,value,data):
        # takes an x1,x2 (x3) value and converts it into an index on the grid
        lower = data.conf["lower"][ax_index]
        N = data.conf["N"][ax_index]
        downsample = data.conf["downsample"]
        size = data.conf["size"][ax_index]

        index = (value-lower)*N/(size*downsample)
        return int(index) # as its an index
    
@register_post
def draw_field_line(data, ax, th_foot=None, r_max=None, color='green', linewidth=1, plot_left=False, **kwargs):
    """
    Draw magnetic field lines as contours of flux at footpoint locations.
    
    Args:
        th_foot: Theta footpoint angle(s) (list or scalar). Default: π/4
        r_max: Radius value(s) to compute footpoints from. Overrides th_foot if provided.
        color: Line color
        linewidth: Line width
        plot_left: If True, plot on left side (negative x)
    """
    if th_foot is None and r_max is None:
        th_foot = np.pi / 4
    
    if not isinstance(th_foot, list) if th_foot is not None else False:
        th_foot = [th_foot] if th_foot is not None else []
    else:
        th_foot = list(th_foot) if th_foot is not None else []
    
    # If r_max provided, compute th_foot from it
    if r_max is not None:
        if not isinstance(r_max, list):
            r_max = [r_max]
        th_foot = []
        for r in r_max:
            th_foot.append(np.arcsin(np.sqrt(1.0 / r)))
    
    # Compute magnetic flux
    flux = np.cumsum(data.B1 * data._rv * data._rv * np.sin(data._thetav) * data._dtheta, axis=0)
    
    # Find the index corresponding to the footpoint (r=1, i.e., surface)
    r_foot_index = convert_to_index(0, np.log(1.0), data)
    
    # Compute footpoint fluxes
    footpoint_fluxes = []
    for th in th_foot:
        theta_foot_index = convert_to_index(1, th, data)
        flux_foot = flux[theta_foot_index, r_foot_index]
        footpoint_fluxes.append(flux_foot)
    
    # Plot contours at footpoint flux levels.
    # Matplotlib requires increasing levels, so we sort levels but keep
    # color mapping consistent with the original input order.
    flux_array = np.asarray(footpoint_fluxes)
    sort_idx = np.argsort(flux_array)
    sorted_levels = flux_array[sort_idx]

    contour_colors = color
    if isinstance(color, (list, tuple, np.ndarray)):
        color_list = list(color)
        if len(color_list) == 1:
            contour_colors = color_list
        elif len(color_list) < len(sorted_levels):
            contour_colors = [color_list[i % len(color_list)] for i in sort_idx]
        else:
            contour_colors = [color_list[i] for i in sort_idx]

    if plot_left:
        contour = ax.contour(-data.x1, data.x2, flux, levels=sorted_levels, colors=contour_colors, linewidths=linewidth)
    else:
        contour = ax.contour(data.x1, data.x2, flux, levels=sorted_levels, colors=contour_colors, linewidths=linewidth)
    
    return contour


################################################################################
# Convenience functions
################################################################################

def close_all():
    """Close all apt_fig figures."""
    global _figure_objects
    for fig in _figure_objects.values():
        if fig.fig is not None:
            plt.close(fig.fig)
    _figure_objects.clear()


def get_figure(name: str) -> apt_fig:
    """Get an apt_fig by name."""
    return _figure_objects.get(name)
