# Aperture_plotting

## current issues to resolve
when calling make_fig again in another cell, the cbar ticks are duplicated weirdly

I'm not sure if the default key values for add_plot works, it seems to keep it blank

## Example use case
```python
import plotting
from plotting import apt_fig
import matplotlib.pyplot as plt
data = plotting.DataSph('[path to data folder]/data/')
```
```python
plt.ioff() # Stops duplicate showing
afig = apt_fig(data,"afig")

afig.add_plot("EdotB")
afig.add_plot("B3")
afig.add_post(["draw_field_lines1", "draw_NS"])

afig.add_parameters("all", xlim = [0,10], ylim = [-5,5]
                    , title_fontsize = 24, tick_fontsize=12
                    , vmin = -0.2, vmax = 0.2)
afig.add_parameters("B3", vmin = -1, vmax = 1)

afig.step = 100

afig.make_fig(fontsize=12)
```
![Example output](README_images/mdexample.png)


## Structure of the code

There are four main classes in the code: apt_fig, apt_plot, apt_post, DataSph/Data


## `apt_fig` Class

The `apt_fig` class is the central class responsible for creating figures. It manages the data object and orchestrates the plotting functions of `apt_plot` and `apt_post` objects.

### Attributes

- `data`: The data object that contains the information to be plotted.
- `plots`: A dictionary of `apt_plot` objects that are part of the figure.
- `post_process`: A dictionary of `apt_post` objects for post-processing steps.

### `apt_fig` Methods

#### `add_plot(plot, pos=None)`:
Adds a plot to the figure. Saves the `apt_plot` object in the `plots` dictionary.

- **Parameters:**
  -  `apt_plot`: an `apt_plot` object, a string corresponding to an `apt_plot` object in `apt_plot_types`, or a string corresponding to a field value in the `data` object (e.g., `"B3"` refers to `data.B3`).
  - `pos` (optional): Specifies the grid position of the plot in the figure as a tuple (e.g., `(0, 1)` for first row, second column). The default is `None`, which adds the plot to a new column.
  - `kwargs`: Additional arguments to override the parameters of the `apt_plot` object.
- **Returns:** None.
- **Example Use**
```python
afig.add_plot("B3", pos=(0, 1))
```

#### `add_post(apt_post_obj)`:
Adds a post-processing step to the figure. Saves the `apt_post` object in the `post_process` dictionary.
- **Parameters:**
  - `apt_post_obj`: a list or a single entry of `apt_post` objects or strings corresponding to an `apt_post` object in `apt_post_types`.
  - `kwargs`: additional arguments to override the parameters of the `apt_post` object.

- **Returns:** None.
- **Example Use**
```python
afig.add_post(["draw_field_lines1", "draw_NS"])
```

#### `add_parameters(plots)`:
Adds parameters to the figure. These parameters overwrite the default parameters of the `apt_plot` objects.
- **Parameters:**
  - `plots`: a string corresponding to an `apt_plot` object in `apt_plot_types`, a list of strings corresponding to multiple `apt_plot` objects, or `"all"` to apply to all plots.
  - `kwargs`: Whatever parameters you want to pass to overwrite the default parameters of the `apt_plot` objects.

#### `make_fig()`:
Creates the figure based on the added plots and post-processing steps.
- **Parameters:** 
  - `step` (optional): The step number to plot, overwrites `afig.step`.
  - `kwargs`: Additional arguments to override (not overwrite) the parameters of the `apt_plot` and `apt_post` objects.
- **Returns:** The matplotlib figure.

- **Example Use**
```python
afig.make_fig(fontsize=12)
```

#### `update_fig(step,set_plot_attr)`:
Updates the figure based on a new step. Does not quite work if you change a bunch of attributes, in those cases just call make_fig again.

- **Parameters:** 
  - `step` (optional): The step number to plot.
  - `set_plot_attr` (optional): True/False to update the plot attributes like vmin, vmax, etc. Default is False so it will not have to redraw these attributes
- **Returns:** The updated matplotlib figure.

- **Example Use**
```python
afig.update_fig(60)
```

#### `make_movie(save_name, start, end, increment)`:
Creates a movie of the figure based on a range of steps and saves it in the movies folder.

- **Parameters:** 
  - `save_name` (default: Untitled): The name of the movie file.
  - `start` (optional): The starting step number. Defaults to 0.
  - `end` (optional): The ending step number. Defaults to the last step in the data.
  - `increment` (optional): The step increment. Defaults to 1.
- **Returns:** None.
- **Creates** a .mp4 file in the movies folder.
- **Example Use**
```python
afig.make_movie("Example_movie")
```

### Other `apt_fig` Methods
#### External Use:
- `del_plot(name)`: Deletes a plot from the figure.
- `move_plot(name, pos)`: Moves a plot to a new position in the figure.
- `rescale_figure(target_size)`: Rescales the figure to the target size, maintaining the aspect ratio.

#### Internal Use (probably should not be used directly):
- `override_params()`: returns the parameters of the figure, overriding (not overwriting) the default parameters of the `apt_plot` objects.
- `check_position_taken(pos)`: Checks if a position in the figure grid is already taken.
- `resize_row_col(pos)`: Resizes the num_col and num_row if the input position is larger than the current size.
- `set_fig_shape(num_rows,num_columns)`: creates a new figure with the specified number of rows and columns, copies over all existing axes, then deletes the old figure.
- `construct_plot_obj(key)`: Constructs an `apt_plot` object based on the key value. key needs to be a key of the data object
- `set_fontsize()`: Sets the fontsizes of everything based on kwargs and the parameters of the individual objects. 

