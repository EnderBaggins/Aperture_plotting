# Aperture_plotting

## current issues to resolve
when calling make_fig again in another cell, the cbar ticks are duplicated weirdly

I will be writing the documentation for current implementation asap
what follows is depracated, might still work but its not the best way to use the code

## Functions to use when plotting
self is your apt_fig object

self.add_plot(apt_plot_obj)
self.del_plot(name)
self.move_plot(position)
self.rescale_figure(size)

## Changing parameters
self.plots[name].parameters[attr] = value
note this will not reload the figure when you change the value, I don't yet have a simple function to do so

you can change things directly by
self.plots[name].ax.set_xlim() for example

## Example Use

from plotting_files.plotting import *

data = DataSph('../data/rthin_high/')

magnetar = apt_fig(data)
magnetar.step = 10
magnetar.add_plot(EdotB("test"))
