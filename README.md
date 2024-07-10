# Aperture_plotting

## Example Use

from plotting_files.plotting import *

data_rthin_high = DataSph('../data/rthin_high/')

magnetar = figure_plotting()
magnetar.add_dataset(data_rwide, name='data_rwide')
magnetar.add_plot(EdotB("test"))
_=plot_fig(magnetar)
