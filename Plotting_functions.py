
import plotting as apl
from Useful_functions import *
from plotting import register_post,register
import matplotlib.pyplot as plt 
import numpy as np
import os
import matplotlib
import scipy.integrate as integrate
import multiprocessing as mp

@register_post
def draw_NS(data, ax, Radius=1, color='black', transparency=0.5, **kwargs):
    """
    Draw a neutron star circle on the plot.
    
    Args:
        Radius: Radius of the neutron star (default: 1)
        color: Color of the neutron star (default: 'black')
        transparency: Alpha transparency (default: 0.5)
    """
    star = ax.add_patch(plt.Circle((0, 0), Radius, fill=True, color=color, alpha=transparency))
    return star


@register_post
def r_th_region(data, ax, rmin=0, rmax=20, thmin=0, thmax=np.pi/2, 
                outline_color='red', fill_color='purple', alpha=0.5, linewidth=0.4, **kwargs):
    """
    Show a specific r, theta region on top of a previous colorplot.
    
    Args:
        rmin: Minimum radius (default: 0)
        rmax: Maximum radius (default: 20)
        thmin: Minimum theta (default: 0)
        thmax: Maximum theta (default: π/2)
        outline_color: Color of the outline (default: 'red')
        fill_color: Color of the fill (default: 'purple')
        alpha: Transparency of the fill (default: 0.5)
        linewidth: Width of the outline (default: 0.4)
    """
    r_grid = data._rv
    theta_grid = data._thetav
    
    mask = (r_grid >= rmin) & (r_grid <= rmax) & (theta_grid >= thmin) & (theta_grid <= thmax)
    alpha_mask = np.where(mask == 1, alpha, 0)
    fill = ax.pcolormesh(data.x1, data.x2, mask, shading='auto', alpha=alpha_mask, color=fill_color)
    numerical_mask = mask.astype(int)
    outline = ax.contour(data.x1, data.x2, numerical_mask, colors=outline_color, linewidths=linewidth)
    return [fill, outline]
@register
def fieldline_value_vs_rmax(data, ax, value_func=None, rmax_list=None, value_func_kwargs=None, ylabel=r"$y$", title=None, **kwargs):
    """
    Generic plotter for a value computed along field lines vs r_max.

    Parameters:
    - data: Data object (passed by plotting framework)
    - ax: matplotlib axis
    - value_func: function to compute the value (e.g., compute_Jpar, compute_voltage)
    - rmax_list: list/array of r_max values
    - value_func_kwargs: dict of additional arguments for value_func
    - ylabel: y-axis label (LaTeX string)
    - title: plot title (LaTeX string)
    - kwargs: passed to ax.plot (including style args)
    """
    postprocess = kwargs.pop('postprocess', None)
    if value_func is None or rmax_list is None:
        raise ValueError("Must provide value_func and rmax_list")
    if value_func_kwargs is None:
        value_func_kwargs = kwargs
    values = value_func(data, rmax_list, **value_func_kwargs)
    if isinstance(values, tuple):  # If function returns (values, stats)
        values = values[0]
    if postprocess is not None:
        values = postprocess(values)

    # Style argument normalization (like other plotters)
    def to_list(val):
        if val is None:
            return None
        if isinstance(val, (list, tuple, np.ndarray)):
            return list(val)
        return [val]

    colors = to_list(kwargs.pop('colors', kwargs.pop('color', None)))
    labels = to_list(kwargs.pop('labels', kwargs.pop('label', None)))
    linestyles = to_list(kwargs.pop('linestyles', kwargs.pop('linestyle', None)))
    linewidths = to_list(kwargs.pop('linewidths', kwargs.pop('linewidth', None)))
    markersize = kwargs.pop('markersize', kwargs.pop('ms', None))

    def _style_value(style_values, idx):
        if style_values is None or len(style_values) == 0:
            return None
        if len(style_values) == 1:
            return style_values[0]
        if idx < len(style_values):
            return style_values[idx]
        return None

    # Support for multiple lines (if rmax_list is 2D, e.g. for grouped fieldlines)
    if hasattr(rmax_list, 'ndim') and getattr(rmax_list, 'ndim', 1) > 1:
        for i, (rmaxs, vals) in enumerate(zip(rmax_list, values)):
            line_kw = {}
            color_val = _style_value(colors, i)
            label_val = _style_value(labels, i)
            linestyle_val = _style_value(linestyles, i)
            linewidth_val = _style_value(linewidths, i)
            if color_val is not None:
                line_kw['color'] = color_val
            if label_val is not None:
                line_kw['label'] = label_val
            if linestyle_val is not None:
                line_kw['linestyle'] = linestyle_val
            if linewidth_val is not None:
                line_kw['linewidth'] = linewidth_val
            if markersize is not None:
                line_kw['markersize'] = markersize
            ax.plot(rmaxs, vals, marker='o', **line_kw)
    else:
        line_kw = {}
        color_val = _style_value(colors, 0)
        label_val = _style_value(labels, 0)
        linestyle_val = _style_value(linestyles, 0)
        linewidth_val = _style_value(linewidths, 0)
        if color_val is not None:
            line_kw['color'] = color_val
        if label_val is not None:
            line_kw['label'] = label_val
        if linestyle_val is not None:
            line_kw['linestyle'] = linestyle_val
        if linewidth_val is not None:
            line_kw['linewidth'] = linewidth_val
        if markersize is not None:
            line_kw['markersize'] = markersize
        ax.plot(rmax_list, values, marker='o', **line_kw)

    ax.set_xlabel(r"$r_{\rm max}$ (R$_{\rm ns}$)")
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    return ax

@register_post
def draw_ptc_paths(data, ax, step=0, num=2, final_step=None, ids=None, **kwargs):
    """
    Draw particle paths for selected particles.
    
    Args:
        step: The step at which to select particles (default: 0)
        num: Number of particles to track (default: 2)
        final_step: Final step to plot (optional)
        ids: Specific particle IDs to track (optional)
        rmin, rmax: Radial bounds for particle selection (optional)
        thmin, thmax: Angular bounds for particle selection (optional)
    """
    ids_to_use = ids if ids is not None else get_ptc_id(data, step, num, **kwargs)
    df = particle_series(data, ids_to_use, ['tracked_ptc_x1', 'tracked_ptc_x2'])
    
    scatter_points = []
    for ptc_id in ids_to_use:
        ptc_data = df[df['ptc_id'] == ptc_id]
        
        # Extract x1 (log r) and x2 (theta) data
        x1_data = ptc_data[ptc_data['key'] == 'tracked_ptc_x1'].sort_values('step')['value'].values
        x2_data = ptc_data[ptc_data['key'] == 'tracked_ptc_x2'].sort_values('step')['value'].values
        
        # Apply final_step limit
        if final_step is not None:
            x1_data = x1_data[:final_step]
            x2_data = x2_data[:final_step]
        
        # Convert to Cartesian coordinates
        r = np.exp(x1_data)
        th = x2_data
        x = r * np.sin(th)
        y = r * np.cos(th)
        
        points = ax.scatter(x, y, label=f'ID: {ptc_id}', s=1)
        scatter_points.append(points)
    
    return scatter_points


@register
def draw_ptc_energy(data, ax,ids, step=None, num=None, final_step=None,  **kwargs):
    """
    Plot particle energy as a function of time for selected particles.
    
    Args:
        step: The step at which to select particles (default: 0)
        num: Number of particles to track (default: 2)
        final_step: Final step to plot (optional)
        ids: Specific particle IDs to track (optional)
    """
    ids_to_use = ids if ids is not None else get_ptc_id(data, step, num, **kwargs)
    df = particle_series(data, ids_to_use, ['tracked_ptc_E'])
    labels = kwargs.get('labels', None)
    lines = []
    for i,ptc_id in enumerate(ids_to_use):
        ptc_data = df[df['ptc_id'] == ptc_id]
        energy_data = ptc_data[ptc_data['key'] == 'tracked_ptc_E'].sort_values('step')

        steps = energy_data['step'].values
        energies = energy_data['value'].values
        times = energy_data['time'].values
        if final_step is not None:
            steps = steps[:final_step]
            energies = energies[:final_step]
        label = labels[i] if labels and i < len(labels) else f'ID: {ptc_id}'
        line, = ax.plot(times, energies, label=label)

        lines.append(line)

    if kwargs.get('set_labels', True):
        ax.set_xlabel(r"$\\mathrm{step}$")
        ax.set_ylabel(r"$E$")

    return lines


@register_post
def draw_radial_lines(data, ax, thetas=None, r_max=10, colors=None, linestyle='--', linewidth=1, **kwargs):
    """
    Draw radial lines from the center.
    
    Args:
        thetas: Angles in degrees for the radial lines
        r_max: Length of radial lines (default: 10)
        colors: Colors for the lines (default: ["red", "blue", "green"])
        linestyle: Line style (default: '--')
        linewidth: Line width (default: 1)
    """
    if thetas is None:
        raise ValueError("Must provide 'thetas' in kwargs.")
    
    thetas = np.radians(thetas) if not isinstance(thetas, str) else thetas
    if colors is None:
    
        colors = ["red", "blue", "green"]
    
    if not isinstance(thetas, (list, np.ndarray)):
        thetas = [thetas]
    if not isinstance(colors, (list, np.ndarray)):
        colors = [colors]
    
    lines = []
    for i, theta in enumerate(thetas):
        r = r_max if np.isscalar(r_max) else r_max[i % len(r_max)]
        color = colors[i % len(colors)]
        x_line = [0, r * np.sin(theta)]
        y_line = [0, r * np.cos(theta)]
        line, = ax.plot(x_line, y_line, color=color, linestyle=linestyle, linewidth=linewidth)
        lines.append(line)
    
    return lines


@register_post
def draw_vertical_lines(data, ax, x_values=None, colors=None, linestyle='-', linewidth=1, labels=None, **kwargs):
    """
    Draw vertical lines at specified x-values.
    
    Args:
        x_values: X-values where vertical lines should be drawn
        colors: Colors for the lines (default: ["red", "blue", "green"])
        linestyle: Line style (default: '-')
        linewidth: Line width (default: 1)
        labels: Labels for the lines (optional)
    """
    if x_values is None:
        raise ValueError("Must provide 'x_values' in kwargs.")
    
    if colors is None:
        colors = ["red", "blue", "green"]
    
    if not isinstance(x_values, (list, np.ndarray)):
        x_values = [x_values]
    if not isinstance(colors, (list, np.ndarray)):
        colors = [colors]
    
    lines = []
    for i, x_val in enumerate(x_values):
        color = colors[i % len(colors)]
        label = labels[i] if labels and isinstance(labels, (list, tuple)) and i < len(labels) else None
        line = ax.axvline(x_val, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
        lines.append(line)
    
    return lines
@register_post
def draw_horizontal_lines(data, ax, y_values=None, colors=None, linestyle='-', linewidth=1, labels=None, **kwargs):
    """
    Draw vertical lines at specified x-values.
    
    Args:
        y_values: y-values where horizontal lines should be drawn
        colors: Colors for the lines (default: ["red", "blue", "green"])
        linestyle: Line style (default: '-')
        linewidth: Line width (default: 1)
        labels: Labels for the lines (optional)
    """
    if y_values is None:
        raise ValueError("Must provide 'y_values' in kwargs.")
    
    if colors is None:
        colors = ["red", "blue", "green"]
    
    if not isinstance(y_values, (list, np.ndarray)):
        y_values = [y_values]
    if not isinstance(colors, (list, np.ndarray)):
        colors = [colors]
    
    lines = []
    for i, y_val in enumerate(y_values):
        color = colors[i % len(colors)]
        label = labels[i] if labels and isinstance(labels, (list, tuple)) and i < len(labels) else None
        # line = ax.axvline(x_val, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
        line = ax.axhline(y_val, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
        lines.append(line)
    
    return lines




@register_post
def legend(data, ax,
           color_labels=None, colors=None,
           linestyle_labels=None, linestyles=None,
           color_title=None, linestyle_title=None,
           color_loc='upper right', linestyle_loc='upper left',
           linestyle_color='black', clear_existing=True,
           linewidth=2.0, frameon=True, separate_boxes=False, **kwargs):
    """
    Add clean custom legend(s), with color and linestyle entries.

    Examples:
        afig.legend(
            "plotname",
            color_labels={r"$r_{\\rm max}=10$": "blue", r"$r_{\\rm max}=15$": "gold"},
            linestyle_labels={r"$E_\\parallel$": "-", r"$(\\nabla\\times B)_\\parallel$": "--"}
        )

    Notes:
        - `color_labels` can be a dict, list of (label, color), or list of labels with `colors=[...]`.
        - `linestyle_labels` can be a dict, list of (label, linestyle), or list of labels with `linestyles=[...]`.
        - By default both sets are merged into one legend box. Set `separate_boxes=True` to split.
    """
    from matplotlib.lines import Line2D

    def _to_entries(labels, values, entry_name):
        if labels is None:
            return []

        if isinstance(labels, dict):
            return list(labels.items())

        label_list = to_list(labels)
        if len(label_list) == 0:
            return []

        first = label_list[0]
        if isinstance(first, (tuple, list)) and len(first) == 2:
            return [(item[0], item[1]) for item in label_list]

        if values is None:
            raise ValueError(f"{entry_name}: provide values or pass list of (label, value) pairs")

        value_list = to_list(values)
        if len(value_list) == 1 and len(label_list) > 1:
            value_list = value_list * len(label_list)
        if len(value_list) != len(label_list):
            raise ValueError(f"{entry_name}: labels and values must have the same length")

        return list(zip(label_list, value_list))

    color_entries = _to_entries(color_labels, colors, 'color_labels')
    linestyle_entries = _to_entries(linestyle_labels, linestyles, 'linestyle_labels')

    if clear_existing:
        old_legend = ax.get_legend()
        if old_legend is not None:
            old_legend.remove()

    color_handles = [
        Line2D([0], [0], color=color, linestyle='-', linewidth=linewidth, label=label)
        for label, color in color_entries
    ]
    linestyle_handles = [
        Line2D([0], [0], color=linestyle_color, linestyle=linestyle, linewidth=linewidth, label=label)
        for label, linestyle in linestyle_entries
    ]

    fontsize = kwargs.get('legend_fontsize', kwargs.get('fontsize', None))
    title_fontsize = kwargs.get('legend_title_fontsize', fontsize)
    legend_loc = kwargs.get('legend_loc', color_loc)
    legend_ncol = kwargs.get('legend_ncol', 1)
    combined_title = kwargs.get('legend_title', None)

    combined_handles = color_handles + linestyle_handles
    if color_title and color_handles and linestyle_title and linestyle_handles and combined_title is None:
        combined_title = f"{color_title}; {linestyle_title}"
    elif combined_title is None:
        combined_title = color_title or linestyle_title

    if not separate_boxes and combined_handles:
        return ax.legend(
            handles=combined_handles,
            title=combined_title,
            loc=legend_loc,
            ncol=legend_ncol,
            frameon=frameon,
            fontsize=fontsize,
            title_fontsize=title_fontsize,
        )

    if color_handles and linestyle_handles:
        legend_color = ax.legend(
            handles=color_handles,
            title=color_title,
            loc=color_loc,
            frameon=frameon,
            fontsize=fontsize,
            title_fontsize=title_fontsize,
        )
        ax.add_artist(legend_color)
        legend_style = ax.legend(
            handles=linestyle_handles,
            title=linestyle_title,
            loc=linestyle_loc,
            frameon=frameon,
            fontsize=fontsize,
            title_fontsize=title_fontsize,
        )
        return [legend_color, legend_style]

    if color_handles:
        return ax.legend(
            handles=color_handles,
            title=color_title,
            loc=color_loc,
            frameon=frameon,
            fontsize=fontsize,
            title_fontsize=title_fontsize,
        )

    if linestyle_handles:
        return ax.legend(
            handles=linestyle_handles,
            title=linestyle_title,
            loc=linestyle_loc,
            frameon=frameon,
            fontsize=fontsize,
            title_fontsize=title_fontsize,
        )

    return None


@register
def lineout(data, ax, fld_val=None, axis='theta', values=None, split_hemisphere=False, **kwargs):
    """
    Plot field value as a lineout along one axis at fixed values on the other.
    
    Args:
        fld_val: Field to plot (string, lambda, or callable)
        axis: Axis to vary ('theta' or 'r')
        values: List of values or single value on the fixed axis
        split_hemisphere: If True and axis='theta', plot north/south separately
        colors, labels, linestyles, linewidths: Can be single values or lists
    """
    # Resolve fld_val
    if callable(fld_val):
        fld = fld_val(data)
    elif isinstance(fld_val, str) and hasattr(data, fld_val):
        fld = getattr(data, fld_val)
    elif fld_val is None:
        raise ValueError("fld_val must be specified for lineout")
    else:
        raise ValueError(f"Cannot resolve fld_val: {fld_val}")
    
    # Ensure fixed_values is a list
    values = to_list(values)
    
    # Extract and normalize style arguments once
    colors_raw = kwargs.pop('colors', kwargs.pop('color', None))
    labels_raw = kwargs.pop('labels', kwargs.pop('label', None))
    linestyle_raw = kwargs.pop('linestyle', kwargs.pop('linestyles', kwargs.pop('linestyes', None)))
    linewidths_raw = kwargs.pop('linewidths', kwargs.pop('linewidth', None))

    colors = to_list(colors_raw) if colors_raw is not None else None
    labels = to_list(labels_raw) if labels_raw is not None else None
    linestyles = to_list(linestyle_raw) if linestyle_raw is not None else None
    linewidths = to_list(linewidths_raw) if linewidths_raw is not None else None

    def _style_value(style_values, idx):
        if style_values is None or len(style_values) == 0:
            return None
        if len(style_values) == 1:
            return style_values[0]
        if idx < len(style_values):
            return style_values[idx]
        return None
    
    thetas = data._theta
    lines = []
    
    # Plot lineout for each fixed value
    for i, fixed_val in enumerate(values):
        # Build line kwargs for this iteration
        line_kw = {}
        color_val = _style_value(colors, i)
        label_val = _style_value(labels, i)
        linestyle_val = _style_value(linestyles, i)
        linewidth_val = _style_value(linewidths, i)

        if color_val is not None:
            line_kw['color'] = color_val
        if label_val is not None:
            line_kw['label'] = label_val
        if linestyle_val is not None:
            line_kw['linestyle'] = linestyle_val
        if linewidth_val is not None:
            line_kw['linewidth'] = linewidth_val
        
        if axis.lower() == 'theta':
            r_ind = convert_to_index(0, np.log(fixed_val), data) if fixed_val is not None else 0
            lineout_values = fld[:, r_ind]
            
            if split_hemisphere:
                north_mask = thetas <= np.pi/2
                south_mask = thetas > np.pi/2
                
                line1 = ax.plot(thetas[north_mask], lineout_values[north_mask], 
                               **(line_kw | {'label': line_kw.get('label', '') + ' (N)' if line_kw.get('label') else None}))
                line2 = ax.plot(thetas[south_mask], lineout_values[south_mask],
                               **(line_kw | {'label': line_kw.get('label', '') + ' (S)' if line_kw.get('label') else None}))
                lines.extend([line1, line2])
            else:
                line = ax.plot(thetas, lineout_values, **line_kw)
                lines.append(line)
        
        elif axis.lower() == 'r':
            theta_ind = convert_to_index(1, fixed_val, data) if fixed_val is not None else 0
            lineout_values = fld[theta_ind, :]
            line = ax.plot(data._r, lineout_values, **line_kw)
            lines.append(line)
        
        else:
            raise ValueError(f"axis must be 'theta' or 'r', got '{axis}'")
    
    return lines


@register
def lineout_fieldline(data, ax, fld_val=None, th_foot=None, r_max=None, distance_norm=False, **kwargs):
    """
    Plot a field sampled along one or more magnetic field lines.

    Args:
        fld_val: Field to sample (string name on data, or callable(data) -> 2D array)
        th_foot: Footpoint theta (rad), scalar or list
        r_max: Alternative to th_foot, converted via th = arcsin(1/sqrt(r_max)); scalar or list
        npts: Optional number of points to resample each line to uniform arc-length spacing
        distance_norm: If True, x-axis is normalized as s/s_max
        y_shifts: Optional scalar or list of vertical offsets, one per field line

    Returns:
        List of matplotlib Line2D objects
    """
    if fld_val is None:
        raise ValueError("fld_val must be specified for draw_fieldline")

    if th_foot is None and r_max is None:
        raise ValueError("Provide either th_foot or r_max")

    # Backward-compatible aliases: axis='normalized' or x_axis='normalized'
    if th_foot is None:
        r_vals = np.atleast_1d(r_max)
        th_vals = r_max_to_th_foot(r_vals)
    else:
        th_vals = np.atleast_1d(th_foot)

    colors_raw = kwargs.pop('colors', kwargs.pop('color', None))
    labels_raw = kwargs.pop('labels', kwargs.pop('label', None))
    linestyles_raw = kwargs.pop('linestyles', kwargs.pop('linestyle', None))
    linewidths_raw = kwargs.pop('linewidths', kwargs.pop('linewidth', None))
    y_shifts_raw = kwargs.pop('y_shifts', 0.0)

    colors = to_list(colors_raw) if colors_raw is not None else None
    labels = to_list(labels_raw) if labels_raw is not None else None
    linestyles = to_list(linestyles_raw) if linestyles_raw is not None else None
    linewidths = to_list(linewidths_raw) if linewidths_raw is not None else None
    y_shifts = to_list(y_shifts_raw)

    def _style_value(style_values, idx):
        if style_values is None or len(style_values) == 0:
            return None
        if len(style_values) == 1:
            return style_values[0]
        if idx < len(style_values):
            return style_values[idx]
        return None

    lines = []
    for i, th in enumerate(th_vals):
        result = values_on_fieldline(data, th_foot=float(th), fld_val=fld_val)

        line_kw = {}
        color_val = _style_value(colors, i)
        label_val = _style_value(labels, i)
        linestyle_val = _style_value(linestyles, i)
        linewidth_val = _style_value(linewidths, i)
        y_shift_val = _style_value(y_shifts, i)

        if color_val is not None:
            line_kw['color'] = color_val
        if label_val is not None:
            line_kw['label'] = label_val
        if linestyle_val is not None:
            line_kw['linestyle'] = linestyle_val
        if linewidth_val is not None:
            line_kw['linewidth'] = linewidth_val

        if line_kw.get('label') is None:
            line_kw['label'] = rf"$\theta_{{\rm foot}}={th:.3f}$"

        s_vals = result['s']
        if distance_norm:
            smax = s_vals[-1] if len(s_vals) > 0 else 0.0
            x_vals = s_vals / (smax + 1e-30)
        else:
            x_vals = s_vals

        y_vals = result['values'] + (0.0 if y_shift_val is None else y_shift_val)
        line = ax.plot(x_vals, y_vals, **line_kw)[0]
        lines.append(line)

    return lines

@register
def compare_energies(data, ax, datas=None, csv_file='twist_Energy.csv', labels=None, **kwargs):
    """
    Compare twist energies across multiple data runs.
    
    Parameters:
    - datas: list of data objects to compare
    - csv_file: CSV filename with twist energy data
    - data_dict: dict mapping names to data objects
    - labels: optional list of labels for legend
    - afig: apt_fig object (passed automatically from framework)
    """
    if datas is None:
        return
    
    # Extract afig object from kwargs to access plot_folder
    afig = kwargs.pop('afig', None)
    plot_folder = afig.plot_folder if afig else os.getcwd()
    
    df = pd.read_csv(plot_folder + '/' + csv_file)
    
    # ax2 = ax.twinx() if len(ax.figure.axes) == 1 else ax.figure.axes[1]
    
    for i, d in enumerate(datas):
        name = get_variable_name(d) #[n for n, obj in data_dict.items() if obj == d][0]
        if labels is not None and i < len(labels):
            label = labels[i]
        else:
            label = name
    
        ec = f"twist_energy_{name}"
        tc = f"times_{name}"

        if ec not in df.columns or tc not in df.columns:
            print(f"Skipping {name}: missing columns {ec}, {tc}")
            continue

        twist_time = d._conf.get("twist_time", None)
        if twist_time is None:
            print(f"Skipping {name}: no twist_time")
            continue

        t = df[tc].to_numpy()
        y = df[ec].to_numpy()

        if twist_time < t.min() or twist_time > t.max():
            print(f"Skipping {name}: twist_time out of range")
            continue

        E_twist = np.interp(twist_time, t, y)
        if E_twist == 0:
            continue
        normalize = kwargs.get("normalize",True)
        x_axis = kwargs.get("x_axis","time")
        if normalize:
            y = y / E_twist

        if x_axis == "end_twist":
            t = t - twist_time

        ax.plot(t, y, label=label)
    
    # ax.legend()
    # ax.grid(True, alpha=0.3)
    
    return {'ax': ax}