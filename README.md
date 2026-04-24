# Aperture Plotting — How to Use

A matplotlib wrapper for visualizing spherical-coordinate PIC magnetar simulation data.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Quick Start Examples](#2-quick-start-examples)
3. [Core API Reference](#3-core-api-reference)
4. [Built-in Plotters and Post-Processors](#4-built-in-plotters-and-post-processors)
5. [Writing Your Own Plotting Functions](#5-writing-your-own-plotting-functions)
6. [Useful Functions Reference](#6-useful-functions-reference)
7. [Parameters and Customization](#7-parameters-and-customization)

---

## 1. Architecture Overview

### Files

| File | Role |
|------|------|
| `plotting.py` | Core framework — `apt_fig`, `apt_plot`, `register`/`register_post`, built-in plotters |
| `Plotting_functions.py` | Additional registered plotters and post-processors (import to activate) |
| `Useful_functions.py` | Pure utility functions for physics calculations, not registered plotters |
| `additional_plotting.py` | Legacy/supplemental plotters |

You always import `plotting` as the base. `Plotting_functions` is optional but contains most domain-specific plotters. `Useful_functions` are standalone helpers used directly in notebooks.

```python
import plotting as apl
import Plotting_functions          # registers draw_NS, lineout, legend, etc.
from Useful_functions import *     # brings in compute_twist_energy, etc.
```

### Two Key Objects

**`apt_fig`** — the figure. Holds one or more subplots, knows which data object to use, and renders everything when you call `make_fig()`.

**`apt_plot`** — a single subplot inside a figure. You rarely construct these directly; `apt_fig` creates them for you when you call a plotter.

### How Rendering Works: Two Phases

Every subplot has two queues:
- **Plotters** — produce the main visual content (colormaps, line plots, etc.)
- **Post-processors** — run after plotters and overlay annotations (NS circle, field lines, time labels, legends, etc.)

When you call `afig.colorplot("myplot", ...)` you are *queuing* a plotter, not rendering it. Rendering happens only when `make_fig()` is called. This means the same `apt_fig` can be cheaply re-rendered at a different timestep, and movies are just repeated renders.

### The `fld_val` Pattern

Everywhere a field quantity is needed, you pass either:
- A **string** attribute name on the data object: `"B3"`, `"Rho_e"`
- A **callable** `(data) -> 2D array` for derived quantities: `lambda d: d.E2*d.B3 - d.E3*d.B2` (This is the preferred method for custom plotters since it avoids hardcoding field names.)

This avoids pre-computing arrays and lets the framework re-evaluate fields at each timestep during movies.

### Registration and Dynamic Dispatch (Defining and running new plotters)

Plotters and post-processors are stored in a global registry. When you call `afig.colorplot(...)`, the framework intercepts that via `__getattr__`, looks up `"colorplot"` in the registry, and queues it. Any function you decorate with `@register` or `@register_post` becomes instantly callable the same way on any `apt_fig`.

---

## 2. Quick Start Examples

### Single colorplot

```python
import plotting as apl
import Plotting_functions

data = apl.DataSph('/path/to/simulation/')
afig = apl.apt_fig(data, "fig1")

afig.colorplot("main",
               fld_val=lambda d: d.Rho_e + d.Rho_p,
               vmin=1e0, vmax=1e6, logscale=True,
               cmap="plasma",
               title=r"$n_{e^\pm}$",
               xlim=(0, 30), ylim=(-10, 10))
afig.draw_NS("all")
afig.draw_field_line("all", r_max=[10, 20], color='white')
afig.draw_time("all")
afig.set_size(8, 8)
afig.make_fig(step=40)
```

### Multi-panel figure

Panels are placed on a grid using `pos=(row, col)`. Each panel needs a unique name string.

```python
afig = apl.apt_fig(data, "comparison")

afig.colorplot("Bphi",  fld_val="B3", pos=(0,0), vmin=-1, vmax=1, title=r"$B_\phi$")
afig.colorplot("Ephi",  fld_val="E3", pos=(0,1), vmin=-1, vmax=1, title=r"$E_\phi$")
afig.colorplot("density", fld_val=lambda d: d.Rho_e + d.Rho_p,
               pos=(1,0), logscale=True, vmin=1, vmax=1e6, title=r"$n$")

afig.draw_NS("all")
afig.draw_time("all")
afig.add_parameters("all", xlim=(0, 20), ylim=(-10, 10), fontsize=20)
afig.set_size(12, 10)
afig.make_fig(step=40)
```

### Different data per panel

Override the data source for individual panels using the `data=` keyword.

```python
afig = apl.apt_fig(data_A, "multi_run")

for i, data in enumerate([data_A, data_B, data_C]):
    afig.colorplot(f"B3_{i}", fld_val="B3", data=data, pos=(i, 0),
                   vmin=-1, vmax=1, title=f"Run {i}")

afig.draw_NS("all")
afig.make_fig(step=40)
```

### Making a movie

```python
afig.make_fig(step=0)                          # renders the first frame
afig.make_movie(save_name='my_movie.mp4', fps=10)
```

`make_movie` iterates over all available field steps. Use `start`, `end`, `increment` to control the range:

```python
afig.make_movie(save_name='movie.mp4', fps=5, start=10, end=100, increment=2)
```

### Updating to a new step (without re-building the figure)

```python
afig.update_fig(step=50)
```

This clears artists and re-renders without re-queuing plotters, which is useful for interactive exploration.

### 1D lineouts

```python
# Plot E_par vs theta at fixed radii
afig.lineout("surface", fld_val=lambda d: d.E3,
             axis='theta', values=[1.0, 1.1, 1.5],
             labels=[r"$r=1.0$", r"$r=1.1$", r"$r=1.5$"],
             legend=True, xlabel=r"$\theta$", ylabel=r"$E_\parallel$")
```

```python
# Plot B_phi vs r at fixed polar angles
afig.lineout("radial", fld_val="B3",
             axis='r', values=[0.3, 0.6, 0.9],
             colors=['blue', 'red', 'green'])
```

### Profile along a magnetic field line

```python
afig.lineout_fieldline("efield_line",
                       fld_val=lambda d: d.E3,
                       r_max=[10, 20, 30],
                       distance_norm=True,
                       labels=["rmax=10", "rmax=20", "rmax=30"],
                       legend=True)
```

### Integrated energies with CSV caching

```python
from Useful_functions import *

times, twist_E = compute_twist_energy(data, label='run_A', csv_loc='energies.csv')
times, JdotE   = integrate_JdotE(data, label='run_A', csv_loc='energies.csv')
```

Results are written to CSV so re-running the notebook skips recomputation. Pass `label=` to distinguish between multiple runs in the same CSV.

---

## 3. Core API Reference

### `apl.DataSph(path)`

Loads simulation data from `path`.

```python
data = apl.DataSph('/path/to/simulation/')
data.load(step)          # load fields for a given step index
```

Key attributes after `data.load(step)`:

| Attribute | Description |
|-----------|-------------|
| `data.B1`, `data.B2`, `data.B3` | B field components (r, θ, φ), shape `(Nθ, Nr)` |
| `data.E1`, `data.E2`, `data.E3` | E field components |
| `data.J1`, `data.J2`, `data.J3` | Current components |
| `data.Rho_e`, `data.Rho_p` | Electron/positron charge densities |
| `data._r`, `data._theta` | 1D coordinate arrays |
| `data._rv`, `data._thetav` | 2D meshgrid arrays |
| `data.x1`, `data.x2` | Cartesian plot coordinates (x = r sinθ, y = r cosθ) |
| `data.time` | Simulation time at current step |
| `data.fld_steps` | List of all available field step indices |
| `data._conf` | Config dict (`'N'`, `'dt'`, `'Bp'`, `'twist_rmax_1'`, etc.) |

### `apl.apt_fig(data, name, **kwargs)`

Creates a figure object.

```python
afig = apl.apt_fig(data, "unique_name", plot_folder='output/')
```

The `name` string must be unique across figures in the same session (used internally for figure tracking).

**Key methods:**

```python
afig.add_plot("name", pos=(row,col))     # add a blank subplot
afig.del_plot("name")                    # remove a subplot

afig.set_step(step, plots="all")         # set step for all or specific plots
afig.add_parameters("all", xlim=..., fontsize=...)  # set parameters on subplots

afig.make_fig(step=None)                 # render the figure
afig.update_fig(step=None)               # re-render at a new step
afig.make_movie(save_name, fps, ...)     # render all steps as a movie

afig.set_size(xsize, ysize)             # figure dimensions in inches
afig.set_fontsize(fontsize=20)          # global font size override
afig.print_info()                        # debug: show all subplots and their queued plotters
```

**Calling a plotter or post-processor:**

```python
afig.<plotter_name>(plots, **kwargs)
```

`plots` is either `"all"`, a single subplot name string, or a list of subplot names.

On first call to a new plotter name, `apt_fig` creates a new subplot automatically. On subsequent calls with the same name, it queues additional plotters to the existing subplot.

---

## 4. Built-in Plotters and Post-Processors

### From `plotting.py`

#### `colorplot`

2D pseudocolor plot of a field.

```python
afig.colorplot("plot_name",
    fld_val=None,           # field: string attr or callable (data) -> 2D array
    cmap="hot_and_cold",    # colormap name
    vmin=None, vmax=None,   # color scale limits
    logscale=False,         # log color scale
    symlogscale=False,      # symmetric log scale
    linthresh=1.0,          # linear threshold for symlog
    include_colorbar=True,
    cbar_pos='right',       # 'right', 'left', 'top', 'bottom'
    cbar_label=None,
    x_grid=None,            # override x grid (2D array)
    y_grid=None,            # override y grid (2D array)
)
```

The built-in `hot_and_cold` colormap is diverging blue-white-red.

#### `lineplot`

1D line plot.

```python
afig.lineplot("plot_name",
    x_val=None,     # string attr, callable, or np.ndarray
    y_val=None,     # same
    label=None,
)
```

#### `draw_time` (post-processor)

Overlays the current simulation time as a text label.

```python
afig.draw_time("all",
    fmt=r"$t = {:.2f}$",
    pos=(0.05, 0.95),       # axes-fraction coordinates
    fontsize=15,
)
```

#### `draw_field_line` (post-processor)

Draws magnetic field line contours.

```python
afig.draw_field_line("all",
    th_foot=None,       # footpoint angles in radians (scalar or list)
    r_max=None,         # or specify by max radius (scalar or list)
    color='green',
    linewidth=1,
    plot_left=False,    # also draw mirror hemisphere
)
```

`th_foot` and `r_max` are interchangeable: `r_max` is converted via `arcsin(1/sqrt(r_max))`.

---

### From `Plotting_functions.py`

Import this file to register all of the following:

```python
import Plotting_functions
```

#### `draw_NS` (post-processor)

Draws the neutron star as a filled circle.

```python
afig.draw_NS("all", Radius=1, color='black', transparency=0.5)
```

#### `lineout`

1D slices of a 2D field.

```python
afig.lineout("name",
    fld_val=None,           # field callable or string
    axis='theta',           # 'theta' (vary θ at fixed r) or 'r' (vary r at fixed θ)
    values=None,            # list of fixed-axis values (physical units)
    split_hemisphere=False,
    colors=None, labels=None, linestyles=None, linewidths=None,
)
```

#### `lineout_fieldline`

Sample a field quantity along a magnetic field line, plotted vs arc-length.

```python
afig.lineout_fieldline("name",
    fld_val=None,
    th_foot=None,           # footpoint angles (radians)
    r_max=None,             # or specify by r_max
    distance_norm=False,    # normalize x-axis to [0, 1]
    colors=None, labels=None, linestyles=None, linewidths=None,
)
```

#### `fieldline_value_vs_rmax`

Plot a scalar quantity (from a `value_func`) vs field line r_max.

```python
afig.fieldline_value_vs_rmax("name",
    value_func=my_func,             # callable(data, rmax_list, **kwargs) -> array
    rmax_list=[5, 10, 15, 20],
    ylabel=r"$y$",
    title=None,
)
```

#### `compare_energies`

Plot twist energy vs time from CSV files, for comparing multiple runs.

```python
afig.compare_energies("name",
    datas=[data_A, data_B],
    csv_file='twist_Energy.csv',
    labels=['Run A', 'Run B'],
    normalize=True,
    x_axis='time',          # 'time' or 'end_twist'
)
```

#### `legend` (post-processor)

Flexible legend overlay. Supports color-only legends, linestyle-only legends, or both, optionally in separate boxes.

```python
afig.legend("all",
    color_labels={'label A': 'red', 'label B': 'blue'},   # dict or list of (label, color) tuples
    linestyle_labels=['solid', 'dashed'],
    linestyles=['-', '--'],
    color_loc='upper right',
    linestyle_loc='upper left',
    separate_boxes=False,   # True: draw two separate legend boxes
)
```

#### `r_th_region` (post-processor)

Shade a rectangular region in (r, θ) space.

```python
afig.r_th_region("all",
    rmin=0, rmax=20, thmin=0, thmax=1.57,
    outline_color='red', fill_color='purple', alpha=0.5,
)
```

#### `draw_radial_lines` / `draw_vertical_lines` / `draw_horizontal_lines` (post-processors)

```python
afig.draw_radial_lines("all", thetas=[0.3, 0.6], r_max=10, colors=['red', 'blue'])
afig.draw_vertical_lines("all", x_values=[10, 20], colors=['gray'], labels=['r=10', 'r=20'])
afig.draw_horizontal_lines("all", y_values=[0], linestyle='--')
```

---

## 5. Writing Your Own Plotting Functions

### The Plotter Signature

Any function you register as a plotter must have this signature:

```python
def my_plotter(data, ax, **kwargs):
    ...
    return ax  # optional but conventional
```

- `data` — the simulation data object (already loaded to the correct step)
- `ax` — the matplotlib `Axes` object for this subplot
- `**kwargs` — any keyword arguments the user passed when calling `afig.my_plotter(..., kwarg=value)`

### Registering as a Plotter

```python
from plotting import register

@register
def energy_ratio(data, ax, **kwargs):
    ratio = data.Rho_e / (data.Rho_e + data.Rho_p + 1e-10)
    ax.pcolormesh(data.x1, data.x2, ratio, vmin=0, vmax=1)
    return ax
```

After this, `afig.energy_ratio("myplot")` works on any `apt_fig`.

You can also register with a custom name:

```python
register(my_func, name="custom_name")
afig.custom_name("plot1")
```

### Registering as a Post-Processor

Post-processors run after plotters and are meant for overlays.

```python
from plotting import register_post

@register_post
def draw_equator(data, ax, color='white', **kwargs):
    ax.axhline(0, color=color, linewidth=0.5, linestyle='--')
```

Then:

```python
afig.draw_equator("all", color='gray')
```

### Using the `fld_val` Pattern in Custom Plotters

If your plotter accepts a field quantity, follow the library convention so users can pass either strings or callables:

```python
@register
def my_contour(data, ax, fld_val=None, levels=10, **kwargs):
    if callable(fld_val):
        field = fld_val(data)
    else:
        field = getattr(data, fld_val)
    ax.contour(data.x1, data.x2, field, levels=levels)
    return ax
```

### Supporting Multiple Lines / Colors

If your plotter loops over a list of values, accept lists of style kwargs:

```python
@register
def my_lineout(data, ax, values=None, colors=None, labels=None, **kwargs):
    for i, val in enumerate(values or []):
        color = colors[i] if colors and i < len(colors) else None
        label = labels[i] if labels and i < len(labels) else None
        idx = ...  # find grid index for val
        ax.plot(data._theta, data.B3[:, idx], color=color, label=label)
    return ax
```

### Tips

- Access the simulation grid via `data._r`, `data._theta`, `data._rv`, `data._thetav`
- Use `data.x1`, `data.x2` as the x/y coordinates for 2D `pcolormesh` calls (they are already in Cartesian plot space)
- The data is already loaded to the right step when your function is called — you don't need to call `data.load()` yourself
- If you need to override data per-call, check for a `data` kwarg in `**kwargs` — the framework passes it through

---

## 6. Useful Functions Reference

These are standalone functions imported with `from Useful_functions import *`. They are not registered plotters — call them directly.

### Field Line Geometry

```python
vertices = compute_vertices(data, th_foot, tol=1e-2)
# Returns list of (x, y) coordinate arrays for the field line contour

result = values_on_fieldline(data, th_foot, fld_val, npts=None)
# Returns dict with keys: 'theta', 'r', 'x1', 's', 'values'
# th_foot can be scalar or list; returns list of dicts if list
```

### Coordinate Conversion

```python
th = r_max_to_th_foot(r_max)   # arcsin(1/sqrt(r_max))
r  = th_foot_to_r_max(th)      # 1/sin(th)^2
step = get_closest_step(data, target_time)
t = get_time(data, step)
```

### Physics Along Field Lines

```python
voltages = compute_voltages(data, r_maxs, npts=100)
# Parallel electric potential drop along field lines

jpar = compute_Jpar(data, r_maxs, npts=None, B_normalize=True, use='median')
# Median parallel current along field lines

L = length_of_fieldline(rmax, r_star=1)            # analytic
L = length_of_fieldline_numerical(r_max)           # numerical
```

### Volume / Surface Integrals with CSV Caching

These functions compute the integral over all available timesteps and cache to CSV. Re-running returns the cached result instantly.

```python
times, values = total_integral(data, fld_val, name="integral",
                                label="run_A", csv_loc="cache.csv",
                                rmin=None, rmax=None, thmin=None, thmax=None)

times, flux   = surface_flux(data, fld_val, r_surf,
                              name="flux", label="run_A", csv_loc="cache.csv")

times, E_twist = compute_twist_energy(data, label="run_A", csv_loc="energies.csv")

times, JdotE   = integrate_JdotE(data, label="run_A", csv_loc="energies.csv")
```

`label` distinguishes different runs in the same CSV file. Use `remove_data_from_csv(label, csv_loc)` to invalidate a cached entry.

### Differential Operators

Finite-volume curl operators on the spherical grid. Return arrays with the same shape as the inputs by default (edge-padded).

```python
curlphi   = curl_phi(data, Br=data.B1, Btheta=data.B2, Bphi=data.B3)
curlr     = curl_r(data, Br=data.B1, Btheta=data.B2, Bphi=data.B3)
curltheta = curl_theta(data, Br=data.B1, Btheta=data.B2, Bphi=data.B3)
curlpar   = curl_par(data, Br=data.B1, Btheta=data.B2, Bphi=data.B3)
```

### Time Derivatives

```python
dBphi = ddt("B3", data)           # finite difference at neighboring timestep
dEphi = ddt(lambda d: d.E3, data)
```

### Particle Tracking

```python
ids   = get_ptc_id(data, step, num=10, rmin=1, rmax=20, species='e', direction='outflow')
df    = particle_series(data, ptc_ids=ids, keys=['x1', 'p1', 'E'])
# Returns DataFrame with columns: ptc_id, key, step, time, value
```

---

## 7. Parameters and Customization

### Setting Parameters

Parameters can be set at three levels (higher priority overrides lower):

1. **Figure-level**: `afig.add_parameters("all", fontsize=20, xlim=(0,30))`
2. **Per-subplot**: `afig.add_parameters("myplot", title="B field")`
3. **Per-call**: `afig.colorplot("myplot", vmin=-1, vmax=1)`

### Common Parameters

| Parameter | Effect |
|-----------|--------|
| `xlim`, `ylim` | Axis limits: `(min, max)` |
| `xlabel`, `ylabel` | Axis labels |
| `title` | Subplot title |
| `xscale`, `yscale` | `'linear'`, `'log'`, etc. |
| `aspect` | `'equal'` or numeric |
| `fontsize` | Master font size for all text in subplot |
| `label_fontsize` | Override for axis labels only |
| `title_fontsize` | Override for title only |
| `tick_fontsize` | Override for tick labels only |
| `legend_fontsize` | Override for legend text |
| `legend`, `legend_loc` | Show/place matplotlib legend |
| `hide_xaxis`, `hide_yaxis` | Boolean, removes axis labels/ticks |
| `wspace`, `hspace` | Subplot spacing passed to `subplots_adjust` |

### Layout: Spanning Multiple Columns or Rows

```python
afig.add_plot("wide_panel", pos=(0, 0), colspan=2)
afig.add_plot("tall_panel", pos=(0, 2), rowspan=2)
```

Then queue plotters to those subplots by name:

```python
afig.colorplot("wide_panel", fld_val="B3")
```

### Saving Figures

```python
afig.make_fig(step=40, save_fig=True, save_name='my_plot.png')
```

Or set it as a parameter:

```python
afig.add_parameters("all", save_fig=True, save_name='output.pdf')
afig.make_fig(step=40)
```

### Custom Colormaps

The `hot_and_cold` diverging colormap is registered at import. Pass any matplotlib colormap name to `cmap=`.

For symmetric data (positive and negative), `hot_and_cold` with `vmin=-v, vmax=v` is a natural choice.

### LaTeX in Labels

LaTeX rendering is enabled globally at import. Use raw strings:

```python
title=r"$B_\phi / B_0$"
xlabel=r"$r / R_\mathrm{NS}$"
```

### Listing Available Plotters

```python
apl.list_plotters()
apl.list_post_processors()
```
