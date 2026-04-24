import plotting as apl
import matplotlib.pyplot as plt 
import numpy as np
import os
import matplotlib
import scipy.integrate as integrate
from scipy.interpolate import RegularGridInterpolator
import multiprocessing as mp
import pandas as pd
from tqdm import tqdm

#Just a compilation of useful functions I make

#
def compare_data_conf(data1, data2):
    conf1 = data1._conf
    conf2= data2._conf
    keys1 = set(conf1.keys())
    keys2 = set(conf2.keys())
    name1 = get_variable_name(data1) or "Data 1"
    name2 = get_variable_name(data2) or "Data 2"
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    common = keys1 & keys2

    if only_in_1:
        print(f"\nKeys only in {name1}:")
        for k in sorted(only_in_1):
            print(f"  {k}: {conf1[k]}")

    if only_in_2:
        print(f"\nKeys only in {name2}:")
        for k in sorted(only_in_2):
            print(f"  {k}: {conf2[k]}")

    print("\nDifferent values:")
    for k in sorted(common):
        v1 = conf1[k]
        v2 = conf2[k]
        if v1 != v2:
            print(f"  {k}:")
            print(f"    {name1}: {v1}")
            print(f"    {name2}: {v2}")

def convert_to_index(ax_index,value,data):
    # takes an x1,x2 (x3) value and converts it into an index on the grid
    lower = data.conf["lower"][ax_index]
    N = data.conf["N"][ax_index]
    downsample = data.conf["downsample"]
    size = data.conf["size"][ax_index]

    index = (value-lower)*N/(size*downsample)
    return int(index) # as its an index

    
def particle_series(data, ptc_ids, keys):
    '''
    Extracts time series data for specified particle IDs and keys.'''
    ptc_ids = np.atleast_1d(ptc_ids).tolist()
    keys = np.atleast_1d(keys).tolist()

    rows = []
    for step_idx, n in enumerate(data.ptc_steps):
        data.load(n)
        tracked_id = data.tracked_ptc_id
        if tracked_id is None:
            continue
        for pid in ptc_ids:
            mask = tracked_id == pid
            if not np.any(mask):
                continue
            for k in keys:
                if k =="r":
                    vals = getattr(data, "tracked_ptc_x1")[mask]
                    val = np.exp(vals[0]) if len(vals) == 1 else np.nan
                else:
                    vals = getattr(data, k)[mask]
                    val = vals[0] if len(vals) == 1 else np.nan
                rows.append({
                    "ptc_id": pid,
                    "key": k,
                    "step": n,
                    "time": data.time,
                    "value": val,
                })
    return pd.DataFrame(rows)

def compute_vertices(data, th_foot, tol=1e-2,):
    ## Make this work with a list of th_foot values
    # from contourpy import contour_generator
    '''
      th_foot must be in radians
    '''
    # compute flux which will be used to define the field line
    flux    = np.cumsum(data.B1 * data._rv * data._rv * np.sin(data._thetav) * data._dtheta, axis=0)

    # Find the index corresponding to the footpoint
    r = np.log(1) # we always start on the surface
    th_foots = np.atleast_1d(th_foot)  # Ensure th_foot is an array
    flux_list = []
    
    r_foot_index = convert_to_index(0,r,data)
    theta_foot_indices = np.array([convert_to_index(1, th, data) for th in th_foots])
    flux_list = [flux[th_ind,r_foot_index] for th_ind in theta_foot_indices]

    # Sort flux levels (required by matplotlib) and track original order
    flux_array = np.array(flux_list)
    sorted_indices = np.argsort(flux_array)
    sorted_flux_list = flux_array[sorted_indices].tolist()
    
    # Remove near-duplicate flux values (matplotlib requires strictly increasing levels)
    unique_flux = []
    unique_indices_map = {}  # Maps original index to unique flux index
    tol_flux = 1e-10  # Tolerance for considering flux values as duplicates
    
    for i, (flux_val, orig_idx) in enumerate(zip(sorted_flux_list, sorted_indices)):
        # Check if this flux value is sufficiently different from the last one
        if not unique_flux or abs(flux_val - unique_flux[-1]) > tol_flux:
            unique_indices_map[orig_idx] = len(unique_flux)
            unique_flux.append(flux_val)
        else:
            # Map to the previous unique flux
            unique_indices_map[orig_idx] = len(unique_flux) - 1
    
    # define the field line as a contour of that flux,
    # then we extract the path of the field line as a series of vertices
    with matplotlib.rc_context({'interactive':False}): # as plt.contour creates it's own figure
        fig, ax = plt.subplots()
        field_line = plt.contour(data.x1, data.x2, flux, levels=unique_flux)
        plt.close(fig)

    # Restore original order to match input theta_values
    # For each level choose the contour segment closest to the requested footpoint.
    # This avoids selecting the wrong branch when a level has multiple segments.
    all_level_segments = field_line.allsegs
    vertices_list = [None] * len(th_foots)
    for orig_idx in range(len(th_foots)):
        unique_idx = unique_indices_map[orig_idx]
        if unique_idx >= len(all_level_segments):
            continue

        segments = all_level_segments[unique_idx]
        if not segments:
            continue

        th = th_foots[orig_idx]
        x_fp = np.sin(th)
        y_fp = np.cos(th)

        best_seg = None
        best_dist = np.inf
        for seg in segments:
            if seg is None or len(seg) < 2:
                continue
            dmin = np.min(np.hypot(seg[:, 0] - x_fp, seg[:, 1] - y_fp))
            if dmin < best_dist:
                best_dist = dmin
                best_seg = seg

        vertices_list[orig_idx] = best_seg
    
    return vertices_list

def compute_voltages(data, r_maxs, npts=100,**kwargs):
    """
    Compute voltage (integral of E_parallel) along each field line for a list of r_maxs,
    using values_on_fieldline for sampling. Returns the voltage from footpoint to equator (halfway point).
    Returns: np.ndarray of voltages, one per r_max.
    """
    #access cumulative trapezoidal integration function
    from scipy.integrate import cumulative_trapezoid
    th_foots = [r_max_to_th_foot(r) for r in r_maxs]
    voltages = []
    for th in th_foots:
        # Use values_on_fieldline to get E_parallel along the field line
        line = values_on_fieldline(
            data, th,
            lambda d: (d.E1*d.B1 + d.E2*d.B2 + d.E3*d.B3) / (np.sqrt(d.B1**2 + d.B2**2 + d.B3**2) + 1e-15),
            npts=npts
        )
        vals = line['values']
        s = line['s']
        if len(vals) < 2 or np.all(np.isnan(vals)):
            voltages.append(np.nan)
            continue
        # Integrate E_parallel along the field line up to the equator (halfway)
        integral_cumulative = np.concatenate(([0], cumulative_trapezoid(vals, s)))
        voltage = integral_cumulative[len(integral_cumulative)//2]
        voltages.append(voltage)
    return np.array(voltages)
def compute_Jpar(data, r_maxs, npts=None, B_normalize=True, use='median',**kwargs):
    """
    For each r_max, compute mean or median of the field along the field line.
    use: 'mean' or 'median' (default: 'median')
    Returns: array of Jpar values (one per r_max), plus dict of stats for further analysis.
    """
    B_mag = lambda d: np.sqrt(d.B1**2 + d.B2**2 + d.B3**2)
    Jpar = lambda d: (d.J1*d.B1 + d.J2*d.B2 + d.J3*d.B3)/(B_mag(d) + 1e-15)
    if B_normalize:
        fld_val = lambda d: Jpar(d) / (B_mag(d) + 1e-15)
    else:
        fld_val = Jpar
    th_foots = [r_max_to_th_foot(r) for r in r_maxs]
    lines = values_on_fieldline(data, th_foots, fld_val, npts=npts)
    if isinstance(lines, dict):
        lines = [lines]
    means = []
    medians = []
    for line in lines:
        vals = line['values']
        if len(vals) == 0 or np.all(np.isnan(vals)):
            means.append(np.nan)
            medians.append(np.nan)
            continue
        means.append(np.nanmean(vals))
        medians.append(np.nanmedian(vals))
    means = np.array(means)
    medians = np.array(medians)
    if use == 'mean':
        Jpar_vals = means
    else:
        Jpar_vals = medians
    return Jpar_vals
def make_rmax_color_labels(rmaxs, colors, label_fmt=r"$r_{{\rm max}}={}$"):
    """
    Build a color-label mapping for legends from rmaxs and colors.

    If `colors` is longer than `rmaxs`, extra colors are ignored.
    """
    rmax_list = to_list(rmaxs)
    color_list = to_list(colors)

    def _safe_format(fmt, value):
        token = "__RMAX_VALUE__"
        escaped = fmt.replace("{}", token).replace("{", "{{").replace("}", "}}")
        escaped = escaped.replace(token, "{}")
        return escaped.format(value)

    n = min(len(rmax_list), len(color_list))
    return {_safe_format(label_fmt, rmax_list[i]): color_list[i] for i in range(n)}


def values_on_fieldline(data, th_foot, fld_val, npts=None, fill_value=np.nan, vertices_coords='xy'):
    """
    Sample a field along one or more magnetic field lines.

    Args:
        th_foot: Footpoint theta in radians (scalar or list).
        fld_val: Field to sample (string name on data, or callable(data) -> 2D array).
        npts: Optional number of points to resample each line to a uniform arc-length grid.
        fill_value: Value used when interpolation queries go out of bounds.
        vertices_coords: Coordinate meaning of compute_vertices output:
            'xy' for Cartesian (x,y), 'logr_theta' for (log(r), theta).

    Returns:
        If th_foot is scalar: dict with keys
            ['th_foot', 'theta', 'r', 'x1', 's', 'values', 'vertices']
        If th_foot is list/array: list of the same dict structure.
    """
    if callable(fld_val):
        fld = fld_val(data)
    elif isinstance(fld_val, str) and hasattr(data, fld_val):
        fld = getattr(data, fld_val)
    else:
        raise ValueError("fld_val must be a callable or a valid field name on data")

    full_nt = len(data._theta)
    full_nr = len(data._r)

    if fld.shape == (full_nt, full_nr):
        theta_grid = data._theta
        r_grid = data._r
    else:
        # Support common interior-trimmed arrays (e.g. Nt-2, Nr-2 from derivatives/curl).
        nt, nr = fld.shape
        dth = full_nt - nt
        dr = full_nr - nr

        if dth < 0 or dr < 0:
            raise ValueError(
                f"Field shape {fld.shape} is larger than grid {(full_nt, full_nr)}"
            )

        if dth % 2 != 0 or dr % 2 != 0:
            raise ValueError(
                f"Field shape {fld.shape} does not align with centered trimming from {(full_nt, full_nr)}"
            )

        th0 = dth // 2
        r0 = dr // 2
        theta_grid = data._theta[th0:th0 + nt]
        r_grid = data._r[r0:r0 + nr]

        if len(theta_grid) != nt or len(r_grid) != nr:
            raise ValueError(
                f"Failed to build coordinate grids for field shape {fld.shape}"
            )

    interpolator = RegularGridInterpolator(
        (theta_grid, r_grid),
        fld,
        method='linear',
        bounds_error=False,
        fill_value=fill_value,
    )

    th_foots = np.atleast_1d(th_foot)
    vertices_list = compute_vertices(data, th_foots)
    results = []

    for th0, vertices in zip(th_foots, vertices_list):
        if vertices is None or len(vertices) < 2:
            results.append({
                "th_foot": th0,
                "theta": np.array([]),
                "r": np.array([]),
                "x1": np.array([]),
                "s": np.array([]),
                "values": np.array([]),
                "vertices": vertices,
            })
            continue

        if vertices_coords == 'xy':
            x = vertices[:, 0]
            y = vertices[:, 1]
            r = np.sqrt(x**2 + y**2)
            theta = np.pi/2 - np.arctan2(y, x)
            s = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))))
            x1 = np.log(r + 1e-30)
        elif vertices_coords == 'logr_theta':
            x1 = vertices[:, 0]
            theta = vertices[:, 1]
            r = np.exp(x1)
            dr = np.diff(r)
            dtheta = np.diff(theta)
            r_mid = 0.5 * (r[1:] + r[:-1])
            ds = np.sqrt(dr**2 + (r_mid * dtheta)**2)
            s = np.concatenate(([0.0], np.cumsum(ds)))
        else:
            raise ValueError("vertices_coords must be 'xy' or 'logr_theta'")

        if npts is not None and npts >= 2:
            s_target = np.linspace(0.0, s[-1], int(npts))
            x1 = np.interp(s_target, s, x1)
            theta = np.interp(s_target, s, theta)
            r = np.exp(x1)
            s = s_target

        sample_points = np.column_stack((theta, r))
        values = interpolator(sample_points)

        results.append({
            "th_foot": th0,
            "theta": theta,
            "r": r,
            "x1": x1,
            "s": s,
            "values": values,
            "vertices": vertices,
        })

    return results[0] if np.ndim(th_foot) == 0 else results

## Functions to convert between rmax and th_foot
def r_max_to_th_foot(rmax):
    # sintheta^2 = 1/rmax
    rmax = np.asarray(rmax)
    return np.arcsin(1/np.sqrt(rmax))

def th_foot_to_r_max(th):
    # rmax = 1/sintheta^2
    th = np.asarray(th)
    return 1/np.sin(th)**2

def get_time(data,step):
    conf = data._conf
    return conf["dt"] * conf["fld_output_interval"] * step
def get_closest_step(data, target_time):


    """Return the step in data.fld_steps whose time is closest to target_time."""
    conf = data._conf
    times = [conf["dt"] * conf["fld_output_interval"] * s for s in data.fld_steps]
    idx = np.argmin(np.abs(np.array(times) - target_time))
    return data.fld_steps[idx]

def checknan(data,step):
    data.load(step)
    isnan = False
    for fld in ["B1","B2","B3","E1","E2","E3","Rho_e","Rho_p","J1","J2","J3"]:
        dat = getattr(data,fld)
        if np.isnan(dat).any():
            print("NaN values found in field",fld)
            isnan = True
    if not isnan:
        print("No NaN values found in field")
        data_ref = data

def to_list(x):
    if not isinstance(x,(list,tuple,np.ndarray)):
        return [x]
    return x


def stack_videos(video_files, output_file='stacked_video.mp4', labels=None, save_folder=None, fontsize=48, fontcolor='white'):

    """
    Stack multiple .mp4 videos vertically into a single video with optional labels.
    
    Args:
        video_files: List of paths to .mp4 files to stack
        output_file: Name of output file (default: 'stacked_video.mp4')
        labels: List of text labels to overlay on each video (optional)
        save_folder: Folder to save output (default: plot_folder/movies)
        fontsize: Font size for labels (default: 48)
        fontcolor: Color of label text (default: 'white')
    """
    if save_folder is None:
        save_folder = f'{plot_folder}/movies'
    
    os.makedirs(save_folder, exist_ok=True)
    
    # Strip .mp4 extension if already present
    if output_file.endswith('.mp4'):
        output_file = output_file[:-4]
    
    output_path = f'{save_folder}/{output_file}.mp4'
    
    if len(video_files) < 2:
        print('Need at least 2 videos to stack')
        return
    
    # Build ffmpeg command
    inputs = ' '.join([f'-i "{f}"' for f in video_files])
    
    # Add labels if provided
    if labels is not None:
        if len(labels) != len(video_files):
            print(f'Warning: {len(labels)} labels provided for {len(video_files)} videos')
        
        # Add text overlay to each video
        labeled_streams = []
        for i, label in enumerate(labels):
            if i < len(video_files):
                # Escape special characters for ffmpeg
                label_escaped = label.replace(':', '\\:').replace("'", "'\\\\\\''")
                labeled_streams.append(f"[{i}:v]drawtext=text='{label_escaped}':x=10:y=10:fontsize={fontsize}:fontcolor={fontcolor}:box=1:boxcolor=black@0.5:boxborderw=5[v{i}]")
        
        # Stack the labeled streams
        filter_complex = ';'.join(labeled_streams) + ';'
        
        if len(video_files) == 2:
            filter_complex += '[v0][v1]vstack=inputs=2[v]'
        else:
            filter_complex += '[v0][v1]vstack=inputs=2[vs0];'
            for i in range(2, len(video_files) - 1):
                filter_complex += f'[vs{i-2}][v{i}]vstack=inputs=2[vs{i-1}];'
            filter_complex += f'[vs{len(video_files)-3}][v{len(video_files)-1}]vstack=inputs=2[v]'
    else:
        # Build filter_complex for vertical stacking without labels
        if len(video_files) == 2:
            filter_complex = '[0:v][1:v]vstack=inputs=2[v]'
        else:
            filter_complex = '[0:v][1:v]vstack=inputs=2[v0];'
            for i in range(2, len(video_files) - 1):
                filter_complex += f'[v{i-2}][{i}:v]vstack=inputs=2[v{i-1}];'
            filter_complex += f'[v{len(video_files)-3}][{len(video_files)-1}:v]vstack=inputs=2[v]'
    
    cmd = f'ffmpeg -y {inputs} -filter_complex "{filter_complex}" -map "[v]" -c:v libx264 -crf 23 -preset fast "{output_path}"'
    
    print(f'Stacking {len(video_files)} videos...')
    os.system(cmd)
    print(f'Stacked video saved to {output_path}')


def total_integral(data, fld_val, name="integral", label=None, csv_loc=None,
                   rmin=None, rmax=None, thmin=None, thmax=None):
    r"""
    Volume-integrate an arbitrary per-cell quantity over the grid at every
    field-output step, with incremental CSV caching.

    Parameters
    ----------
    data : DataSph
        Simulation data object.
    fld_val : callable(data) -> 2D array
        Function that, given a loaded ``data`` object, returns the 2D integrand
        (same shape as ``data.B1``).  Examples::

            lambda d: d.J1*d.E1 + d.J2*d.E2 + d.J3*d.E3   # J·E
            lambda d: d.B3**2 / (8*np.pi)                   # B_phi energy density
            lambda d: (d.B**2 - B_dip**2) / (8*np.pi)       # twist energy density

    name : str
        Short identifier used for the CSV column (e.g. ``"JdotE"``, ``"twist_energy"``).
    label : str, optional
        Run label for the CSV columns.  Defaults to the variable name of *data*.
    csv_loc : str, optional
        CSV file path.  Defaults to ``<plot_folder><name>.csv``.
    rmin, rmax : float, optional
        Radial bounds for the integration region (physical r, not log r).
    thmin, thmax : float, optional
        Polar-angle bounds (radians) for the integration region.

    Returns
    -------
    times : ndarray
        Simulation times at each step.
    values : ndarray
        Volume-integrated values at each step.
    """
    file_loc = f'{plot_folder}{name}.csv' if csv_loc is None else csv_loc

    # ---- CHECK IF DATA ALREADY EXISTS FOR THIS LABEL ----
    start_step_idx = 0
    existing_df = None
    if label is None:
        label = get_variable_name(data)
    time_col = f'times_{label}'
    val_col  = f'{name}_{label}'
    if label is not None and os.path.exists(file_loc) and os.path.getsize(file_loc) > 0:
        existing_df = pd.read_csv(file_loc)

        if val_col in existing_df.columns:
            n_existing = existing_df[val_col].notna().sum()
            n_expected = len(data.fld_steps)

            if n_existing == n_expected:
                times_out = existing_df[time_col].dropna().values if time_col in existing_df.columns else None
                return times_out, existing_df[val_col].dropna().values
            elif n_existing < n_expected:
                start_step_idx = n_existing

    # ---- VOLUME ELEMENT ----
    data.load(0)
    dlogr = np.log(data._r[1]) - np.log(data._r[0])
    dtheta = data._theta[1] - data._theta[0]
    cell_vol = data._rv**3 * dlogr * dtheta * np.sin(data._thetav)
    dphi = 2 * np.pi

    # ---- REGION MASK ----
    mask = np.ones_like(data._rv, dtype=bool)
    if rmin is not None:
        mask &= (data._rv >= rmin)
    if rmax is not None:
        mask &= (data._rv <= rmax)
    if thmin is not None:
        mask &= (data._thetav >= thmin)
    if thmax is not None:
        mask &= (data._thetav <= thmax)

    # ---- INTEGRATE ----
    values = []
    times  = []

    steps_to_compute = data.fld_steps[start_step_idx:]
    for step in tqdm(steps_to_compute, desc=f"Integrating {name} for {label}"):
        data.load(step)
        integrand = fld_val(data)
        values.append(dphi * np.sum(integrand[mask] * cell_vol[mask]))
        times.append(data.time)

    times  = np.array(times)
    values = np.array(values)

    # ---- SAVE ----
    if label is not None:
        if existing_df is not None and len(existing_df) > 0:
            max_rows_needed = start_step_idx + len(times)
            if len(existing_df) < max_rows_needed:
                empty_rows = pd.DataFrame(index=range(len(existing_df), max_rows_needed))
                existing_df = pd.concat([existing_df, empty_rows], ignore_index=True)

            # Write times only if column doesn't exist yet (shared across quantities)
            if time_col not in existing_df.columns:
                existing_df.loc[start_step_idx:start_step_idx + len(times) - 1, time_col] = times
            existing_df.loc[start_step_idx:start_step_idx + len(times) - 1, val_col]  = values
            result_df = existing_df
        else:
            result_df = pd.DataFrame({
                time_col: times,
                val_col: values,
            })

        result_df.to_csv(file_loc, index=False)

    return times, values


def compute_twist_energy(data, label=None, csv_loc=None, rmin=None, rmax=None, thmin=None, thmax=None):
    """Convenience wrapper: volume-integrated twist energy (B^2 - B_dip^2)/(8pi)."""
    data.load(0)
    B_dip = data.B.copy()
    return total_integral(
        data,
        fld_val=lambda d: (d.B**2 - B_dip**2) / (8 * np.pi),
        name="twist_energy",
        label=label,
        csv_loc=csv_loc,
        rmin=rmin, rmax=rmax, thmin=thmin, thmax=thmax,
    )
def integrate_JdotE(data, label=None, csv_loc=None, rmin=None, rmax=None, thmin=None, thmax=None):
    """Convenience wrapper: volume-integrated J·E."""
    return total_integral(
        data,
        fld_val=lambda d: d.J1*d.E1 + d.J2*d.E2 + d.J3*d.E3,
        name="JdotE",
        label=label,
        csv_loc=csv_loc,
        rmin=rmin, rmax=rmax, thmin=thmin, thmax=thmax,
    )


def surface_flux(data, fld_val, r_surf, name="flux", label=None, csv_loc=None,
                 thmin=None, thmax=None):
    r"""
    Compute the radial flux of a vector quantity through a spherical surface at
    ``r = r_surf`` at each field-output step, with CSV caching.

    The surface integral is:

    .. math::
        \Phi = \oint F_r \, r^2 \sin\theta \, d\theta \, d\phi
             = 2\pi \, r_{\rm surf}^2 \sum_i F_{r,i} \sin\theta_i \, \Delta\theta

    Parameters
    ----------
    data : DataSph
        Simulation data object.
    fld_val : callable(data) -> 2D array
        Function returning the **radial component** of the quantity whose flux
        you want.  Must have shape ``(Ntheta, Nr)`` matching the grid.
        Examples::

            lambda d: d.E2*d.B3 - d.E3*d.B2   # Poynting flux S_r (Heaviside, c=1)
            lambda d: d.B1                       # magnetic flux
            lambda d: d.J1                       # current flux

    r_surf : float
        Radius of the spherical surface (physical r, not log r).
    name : str
        Column-name identifier in the CSV.
    label : str, optional
        Run label.  Defaults to variable name of *data*.
    csv_loc : str, optional
        CSV path.  Defaults to ``<plot_folder><name>.csv``.
    thmin, thmax : float, optional
        Polar-angle bounds (radians) to restrict the surface patch.

    Returns
    -------
    times : ndarray
    flux_values : ndarray
    """
    file_loc = f'{plot_folder}{name}.csv' if csv_loc is None else csv_loc

    # ---- CHECK CACHE ----
    start_step_idx = 0
    existing_df = None
    if label is None:
        label = get_variable_name(data)
    time_col = f'times_{label}'
    val_col  = f'{name}_{label}'
    if label is not None and os.path.exists(file_loc) and os.path.getsize(file_loc) > 0:
        existing_df = pd.read_csv(file_loc)
        if val_col in existing_df.columns:
            n_existing = existing_df[val_col].notna().sum()
            n_expected = len(data.fld_steps)
            if n_existing == n_expected:
                times_out = existing_df[time_col].dropna().values if time_col in existing_df.columns else None
                return times_out, existing_df[val_col].dropna().values
            elif n_existing < n_expected:
                start_step_idx = n_existing

    # ---- FIND RADIAL INDEX CLOSEST TO r_surf ----
    data.load(0)
    r_idx = np.argmin(np.abs(data._r - r_surf))
    r_actual = data._r[r_idx]

    # Surface area element for each theta cell: r^2 sin(theta) dtheta * 2pi
    dtheta = data._theta[1] - data._theta[0]
    dphi = 2 * np.pi
    theta_1d = data._theta
    dA = r_actual**2 * np.sin(theta_1d) * dtheta  # shape (Ntheta,)

    # Theta mask
    th_mask = np.ones(len(theta_1d), dtype=bool)
    if thmin is not None:
        th_mask &= (theta_1d >= thmin)
    if thmax is not None:
        th_mask &= (theta_1d <= thmax)

    # ---- INTEGRATE ----
    flux_vals = []
    times = []

    steps_to_compute = data.fld_steps[start_step_idx:]
    for step in tqdm(steps_to_compute, desc=f"Surface flux {name} for {label}"):
        data.load(step)
        Fr = fld_val(data)
        # Extract radial slice and sum over theta
        Fr_slice = Fr[:, r_idx]  # shape (Ntheta,)
        flux_vals.append(dphi * np.sum(Fr_slice[th_mask] * dA[th_mask]))
        times.append(data.time)

    times = np.array(times)
    flux_vals = np.array(flux_vals)

    # ---- SAVE ----
    if label is not None:
        if existing_df is not None and len(existing_df) > 0:
            max_rows_needed = start_step_idx + len(times)
            if len(existing_df) < max_rows_needed:
                empty_rows = pd.DataFrame(index=range(len(existing_df), max_rows_needed))
                existing_df = pd.concat([existing_df, empty_rows], ignore_index=True)
            if time_col not in existing_df.columns:
                existing_df.loc[start_step_idx:start_step_idx + len(times) - 1, time_col] = times
            existing_df.loc[start_step_idx:start_step_idx + len(times) - 1, val_col] = flux_vals
            result_df = existing_df
        else:
            result_df = pd.DataFrame({time_col: times, val_col: flux_vals})
        result_df.to_csv(file_loc, index=False)

    return times, flux_vals


def remove_data_from_csv(label, csv_loc=None):
    """Remove a label's columns from the twist energy CSV."""
    file_loc = f'{plot_folder}twist_energy.csv' if csv_loc is None else csv_loc
    
    if not os.path.exists(file_loc):
        print(f"File {file_loc} does not exist")
        return
    
    df = pd.read_csv(file_loc)
    time_col = f'times_{label}'
    energy_col = f'twist_energy_{label}'
    
    cols_to_drop = [col for col in [time_col, energy_col] if col in df.columns]
    
    if not cols_to_drop:
        print(f"No columns found for label '{label}'")
        return
    
    df = df.drop(columns=cols_to_drop)
    df.to_csv(file_loc, index=False)
    print(f"Removed columns for '{label}': {cols_to_drop}")


def get_ptc_id(data, step, num, **kwargs):
    """
    At a step, looks within an (r, th) rectangle and returns the ids of interior particles.
    """
    rmin = kwargs.get('rmin', 0)
    rmax = kwargs.get('rmax', 15)
    thmin = kwargs.get('thmin', 0)
    thmax = kwargs.get('thmax', np.pi / 2)
    species = kwargs.get('species', False)
    consider_radiation = kwargs.get('consider_radiation', False)
    direction = kwargs.get('direction', None) # 'outflow or inflow
    data.load(step)
    r = np.exp(data.tracked_ptc_x1)
    th = data.tracked_ptc_x2
    pr = data.tracked_ptc_p1
    if direction is None:
        mask = (r >= rmin) & (r <= rmax) & (th >= thmin) & (th <= thmax)
    elif direction == 'outflow':
        mask = (r >= rmin) & (r <= rmax) & (th >= thmin) & (th <= thmax) & (pr > 0)
    elif direction == 'inflow':
        mask = (r >= rmin) & (r <= rmax) & (th >= thmin) & (th <= thmax) & (pr < 0)
    else:
        raise ValueError("Invalid direction specified. Use None, 'outflow', or 'inflow'.")
    if species is not False:
        flag_val = apl.flag_to_species(data.tracked_ptc_flag)
        mask = mask & (flag_val == species)
    current_ids = data.tracked_ptc_id.copy()
    current_mask = mask.copy()
    valid_ids = current_ids[current_mask]
    if consider_radiation is not False:
        data.load(0)
        zero_ids = data.tracked_ptc_id
        ids_not_in_step0 = np.setdiff1d(valid_ids, zero_ids)
        mask = mask & np.isin(current_ids, ids_not_in_step0)
    valid_ids = current_ids[mask]
    if len(valid_ids) > num:
        return np.random.choice(valid_ids, num)
    else:
        return valid_ids


def get_footpoint_ptc_id(data, step, rmax, dtheta=0.01, dr=0.2,**kwargs):
    """Access the id of footpoint particles at a given step for specified rmax values."""
    rmax = to_list(rmax)
    ids = []
    for r in rmax:
        th_foot = r_max_to_th_foot(r)
        found = get_ptc_id(data, step,1,rmin=1,rmax=1 + dr,thmin=th_foot - dtheta,thmax=th_foot + dtheta,**kwargs)
        ids.append(found[0]) if len(found) > 0 else ids.append(None)
    return ids


def _pad_interior_to_shape(arr_interior, target_shape, mode='edge', fill_value=np.nan):
    """Pad an interior-centered 2D array back to a target full-grid shape."""
    if arr_interior.shape == target_shape:
        return arr_interior

    d0 = target_shape[0] - arr_interior.shape[0]
    d1 = target_shape[1] - arr_interior.shape[1]

    if d0 < 0 or d1 < 0 or (d0 % 2 != 0) or (d1 % 2 != 0):
        raise ValueError(
            f"Cannot symmetrically pad array of shape {arr_interior.shape} to {target_shape}"
        )

    pad0 = (d0 // 2, d0 // 2)
    pad1 = (d1 // 2, d1 // 2)

    if mode == 'edge':
        return np.pad(arr_interior, (pad0, pad1), mode='edge')
    if mode in ('nan', 'constant'):
        cval = np.nan if mode == 'nan' else fill_value
        return np.pad(arr_interior, (pad0, pad1), mode='constant', constant_values=cval)

    raise ValueError("mode must be one of {'edge', 'nan', 'constant'}")


def _center_crop_to_shape(arr, target_shape):
    """Center-crop a 2D array to target_shape, requiring symmetric trimming."""
    if arr.shape == target_shape:
        return arr

    d0 = arr.shape[0] - target_shape[0]
    d1 = arr.shape[1] - target_shape[1]

    if d0 < 0 or d1 < 0 or (d0 % 2 != 0) or (d1 % 2 != 0):
        raise ValueError(
            f"Cannot symmetrically crop array of shape {arr.shape} to {target_shape}"
        )

    s0 = d0 // 2
    s1 = d1 // 2
    return arr[s0:s0 + target_shape[0], s1:s1 + target_shape[1]]

def curl_phi(data,Br=None,Btheta=None, Bphi=None, preserve_shape=True, pad_mode='edge', fill_value=np.nan):
    """
    Finite-volume curl_phi on an (theta, r) grid.
    Assumes Br, Btheta are cell-centered with axis 0 = theta, axis 1 = r
    """
    #basically allow for a different vector
    if Br is None and Btheta is None:
        Br = data.B1
        Btheta = data.B2
    r = data._r
    theta = data._theta
    Nt, Nr = Br.shape  # Swapped shape extraction

    # Cell spacings (at faces)
    dr = r[1:] - r[:-1]            # shape (Nr-1,)
    dtheta = theta[1:] - theta[:-1]  # shape (Nt-1,)

    # Face-centered radii
    r_face = 0.5 * (r[1:] + r[:-1])   # shape (Nr-1,)

    # Face-centered fields (average along the appropriate axis)
    Btheta_rface = 0.5 * (Btheta[:, 1:] + Btheta[:, :-1])   # r-faces, shape (Nt, Nr-1)
    Br_thetaface = 0.5 * (Br[1:, :] + Br[:-1, :])           # theta-faces, shape (Nt-1, Nr)

    # Interior cell centers and spacings
    r_c = r[1:-1][None, :]           # shape (1, Nr-2)
    # Centered cell width: distance between adjacent face midpoints
    dr_c = 0.5 * (dr[1:] + dr[:-1])[None, :]   # shape (1, Nr-2)
    dtheta_c = 0.5 * (dtheta[1:] + dtheta[:-1])[:, None]  # shape (Nt-2, 1)

    # Radial circulation: (r Bθ)_{j+1/2} - (r Bθ)_{j-1/2}
    term_r = (
        r_face[1:][None, :] * Btheta_rface[1:-1, 1:]
        - r_face[:-1][None, :] * Btheta_rface[1:-1, :-1]
    ) / dr_c

    # Angular circulation: Br_{i+1/2} - Br_{i-1/2}
    term_theta = (
        Br_thetaface[1:, 1:-1]
        - Br_thetaface[:-1, 1:-1]
    ) / dtheta_c

    curl_phi_val = (term_r - term_theta) / r_c

    if preserve_shape:
        return _pad_interior_to_shape(curl_phi_val, Br.shape, mode=pad_mode, fill_value=fill_value)
    return curl_phi_val
def curl_r(data,Br=None,Btheta=None, Bphi=None, preserve_shape=True, pad_mode='edge', fill_value=np.nan):
    r"""
    Finite-volume $(\nabla \times \mathbf{B})_r$ on an (theta, r) grid.
    
    Computes the radial component of curl B using the finite-volume method.
    Returns shape (Nt-2, Nr-2) after interior truncation.
    
    Assumes Br, Btheta, Bphi are cell-centered with axis 0 = theta, axis 1 = r.
    
    Formula: $(\nabla \times B)_r = \frac{1}{r\sin\theta} \frac{\partial (\sin\theta\, B_\phi)}{\partial \theta}$
    
    For axisymmetric fields (∂/∂φ = 0).
    """
    if Bphi is None and Btheta is None:
        Bphi = data.B3
        Btheta = data.B2
    r = data._r
    theta = data._theta
    Nt, Nr = Bphi.shape

    # Cell spacings (at faces)
    dtheta = theta[1:] - theta[:-1]  # shape (Nt-1,)
    
    # Face-centered theta values and their sin
    theta_face = 0.5 * (theta[1:] + theta[:-1])  # shape (Nt-1,)
    sin_theta_face = np.sin(theta_face)            # shape (Nt-1,)
    
    # Face-centered fields
    Bphi_thetaface = 0.5 * (Bphi[1:, :] + Bphi[:-1, :])  # theta-faces, shape (Nt-1, Nr)
    
    # Interior cell centers and spacings
    r_c = r[1:-1][None, :]           # shape (1, Nr-2)
    theta_c = theta[1:-1][:, None]   # shape (Nt-2, 1)
    sin_theta_c = np.sin(theta_c)    # shape (Nt-2, 1)
    # Centered cell width: distance between adjacent face midpoints
    dtheta_c = 0.5 * (dtheta[1:] + dtheta[:-1])[:, None]  # shape (Nt-2, 1)
    
    # Angular derivative: ∂(sinθ Bφ)/∂θ at interior points
    sinth_Bphi_hi = sin_theta_face[1:, None]  * Bphi_thetaface[1:, 1:-1]   # (Nt-2, Nr-2)
    sinth_Bphi_lo = sin_theta_face[:-1, None] * Bphi_thetaface[:-1, 1:-1]  # (Nt-2, Nr-2)
    
    curl_r_val = (sinth_Bphi_hi - sinth_Bphi_lo) / dtheta_c / (r_c * sin_theta_c)

    if preserve_shape:
        return _pad_interior_to_shape(curl_r_val, Bphi.shape, mode=pad_mode, fill_value=fill_value)
    return curl_r_val


def curl_theta(data,Br=None,Btheta=None, Bphi=None, preserve_shape=True, pad_mode='edge', fill_value=np.nan):
    r"""
    Finite-volume $(\nabla \times \mathbf{B})_\theta$ on an (theta, r) grid.
    
    Computes the polar (theta) component of curl B using the finite-volume method.
    Returns shape (Nt-2, Nr-2) after interior truncation.
    
    Assumes Br, Btheta, Bphi are cell-centered with axis 0 = theta, axis 1 = r.
    
    Formula: $(\nabla \times B)_\theta = -\frac{1}{r} \frac{\partial (r B_\phi)}{\partial r}$
    
    For axisymmetric fields (∂/∂φ = 0).
    """
    if Br is None and Bphi is None:
        Br = data.B1
        Bphi = data.B3
    r = data._r
    theta = data._theta
    Nt, Nr = Br.shape

    # Cell spacings (at faces)
    dr = r[1:] - r[:-1]  # shape (Nr-1,)
    
    # Face-centered fields
    Bphi_rface = 0.5 * (Bphi[:, 1:] + Bphi[:, :-1])  # r-faces, shape (Nt, Nr-1)
    
    # Interior cell centers and spacings
    r_c = r[1:-1][None, :]           # shape (1, Nr-2)
    # Centered cell width: distance between adjacent face midpoints
    dr_c = 0.5 * (dr[1:] + dr[:-1])[None, :]  # shape (1, Nr-2)
    
    # Radial derivative: ∂(r Bφ)/∂r at interior points
    r_face = 0.5 * (r[1:] + r[:-1])  # shape (Nr-1,)
    r_Bphi_r = r_face[1:][None, :] * Bphi_rface[1:-1, 1:] - r_face[:-1][None, :] * Bphi_rface[1:-1, :-1]
    
    curl_theta_val = -r_Bphi_r / dr_c / r_c

    if preserve_shape:
        return _pad_interior_to_shape(curl_theta_val, Br.shape, mode=pad_mode, fill_value=fill_value)
    return curl_theta_val
def make_curl_grid(data):
    r = data._rv[1:-1, 1:-1]
    theta = data._thetav[1:-1, 1:-1]
    X = r * np.sin(theta)
    Y = r * np.cos(theta)
    return X, Y

def make_data_for_curl(grid):
    return grid[1:-1, 1:-1]
def curl_par(data,Br=None,Btheta=None,Bphi=None, preserve_shape=True, pad_mode='edge', fill_value=np.nan):
    r"""
    Compute $(\nabla \times \mathbf{B})_\parallel$, the component of curl B parallel to the magnetic field.
    
    Returns shape (Nt-2, Nr-2) after interior truncation, matching curl_r, curl_theta, curl_phi.
    Note I say B but can replace Br,Btheta,Bphi with any vector field of the same shape.
    Formula:
    $(\nabla \times \mathbf{B})_\parallel = \frac{(\nabla \times \mathbf{B}) \cdot \mathbf{B}}{|\mathbf{B}|}$
    
    This represents the component of the curl aligned with the local magnetic field direction.
    """
    if Br is None and Btheta is None and Bphi is None:
        Br = data.B1
        Btheta = data.B2
        Bphi = data.B3
    elif Br is not None and Btheta is not None and Bphi is not None:
        pass  # Use provided components
    else:
        raise ValueError("Must provide either all of Br, Btheta, Bphi or none of them.")
    # Compute all three components of curl B with consistent shape mode.
    curl_r_val = curl_r(
        data, Br, Btheta, Bphi,
        preserve_shape=preserve_shape,
        pad_mode=pad_mode,
        fill_value=fill_value,
    )
    curl_theta_val = curl_theta(
        data, Br, Btheta, Bphi,
        preserve_shape=preserve_shape,
        pad_mode=pad_mode,
        fill_value=fill_value,
    )
    curl_phi_val = curl_phi(
        data, Br, Btheta, Bphi,
        preserve_shape=preserve_shape,
        pad_mode=pad_mode,
        fill_value=fill_value,
    )

    # Align all arrays to a common centered shape for robust dot products.
    common_shape = (
        min(curl_r_val.shape[0], curl_theta_val.shape[0], curl_phi_val.shape[0], Br.shape[0]),
        min(curl_r_val.shape[1], curl_theta_val.shape[1], curl_phi_val.shape[1], Br.shape[1]),
    )

    curl_r_use = _center_crop_to_shape(curl_r_val, common_shape)
    curl_theta_use = _center_crop_to_shape(curl_theta_val, common_shape)
    curl_phi_use = _center_crop_to_shape(curl_phi_val, common_shape)
    B1_use = _center_crop_to_shape(Br, common_shape)
    B2_use = _center_crop_to_shape(Btheta, common_shape)
    B3_use = _center_crop_to_shape(Bphi, common_shape)

    B_mag = np.sqrt(B1_use**2 + B2_use**2 + B3_use**2 + 1e-30)
    curl_dot_B = curl_r_use * B1_use + curl_theta_use * B2_use + curl_phi_use * B3_use
    curl_par_val = curl_dot_B / B_mag
    
    if preserve_shape:
        return _pad_interior_to_shape(curl_par_val, Br.shape, mode=pad_mode, fill_value=fill_value)
    return curl_par_val

def vec_cart_to_sph(r,theta, Bx,By,Bz):
    """
    Convert Cartesian magnetic field components (Bx, By, Bz) to spherical components (Br, Btheta, Bphi).
    
    Parameters:
    - r: radial distance
    - theta: polar angle (colatitude)
    - Bx, By, Bz: Cartesian magnetic field components
    
    Returns:
    - Br, Btheta, Bphi: Spherical magnetic field components
    """
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    
    Br = (Bx * sin_theta * np.cos(0) + By * sin_theta * np.sin(0) + Bz * cos_theta)
    Btheta = (Bx * cos_theta * np.cos(0) + By * cos_theta * np.sin(0) - Bz * sin_theta)
    Bphi = (-Bx * np.sin(0) + By * np.cos(0))
    
    return Br, Btheta, Bphi
def vec_sph_to_cart(r, theta, Br, Btheta, Bphi):
    """
    Convert spherical magnetic field components (Br, Btheta, Bphi) to Cartesian components (Bx, By, Bz).
    
    Parameters:
    - r: radial distance
    - theta: polar angle (colatitude)
    - Br, Btheta, Bphi: Spherical magnetic field components
    
    Returns:
    - Bx, By, Bz: Cartesian magnetic field components
    """
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    
    Bx = (Br * sin_theta * np.cos(0) + Btheta * cos_theta * np.cos(0) - Bphi * np.sin(0))
    By = (Br * sin_theta * np.sin(0) + Btheta * cos_theta * np.sin(0) + Bphi * np.cos(0))
    Bz = (Br * cos_theta - Btheta * sin_theta)
    
    return Bx, By, Bz



def varname(var, scope=None):
    """
    Get the name of a variable as a string.
    
    Parameters:
    - var: the variable object
    - scope: dict to search in (default: caller's locals())
    
    Returns:
    - str: variable name, or None if not found
    
    Example:
    >>> my_data = DataSph(...)
    >>> varname(my_data)
    'my_data'
    """
    import inspect
    
    if scope is None:
        # Get caller's local variables
        frame = inspect.currentframe()
        scope = frame.f_back.f_locals
    
    for name, value in scope.items():
        if value is var:
            return name
    
    return None
def get_variable_name(obj, skip_names=['data', 'self']):

    """
    Find the variable name by matching the object identity across stack frames.
    More reliable than parsing code strings.
    """
    import inspect
    
    frame = inspect.currentframe()
    try:
        # Search through stack frames
        for depth in range(2, 15):
            caller_frame = frame
            for _ in range(depth):
                caller_frame = caller_frame.f_back
                if caller_frame is None:
                    break
            
            if caller_frame is None:
                continue
            
            # Look for a variable that is the same object
            for var_name, var_value in caller_frame.f_locals.items():
                if var_value is obj:
                    # Skip internal names and common parameter names
                    if not var_name.startswith('_') and var_name not in skip_names:
                        return var_name
    finally:
        del frame
    
    return None
def ensure_timeseries_keys(data, ptc_ids, base_df, required_keys):
    """
    Ensure the dataframe has all required keys; if not, fetch missing via particle_series
    and return a unified dataframe.
    """
    missing = []
    if base_df is None or len(base_df) == 0:
        missing = required_keys
        base_df = pd.DataFrame()
    else:
        have = set(base_df['key'].unique())
        missing = [k for k in required_keys if k not in have]

    if missing:
        fetched = particle_series(data, ptc_ids, missing)
        if len(base_df):
            base_df = pd.concat([base_df, fetched], ignore_index=True)
        else:
            base_df = fetched
    return base_df


def _eval_source(source, data):
    if isinstance(source, str):
        if not hasattr(data, source):
            raise ValueError(f"data has no field '{source}'")
        return np.array(getattr(data, source), copy=True)
    if callable(source):
        return np.array(source(data), copy=True)
    raise TypeError("source must be a field name string or a callable")


def _get_current_and_neighbor_fld_steps(data):
    fld_steps = np.asarray(data.fld_steps, dtype=int)
    if fld_steps.size < 2:
        raise ValueError("Need at least two field steps")

    current_fld_step = getattr(data, "_current_fld_step", None)
    if current_fld_step is None or int(current_fld_step) not in fld_steps:
        current_fld_step = int(fld_steps[0])
        data.load_fld(current_fld_step)
    else:
        current_fld_step = int(current_fld_step)

    idx = int(np.where(fld_steps == current_fld_step)[0][0])
    if idx < fld_steps.size - 1:
        step_other = int(fld_steps[idx + 1])
    elif idx > 0:
        step_other = int(fld_steps[idx - 1])
    else:
        raise ValueError("No neighboring step found")

    return current_fld_step, step_other

def ddt(source, data, eps=1e-10):
    """
    Compute time derivative of a field/component using saved field outputs only.

    Parameters
    ----------
    source : str or callable
        - str: data attribute name, e.g. "B3" or "E3".
        - callable: function/lambda taking ``data`` and returning an array, e.g. ``Epar``.
    data : DataSph-like
        Loaded simulation object with ``fld_steps``, ``load_fld(step)``, and ``time``.

    Returns
    -------
    ndarray
        Time derivative of the requested quantity.
    """
    current_fld_step, step_other = _get_current_and_neighbor_fld_steps(data)

    val_now = _eval_source(source, data)
    t_now = float(data.time)

    try:
        data.load_fld(step_other)
        val_other = _eval_source(source, data)
        t_other = float(data.time)
    finally:
        data.load_fld(current_fld_step)

    dt = t_other - t_now
    if np.isclose(dt, 0.0):
        raise ValueError("Time difference is zero; cannot compute derivative")

    return (val_other - val_now) / (dt + eps)
def dBphi_dt(data, eps=1e-10):
    """Convenience wrapper for d(B_phi)/dt using B3."""
    return ddt("B3", data, eps=eps)

def dEphi_dt(data, eps=1e-10):
    """Convenience wrapper for d(E_phi)/dt using E3."""
    return ddt("E3", data, eps=eps)


def length_of_fieldline_numerical(r_max, to_equator=True):
    r_maxs = to_list(r_max)
    lengths = []
    # invert dipole field line: theta = arcsin(sqrt(r/r_max))
    def r(theta, r_max):
    # Dipole field line: r = r_max * sin^2(theta)
        return r_max * np.sin(theta)**2

    def dr_dtheta(theta, r_max):
        return 2 * r_max * np.sin(theta) * np.cos(theta)

    def length_integral(theta, r_max=10):
        r_val = r(theta, r_max)
        dr = dr_dtheta(theta, r_max)
        return np.sqrt(dr**2 + r_val**2)
    if to_equator:
        theta_val = np.pi/2
    else:
        theta_val = np.pi
    for rmax in r_maxs:
        
        th_star = np.arcsin(np.sqrt(1/rmax))
        thmin = th_star
        thmax = np.pi/2 if to_equator else np.pi -th_star
        length = integrate.quad(length_integral, thmin, thmax, args=(rmax,))[0]

        lengths.append(length)
    return lengths if len(lengths) > 1 else lengths[0]

def length_of_fieldline(rmax, r_star=1, to_equator=True):
    #Analytic expression for length from rstar to rstar along theta for a given rmax field line, derived from the dipole geometry. Matches numerical integration.
    u = np.sqrt(1 - r_star / rmax)
    length = rmax * (
        u * np.sqrt(1 + 3*u**2)
        + (1/np.sqrt(3)) * np.arcsinh(np.sqrt(3)*u)
    )
    return length / 2 if to_equator else length


