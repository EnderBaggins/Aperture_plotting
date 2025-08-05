import plotting as apl
import matplotlib.pyplot as plt 
import numpy as np
import os
import matplotlib
import scipy.integrate as integrate
import multiprocessing as mp
def r_th_region(name="r_th_region", **kwargs):
    '''Apt_post object to show a specific r, theta region on top of a previous colorplot
        Arguments:
            rmin: minimum radius
            rmax: maximum radius
            thmin: minimum theta
            thmax: maximum theta
            outline_color: color of the outline
            fill_color: color of the fill
            alpha: transparency of the fill
            linewidth: width of the outline
        '''
    def func(self, data,ax, **kwargs):
        from scipy.ndimage import convolve
        # Extract parameters
        # data = apt_fig.data
        r_grid = data._rv
        theta_grid = data._thetav
        # y_grid = data.x2

        #Defining the region of interest
        rmin = kwargs.get('rmin', 0)  # Minimum radius
        rmax = kwargs.get('rmax', 20)  # Maximum radius
        thmin = kwargs.get('thmin', 0)  # Minimum theta
        thmax = kwargs.get('thmax', np.pi/2)  # Maximum theta

        outline_color = kwargs.get('outline_color', 'red')
        fill_color = kwargs.get('fill_color', 'purple')
        alpha = kwargs.get('alpha', 0.5)
        linewidth = kwargs.get('linewidth', 0.4)

        # for plot in self.post_plots:
        # ax = plotax
        mask = (r_grid >= rmin) & (r_grid <= rmax) & (theta_grid >= thmin) & (theta_grid <= thmax)
        # print(mask)
        alpha_mask = np.where(mask == 1, alpha, 0)  # Set alpha to 0 where mask is 0
        fill = ax.pcolormesh(data.x1,data.x2,mask,shading='auto',alpha = alpha_mask,color=fill_color)
        numerical_mask = mask.astype(int)
        outline =ax.contour(data.x1,data.x2,numerical_mask,colors=outline_color,linewidths=linewidth)
        return [fill,outline]
    return apl.apt_post(name, func, **kwargs)
# Register the new post-processing type
apl.apt_post_types["r_th_region"] = r_th_region



def e_dir_quiver(name="e_dir_quiver",**kwargs):
    '''Apt_post object to plot the electron momentum as a quiver plot
    requires moments in the data object and spherical coordinates'''
    def func(self,data,ax,**kwargs):
        linescale = kwargs.get('linescale',1)
        pr = data.stress_e01
        pth = data.stress_e02
        pphi = data.stress_e03

        # on sreen
        px = pr*np.sin(data._thetav) + pth*np.cos(data._thetav)
        py = pr*np.cos(data._thetav) - pth*np.sin(data._thetav)
        n2 = kwargs.get('r_spacing',10)
        n1 = kwargs.get('theta_spacing',10)# python flips things I forget
        
        px = px[::n1, ::n2]
        py = py[::n1, ::n2]
        x1 = data.x1[::n1, ::n2]
        x2 = data.x2[::n1, ::n2]
        magnitude = 0.001+np.sqrt(px**2 + py**2)
        log_mag = np.log(magnitude)
        px = px *log_mag/ magnitude
        py = py *log_mag/ magnitude
        from matplotlib.colors import LogNorm
        c=ax.quiver(x1,x2,px,py,scale=linescale, scale_units='xy')
            # setattr(plot,name,c)
        return c
    return apl.apt_post(name,func,**kwargs)
apl.apt_post_types["e_dir_quiver"] = e_dir_quiver

def res_ph_lineout_plot(apt_plot_object, **kwargs):
    '''
    currently needs to be made with apt_plot and then _add_plot
    Any parameters need to be added to the apt_plot_object.parameters
    The res_ph_lineout_plot will show the emitted photon spectra from the 
    resonant scattering scheme along a specific theta angle.
    Will require photon_flux_lower/upper in the config
    and resonant_ph_flux in the data keys
    '''
    ap = apt_plot_object
    data = ap.data
    ax = ap.ax
    #kwargs.setdefault (rasterized???)
    #accessing the ranges of the output photon flux structured as (energy,theta)
    lower_e = data._conf['ph_flux_lower'][0]
    upper_e = data._conf['ph_flux_upper'][0]
    lower_theta = data._conf['ph_flux_lower'][1]
    upper_theta = data._conf['ph_flux_upper'][1]
    
    #creating the meshgrid for the energy and theta values
    X1 = np.exp(np.linspace(np.log10(lower_e), np.log10(upper_e), num=data.resonant_ph_flux.shape[3])*np.log(10))
    Y1 = np.linspace(lower_theta*180/np.pi, upper_theta*180/np.pi, num=data.resonant_ph_flux.shape[2])
    X1, Y1 = np.meshgrid(X1, Y1)

    N = data._conf["N"]                # [1024, 1024]
    lower = data._conf["lower"]        # [0, 0]
    size = data._conf["size"]          # [3.4, 3.14]
    downsample = data._conf["ph_flux_downsample"]  # 16

    N_r = N[0] // downsample           # 64
    N_th = N[1] // downsample          # 64
    r = np.exp(np.linspace(lower[0], lower[0] + size[0], N_r))
    theta = np.linspace(lower[1], lower[1] + size[1], N_th)
    R_max_range = kwargs.get('R_max_range', [1, 100])  # Default range for R_max
    R_max_max = np.max(R_max_range)  # Maximum R_max value
    R_max_min = np.min(R_max_range)  # Minimum R_max value
    def mask_function(r, theta):
        R_max = r/ (np.sin(theta)**2+1e-15)  # r_max = 1/sin(theta)^2
        return (R_max >= R_max_min) & (R_max <= R_max_max)
    r2d, th2d = np.meshgrid(r, theta, indexing='ij')  # shape (N_r, N_th)
    # Create a mask for the r, theta values
    mask = mask_function(r2d, th2d)
    # Apply the mask to the resonant photon flux data which is a 4d array
    # Expand mask to match resonant_ph_flux shape
    mask_expanded = mask[:, :, np.newaxis, np.newaxis]  # shape (N_r, N_th, 1, 1)

    masked_flux = np.where(mask_expanded, data.resonant_ph_flux, 0) 
    flux_data = (1e-13+np.sum(masked_flux, axis=(0,1)))
    theta = kwargs.get('theta', [90]) # degrees
    
    if type(theta) is not list:
        theta = [theta]
    th_index = []
    for th in theta:
        th_index.append(np.argmin(np.abs(Y1[:,0]-th)))
    E_vals = X1[0,:]
    # considering changing to keV
    # E_vals is currently in electron mass units (m_e*c^2) which is about 0.511 MeV
    units = kwargs.get('units', 'm_e*c^2') # default is electron mass units
    if units == 'keV':
        E_vals = E_vals * 511 # convert to keV
    flux_at_th = []
    for i in th_index:
        flux_at_th.append(flux_data[i,:])
    # th_slice = flux_data[th_index,:]

    ap.parameters.setdefault('xscale',"log")
    ap.parameters.setdefault('yscale',"log")
    ylim = ap.parameters.get('ylim', None)
    labels = kwargs.get('labels', None)
    colors = kwargs.get('colors', ["red","blue","green","purple","orange","brown","pink","gray","olive","cyan"])
    show_legend = kwargs.get('show_legend', False)
    if hasattr(ap, 'linemade'):
        #for when I allow multiple lineouts on the same plot
        for i, line in enumerate(ap.plot_object):
            
        #     # line.set_data(E_vals, th_slice)
            if isinstance(line, list):
                if len(line) == 1:
                    line = line[0]
            line.remove()
        #     line.set_data(np.array([E_vals, flux_at_th[i]]).T)
        # return ap.plot_object
    # else:
    plot_object = []
    for i in range(len(theta)):
        label = labels[i] if labels is not None else f"theta = {theta[i]:.2f}\u00B0"
        line = ax.plot(E_vals,E_vals*flux_at_th[i],label=label,color=colors[i])
        plot_object.append(line)
    if ylim is not None:
        ax.set_ylim(ylim)
    if show_legend:
        ax.legend()
    ap.linemade = True
    # ax.set_aspect('equal')
    # print(type(line))
    return plot_object

def convert_to_index(ax_index,value,data):
        # takes an x1,x2 (x3) value and converts it into an index on the grid
        lower = data.conf["lower"][ax_index]
        N = data.conf["N"][ax_index]
        downsample = data.conf["downsample"]
        size = data.conf["size"][ax_index]

        index = (value-lower)*N/(size*downsample)
        return int(index) # as its an index
    
# def index_from_pos(data,x,y, tol=1e-2):
#     ''' no longer necessary since we have data._rv and data._thetav
#     *sigh I lie because i still need to convert x,y to index'''
#     abs_diff_x = np.abs(data.x1 - x)
#     abs_diff_y = np.abs(data.x2 - y)

#     # Find the indices of the minimum values in the absolute difference arrays
#     x_index = np.min(abs_diff_x)
#     y_index = np.min(abs_diff_y)

    
#     mask_x = abs_diff_x < tol
#     mask_y = abs_diff_y <tol
#     comb_mask = mask_x | mask_y
#     common_indices = np.argwhere(mask_x & mask_y)
#     print(common_indices)

#     average_x_index = int(np.median(common_indices[:,0]))
#     average_y_index = int(np.median(common_indices[:,1]))
#     return average_x_index,average_y_index

def data_field_line(data, th_foot, fld_val, tol=1e-2,):
    from contourpy import contour_generator
    '''
      th_foot must be in radians
    '''
    assert callable(fld_val), "fld_val must be a lambda data: function"
    
    # compute flux which will be used to define the field line
    flux    = np.cumsum(data.B1 * data._rv * data._rv * np.sin(data._thetav) * data._dtheta, axis=0)

    # Find the index corresponding to the footpoint
    r = np.log(1) # we always start on the surface
    r_foot_index = convert_to_index(0,r,data)
    theta_foot_index = convert_to_index(1,th_foot,data)
    #Then we compute the flux at that point
    flux_foot = flux[theta_foot_index,r_foot_index] # yes these seem flipped but its the way it is

    # define the field line as a contour of that flux,
    # then we extract the path of the field line as a series of vertices
    with matplotlib.rc_context({'interactive':False}): # as plt.contour creates it's own figure
        
        fig, ax = plt.subplots()
        field_line = plt.contour(data.x1, data.x2, flux, levels=[flux_foot])
        plt.close(fig)
    
    # paths = field_line.collections[0].get_paths() #With one level we get one path
    paths = field_line.get_paths()
    vertices = paths[0].vertices
    
    data_on_fld_line = []
    for (x,y) in vertices:
        r =np.sqrt(x**2+y**2)
        theta = np.pi/2-np.arctan2(y,x)
        # print(r,theta,x,y)
        x_index = convert_to_index(0,np.log(r),data)
        y_index = convert_to_index(1,theta,data)
        # x_index,y_index = index_from_pos(data,x,y,tol)
        data_on_fld_line.append(fld_val(data)[y_index,x_index])
    return data_on_fld_line,vertices
def plot_data_fld_line(data, th_foot, fld_val, tol=1e-2):
    fld_line,vertices = data_field_line(data, th_foot, fld_val, tol)
    # distances = np.concatenate(([0], np.linalg.norm(vertices[1:]-vertices[:-1],axis=1)))
    def func(ax):
        ax.scatter(vertices[:,0],vertices[:,1],c=fld_line)
        ax.set_title(f"th_foot = {th_foot*180/np.pi:.2f}\u00B0")
    # c = ax.scatter(vertices[:,0],vertices[:,1],c=fld_line)
    return func
def integrate_data_fld_line(data, th_foot, fld_val,val_label = None,int_label=None, tol=1e-2,**kwargs):
    def func(ax):
        fld_line,vertices = data_field_line(data, th_foot, fld_val, tol)
        # print(np.shape(vertices))
        # c = ax.scatter(vertices[:,0],vertices[:,1],c=fld_line)
        
        # now I integreat the field_line defined at vertices
        # to get the total charge along the field line
        distances = np.concatenate(([0], np.linalg.norm(vertices[1:]-vertices[:-1],axis=1)))
        tot_distance = np.cumsum(distances)
        integral_cumulative = np.concatenate(([0],integrate.cumulative_trapezoid(y=fld_line,x=tot_distance)))
        xs = np.linspace(0,len(fld_line),len(fld_line))
        # params = match_param(kwargs,ax.scatter)
        ax.scatter(tot_distance,fld_line,c =xs ,label = val_label,s=1)
        ax.legend(loc='center',bbox_to_anchor=(0.5, 0.6))
        # ax.set_ylim(-1500,10)
        # ax.scatter(vertices[:,0],vertices[:,1],c=fld_line)
        # ax.set_xlabel("Distance along field line")
        # ax.set_ylabel("$E \\cdot B$")
        # ax.set_title(f"$E \\cdot B$, th_foot = {th_foot*180/np.pi:.2f}\u00B0")
        if hasattr(ax, 'second_ax'):
             ax.second_ax.remove() # to work with updating
        ax.second_ax = ax.twinx()
        ax.second_ax.plot(tot_distance,integral_cumulative,label=int_label)
        ax.second_ax.set_ylabel(int_label)
        # ax.second_ax.set_ylim(-375,2.5)
        # ax.set_ylim(-1,1)
        ax.second_ax.legend(loc='center',bbox_to_anchor=(0.5, 0.4))
    
    return func
def plot_data_on_fld_line(data, th_foot, fld_val, tol=1e-2,**kwargs):
    def func(ax):
        fld_line,vertices = data_field_line(data, th_foot, fld_val, tol)
        distances = np.concatenate(([0], np.linalg.norm(vertices[1:]-vertices[:-1],axis=1)))
        tot_distance = np.cumsum(distances)
        xs = np.linspace(0,len(fld_line),len(fld_line))
        val_label = kwargs.get('val_label',None)
        ax.scatter(tot_distance,fld_line,c =xs ,label = val_label,s=1)
        if val_label is not None:
            ax.legend(loc='center',bbox_to_anchor=(0.5, 0.6))
    return func


def get_ptc_id(data,step,num,**kwargs):
    '''
    at a step, looks within an (r,th) rectangle and returns the ids of interior particles'''
    rmin = kwargs.get('rmin',0)
    rmax = kwargs.get('rmax',15)
    thmin = kwargs.get('thmin',0)
    thmax = kwargs.get('thmax',np.pi/2)
    species = kwargs.get('species',False)
    data.load(step)
    r = np.exp(data.tracked_ptc_x1)
    th = data.tracked_ptc_x2
    mask = (r >= rmin) & (r <= rmax) & (th >= thmin) & (th <= thmax)
    # valid_ids = data.tracked_ptc_id[mask]
    if species is not False:
        flag_val = apl.flag_to_species(data.tracked_ptc_flag)
        mask = mask & (flag_val == species)
    valid_ids = data.tracked_ptc_id[mask]
    if len(valid_ids) > num:
        return np.random.choice(valid_ids,num)
    else:
        return valid_ids
    
def draw_ptc_paths(name='draw_ptc_paths', **kwargs):
    '''
    Draw particle paths for selected particles.
    
    Parameters:
    -----------
    step : int
        The step at which to select particles
    num : int
        Number of particles to track
    rmin, rmax : float
        Radial bounds for particle selection (optional)
    thmin, thmax : float
        Angular bounds for particle selection (optional)
    '''
    
    def func(self, data, ax, **kwargs):
        step = kwargs.pop('step', 0)
        num = kwargs.pop('num', 2)
        if not hasattr(self, 'ptc_id_found'):
            ids = get_ptc_id(data, step, num, **kwargs)
            paths = particle_series1(data, ids, ['tracked_ptc_x1', 'tracked_ptc_x2'])
            self.ptc_id_found = ids
            self.ptc_id_paths = paths
        paths = self.ptc_id_paths
        ids = self.ptc_id_found
        scatter_points = []
        for i, id in enumerate(ids):
            final_step = kwargs.get('final_step',len(paths[i,0,:]))
            r = np.exp(paths[i, 0, :final_step])
            th = paths[i, 1, :final_step]
            x = r * np.sin(th)
            y = r * np.cos(th)
            points = ax.scatter(x, y, label=f'ID: {id}', s=1)
            scatter_points.append(points)
            
        return scatter_points
    
    return apl.apt_post(name, func, **kwargs)

# Register the new post-processor
apl.apt_post_types['draw_ptc_paths'] = draw_ptc_paths

def particle_series1(data, ptc_id, key):
    '''
    returns a 3D array with the particle ids as the first index,
    the keys as the second index, and the steps as the third index
    '''
    # if isinstance(ptc_id, list) or isinstance(ptc_id, np.ndarray):
    
    if not isinstance(ptc_id, list) and not isinstance(ptc_id, np.ndarray):
        ptc_id = [ptc_id]
    elif isinstance(ptc_id, np.ndarray):
        ptc_id = ptc_id.tolist()
    if not isinstance(key, list):
        key = [key]
    # we make an array with ptc_id as first index, key as second index, 
    # and steps as third index
    result = np.zeros((len(ptc_id), len(key), len(data.ptc_steps)))
    for n in data.ptc_steps:
        data.load_ptc(n)
        tracked_id = data.tracked_ptc_id
        if tracked_id is None:
            # print(f"No tracked particles found in step {n}. Skipping.")
            continue
    
        for i, pid in enumerate(ptc_id):
            if pid not in tracked_id:
                print(f"Particle ID {pid} not found in step {n}. Skipping.")
                continue
            mask = tracked_id == pid
            for j, k in enumerate(key):
                res = data.__getattr__(k)[mask]
                # result[i,j,n] = res[0]
                if len(res) == 1:
                    result[i,j,n] = res[0]
                else:
                    result[i,j,n] = np.nan
    
    return result

    ## Voltages section
    def convert_to_index(ax_index,value,data):
        # takes an x1,x2 (x3) value and converts it into an index on the grid
        lower = data.conf["lower"][ax_index]
        N = data.conf["N"][ax_index]
        downsample = data.conf["downsample"]
        size = data.conf["size"][ax_index]

        index = (value-lower)*N/(size*downsample)
        return int(index) # as its an index
def compute_vertices(data, th_foot, tol=1e-2,):
    ## Make this work with a list of th_foot values
    from contourpy import contour_generator
    '''
      th_foot must be in radians
    '''
    # compute flux which will be used to define the field line
    flux    = np.cumsum(data.B1 * data._rv * data._rv * np.sin(data._thetav) * data._dtheta, axis=0)

    # Find the index corresponding to the footpoint
    r = np.log(1) # we always start on the surface
    th_foots = np.atleast_1d(th_foot)  # Ensure th_foot is an array
    flux_list = []
    for th_foot in th_foots:    
        r_foot_index = convert_to_index(0,r,data)
        theta_foot_index = convert_to_index(1,th_foot,data)
        #Then we compute the flux at that point
        flux_list.append(flux[theta_foot_index,r_foot_index]) # yes these seem flipped but its the way it is

    # define the field line as a contour of that flux,
    # then we extract the path of the field line as a series of vertices
    with matplotlib.rc_context({'interactive':False}): # as plt.contour creates it's own figure
        fig, ax = plt.subplots()
        field_line = plt.contour(data.x1, data.x2, flux, levels=flux_list)
        plt.close(fig)
    # paths = field_line.collections[0].get_paths() #With one level we get one path
    vertices_list = []
    vertices_list = [path.vertices for path in field_line.get_paths()]
    return vertices_list
def data_on_vertices(data, vertices, fld_val):
    '''
    Extracts data values at specified vertices.
    
    Parameters:
    -----------
    data : Data object
        The data object containing the field values.
    vertices : array-like
        An array of shape (N, 2) where N is the number of vertices,
        and each vertex is represented by (x, y) coordinates.
    fld_val : callable
        A function that takes the data object and returns the field values.
    
    Returns:
    --------
    data_on_fld_line : array
        The field values at the specified vertices.
    '''
    x = vertices[:,0]
    y = vertices[:,1]
    r = np.sqrt(x**2+y**2)
    theta = np.pi/2-np.arctan2(y,x)

    # Vectorized index conversion
    x_indices = np.array([convert_to_index(0, np.log(r_val), data) for r_val in r])
    y_indices = np.array([convert_to_index(1, theta_val, data) for theta_val in theta])

    # data extraction   
    fld = fld_val(data)
    # print(fld)
    #clipping
    x_indices = np.clip(x_indices, 0, fld.shape[1] - 1)
    y_indices = np.clip(y_indices, 0, fld.shape[0] - 1)
    
    return fld[y_indices, x_indices]
def plot_data_fld_line(data, th_foot, fld_val,labels=None,tol=1e-2):
    '''
    Plots the field line for a given data object and footpoint angle.
    Parameters:
    -----------
    data : Data object
        The data object containing the field values.
    th_foot : float or list of floats
        The footpoint angle in radians. Can be a single value or a list of values.
    fld_val : callable
        A function that takes the data object and returns the field values.
    labels : str or list of str, optional
        Labels for the field lines. If a single string is provided, it will be used for all lines.
    tol : float, optional
        Tolerance for the contour extraction. Default is 1e-2.
    Returns:
    --------
    func : callable
        A function that takes an axis object and plots the field lines on it.
    '''
    assert callable(fld_val), "fld_val must be a lambda data: function"
    if labels is not None and np.istype(labels, str):
        labels = [labels]
    vertices = compute_vertices(data, th_foot, tol)
    fld_lines = [data_on_vertices(data, v, fld_val) for v in vertices]
    # distances = np.concatenate(([0], np.linalg.norm(vertices[1:]-vertices[:-1],axis=1)))
    def func(ax):
        for i in range(len(vertices)):
            fld_line = fld_lines[i]
            vertex= vertices[i]
            label = labels[i] if labels is not None else None
            ax.scatter(vertex[:,0],vertex[:,1],c=fld_line)
        # ax.set_title(f"th_foot = {th_foot*180/np.pi:.2f}\u00B0")
    # c = ax.scatter(vertices[:,0],vertices[:,1],c=fld_line)
    return func

def Voltage_of_line(data, th_foot,vertices=False):
    Epar = lambda data: (data.E1*data.B1 + data.E2*data.B2 + data.E3*data.B3) / np.sqrt(data.B1**2 + data.B2**2 + data.B3**2)
    #define the field line as a contour of that flux,
    th_foots = np.atleast_1d(th_foot)  # Ensure th_foot is an array
    if vertices is False or vertices is None:
        vertices = compute_vertices(data, th_foots , tol=1e-2)
    # print("Vertices shape:", np.shape(vertices))
    # setting up the list
    voltages = []
    for v in vertices:
        # Extract the field line data at the vertices
        fld_line = data_on_vertices(data, v, Epar)

    # converting 2d positions to 1d distances
        distances = np.concatenate(([0], np.linalg.norm(v[1:]-v[:-1],axis=1)))
        tot_distance = np.cumsum(distances)
    
        integral_cumulative = np.concatenate(([0],integrate.cumulative_trapezoid(y=fld_line,x=tot_distance)))
    # Now since we want to compute voltage up to equatorial plane we merely need to pick the first half of this
        voltages.append(integral_cumulative[len(integral_cumulative)//2])
        
    return voltages
def compute_voltage(args):
    data, th = args
    volt = Voltage_of_line(data, th)
    print(volt)
    return [volt]
def plot_Voltage_th(data,th_min,th_max,num_points,vert_lines = None,use_rmax=False):
    def func(ax):
        voltages = []
        thetas = []
        r_maxs = []
        if vert_lines is not None:
            assert type(vert_lines) is list, "vert_lines must be a list of values"
            for line in vert_lines:
                ax.axvline(line, color='green', linestyle='--')

        thetas = np.linspace(th_min,th_max,num_points)
        # if use_dipole:
        #     #then we just use the first time step
        voltages = Voltage_of_line(data, thetas)
        if use_rmax:
            r_maxs = [1/(np.sin(th)**2) for th in thetas]
            ax.plot(r_maxs,voltages)
        else:
            ax.plot(thetas,voltages)
    
    return func


        

def res_ph_colorplot(apt_plot_object, **kwargs):
    '''
    Creates a colorplot of photon spectra across all theta values versus energy.
    Similar to res_ph_lineout_plot but shows the full theta-energy distribution.
    
    Parameters:
    -----------
    apt_plot_object : AperturePlot object
        The plot object containing the data and axes
    **kwargs:
        logscale : bool, optional
            Whether to use log scale for the colorbar (default: True)
        vmin, vmax : float, optional
            Minimum and maximum values for the colorbar
    '''
    ap = apt_plot_object
    data = ap.data
    ax = ap.ax

    # Get energy and theta ranges from config
    lower_e = data._conf['ph_flux_lower'][0]
    upper_e = data._conf['ph_flux_upper'][0]
    lower_theta = data._conf['ph_flux_lower'][1]
    upper_theta = data._conf['ph_flux_upper'][1]

    # Create meshgrid for energy and theta
    energies = np.exp(np.linspace(np.log10(lower_e), np.log10(upper_e), 
                                 num=data.resonant_ph_flux.shape[3])*np.log(10))
    thetas = np.linspace(lower_theta*180/np.pi, upper_theta*180/np.pi, 
                        num=data.resonant_ph_flux.shape[2])
    E_mesh, theta_mesh = np.meshgrid(energies, thetas)

    # Calculate flux data
    flux_data = 1e-13 + np.sum(data.resonant_ph_flux, axis=(0,1))
    
    # Multiply by energy to get E*dN/dE
    flux_data = flux_data

    # Set up the plot
    logscale = kwargs.get('logscale', True)
    vmin = kwargs.get('vmin', None)
    vmax = kwargs.get('vmax', None)

    if hasattr(ap, 'colormade'):
        ap.plot_object.remove()

    if logscale:
        plot_object = ax.pcolormesh(E_mesh, theta_mesh, flux_data,
                                  norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
                                  shading='auto')
    else:
        plot_object = ax.pcolormesh(E_mesh, theta_mesh, flux_data,
                                  vmin=vmin, vmax=vmax,
                                  shading='auto')

    # Set scales and labels
    ax.set_xscale('log')
    # ax.set_xlabel('Energy [MeV]')
    # ax.set_ylabel('θ [degrees]')

    # Add colorbar
    if not hasattr(ap, 'colormade'):
        plt.colorbar(plot_object, ax=ax, label='E*dN/dE')

    ap.colormade = True
    return plot_object
def r_max_to_th_foot(rmax):
    # sintheta^2 = 1/rmax
    rmax = np.asarray(rmax)
    return np.arcsin(1/np.sqrt(rmax))

def th_foot_to_r_max(th):
    # rmax = 1/sintheta^2
    th = np.asarray(th)
    return 1/np.sin(th)**2
def max_voltage(data,thmin,thmax,num=50,vertices=None):
    th = np.linspace(thmin,thmax,num)
    voltages = Voltage_of_line(data,th,vertices=vertices)
    abs_voltages = np.abs(voltages)
    max_index = np.argmax(abs_voltages)
    max_voltage = voltages[max_index]
    max_voltage_th = th[max_index]
    return max_voltage, max_voltage_th
def voltage_over_time(data,thmin,thmax,num=50,use_dipole=False):
    from tqdm import tqdm
    steps = data.fld_steps

    voltages = []
    times = []
    max_thetas = []
    if use_dipole:
        data.load(0)
        ths = np.linspace(thmin,thmax,num)
        vertices = compute_vertices(data,ths)
    else:
        vertices = None
    for step in tqdm(steps, desc="Calculating voltages over time"):
        data.load(step)
        max,th = max_voltage(data,thmin,thmax,num,vertices=vertices)
        # max = Voltage_of_line(data,0.56) # seems very similar to max voltage so we will stick with the max
        voltages.append(max)
        times.append(data.time)
        max_thetas.append(th)

    return times,voltages,max_thetas
def plot_field_along_curve_all_post(name="plot_field_along_curve_all_post", **kwargs):
    '''
    Apt_post object to plot field values along field lines for all data in a data_dict.

    Parameters:
    -----------
    data_dict : dict of {name: Data object}
    th_foot : float or list of floats (radians)
    fld_val : callable, returns field array from data
    step : int, time step to load
    tol : float, optional
    **kwargs : passed to plt.plot
    '''
    def func(self, data, ax, **kwargs):
        data_dict = kwargs.pop('data_dict', None)
        if data_dict is None:
            raise ValueError("data_dict must be provided")
        th_foot = kwargs.pop('th_foot', 0.56)
        fld_val = kwargs.pop('fld_val', None)
        step = ax.step
        tol = kwargs.pop('tol', 1e-2)
        plot_objects = []
        colors = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray"]
        for name, dataf in data_dict.items():
            dataf.load(step)
            color = colors.pop(0) if colors else 'black'
            vertices_list = compute_vertices(dataf, th_foot, tol)
            for i, vertices in enumerate(vertices_list):
                values = data_on_vertices(dataf, vertices, fld_val)
                distances = np.concatenate(([0], np.cumsum(np.linalg.norm(vertices[1:] - vertices[:-1], axis=1))))
                if isinstance(th_foot, (list, np.ndarray)):
                    label = f"{name}, th_foot={np.round(np.atleast_1d(th_foot)[i]*180/np.pi,2)}°"
                else:
                    label = name
                line, = ax.plot(distances, values, marker='o', label=label,color =color)
                plot_objects.append(line)
        ax.set_xlabel("Distance along field line")
        ax.set_ylabel("Epar")
        # ax.set_xlim(0, 0.2)
        ax.legend(loc='upper right')
        return plot_objects

    return apl.apt_post(name, func, **kwargs)

# Register the new post-processor
apl.apt_post_types['plot_field_along_curve_all_post'] = plot_field_along_curve_all_post
def plot_max_voltage_over_time_post(name="plot_max_voltage_over_time_post", **kwargs):
    '''
    Apt_post object to plot the maximum voltage as a function of time for a single data object.

    Parameters:
    -----------
    data : Data object
    thmin : float, minimum theta value (radians)
    thmax : float, maximum theta value (radians)
    num : int, number of theta points to sample
    use_dipole : bool, optional
    **kwargs : passed to plt.plot
    '''
    def func(self, data, ax, **kwargs):
        import pandas as pd
        thmin = r_max_to_th_foot(data._conf['twist_rmax_1'])
        thmax = r_max_to_th_foot(data._conf['twist_rmax_2'])
        thmin = kwargs.pop('thmin',thmin )
        thmax = kwargs.pop('thmax', thmax)
        num = kwargs.pop('num', 50)
        use_dipole = kwargs.pop('use_dipole', False)
        dataname = kwargs.pop("dataname",None)
        if dataname is None:
            print("To run faster input dataname to save intermediate results")
        # dataset_name = data.name if hasattr(data, 'name') else "dataset"
        csv_filename = f"{dataname}_max_voltage.csv"

         # Check if the CSV file already exists
        if os.path.exists(csv_filename):
            # Read data from the CSV file
            print(f"Loading max voltage data from {csv_filename}")
            df = pd.read_csv(csv_filename)
            times = df['time'].to_numpy()
            voltages = df['voltage'].to_numpy()
        else:
            # Calculate voltages and times
            times, voltages, _ = voltage_over_time(data, thmin, thmax, num=num, use_dipole=use_dipole)

            # Save the results to a CSV file
            df = pd.DataFrame({'time': times, 'voltage': voltages})
            df.to_csv(csv_filename, index=False)
            print(f"Saved max voltage data to {csv_filename}")

        # Plot the data
        params = apl.match_param(kwargs, ax.scatter)
        marker = kwargs.get('marker',"o")
        s = kwargs.get('s', 1)
        label = kwargs.pop('label', dataname)
        color = kwargs.get('color', 'blue')
        # line =  apl.run_function_safely(ax.scatter, times, np.abs(voltages),label=label, **params)
        line, = ax.plot(times, np.abs(voltages), marker=marker, markersize=s,label=label,color=color)
        # ax.set_xlabel("Time")
        # ax.set_ylabel("Max Voltage")
        # ax.legend(loc='upper right')
        return [line]
    # print(kwargs)
    return apl.apt_post(name, func, **kwargs)

# Register the new post-processor
apl.apt_post_types['plot_max_voltage_over_time'] = plot_max_voltage_over_time_post

def plot_field_along_curve_post(name="plot_field_along_curve_post", **kwargs):
    '''
    Apt_post object to plot field values along field lines for a single data object.

    Parameters:
    -----------
    data : Data object
    th_foot : float or list of floats (radians)
    fld_val : callable, returns field array from data
    step : int, time step to load
    tol : float, optional
    **kwargs : passed to plt.plot
    '''
    def func(self, data, ax, **kwargs):
        
        th_foot = kwargs.pop('th_foot', 0.56)
        fld_val = kwargs.pop('fld_val', None)
        step = kwargs.pop('step', 0)
        tol = kwargs.pop('tol', 1e-2)
        color = kwargs.pop('color', 'blue')
        label = kwargs.pop('label', None)

        # Load the data at the specified step
        data.load(step)

        # Compute vertices and field values
        vertices_list = compute_vertices(data, th_foot, tol)
        plot_objects = []
        for i, vertices in enumerate(vertices_list):
            values = data_on_vertices(data, vertices, fld_val)
            distances = np.concatenate(([0], np.cumsum(np.linalg.norm(vertices[1:] - vertices[:-1], axis=1))))
            if isinstance(th_foot, (list, np.ndarray)):
                line_label = label if label else f"th_foot={np.round(np.atleast_1d(th_foot)[i]*180/np.pi,2)}°"
            else:
                line_label = label if label else f"th_foot={th_foot*180/np.pi:.2f}°"
            line, = ax.plot(distances, values, marker='o', label=line_label, color=color)
            plot_objects.append(line)

        # Set axis labels
        ax.set_xlabel("Distance along field line")
        ax.set_ylabel("Epar")
        ax.legend(loc='upper right')
        return plot_objects

    return apl.apt_post(name, func, **kwargs)

# Register the new post-processor
apl.apt_post_types['plot_field_along_curve_post'] = plot_field_along_curve_post

import pandas as pd
# ...existing code...

def plot_twist_energy_post(name="plot_twist_energy_post", csv_path=None, plot_these=None, data_dict=None, **kwargs):
    def func(self, data, ax, **kwargs):
        
        # Use your existing function, but only plot the selected dataset
        plot_twist_energies(
            csv_path,
            ax=ax,
            plot_these=plot_these,
            data_dict=data_dict
        )
        # Draw vertical line at the current time for each dataset in plot_these
        if plot_these is not None and data_dict is not None:
            for d in plot_these:
                # Resolve to data object if string
                if isinstance(d, str):
                    d_obj = data_dict[d]
                else:
                    d_obj = d
                if hasattr(d_obj, "time"):
                    ax.axvline(d_obj.time, color="red", linestyle="--", linewidth=1, alpha=0.7)
    return apl.apt_post(name, func, **kwargs)

def plot_voltage_post(name="plot_voltage_post", csv_path=None, plot_these=None, data_dict=None, **kwargs):
    def func(self, data, ax, **kwargs):
        plot_max_voltages(
            csv_path,
            ax=ax,
            plot_these=plot_these,
            data_dict=data_dict
        )
        if plot_these is not None and data_dict is not None:
            for d in plot_these:
                if isinstance(d, str):
                    d_obj = data_dict[d]
                else:
                    d_obj = d
                if hasattr(d_obj, "time"):
                    ax.axvline(d_obj.time, color="red", linestyle="--", linewidth=1, alpha=0.7)

    return apl.apt_post(name, func, **kwargs)

# Register posts if needed
apl.apt_post_types["plot_twist_energy_post"] = plot_twist_energy_post
apl.apt_post_types["plot_voltage_post"] = plot_voltage_post

def plot_twist_energies(csv_path, ax=None,plot_these=None,data_dict=None):
    df = pd.read_csv(csv_path)
    # Find all datasets by looking for columns ending with '_time'
    dataset_names = [col[:-5] for col in df.columns if col.endswith('_time')]
    if data_dict is None:
        print("Warning: data_dict is None, plot_these will not be resolved.")

    if plot_these is not None:
        plot_these_names = _resolve_plot_these(plot_these, data_dict)
        
        dataset_names = [name for name in dataset_names if name in plot_these_names]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    for name in dataset_names:
        # print("Plotting twist energy for dataset:", name)
        time_col = f"{name}_time"
        energy_col = name
        # Drop NaNs for plotting
        times = df[time_col].dropna()
        energies = df[energy_col].dropna()
        if len(times) != len(energies):
            print(f"Warning: Length mismatch for {name}. Times: {len(times)}, Energies: {len(energies)}")
            min_len = min(len(times), len(energies))
            times = times[:min_len]
            energies = energies[:min_len]
        ax.plot(times, energies, label=name)
    ax.axvline(20, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('Twist Energy')
    # ax.set_title('Twist Energy vs Time')
    ax.set_yscale('log')
    ax.legend()
    # ax.set_xlim(0, 60)
    if ax is None:
        fig.tight_layout()
        plt.show()
    # fig.tight_layout()
    # plt.show()
def _resolve_plot_these(plot_these, data_dict):
    """Helper to convert plot_these (names or data objects) to a list of names."""
    if plot_these is None:
        return None
    resolved = []
    for item in plot_these:
        if isinstance(item, str):
            resolved.append(item)
        else:
            # Assume it's a data object, find its key in data_dict
            for k, v in data_dict.items():
                if v is item:
                    resolved.append(k)
                    break
    return resolved
def plot_max_voltages(csv_path,ax=None, plot_these=None,data_dict=None):
    import pandas as pd
    import matplotlib.pyplot as plt
    if data_dict is None:
        raise
    df = pd.read_csv(csv_path)
    # Find all datasets by looking for columns ending with '_max_voltage'
    dataset_names = [col[:-12] for col in df.columns if col.endswith('_max_voltage')]

    if plot_these is not None:
        plot_these_names = _resolve_plot_these(plot_these, data_dict)
        
        dataset_names = [name for name in dataset_names if name in plot_these_names]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    for name in dataset_names:
        time_col = f"{name}_time"
        voltage_col = f"{name}_max_voltage"
        if time_col not in df.columns or voltage_col not in df.columns:
            continue
        times = df[time_col].dropna()
        voltages = np.abs(df[voltage_col].dropna())
        ax.plot(times, voltages, label=name)
    ax.axvline(20, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('Max Voltage')
    ax.set_yscale('log')
    # ax.set_title('Max Voltage vs Time')
    ax.legend()
    ax.set_xlim(0, 60)
    # ax.set_yscale('log')
    if ax is None:
        fig.tight_layout()
        plt.show()


def spectra(name="spectral_lineout", **kwargs):
    '''
    Apt_post object to plot the resonant photon flux as a function of energy and theta.
    Args:
    theta : float or list of floats
        The observer theta angle(s) in degrees.
    R_max_range : the (small,large)  range of R_max for the flux tube
    units: 'keV' or default 'm_e*x^2'
    labels: labels for theta values (defaults to the theta values in degrees)
        '''
    def func(self, data, ax, **kwargs):
        #accessing the ranges of the output photon flux structured as (energy,theta)
        lower_e = data._conf['ph_flux_lower'][0]
        upper_e = data._conf['ph_flux_upper'][0]
        lower_theta = data._conf['ph_flux_lower'][1]
        upper_theta = data._conf['ph_flux_upper'][1]
        #creating the meshgrid for the energy and theta values
        X1 = np.exp(np.linspace(np.log10(lower_e), np.log10(upper_e), num=data.resonant_ph_flux.shape[3])*np.log(10))
        Y1 = np.linspace(lower_theta*180/np.pi, upper_theta*180/np.pi, num=data.resonant_ph_flux.shape[2])
        X1, Y1 = np.meshgrid(X1, Y1)
        # spatial grid including downsampling
        N = data._conf["N"]                # [1024, 1024]
        lower = data._conf["lower"]        # [0, 0]
        size = data._conf["size"]          # [3.4, 3.14]
        downsample = data._conf["ph_flux_downsample"]  # 16
        
        N_r = N[0] // downsample           # 64
        N_th = N[1] // downsample          # 64
        r = np.exp(np.linspace(lower[0], lower[0] + size[0], N_r))
        theta = np.linspace(lower[1], lower[1] + size[1], N_th)
        # accessing the specific flux tube to see the specta for
        R_max_range = kwargs.get('R_max_range', None)  # Default range for R_max
        if R_max_range is not None:
            R_max_max = np.max(R_max_range)  # Maximum R_max value
            R_max_min = np.min(R_max_range)  # Minimum R_max value
            def flux_tube_mask(r,theta):
            # mask for the r,theta values
                R_max = r/ (np.sin(theta)**2+1e-15)  # r_max = 1/sin(theta)^2  #TODO: Double check
                return (R_max >= R_max_min) & (R_max <= R_max_max)
            
            r2d, th2d = np.meshgrid(r, theta, indexing='ij')  # shape (N_r, N_th)
            mask = flux_tube_mask(r2d, th2d)
            mask_reshape = mask[:,:, np.newaxis, np.newaxis]  # shape (N_r, N_th, 1, 1)
            masked_flux = np.where(mask_reshape, data.resonant_ph_flux, 0)
            flux_data = (1e-13+np.sum(masked_flux, axis=(0, 1))) # shape (N_e, N_th)
        else:
            flux_data = (1e-13+np.sum(data.resonant_ph_flux, axis=(0, 1)))
        theta = kwargs.get('theta', [90])
        
        if type(theta) is not list:
            theta = [theta]
        th_index = []
        for th in theta:
            th_index.append(np.argmin(np.abs(Y1[:,0]-th)))
        
        labels = kwargs.get('labels', None)
        colors = kwargs.get('colors', ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray'])
    # E_vals is currently in electron mass units (m_e*c^2) which is about 0.511 MeV
        E_vals = X1[0,:]
        units = kwargs.get('units', 'm_e*c^2') # default is electron mass units
        if units == 'keV':
            E_vals = E_vals * 511 # convert to keV
        flux_at_th = []
        for i in th_index:
            flux_at_th.append(flux_data[i,:])
        
        for i in range(len(theta)):
            label = labels[i] if labels is not None else f"theta = {theta[i]:.2f}\u00B0"
            ax.plot(E_vals,E_vals*flux_at_th[i],label=label,color=colors[i])
        ax.set_yscale('log')
        ax.set_xscale('log')


    return apl.apt_post(name, func, **kwargs)

apl.apt_post_types["spectra"] = spectra
