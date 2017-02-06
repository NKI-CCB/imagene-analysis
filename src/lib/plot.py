import contextlib
from functools import wraps

import matplotlib
import matplotlib.pyplot
import matplotlib.cm
import IPython.display

import numpy as np
import xarray as xr


def display(fig):
    IPython.display.display(fig)


def _autoplot(f):

    @wraps(f)
    def autoplot_wrapper(*args, ax=None, disp=display, **kwargs):
        fig = None
        if ax is None:
            fig, ax = matplotlib.pyplot.subplots()
        try:
            f(*args, ax=ax, **kwargs)
            if fig is not None:
                disp(fig)
        finally:
            if fig is not None:
                matplotlib.pyplot.close(fig)

    return autoplot_wrapper


@contextlib.contextmanager
def figure(*args, **kwargs):
    fig = matplotlib.pyplot.figure(*args, **kwargs)
    yield(fig)
    try:
        display(fig)
    finally:
        matplotlib.pyplot.close(fig)


@contextlib.contextmanager
def subplots(*args, **kwargs):
    fig, axs = matplotlib.pyplot.subplots(*args, **kwargs)
    yield fig, axs
    display(fig)
    matplotlib.pyplot.close(fig)


def _infer_set_ticklabels(ticklabels):
    if isinstance(ticklabels, str) and ticklabels == 'index':
        set_ticklabels = False  # Matplotlib sets them automatically
        ticklabels = None
    else:
        set_ticklabels = ticklabels is not None

    return (ticklabels, set_ticklabels)


@_autoplot
def heatmap(x, y=None, z=None, mask=None, *, xticklabels=None,
            yticklabels=None,
            xlabel=None, ylabel=None, zlabel=None, zlim=None,
            row_dendrogram=False, row_dist_metric='euclidean',
            row_cluster_method='average', col_dendrogram=False,
            col_dist_metric='euclidean', col_cluster_method='average',
            ax, cmap=None):
    if row_dendrogram or col_dendrogram:
        import scipy.cluster.hierarchy as scipy_ch
        import scipy.spatial.distance as scipy_sd

    if hasattr(x, 'shape') and len(x.shape) == 2:
        # x is a matrix
        if y is not None or z is not None:
            raise TypeError("y and z should not be specified if x is matrix.")
        if isinstance(x, xr.DataArray):
            Z = x.values
            if xlabel is None:
                xlabel = x.dims[1]
            if ylabel is None:
                ylabel = x.dims[0]
            if zlabel is None:
                zlabel = x.name
            if xticklabels is None:
                xticklabels = x.coords[x.dims[1]].values
            if yticklabels is None:
                yticklabels = x.coords[x.dims[0]].values
        else:
            Z = np.array(x)
        if mask is not None:
            Zmask = np.array(mask)
        else:
            Zmask = np.zeros(Z.shape, dtype='bool')
    else:
        # Only matrix is supported for now. We should also support DataFrames
        # in tidy format, and separate x, y and z values.
        return NotImplemented()

    Zo = Z.copy()
    Zo[Zmask] = np.nan
    row_order = None
    if row_dendrogram:
        row_dist = scipy_sd.pdist(Z, row_dist_metric)
        row_linkage = scipy_ch.linkage(row_dist, row_cluster_method)
        row_order = scipy_ch.leaves_list(row_linkage)
        Zo = Zo[row_order, :]
    col_order = None
    if col_dendrogram:
        col_dist = scipy_sd.pdist(Z.T, col_dist_metric)
        col_linkage = scipy_ch.linkage(col_dist, col_cluster_method)
        col_order = scipy_ch.leaves_list(col_linkage)
        Zo = Zo[:, col_order]

    if zlim is None:
        zlim = (np.min(Z), np.max(Z))

    if zlim[0] < 0.0 and zlim[1] > 1:
        # Automatic symmetric colormap
        zlim_p = max((-zlim[0], zlim[1]))
        zlim = (-zlim_p, zlim_p)
        if cmap is None:
            cmap = 'coolwarm'

    c = ax.imshow(Zo, cmap=cmap, origin='lower', interpolation='none',
                  aspect='auto', vmin=zlim[0], vmax=zlim[1])
    cbar = ax.figure.colorbar(c)

    xticklabels, set_xticklabels = _infer_set_ticklabels(xticklabels)
    if set_xticklabels:
        if col_order is None:
            ax.set_xticks(range(len(xticklabels)))
        else:
            ax.set_xticks(col_order)
        ax.set_xticklabels(xticklabels)
    yticklabels, set_yticklabels = _infer_set_ticklabels(yticklabels)
    if set_yticklabels:
        if row_order is None:
            ax.set_yticks(range(len(yticklabels)))
        else:
            ax.set_yticks(row_order)
        ax.set_yticklabels(yticklabels)

    if zlabel is not None:
        cbar.ax.set_ylabel(zlabel)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)


# From an example at:
# http://matplotlib.org/examples/api/histogram_path_demo.html
@_autoplot
def hist(x, bins='sturges', *, range=None, weights=None, density=False,
         facecolor='white', edgecolor='black',
         xlabel=None, ylabel=None,
         ax):
    n, bins = np.histogram(x, bins=bins, range=range, weights=weights,
                           density=density)

    # get the corners of the rectangles for the histogram
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n

    # we need a (numrects x numsides x 2) numpy array for the path helper
    # function to build a compound path
    XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T

    # get the Path object
    barpath = matplotlib.path.Path.make_compound_path_from_polys(XY)

    # make a patch out of it
    patch = matplotlib.patches.PathPatch(
        barpath, facecolor=facecolor, edgecolor=edgecolor)
    ax.add_patch(patch)

    # update the view limits
    view_x_eps = (right[-1] - left[0])/100
    ax.set_xlim(left[0] - view_x_eps, right[-1] + view_x_eps)
    ax.set_ylim(bottom.min(), top.max())

    if xlabel is None:
        xlabel = False
    if ylabel is None:
        if density:
            ylabel = 'Density'
        else:
            ylabel = 'Frequency'
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


def interpolate_quantiles(x, n):
    idx = np.array(np.floor(np.linspace(0, len(x)-1, num=n)), dtype='int')
    return x[idx]


@_autoplot
def qqplot(x, y, *, ax, diagonal=False, xlabel=None, ylabel=None):

    # Arguments
    if xlabel is None:
        xlabel = False
    if ylabel is None:
        ylabel = False

    # Calculations
    x = np.sort(np.array(x))
    y = np.sort(np.array(y))

    if len(x) > len(y):
        x = interpolate_quantiles(x, len(y))
    elif len(y) > len(x):
        y = interpolate_quantiles(y, len(x))

    # Plot
    if diagonal:
        ax.axis('equal')
        min_val = min(x[0], y[0])
        max_val = max(x[-1], y[-1])
        ax.plot([min_val, max_val], [min_val, max_val], ls="--")

    ax.scatter(x, y, c='none', edgecolors='k')

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
