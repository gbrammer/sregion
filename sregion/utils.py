"""
Utilities
"""

import numpy as np
from .sregion import SRegion


def concave_hull(
    points, smooth=0.08, sky=True, ref_density=250, scale_power=0.5, **kwargs
):
    """
    Compute the approximate concave hull of a set of points.

    Adapted from
    https://gist.github.com/jclosure/d93f39a6c7b1f24f8b92252800182889
    by GitHub user @jclosure

    Parameters
    ----------
    points : (N, 2) array
        Point list

    smooth : float
        Smoothing parameter

    sky : bool
        Treat ``points`` as sky coordinates in decimal degrees

    ref_density : float, None
        Reference sky surface density (N per square arcmin).
        If provided, scale the input ``smooth`` by the surface density of ``points``
        relative to ``ref_density``.

    scale_power : float
        Power on the relative surface density

    Returns
    -------
    output_hull : `~sregion.sregion.SRegion`
        Concave hull of the input ``points`` array

    scale_smooth : float
        Smoothing parameter used, perhaps rescaled

    """
    from shapely.ops import unary_union, polygonize
    from shapely import geometry
    from scipy.spatial import Delaunay, ConvexHull

    if points.shape[0] == 2:
        points = points.T

    if len(points) < 4:
        # When you have a triangle, there is no sense in computing an smooth
        # shape.
        return geometry.MultiPoint(list(points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    center = np.median(points, axis=0)
    cosd = np.array([np.cos(center[1] / 180 * np.pi), 1.0]) if sky else np.array([1, 1])

    if (ref_density is not None) & (sky):
        hull = ConvexHull(points)
        sr = SRegion(points[hull.vertices, :])
        scale_smooth = (
            smooth
            * (ref_density / (len(points) / sr.sky_area()[0].value)) ** scale_power
        )
    else:
        scale_smooth = smooth

    if sky:
        cosd *= 60

    coords = points - center
    coords *= cosd

    tri = Delaunay(coords)
    edges = set()
    edge_points = []

    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]

        # Lengths of sides of triangle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        # Semiperimeter of triangle
        s = (a + b + c) / 2.0

        # Area of triangle by Heron's formula
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        # radius = np.sqrt(area / np.pi)

        # Here's the radius filter.
        # print circum_r
        # if circum_r < 1.0/smooth:
        if circum_r < scale_smooth:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)

    edge_points = [p / cosd + center for p in edge_points]

    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    output_hull = SRegion(unary_union(triangles))

    return output_hull, scale_smooth


def get_sky_footprint(ra, dec, make_figure=False, **kwargs):
    """
    Get the footprint of a list of sky coordinates

    Parameters
    ----------
    ra, dec : array-like
        Sky coordinates, decimal degrees

    make_figure : bool
        Make a diagnostic figure

    Returns
    -------
    catalog_hull : `~sregion.sregion.SRegion`
        Footprint

    fig : `~matplotlib.Figure`
        Figure object

    """
    import matplotlib.pyplot as plt

    points = np.array([ra, dec]).T
    catalog_hull, scale_smooth = concave_hull(points, **kwargs)

    if make_figure:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.scatter(*points.T, color="0.5", alpha=0.05, marker=".", zorder=100)
        label = f"concave hull, area={catalog_hull.sky_area()[0]:.1f}"
        catalog_hull.add_patch_to_axis(ax, fc="tomato", ec="None", alpha=0.2)
        catalog_hull.add_patch_to_axis(
            ax, fc="None", ec="tomato", alpha=0.5, label=label
        )
        ax.set_xlim(*ax.get_xlim()[::-1])
    else:
        fig = None

    return catalog_hull, fig


def moc_patch_union(moc, wcs, buffer=0.001, optimize=True, sky=False):
    """
    Compute SRegion from union of healpix pixels in a Multi-object Coverage
    (MOC) map
    
    Parameters
    ----------
    moc : `mocpy.moc.MOC`
        Multi-object Coverage object

    wcs : `astropy.wcs.WCS`
        WCS to evaluate from the healpix pixels

    buffer : float
        Small buffer to avoid numerical issues with exact healpix region
        overlaps

    optimize : bool
        Passed to `mocpy.moc.plot.fill.build_plotting_moc`.

    sky : bool
        If True, return region(s) in celestial coordinates.  Otherwise,
        return regions in the pixels of ``wcs``

    Returns
    -------
    sr : `~sregion.sregion.SRegion`
        Region object

    """
    from matplotlib.path import Path
    from shapely import Polygon

    from mocpy.moc.plot.fill import compute_the_patches, build_plotting_moc

    moc_to_plot = build_plotting_moc(moc=moc, wcs=wcs, optimize=optimize)

    patches = compute_the_patches(moc=moc_to_plot, wcs=wcs)

    poly = None
    for p in patches:
            
        last = np.append(-1, np.where(p[1] == Path.CLOSEPOLY)[0])

        if sky:
            xy = np.array(wcs.all_pix2world(*p[0].T, 0)).T
        else:
            xy = p[0]
        
        for j in range(len(last)-1):
            sl = slice(last[j]+1, last[j+1])
            poly_i = Polygon(xy[sl,:]).buffer(buffer)
            if poly is None:
                poly = poly_i
            else:
                poly = poly.union(poly_i)

    sr = SRegion(poly, wrap=(sky))

    return sr

