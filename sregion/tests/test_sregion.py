import numpy as np
import astropy.units as u

from .. import SRegion, patch_from_polygon


def test_sregion():
    """
    Test SRegion object
    """
    # From arrays
    x = np.array([0, 0, 1, 1])
    y = np.array([0, 1, 1, 0])
    sr = SRegion(np.array([x, y]).T)

    assert(sr.N == 1)

    assert(np.allclose(sr.centroid[0], 0.5, rtol=1.e-3))

    assert(sr.area[0] == 1.0)

    # Converters
    assert(hasattr(sr.path[0], 'contains_point'))
    assert(hasattr(sr.shapely[0], 'boundary'))
    assert(hasattr(sr.patch()[0], 'get_fc'))

    regstr = ('polygon(0.000000,0.000000,0.000000,1.000000,' +
              '1.000000,1.000000,1.000000,0.000000)')
    assert(sr.region[0] == regstr)

    sr.label = 'test'
    assert(sr.region[0] == regstr + ' #  text={test}')

    # SQL string
    pstr = '((0.000,0.000),(0.000,1.000),(1.000,1.000),(1.000,0.000))'
    assert(sr.polystr(precision=3)[0] == pstr)
    
    # From s_region string
    pstr = ('POLYGON 0.000000 0.000000 0.000000 1.000000 ' +
            '1.000000 1.000000 1.000000 0.000000')

    assert(pstr == sr.s_region)

    snew = SRegion(pstr)
    assert(snew.area[0] == 1.0)

    # From polygon
    snew = SRegion(sr.shapely[0])
    assert(snew.area[0] == 1.0)
    
    # Compound regions
    x2 = np.array([0, 0, 1, 1]) + 2
    y2 = np.array([0, 1, 1, 0]) + 2
    s2 = SRegion(np.array([x2, y2]).T)

    un = sr.union(s2.shapely[0], as_polygon=True)

    assert(un.area == 2.0)

    # Multiple string
    comp = SRegion(' '.join([sr.s_region, s2.s_region]))
    assert(np.allclose(comp.area, 1.0, rtol=1.e-3))

    un = comp.union(as_polygon=True)
    assert(un.area == 2.0)

    # Intersects
    assert(sr.intersects(s2.shapely[0]) is False)

    x3 = np.array([0, 0, 1, 1]) + 0.5
    y3 = np.array([0, 1, 1, 0]) + 0.5
    s3 = SRegion(np.array([x3, y3]).T)

    assert(sr.intersects(s3.shapely[0]) is True)


def test_wrap():
    """
    """
    x = np.array([0, 0, 1, 1]) - 177
    y = np.array([0, 1, 1, 0])
    sr = SRegion(np.array([x, y]).T, wrap=True)

    assert(sr.N == 1)
    assert(np.allclose(sr.area, 1.))
    assert(np.allclose(sr.centroid, [360-177+0.5, 0.5]))


def test_circles():
    """
    Initialize from ``CIRCLE X Y R``
    """
    # CIRCLE string
    circ = SRegion('CIRCLE 5. 5. 1', ncircle=256)
    assert(np.allclose(circ.area, np.pi, rtol=1.e-3))
    assert(np.allclose(circ.centroid[0], 5., rtol=1.e-3))

    # R = 2, no units
    circ2 = SRegion('CIRCLE 5. 5. 2', ncircle=256)
    assert(np.allclose(circ2.area, 2**2*np.pi, rtol=1.e-3))
    assert(np.allclose(circ2.centroid[0], 5., rtol=1.e-3))

    # R = 1 arcsec
    cosd = np.cos(45./180*np.pi)
    csky = SRegion('CIRCLE 45. 45. 1"', ncircle=256)
    assert(np.allclose(csky.area, np.pi/3600.**2/cosd, rtol=1.e-3))
    assert(np.allclose(csky.sky_area(unit=u.arcsec**2),
                       np.pi*u.arcsec**2, rtol=1.e-3))

    # R = 1 arcmin
    cosd = np.cos(45./180*np.pi)
    csky = SRegion('CIRCLE 45. 45. 1\'', ncircle=256)
    assert(np.allclose(csky.area, np.pi/3600./cosd, rtol=1.e-3))
    assert(np.allclose(csky.sky_area(unit=u.arcsec**2),
                       3600*np.pi*u.arcsec**2, rtol=1.e-3))
    assert(np.allclose(csky.centroid[0], 45., rtol=1.e-3))

    # Sky buffer
    cosd = np.cos(45./180*np.pi)
    csky = SRegion('CIRCLE 45. 45. 0.001"', ncircle=256)
    csky.sky_buffer(1.)

    assert(np.allclose(csky.area, np.pi/cosd, rtol=1.e-3))
    assert(np.allclose(csky.sky_area(unit=u.deg**2),
                       np.pi*u.deg**2, rtol=1.e-3))


def test_stcs():
    """
    Test stripping values from STC-S specified strings
    """
    stcs = """Union ICRS ( Polygon 239.807341 -18.296691 239.803564 -18.300277 239.799786 -18.296691 239.803563
    -18.293105 Polygon 239.797826 -18.295944 239.794049 -18.299530 239.790272 -18.295944 239.794049
    -18.292358)"""
    
    sr = SRegion(stcs, verbose=False)
    
    assert(len(sr.xy) == 2)
    assert(np.allclose(sr.area, 2.709e-5, rtol=0.01))


def test_whitespace():
    """
    """
    pstr = 'CIRCLE 5. 5. 1'
    for str_i in [pstr, pstr.replace(' ', '   '), '  '+pstr, pstr+'\n']:
        sr = SRegion(str_i, ncircle=256)
        assert(sr.N == 1)
        assert(np.allclose(sr.area, np.pi, rtol=1.e-3))
        assert(np.allclose(sr.centroid[0], 5., rtol=1.e-3))

    pstr = ('POLYGON 0.000000 0.000000 0.000000 1.000000 ' +
            '1.000000 1.000000 1.000000 0.000000')

    for str_i in [pstr, pstr.replace(' ', '  '), '  '+pstr, pstr+'\n']:
        sr = SRegion(str_i)
        assert(sr.N == 1)
        assert(np.allclose(sr.area, 1.))
        assert(np.allclose(sr.centroid, 0.5))


def test_fromwcs():
    """
    Initialize from `astropy.wcs.WCS`
    """
    from astropy.io.fits import Header
    import astropy.wcs as pywcs

    # From WCS
    header_string = """
CRPIX1  =                  1.0                                                  
CRPIX2  =                  1.0                                                  
CRVAL1  =                 90.0                                                  
CRVAL2  =                  0.0                                                  
CD1_1   =             -0.00001                                                  
CD1_2   =                  0.0                                                  
CD2_1   =                  0.0                                                  
CD2_2   =              0.00001                                                  
NAXIS1  =                 1001                                                  
NAXIS2  =                 1001                                                  
CTYPE1  = 'RA---TAN'                                                            
CTYPE2  = 'DEC--TAN'                                                            
RADESYS = 'ICRS    '                                                            
EQUINOX =               2000.0                                                  
LATPOLE =                    0                                                  
LONPOLE =                180.0
    """

    head = Header.fromstring(header_string.replace('\n', ''))
    wcs = pywcs.WCS(head)
    sw = SRegion(wcs)
    pixel_area = np.abs((head['NAXIS1']-1)*head['CD1_1'] *
                        (head['NAXIS2']-1)*head['CD2_2'])

    assert(np.allclose(sw.area, pixel_area, atol=1.e-5))

    assert(np.allclose(sw.sky_area(unit=u.deg**2), pixel_area*u.deg**2))
    assert(np.allclose(sw.sky_area(unit=u.arcmin**2),
                       (pixel_area*u.deg**2).to(u.arcmin**2)))


def test_patch():
    """
    Test patch function
    """
    
    circ = SRegion('CIRCLE 5. 5. 1', ncircle=256)
    patch = patch_from_polygon(circ.shapely[0], fc='k')


def test_draw_patch():
    """
    Draw in matplotlib axis
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(1,1)
    ax.plot([-10,10], [-10,10])
    
    circ = SRegion('CIRCLE 5. 5. 1', ncircle=256)
    
    circ.add_patch_to_axis(ax, fc='r', alpha=0.5)
    
    plt.close('all')
    
