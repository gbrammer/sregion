![python package](https://github.com/gbrammer/sregion/actions/workflows/python-package.yml/badge.svg)

# sregion
Parsing of IVOA S_REGION strings

The STS-C formalism is described at http://www.ivoa.net/Documents/latest/STC-S.html, though [it seems](https://github.com/astropy/regions/issues/21) that it was never adopted as an official standard.  Nevertheless, the `s_region` strings do seem to have been adopted as a sort of pseudostandard in [IVOA-compliant](https://wiki.ivoa.net/twiki/bin/view/IVOA/DCPToolsFITS) datasets / databases.

[`astropy-regions`](https://github.com/astropy/regions) would probably be a better place to put this, but I'm not interested in all of the full astropy coordinate compatibility for now.

## Examples

```python
>>> import numpy as np
>>> from sregion import SRegion

#
# Polygon string
#
>>> sr = SRegion('POLYGON 0.0 0.0 0.0 1.0 1.0 1.0 1.0 0.0')
>>> print(sr.area)
[1.0]
>>> print(sr.centroid)
[array([0.5, 0.5])]

#
# Circle string
#
>>> for i in range(4,10):
>>>     sr = SRegion('CIRCLE 10 10 1', ncircle=2**i)
>>>     print(f'ncircle={2**i:>3} {sr.area[0]/np.pi:.5f} {sr.centroid[0]}')
ncircle= 16 0.97450 [10. 10.]
ncircle= 32 0.99359 [10. 10.]
ncircle= 64 0.99839 [10. 10.]
ncircle=128 0.99960 [10. 10.]
ncircle=256 0.99990 [10. 10.]
ncircle=512 0.99997 [10. 10.]

# Circle with radius in angular units
>>> import astropy.units as u
>>> sr = SRegion('CIRCLE 10 10 1"', ncircle=256)
>>> print(f'{sr.sky_area(unit=u.arcsec**2)[0]:.5f}')
3.14128 arcsec2

#
# From WCS objects
#
>>> from astropy.wcs import WCS
>>> wcs = WCS()
>>> wcs.pixel_shape = [601,601]
>>> wcs.wcs.cdelt *= 0.1/3600
>>> wcs.wcs.crpix[1] = 300
>>> wcs.wcs.crval = [0,0]
>>> print(SRegion(wcs).sky_area())
[<Quantity 1. arcmin2>]

#
# From arrays
#
>>> x = np.array([0, 0, 1, 1])
>>> y = np.array([0, 1, 1, 0])
>>> sr = SRegion(np.array([x, y]).T)
>>> print(sr.area)
[1.0]
>>> print(sr.centroid)
[array([0.5, 0.5])]

# 
# To s_region string
#
>>> print(sr.s_region)
POLYGON 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 0.000000

#
# To matplotlib path object(s)
#
>>> print(sr.path[0].contains_point([0.5, 0.5]))
True
>>> print(sr.path[0].contains_points([[0.5, 0.5], [2.0, 2.0]]))
[ True False]

#
# To matplotlib patch(es)
#
>>> import matplotlib.pyplot as plt
>>> fig, ax = plt.subplots(1,1,figsize=(2,2))
>>> for p in sr.patch(alpha=0.5, fc='r'):
>>>     ax.add_patch(p)
>>> ax.set_xlim(-1, 2)
>>> ax.set_ylim(*ax.get_xlim())
>>> ax.grid()

#
# To shapely polygons
# 
>>> sr.shapely
[<shapely.geometry.polygon.Polygon at 0x18055b910>]

#
# To DS9 region(s)
#
>>> for r in sr.region:
>>>    print(r)
polygon(0.000000,0.000000,0.000000,1.000000,1.000000,1.000000,1.000000,0.000000)

>>> sr.ds9_properties = 'color=red width=2'
>>> sr.label = 'my_group'
>>> for r in sr.region:
>>>    print(r)
polygon(0.000000,0.000000,0.000000,1.000000,1.000000,1.000000,1.000000,0.000000) # color=red width=2 text={my_group}

    
```
