import numpy as np
import astropy.units as u
import astropy.wcs as pywcs

class SRegion(object):
    def __init__(self, inp, label=None, **kwargs):
        """
        Helper class for parsing an S_REGION strings and general polygon 
        tools
        
        Parameters
        ----------
        inp : str, (M,2) array, `~astropy.wcs.WCS`, `shapely.geometry.polygon.Polygon`
            
            Can be a "S_REGION" string, an array, or an `astropy.wcs.WCS` 
            object from which the corners will be extracted with 
            `~astropy.wcs.WCS.calc_footprint`.
        
        label : str
            Optional label attached to regions and patches
            
        """
        if isinstance(inp, str):
            self.xy = self._parse_sregion(inp, **kwargs)

        elif hasattr(inp, 'sum'):
            # NDarray
            sh = inp.shape
            if inp.ndim != 2:
                raise ValueError(f'Input shape {sh} is not (M,2)')
            else:
                if inp.shape[1] != 2:
                    if inp.shape[0] != 2:
                        raise ValueError(f'Input shape {sh} is not (M,2)')
                    else:
                        inp = inp.T

            self.xy = [inp]

        elif isinstance(inp, list):
            self.xy = inp

        elif isinstance(inp, pywcs.WCS):
            self.xy = [inp.calc_footprint()]

        elif hasattr(inp, 'buffer'):
            # Shapely polygon
            if hasattr(inp, '__len__'):
                self.xy = [np.array(p.boundary.xy).T for p in inp]
            else:
                self.xy = [np.array(inp.boundary.xy).T]

        else:
            raise IOError('input must be ``str``, ``list``, or ``np.array``')

        self.inp = inp
        self.ds9_properties = ''
        self.label = label


    @staticmethod   
    def _parse_sregion(sregion, ncircle=32, wrap=False, **kwargs):
        """
        Parse an S_REGION string with CIRCLE or POLYGON
        """

        from astropy.coordinates import Angle
        import astropy.units as u

        if hasattr(sregion, 'decode'):
            decoded = sregion.decode('utf-8').strip().upper()
        else:
            decoded = sregion.strip().upper()

        polyspl = decoded.replace('POLYGON','xxx').replace('CIRCLE','xxx')
        polyspl = polyspl.split('xxx')

        poly = []
        for pp in polyspl:
            if not pp:
                continue

            pp = pp.replace('(','').replace(')','')

            if ',' in pp:
                spl = pp.strip().split(',')
            else:
                spl = pp.strip().split()

            for ip, p in enumerate(spl):
                # Find index of first float
                try:
                    pf = float(p)
                    break
                except ValueError:
                    continue
                            
            if len(spl[ip:]) == 3:
                # Circle
                x0, y0 = np.cast[float](spl[ip:-1])
                cosd = np.cos(y0/180*np.pi)

                r0 = spl[-1]
                if r0.endswith('"'):
                    scl = float(r0[:-1])/3600
                elif r0.endswith('\''):
                    scl = float(r0[:-1])/60
                else:
                    try:
                        scl = float(r0)
                    except ValueError:
                        scl = 1.
                        
                    cosd = 1.

                _theta = np.linspace(0, 2*np.pi, ncircle+1)[:-1]
                _xc = np.cos(_theta)
                _yc = np.sin(_theta)

                poly_i = np.array([_xc/cosd*scl+x0, _yc*scl+y0]).T
            else:
                poly_i = np.cast[float](spl[ip:]).reshape((-1,2))

            if wrap:
                ra = Angle(poly_i[:,0]*u.deg).wrap_at(360*u.deg).value
                poly_i[:,0] = ra
                
            if len(poly_i) < 2:
                continue

            poly.append(poly_i)

        return poly


    @property 
    def N(self):
        return len(self.xy)


    @property
    def centroid(self):
        return [np.mean(fp, axis=0) for fp in self.xy]


    def sky_buffer(self, buffer_deg):
        """
        Buffer polygons accounting for cos(dec)
        """
        from shapely.geometry import Polygon
        
        new_xy = []
        for xy, c in zip(self.xy, self.centroid):
            cosd = np.array([1, np.cos(c[1]/180*np.pi)])
            p = Polygon((xy - c)*cosd).convex_hull.buffer(buffer_deg)
            #p = Polygon(np.array(p.boundary.xy).T).buffer(buffer_deg)
            xyp = np.array(p.boundary.xy).T
            new_xy.append((xyp/cosd)+c)
        
        self.xy = new_xy


    @property
    def path(self):
        """
        `~matplotlib.path.Path` object
        """
        import matplotlib.path
        return [matplotlib.path.Path(fp) for fp in self.xy]


    @property
    def shapely(self):
        """
        `~shapely.geometry.Polygon` object.
        """
        from shapely.geometry import Polygon
        return [Polygon(fp).convex_hull for fp in self.xy]


    @property 
    def area(self):
        """
        Area of shapely polygons
        """
        return [sh.area for sh in self.shapely]


    def sky_area(self, unit=u.arcmin**2):
        """
        Assuming coordinates provided are RA/Dec degrees, compute area
        """
        cosd = np.cos(self.centroid[0][1]/180*np.pi)
        return [(sh.area*cosd*u.deg**2).to(unit)
                for sh in self.shapely]


    def patch(self, **kwargs):
        """
        `~descartes.PolygonPatch` object
        """
        from descartes import PolygonPatch
        if 'label' not in kwargs:
            kwargs['label'] = self.label

        return [PolygonPatch(p, **kwargs) for p in self.shapely]


    def get_patch(self, **kwargs):
        """
        `~descartes.PolygonPatch` object
        
        ** Deprecated, use patch **
        """
        from descartes import PolygonPatch
        return self.patch(**kwargs)


    def union(self, shape=None, as_polygon=False):
        """
        Union of self and `shape` object.  If no `shape` provided, then 
        return union of (optional) multiple components of self
        """
        if shape is None:
            un = self.shapely[0]
        else:
            un = shape
            
        for s in self.shapely:
            un = un.union(s)
        
        if as_polygon:
            return un
        else:
            return SRegion(un)


    def intersects(self, shape):
        """
        Union of self and `shape` object
        """
        test = False
        for s in self.shapely:
            test |= s.intersects(shape)

        return test


    @property
    def region(self):
        """
        Polygon string in DS9 region format
        """
        pstr = 'polygon({0})'
        if hasattr(self, 'ds9_properties'):
            tail = '{0}'.format(self.ds9_properties)
        else:
            tail = ''

        if hasattr(self, 'label'):    
            if self.label is not None:
                tail += ' text={xx}'.replace('xx', self.label)

        if tail:
            tail = ' # '+tail

        return [pstr.format(','.join([f'{c:.6f}' for c in fp.flatten()]))+tail
                for fp in self.xy]


    @property
    def s_region(self):
        """
        Polygon as VO s_region
        """
        pstr = 'POLYGON {0}'
        polys = [pstr.format(' '.join([f'{c:.6f}' for c in fp.flatten()]))
                for fp in self.xy]
        return ' '.join(polys)