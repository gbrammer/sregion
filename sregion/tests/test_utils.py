def test_concave_hull():
    """
    Test concave hull function
    """
    import numpy as np

    from ..sregion import SRegion
    from ..utils import concave_hull

    # Demo concave polygon
    xy = np.array(
        [
            [-0.05340819, 0.01697288],
            [0.00633657, 0.08214911],
            [0.07965975, -0.01697297],
            [0.01991492, -0.00746813],
            [0.00090514, -0.03734045],
            [-0.05340819, -0.03734045],
        ]
    )

    nominal_area = 34.26
    xy[:, 0] += 10

    sr = SRegion(xy, wrap=True)

    # Random points within the test polygon
    np.random.seed(1)
    rnd = np.random.rand(2000, 2)

    mi, ma = xy.min(axis=0), xy.max(axis=0)
    si = ma - mi
    rnd = rnd * si + mi

    inside = sr.path[0].contains_points(rnd)
    rnd = rnd[inside, :]

    assert np.allclose(
        [area.value for area in sr.sky_area()], nominal_area, rtol=1.0e-2
    )

    # Test that get the same answer at different locations on the sky
    for ra in [45, 150]:
        for dec in [-80, -10, 0, 30, 50]:
            cosd = np.cos(dec / 180 * np.pi)

            xy_i = xy / np.array([cosd, 1.0]) + np.array([ra, dec])
            sr_i = SRegion(xy_i, wrap=True)

            assert np.allclose(
                [area.value for area in sr_i.sky_area()], nominal_area, rtol=1.0e-2
            )

            rnd_i = rnd / np.array([cosd, 1.0]) + np.array([ra, dec])

            hull_i, _alpha = concave_hull(
                rnd_i,
                alpha=0.08,
                sky=True,
                ref_density=250,
                scale_power=0.5,
            )

            print(ra, dec, hull_i.sky_area()[0], sr_i.sky_area()[0], nominal_area)

            assert hull_i.sky_area()[0] < sr_i.sky_area()[0]

            assert np.allclose(
                [area.value for area in hull_i.sky_area()], nominal_area, rtol=5.0e-2
            )
