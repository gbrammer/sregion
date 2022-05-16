![python package](https://github.com/gbrammer/sregion/actions/workflows/python-package.yml/badge.svg)

# sregion
Parsing of IVOA S_REGION strings

The STS-C formalism is described at http://www.ivoa.net/Documents/latest/STC-S.html, though [it seems](https://github.com/astropy/regions/issues/21) that it was never adopted as an official standard.  Nevertheless, the `s_region` strings do seem to have been adopted as a sort of pseudostandard in [IVOA-compliant](https://wiki.ivoa.net/twiki/bin/view/IVOA/DCPToolsFITS) datasets / databases.

[`astropy-regions`](https://github.com/astropy/regions) would probably be a better place to put this, but I'm not interested in all of the full astropy coordinate compatibility for now.
