from setuptools import setup

version = "0.1"  # First version
version = "1.0" # works well, includes polystr

version_str = """#
__version__ = "{0}"\n""".format(version)

with open('sregion/version.py', 'w') as fp:
    fp.write(version_str)
    fp.close()

setup(name='sregion',
      version=version,
      description='Parsing of IVOA S_REGION strings',
      install_requires=['numpy', 'astropy', 'shapely', 'descartes'],
      author='Gabriel Brammer',
      author_email='gbrammer@gmail.com',
      url='https://github.com/gbrammer/sregion/',
      packages=['sregion', 'sregion/tests'],
      )
