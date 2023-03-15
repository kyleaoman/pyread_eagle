from setuptools import setup
from pyread_eagle.__version__ import __version__

setup(
    name="pyread_eagle",
    version=__version__,
    description="Pure-python port of J. Helly's read_eagle.",
    url="https://github.com/kyleaoman/pyread_eagle",
    author="Kyle Oman",
    author_email="kyle.a.oman@durham.ac.uk",
    license="GNU GPL v3",
    packages=["pyread_eagle"],
    install_requires=["numpy", "h5py"],
    include_package_data=True,
    zip_safe=False,
)
