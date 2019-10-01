from setuptools import setup

setup(
    name='pyread_eagle',
    version='0.1',
    description="Pure-python port of J. Helly's read_eagle.",
    url='',
    author='Kyle Oman',
    author_email='koman@astro.rug.nl',
    license='',
    packages=['pyread_eagle'],
    install_requires=['numpy', 'h5py'],
    include_package_data=True,
    zip_safe=False
)
