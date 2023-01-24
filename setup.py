from setuptools import setup

setup(
    name='molgroups',
    version='0.1.0',
    packages=['molgroups', 'molgroups.infotheory', 'molgroups.support'],
    url='https://github.com/criosx/molgroups',
    license='GNU GENERAL PUBLIC LICENSE, Version 3',
    author='Frank Heinrich, David Hoogerheide, Alyssa Thomas',
    author_email='mail@frank-heinrich.net',
    description='Molecular Modeling for Scattering Data Analysis',
    install_requires=[
        "bumps", "refl1d", "periodictable", "scipy", "numpy", "matplotlib", "pandas", "scikit-learn", "shapely",
        "sasmodels", "sasdata", "dill", "gpcam"
    ]
)
