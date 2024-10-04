from setuptools import setup

setup(
    name='molgroups',
    version='0.2.0',
    packages=['molgroups', 'molgroups.refl1d_interface'],
    url='https://github.com/criosx/molgroups',
    license='MIT License',
    author='Frank Heinrich, David Hoogerheide, Alyssa Thomas',
    author_email='mail@frank-heinrich.net',
    description='Molecular Modeling for Scattering Data Analysis',
    install_requires=[
        "bumps", "refl1d", "periodictable", "scipy", "numpy", "matplotlib", "pandas", "scikit-learn",
        "sasmodels", "sasdata", "dill",
    ]
)
