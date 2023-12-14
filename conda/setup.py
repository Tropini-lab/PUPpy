from setuptools import setup, find_packages

setup(
    name='puppy',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'primer3-py>=0.6.1',
        'seaborn>=0.13.0',
        'colorama>=0.4.1'
    ],
    entry_points={
        'console_scripts': [
            'puppy-align=puppy.puppy_align:main',
            'puppy-primers=puppy.puppy_primers:main',
            'puppy-GUI=puppy.puppy_GUI:main'
        ]
    }
)
