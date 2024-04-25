from setuptools import setup, find_packages

setup(
    name='puppy',
    use_scm_version={
            "local_scheme": "node-and-timestamp",
            "fallback_version": "1.1.2",
        },
    packages=find_packages(),
    install_requires=[
        'primer3-py>=0.6.1',
        'seaborn>=0.13.0',
        'colorama>=0.4.1'
    ],
    setup_requires=['setuptools_scm'],
    entry_points={
        'console_scripts': [
            'puppy-align=puppy.puppy_align:main',
            'puppy-primers=puppy.puppy_primers:main',
            'puppy-GUI=puppy.puppy_GUI:main'
        ]
    }
)
