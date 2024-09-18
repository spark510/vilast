from setuptools import setup, find_packages

setup(
    name="vilast",
    version="0.1.0",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    description="A taxonomy-based Python package",
    author="Your Name",
    author_email="your.email@example.com",
    install_requires=[
        'lxml',  # Add lxml as a required dependency
        'pandas',
        'sqlite3'
    ],
    entry_points={
        'console_scripts': [
            'vilast=vilast.lpb_virus:main',  # Replace this with your entry point
        ],
    },
)