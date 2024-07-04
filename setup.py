from setuptools import setup, find_packages

setup(
    name='tcr_antigen_prediction',
    version='0.0.0',
    author='Benjamin McMaster',
    author_email='benjamin.mcmaster@rdm.ox.ac.uk',
    packages=find_packages(where='src', include=['tcr_antigen_prediction']),
    package_dir={"": "src"},
    install_requires=[],
    extras_require={
        'develop': ['flake8', 'pylint'],
    },
)
