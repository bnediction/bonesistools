from setuptools import setup, find_packages

NAME = 'bonesis_tools'
VERSION = '9999'

setup(name=NAME,
    version=VERSION,
    keywords="computational systems biology",
    packages = find_packages(),
    package_data = {"databases": ['.*.tsv']}
)


