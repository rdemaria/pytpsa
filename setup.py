from setuptools import setup, find_packages

setup(
    name="pytpsa",
    version="0.0.0",
    description="A Python TPSA",
    author="Riccardo De Maria",
    author_email="riccardo.de.maria@cern.ch",
    url='https://github.com/rdemaria/pytpsa',
    packages=find_packages(),
    install_requires=['numpy'],
)
