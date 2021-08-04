from setuptools import setup

setup(
    name="pypoptim",
    version="2.1.0",
    packages=["pypoptim"],
    package_dir={"": "src"},
    url="https://github.com/humanphysiologylab/pypoptim",
    author="Andrey Pikunov",
    author_email="pikunov@phystech.edu",
    install_requires=["numpy", "pandas", "scikit-learn", "numba", "pytest"],
)
