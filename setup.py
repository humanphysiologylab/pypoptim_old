from setuptools import setup

setup(
    name='pypoptim',
    version='2.0.1',
    packages=["pypoptim"],
    package_dir={"": "src"},
    url='https://github.com/humanphysiologylab/pypoptim',
    author='Andrey Pikunov',
    author_email='pikunov@phystech.edu',
    install_required=['numpy',
                      'pandas',
                      'scikit-learn',
                      'numba',
                      'pytest']
)
