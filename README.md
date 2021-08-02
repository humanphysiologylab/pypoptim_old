# pypoptim
Population-based algorithms for global optimization of the cardiac action potential models.

## Getting started
1. Download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
2. Create environment (in this folder) and activate it:
```sh
conda env create -f environment.yml --prefix ./env
conda activate ./env
```
3. Install `pypoptim` package:
```sh
pip install -e .
```

4. Test:
```sh
pytest --pyargs pypoptim
```

5. Run:
```sh
mpirun -n 2 python mpi_scripts/path/to/some/mpi_script.py configs/path/to/some/config.json
```

### ToDo:
Split this repository into three or so:
- pypoptim itself
- models (C and Python variants)
- mpi_scripts/
- add tests for all
