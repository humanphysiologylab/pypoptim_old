# pypoptim
>Population-based algorithms for global optimization of the cardiac action potential models.

## Getting started
### Install
1. Download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
(optional) Install mamba and replace further `conda` calls by `mamba` calls:
    ```sh
    conda install mamba -n base -c conda-forge
    ```
2. Create environment (in this folder) and activate it:
```sh
conda env create -f environment.yml --prefix ./env
conda activate ./env
```
3. Install `pypoptim` package:
```sh
pip install -e .
```

### Test
```sh
pytest --pyargs pypoptim
```

### Usage
See [examples](./examples) folder.

### Pre-commit
```sh
pre-commit install  # only once
git add FILES
pre-commit run
git commit -m "MESSAGE"
```

### ToDo:
Split this repository into three or so:
- pypoptim itself
- models (C and Python variants)
- mpi_scripts/
- add tests for all
