# pypoptim
Population-based algorithms for global optimization

How to do something:
```sh
conda env create -f environment.yml --prefix ./env
conda activate ./env
mpirun -n 2 python mpi_scripts/cardio/mpi_script.py configs/kernik_clancy/config_syn_IK1_x2_biphasic.json
conda deactivate  # if you need
```

How to remove the environment:
```sh
conda remove --prefix ./env --all
```
