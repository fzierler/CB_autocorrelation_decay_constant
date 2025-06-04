Clone the repository including all submodules:

```
git clone --recurse-submodules https://github.com/fzierler/CB_autocorrelation_decay_constant
```

Then copy the directories `raw_data`, `external_data` and `metadata` into the repository. 

Activate the conda environment in the repository:

```
cd CB_autocorrelation_decay_constant
conda env create --name wall_decay_constant --file=environment.yml
conda activate wall_decay_constant
```

The julia programming language is not available in conda-forge. If it is not installed on your system, you can install it with the `juliaup` utility included in the provided conda environment using

```
export JULIA_DEPOT_PATH=$CONDA_PREFIX/.julia
juliaup add 1.11.5
```

The analysis can then be performed by executing

```
bash main.sh
```