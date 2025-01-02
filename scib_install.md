
# Step 1: Install the R package dependencies:

``` r
if(!require(pacman)) install.packages("pacman")

if(!require(RcppPlanc)) 
  install.packages('RcppPlanc', repos = 'https://welch-lab.r-universe.dev')

pacman::p_load(tidyverse, Seurat, SeuratWrappers,
               remotes, reticulate, leidenAlg,
               conos, harmony, rliger)

if(!require(scMC)) install_github("amsszlh/scMC")

if(!require(scib)) install_github("mthimma/scib")
```

# Step 2: Install packages in conda environment for Windows OS

``` r
if( !file.exists(conda_binary()) ) install_miniconda()

system('curl -LJO https://raw.githubusercontent.com/mthimma/scib/refs/heads/main/scib_win.yaml')

conda_create(envname = "scib", environment = "scib_win.yaml")
```

# Step 2: Install packages in conda environment for Linux and MacOS

If conda does not exist, you can install miniforge or miniconda.

``` bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh
```

Then restart the terminal or open a new terminal for the changes to take
effect.

``` bash
curl -LJO https://raw.githubusercontent.com/mthimma/scib/refs/heads/main/scib_linux.yaml

conda env create --name scib --file=scib_linux.yaml
```

# Set python path before usage

R sometimes can source the system

R uses the systems default Python. Therefore, we need to explicitly
specify the location of the python used for `scib`.

``` r
python_binary <- conda_list() %>%
  filter(name == "scib") %>%
  pull(python) %>% 
  normalizePath(winslash = "/")

python_binary
# "C:/Users/mthimma/AppData/Local/r-miniconda/envs/scib/python.exe"

use_python(python_binary, required = TRUE)
```

`IntegrateLayers()` function also requires the path to conda enviroment
which is the folder where the python binary is located in.

``` r
conda_env <- dirname(python_binary)

conda_env
# "C:/Users/mthimma/AppData/Local/r-miniconda/envs/scib/"
```