# Step 1: Install the R package dependencies:

``` r
if(!require(pacman)) install.packages("pacman")

if(!require(RcppPlanc)) 
  install.packages('RcppPlanc', repos = 'https://welch-lab.r-universe.dev')

pacman::p_load(tidyverse, janitor, remotes, reticulate, 
               Seurat, SeuratWrappers, leidenAlg,
               conos, harmony, rliger)

if(!require(scMC)) install_github("amsszlh/scMC")

if(!require(scib)) install_github("mthimma/scib")
```


# Step 2A: Install packages in conda environment for Windows OS

You can execute the following commands in R:

``` r
if( !file.exists(conda_binary()) ) install_miniconda()

system("curl -LJO https://raw.githubusercontent.com/mthimma/scib/refs/heads/main/scib_env_windows.yaml")

conda_create(envname = "scib", environment = "scib_win.yaml")
```


# Step 2B: Install packages in conda environment for Linux and MacOS

If you do not have conda, you can install it via miniforge (see below) or miniconda. 

``` bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh
```

Then restart the terminal or open a new terminal for the changes to take
effect. Next we create a conda environment with the required softwares:

``` bash
curl -LJO https://raw.githubusercontent.com/mthimma/scib/refs/heads/main/scib_env_linux.yaml

conda env create --name scib --file=scib_env_linux.yaml
```


# Step 3: Set python path before each usage

R typically uses the systems default Python. Therefore, we need to explicitly
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


[Example of running and evaluating one data integration method](scib_one_method.md)

[Example of running and evaluating **multiple** data integration methods](scib_multiple_methods.md)

