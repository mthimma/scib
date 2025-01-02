# scib

The goal of scib (single-cell integration and benchmarking) is to
provide a unified way to perform data integration and benchmarking when
the ground truth is available.


# Data integration methods implemented via scib and other packages

| method<sup>*</sup> | Package for wrapper | Publication                                     | DOI                                           |
|---         |---                  |---                                              |---                                            |
| CCA        | Seurat              | Stuart T. et al. Cell (2019)                    | https:://doi.org/10.1016/j.cell.2019.05.031   |
| RPCA       | Seurat              | Stuart T. et al. Cell (2019)                    | https:://doi.org/10.1016/j.cell.2019.05.031   |
| Harmony    | Seurat              | Korsunsky, I. et al. Nat Methods (2019)         | https://doi.org/10.1038/s41592-019-0619-0     |
| FastMNN    | SeuratWrapper       | Haghverdi, L. et al. Nat. Biotechnol. (2018)    | https://doi.org/10.1038/s41421-019-0114-x     |
| SCVI       | SeuratWrapper       | Lopez, R.et al. Nat Methods (2018)              | https://doi.org/10.1038/s41592-018-0229-2     |
| BBKNN      | scib (this package) | Krzysztof Polański et al. Bioinformatics (2020) | https://doi.org/10.1093/bioinformatics/btz625 |
| CONOS      | scib (this package) | Barkas, N et al. Nat Methods (2019)             | https://doi.org/10.1038/s41592-019-0466-z     |
| LIGER      | scib (this package) | Liu, J.et al. Nat. Protocols (2020)             | https://doi.org/10.1038/s41596-020-0391-8     |
| SCANORAMA  | scib (this package) | Hie, B et al. Nat Biotechnology (2019)          | https://doi.org/10.1038/s41587-019-0113-3     |
| SCMC       | scib (this package) | Zhang, L. et al. Genome Biology (2021)          | https://doi.org/10.1186/s13059-020-02238-2    |


<sup>*</sup> The specific method can be called by appending the word "Integration" to the method name (e.g. HarmonyIntegration) for `IntegrateLayers2()` in `scib` package.


# Guides

[Installation of scib and dependencies](scib_install.md)

[Example of running and evaluating one data integration method](scib_one_method.md)

[Example of running and evaluating *multiple* data integration methods](scib_multiple_methods.md)


