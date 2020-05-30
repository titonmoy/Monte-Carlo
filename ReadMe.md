# Monte Carlo Simulation of Light-Tissue Interaction

This is a class project for Contemporary Approaches to Mathematical and Quantitative Biology. It contains the python implementation for Monte Carlo simulation of light-tissue interaction. This is a simple implementaion considering only a single layer of tissue in a surrounding medium. The project mainly focuses on simulating the light absorbing charactericts of melanin in human epidermis.

The algorithm used in this project is based on the work of Wang et al. [1]. The optical properties of epidermis is taken from [2].

The following libraries have been used for the python implementation
- math (for numeric computation)
- numpy (for array variables)
- matplotlib (for graphical representation of simulation results)
- tabulate (for tabular representation of simulation results)

The conda environment used for the project can be recreated using [environment.yml](environment.yml).

A detailed description of the project is given in the [Melanin.ipynb](Melanin.ipynb) notebook.


## References
1. Wang, L., Jacques, S.L. and Zheng, L., 1995. MCML—Monte Carlo modeling of light transport in multi-layered tissues. Computer methods and programs in biomedicine, 47(2), pp.131-146.
2. Bashkatov, A.N., Genina, E.A. and Tuchin, V.V., 2011. Optical properties of skin, subcutaneous, and muscle tissues: a review. Journal of Innovative Optical Health Sciences, 4(01), pp.9-38.
