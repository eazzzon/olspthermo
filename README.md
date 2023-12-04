# Extended calibration of olivine-spinel Al exchange thermometers (OSAT)

Python scripts for the new models of the olivine-spinel Al exchange thermometers (OSAT). 

## Installation
```pip install https://github.com/eazzzon/olspthermo.git```

## Requirement
`numpy`

`pandas`

`scipy`

`periodictable`

## Uninstallation
```pip uninstall olspthermo```

## Example

To quick run, navigate to the `example` folder, call the python script as:

```python
python example_calculation.py
```

or run through the jupiter file.

**DON'T** change the header/column name of the example excel template files. 

## Monte Carlo simulation for error progagation

We also provide a MC approach to progagate the analytical errors to the temperature estimation. an example is given by calling 

```
python example_mc.py
```

or run in the jupyter example.

## Excel Spreadsheet & Jupyter notebook

We also provide an excel spreadsheet (`OSAT.xlsx` ) and example notebook files to calculate temperature with all models. Error progagation will have to go with MC in python

## Results

| String name     | property                                                     |
| --------------- | ------------------------------------------------------------ |
| t_thermo        | Temperature estimates from Eq. 5                             |
| t_kdal          | Temperature estimates from Eq. 7                             |
| t_kdcr          | Temperature estimates from Eq. 8                             |
| t_coogan        | Temperature estimates from Coogan et al. (2014)              |
| t_zThermoAlCr   | Temperature chose from comparison of Eqs. 5, 7 & 8           |
| err_zThermoAlCr | Temperature uncertainty chose from comparison of Eqs. 5, 7 & 8 |
| t_zThermoAl     | Temperature chose from comparison of Eqs. 5 & 7              |
| err_zThermoAl   | Temperature uncertainty chose from comparison of Eqs. 5 & 7  |
| t_zThermoCr     | Temperature chose from comparison of Eqs. 5 & 8              |
| err_zThermoCr   | Temperature uncertainty chose from comparison of Eqs. 5 & 8  |
| z_thermo_cr     | Z-score between Eq. 5 & 8                                    |
| z_thermo_al     | Z-score between Eq. 5 & 7                                    |
| p_thermo_cr     | p-value for z_thermo_cr                                      |
| p_thermo_al     | p-value for p_thermo_al                                      |

## How to cite

Zhang Y, Namur O, Li W, Shorttle O, Gazel E, Jennings E, Thy P, Grove TL, Charlier B (2023). An extended calibration of the olivine–spinel aluminum exchange thermometer: Application to the melting conditions and mantle lithologies of large igneous provinces. *Journal of Petrology*, *64*(11), egad077. https://doi.org/10.1093/petrology/egad077

## Contact

If you find any errors or would like to incorporate the model in your work, don't hesitate to contact me at yishen.zhang@kuleuven.be
